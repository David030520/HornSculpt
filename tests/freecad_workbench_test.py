"""Run inside a separate FreeCAD GUI with temporary user/system configuration.

The launcher supplies HORNSCULPT_TEST_SOURCE and HORNSCULPT_TEST_RESULT.
Only the test's temporary discovery descriptor is created in shared app data.
"""

import importlib.util
import base64
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import traceback
import uuid

import FreeCAD as App
import FreeCADGui as Gui
import Part
from PySide import QtCore, QtNetwork


source = Path(os.environ["HORNSCULPT_TEST_SOURCE"])
result_path = Path(os.environ["HORNSCULPT_TEST_RESULT"])
integration = source / "integrations" / "freecad"
sys.path[:0] = [str(integration / "workbench"), str(integration)]
for module_name in ("bridge", "hornsculpt_commands"):
    sys.modules.pop(module_name, None)
import hornsculpt_commands as commands
import bridge

server = QtNetwork.QTcpServer()
assert server.listen(QtNetwork.QHostAddress.LocalHost, 0)
port = server.serverPort()
token = str(uuid.uuid4())
discovery_directory = Path(QtCore.QStandardPaths.writableLocation(
    QtCore.QStandardPaths.GenericDataLocation)) / "HornSculpt" / "freecad-link"
discovery_directory.mkdir(parents=True, exist_ok=True)
descriptor_path = discovery_directory / ("integration-test-" + token + ".json")
descriptor_path.write_text(json.dumps(dict(protocol=2, host="127.0.0.1", port=port,
                                         token=token, pid=os.getpid(), name="Integration test")),
                           encoding="utf-8")
messages = []
clients = []
finished = False
document = App.newDocument("HornSculptWorkbenchTest")
box = document.addObject("Part::Box", "Box")
box.Length, box.Width, box.Height = 10, 12, 14
document.recompute()


def finish(error=None):
    global finished
    if finished:
        return
    finished = True
    bridge.stop()
    descriptor_path.unlink(missing_ok=True)
    server.close()
    result_path.write_text(json.dumps(dict(passed=error is None, error=error,
                                          freecad_version=App.Version(),
                                          messages=[m["type"] for m in messages],
                                          snapshot_revisions=[m["revision"] for m in messages
                                                              if m["type"] == "snapshot"])),
                           encoding="utf-8")
    App.closeDocument(document.Name)
    QtCore.QTimer.singleShot(0, QtCore.QCoreApplication.quit)


def guarded(action):
    def run():
        try:
            action()
        except Exception:
            finish(traceback.format_exc())
    return run


def receive(socket, buffer):
    buffer.extend(bytes(socket.readAll()))
    while len(buffer) >= 4:
        length = struct.unpack(">I", buffer[:4])[0]
        assert 0 < length <= 64 * 1024 * 1024
        if len(buffer) < length + 4:
            return
        message = json.loads(bytes(buffer[4:length + 4]))
        del buffer[:length + 4]
        messages.append(message)
        if message["type"] == "hello":
            assert message["token"] == token
            payload = json.dumps(dict(type="ready", protocol=2)).encode()
            socket.write(struct.pack(">I", len(payload)) + payload)


def accept():
    socket = server.nextPendingConnection()
    clients.append(socket)
    buffer = bytearray()
    socket.readyRead.connect(guarded(lambda: receive(socket, buffer)))


server.newConnection.connect(guarded(accept))
# Keep actual discovery and command handling, choosing our isolated listener if
# the user has other running HornSculpt instances.
commands.choose_instance = lambda instances, parent: next(
    instance for instance in instances if instance["port"] == port and instance["token"] == token)
commands._report_error = lambda message: finish(message)
if "HornSculptWorkbench" in Gui.listWorkbenches():
    Gui.activateWorkbench("PartWorkbench")
    Gui.removeWorkbench("HornSculptWorkbench")
spec = importlib.util.spec_from_file_location("hornsculpt_test_init_gui",
                                               integration / "workbench" / "InitGui.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
Gui.activateWorkbench("HornSculptWorkbench")
assert "HornSculpt_Connect" in Gui.listCommands()
assert "HornSculpt_Disconnect" in Gui.listCommands()
Gui.runCommand("HornSculpt_Connect")


def no_op_recompute():
    assert [m["revision"] for m in messages if m["type"] == "snapshot"] == ["1"]
    for unused in range(10):
        document.recompute()


def edit():
    assert [m["type"] for m in messages] == ["hello", "geometryPending", "snapshot"]
    box.Length = 20
    document.recompute()


def check_and_disconnect():
    assert [m["revision"] for m in messages if m["type"] == "snapshot"] == ["1", "2"]
    assert [m["type"] for m in messages].count("geometryPending") == 2
    snapshot = messages[-1]
    assert snapshot["protocol"] == 2 and "vertices" not in snapshot and "triangles" not in snapshot
    part = snapshot["brepParts"][0]
    data = base64.b64decode(part["data"], validate=True)
    assert hashlib.sha256(data).hexdigest() == part["sha256"]
    shape = Part.Shape()
    shape.importBrepFromString(data.decode("utf-8"), False)
    assert shape.isValid() and shape.isClosed() and len(shape.Faces) == 6
    assert abs(shape.BoundBox.XMax * snapshot["brepMetresPerUnit"] - .02) < 1e-12
    Gui.runCommand("HornSculpt_Disconnect")
    bridge.stop()
    assert not hasattr(App, "_hornsculpt_live_bridge")
    finish()


QtCore.QTimer.singleShot(1200, guarded(no_op_recompute))
QtCore.QTimer.singleShot(2400, guarded(edit))
QtCore.QTimer.singleShot(4000, guarded(check_and_disconnect))
QtCore.QTimer.singleShot(12000, lambda: finish("GUI integration timed out"))
