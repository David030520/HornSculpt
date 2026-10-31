"""Native FreeCAD commands for HornSculpt's authenticated local link."""

import json
from pathlib import Path

import FreeCAD as App
import FreeCADGui as Gui
from PySide import QtCore, QtWidgets

import bridge


_DESCRIPTOR_KEYS = {"protocol", "host", "port", "token", "pid", "name"}
_MAX_DESCRIPTOR_BYTES = 16 * 1024


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("Doppelter Discovery-Schlüssel.")
        result[key] = value
    return result


def _descriptor(value):
    if not isinstance(value, dict) or set(value) != _DESCRIPTOR_KEYS:
        raise ValueError("Unbekanntes Discovery-Format.")
    if type(value["protocol"]) is not int or value["protocol"] != 2:
        raise ValueError("Unbekanntes Link-Protokoll.")
    if value["host"] != "127.0.0.1":
        raise ValueError("Der FreeCAD-Link erlaubt nur lokale Verbindungen.")
    if type(value["port"]) is not int or not 0 < value["port"] < 65536:
        raise ValueError("Ungültiger Link-Port.")
    if type(value["pid"]) is not int or not 0 < value["pid"] < 2**63:
        raise ValueError("Ungültige Prozess-ID.")
    for key, maximum in (("token", 256), ("name", 128)):
        text = value[key]
        if not isinstance(text, str) or not 0 < len(text) <= maximum:
            raise ValueError("Ungültiger Discovery-Text.")
        if any(ord(character) < 32 for character in text):
            raise ValueError("Ungültiger Discovery-Text.")
    return value


def find_instances():
    """Read descriptors only; probing sockets would interrupt an existing CAD client."""
    location = QtCore.QStandardPaths.writableLocation(
        QtCore.QStandardPaths.GenericDataLocation
    )
    if not location:
        raise RuntimeError("Das lokale Anwendungsverzeichnis ist nicht verfügbar.")
    directory = Path(location) / "HornSculpt" / "freecad-link"
    found = {}
    for path in directory.glob("*.json"):
        try:
            with path.open("rb") as stream:
                data = stream.read(_MAX_DESCRIPTOR_BYTES + 1)
            if len(data) > _MAX_DESCRIPTOR_BYTES:
                continue
            descriptor = _descriptor(json.loads(data, object_pairs_hook=_unique_object))
            key = descriptor["pid"], descriptor["port"]
            timestamp = path.stat().st_mtime_ns
            if key not in found or timestamp > found[key][0]:
                found[key] = timestamp, descriptor
        except (OSError, ValueError, UnicodeError, RecursionError):
            continue
    return [found[key][1] for key in sorted(found)]


def choose_instance(instances, parent):
    if len(instances) == 1:
        return instances[0]
    choices = [
        f'{item["name"]} · PID {item["pid"]} · Port {item["port"]}'
        for item in instances
    ]
    choice, accepted = QtWidgets.QInputDialog.getItem(
        parent, "HornSculpt", "Instanz", choices, 0, False
    )
    return instances[choices.index(choice)] if accepted else None


def _status(message):
    Gui.getMainWindow().statusBar().showMessage("HornSculpt · " + message, 8000)


def _report_error(message):
    App.Console.PrintError("HornSculpt-Link: " + message + "\n")
    QtWidgets.QMessageBox.warning(Gui.getMainWindow(), "HornSculpt", message)


def _connection_started(connection, instance):
    _status(f'Verbindung zu PID {instance["pid"]} wird aufgebaut')

    def check_ready():
        if getattr(App, "_hornsculpt_live_bridge", None) is not connection or connection.stopped:
            return
        if connection.ready:
            _status(f'Verbunden · PID {instance["pid"]} · Port {instance["port"]}')
        else:
            bridge.stop()
            _report_error("Die HornSculpt-Instanz antwortet nicht. Die Verbindung wurde beendet.")

    QtCore.QTimer.singleShot(3000, check_ready)


class ConnectCommand:
    def GetResources(self):
        return {
            "MenuText": "Verbinden",
            "ToolTip": "Aktives Dokument mit einer laufenden HornSculpt-Instanz verbinden. "
                       "Eine Körperauswahl begrenzt den Export; sonst gelten sichtbare Körper.",
        }

    def IsActive(self):
        return True

    def Activated(self):
        try:
            instances = find_instances()
            if not instances:
                _report_error("Keine HornSculpt-Instanz gefunden. Den FreeCAD-Link in HornSculpt starten.")
                return
            instance = choose_instance(instances, Gui.getMainWindow())
            if instance is not None:
                connection = bridge.start(instance["port"], instance["token"])
                _connection_started(connection, instance)
        except (OSError, ValueError, RuntimeError) as error:
            _report_error(str(error))


class DisconnectCommand:
    def GetResources(self):
        return {"MenuText": "Trennen", "ToolTip": "Lokale HornSculpt-Verbindung beenden."}

    def IsActive(self):
        connection = getattr(App, "_hornsculpt_live_bridge", None)
        return connection is not None and not connection.stopped

    def Activated(self):
        bridge.stop()
        _status("Getrennt")


def register_commands():
    Gui.addCommand("HornSculpt_Connect", ConnectCommand())
    Gui.addCommand("HornSculpt_Disconnect", DisconnectCommand())
