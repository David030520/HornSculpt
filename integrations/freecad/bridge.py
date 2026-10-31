"""FreeCAD-only API bridge for HornSculpt's versioned local geometry stream."""

import base64
import hashlib
import json
import math
import struct
import uuid

import FreeCAD as App
from PySide import QtCore, QtNetwork


_MAX_FRAME_BYTES = 64 * 1024 * 1024


def _digest(value):
    payload = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _components(vector, scale=1.0):
    return [float(vector.x) * scale, float(vector.y) * scale, float(vector.z) * scale]


def _object_identity(obj):
    # FreeCAD persists ID through save/reopen; labels and FaceN ordinals are not identities.
    # Reading the built-in identity avoids touching and invalidating the CAD model.
    return str(obj.Name) + ":" + str(obj.ID)


def _face_description(face, object_box, mapped_name, object_id, label):
    box = face.BoundBox
    center = face.CenterOfMass
    bounds = [object_box.XLength, object_box.YLength, object_box.ZLength]
    object_min = [object_box.XMin, object_box.YMin, object_box.ZMin]
    center_values = _components(center)
    normalized_center = [
        round((center_values[i] - object_min[i]) / bounds[i], 7) if bounds[i] > 1e-10 else 0.5
        for i in range(3)
    ]
    u0, u1, v0, v1 = face.ParameterRange
    normal = face.normalAt((u0 + u1) / 2.0, (v0 + v1) / 2.0)
    normal_values = [round(value, 8) for value in _components(normal)]
    surface_type = type(face.Surface).__name__
    structure = [surface_type, len(face.Edges), len(face.Wires)]
    topology_key = _digest([structure, normalized_center, normal_values])
    # Local SI signatures make rigid object placement independent of face matching.
    geometry = [
        structure,
        round(float(face.Area) * 1e-6, 13),
        [round(value, 12) for value in _components(center, 0.001)],
        normal_values,
        [round(float(value) * 0.001, 12) for value in
         (box.XMin, box.YMin, box.ZMin, box.XMax, box.YMax, box.ZMax)],
        sorted(round(float(edge.Length) * 0.001, 12) for edge in face.Edges),
    ]
    return {
        "objectId": object_id,
        "mappedName": mapped_name,
        "label": label,
        "surfaceType": surface_type,
        "geometryKey": _digest(geometry),
        "topologyKey": topology_key,
        "area": float(face.Area) * 1e-6,
        "center": _components(center, 0.001),
        "normal": _components(normal),
        "edgeCount": len(face.Edges),
        "wireCount": len(face.Wires),
    }, structure


def _candidate_objects(document, object_names):
    if object_names:
        objects = [document.getObject(name) for name in object_names]
        if any(obj is None for obj in objects):
            raise ValueError("Ein ausgewaehltes FreeCAD-Objekt wurde geloescht.")
        return objects
    candidates = []
    for obj in document.Objects:
        if not hasattr(obj, "Shape") or obj.Shape.isNull():
            continue
        if App.GuiUp and not obj.ViewObject.Visibility:
            continue
        candidates.append(obj)
    # A visible body and its visible tip describe the same solid. Export the container once.
    included = set(candidates)
    return [obj for obj in candidates if not any(parent in included for parent in obj.InList)]


def _snapshot(document, object_names, origin, axis, session_id, revision):
    faces = []
    objects = []
    brep_parts = []
    for obj in _candidate_objects(document, object_names):
        if not hasattr(obj, "Shape") or obj.Shape.isNull():
            raise ValueError("Ausgewaehltes Objekt besitzt keine Randgeometrie.")
        if "Invalid" in obj.State or "Touched" in obj.State:
            raise ValueError("FreeCAD-Modell muss erfolgreich neu berechnet sein.")
        shape = obj.Shape.copy()
        shape.Placement = App.Placement()
        if not shape.isValid() or not shape.isClosed() or not shape.Solids:
            raise ValueError("FreeCAD-Link benoetigt geschlossene gueltige Volumenkoerper.")
        placement = obj.getGlobalPlacement()
        matrix = placement.toMatrix()
        pose = [float(value) for value in matrix.A]
        for index in (3, 7, 11):
            pose[index] *= 0.001
        object_id = _object_identity(obj)
        object_box = shape.BoundBox
        first_face = len(faces)
        topology = []
        for index, face in enumerate(shape.Faces):
            ordinal = "Face" + str(index + 1)
            mapped_name = shape.getElementMappedName(ordinal)
            description, structure = _face_description(
                face, object_box, str(mapped_name or ""), object_id,
                str(obj.Label) + " / " + ordinal,
            )
            faces.append(description)
            topology.append(structure)
        # Export only the exact shape, without FreeCAD's display triangulations.
        global_shape = shape.cleaned()
        global_shape.Placement = placement
        brep = global_shape.exportBrepToString().encode("utf-8")
        brep_parts.append({
            "objectId": object_id,
            "data": base64.b64encode(brep).decode("ascii"),
            "sha256": hashlib.sha256(brep).hexdigest(),
            "faceIndices": list(range(first_face, len(faces))),
        })
        objects.append({
            "id": object_id,
            "name": str(obj.Name),
            "label": str(obj.Label),
            "placement": pose,
            "topologyFingerprint": _digest(sorted(topology)),
        })
    if not objects:
        raise ValueError("Kein sichtbarer oder ausgewaehlter geschlossener Koerper vorhanden.")
    return {
        "type": "snapshot", "protocol": 2, "units": "m",
        "sessionId": session_id, "documentId": str(document.Uid),
        "documentName": str(document.Label), "revision": str(revision),
        "acousticOrigin": list(origin), "acousticAxis": list(axis),
        "objects": objects, "faces": faces,
        "brepMetresPerUnit": 0.001, "brepParts": brep_parts,
    }


class _Bridge(QtCore.QObject):
    def __init__(self, port, token, object_names, origin, axis):
        super().__init__()
        self.port = int(port)
        self.token = str(token)
        self.object_names = list(object_names)
        self.origin = tuple(float(value) for value in origin)
        norm = math.sqrt(sum(float(value) ** 2 for value in axis))
        if not math.isfinite(norm) or norm <= 0 or len(axis) != 3 or len(self.origin) != 3:
            raise ValueError("Akustische Achse oder Ursprung ungueltig.")
        self.axis = tuple(float(value) / norm for value in axis)
        if not all(math.isfinite(value) for value in self.origin):
            raise ValueError("Akustischer Ursprung muss endlich sein.")
        if not 0 < self.port <= 65535 or not self.token:
            raise ValueError("Verbindung ungueltig.")
        self.session_id = str(uuid.uuid4())
        self.revision = 0
        self.document_id = None
        self.dirty = True
        self.ready = False
        self.stopped = False
        self.received = bytearray()
        self.last_error = None
        self.last_snapshot = None
        self.pending_document = None
        self.socket = QtNetwork.QTcpSocket(self)
        self.socket.connected.connect(self._connected)
        self.socket.disconnected.connect(self._disconnected)
        self.socket.readyRead.connect(self._read)
        self.socket.errorOccurred.connect(self._connection_error)
        self.debounce = QtCore.QTimer(self)
        self.debounce.setSingleShot(True)
        self.debounce.timeout.connect(self._flush)
        self.reconnect = QtCore.QTimer(self)
        self.reconnect.setSingleShot(True)
        self.reconnect.timeout.connect(self._connect)
        self.document_poll = QtCore.QTimer(self)
        self.document_poll.setInterval(700)
        self.document_poll.timeout.connect(self._check_document)
        self.document_poll.start()
        App.addDocumentObserver(self)
        self._connect()

    def _connect(self):
        if not self.stopped and self.socket.state() == QtNetwork.QAbstractSocket.UnconnectedState:
            self.socket.connectToHost("127.0.0.1", self.port)

    def _send(self, message):
        payload = json.dumps(message, separators=(",", ":"), allow_nan=False).encode("utf-8")
        if len(payload) > _MAX_FRAME_BYTES:
            raise ValueError("FreeCAD-Snapshot ueberschreitet das Protokolllimit von 64 MiB.")
        if self.socket.bytesToWrite() + len(payload) > _MAX_FRAME_BYTES:
            raise ValueError("FreeCAD-Link wartet auf Uebertragung des vorherigen Snapshots.")
        self.socket.write(struct.pack(">I", len(payload)) + payload)

    def _connected(self):
        self.received.clear()
        self._send({"type": "hello", "protocol": 2, "token": self.token, "sessionId": self.session_id})

    def _read(self):
        self.received.extend(bytes(self.socket.readAll()))
        if len(self.received) > 4096:
            self.socket.abort()
            return
        while len(self.received) >= 4:
            length = struct.unpack(">I", self.received[:4])[0]
            if length > 4092 or length == 0:
                self.socket.abort()
                return
            if len(self.received) < length + 4:
                return
            payload = bytes(self.received[4:length + 4])
            del self.received[:length + 4]
            message = json.loads(payload)
            if message.get("type") != "ready" or message.get("protocol") != 2:
                self.socket.abort()
                return
            self.ready = True
            # A newly authenticated receiver needs a complete snapshot even if CAD is unchanged.
            self.last_snapshot = None
            self.pending_document = None
            self._mark_dirty()

    def _disconnected(self):
        self.ready = False
        self.dirty = True
        if not self.stopped:
            self.reconnect.start(1000)

    def _connection_error(self, _error):
        if not self.stopped:
            self.reconnect.start(1000)

    def _check_document(self):
        document = App.ActiveDocument
        identity = str(document.Uid) if document else "none"
        if identity != self.document_id:
            self.document_id = identity
            self._mark_dirty()

    def _mark_dirty(self):
        if self.stopped:
            return
        self.dirty = True
        if self.ready:
            self.debounce.start(150)

    def _publish_pending(self, document):
        identity = str(document.Uid) if document else "none"
        if self.pending_document != identity:
            self._send({
                "type": "geometryPending", "protocol": 2, "sessionId": self.session_id,
                "documentId": identity, "revision": str(self.revision + 1),
            })
            self.pending_document = identity

    def _flush(self):
        if not self.ready or not self.dirty:
            return
        document = App.ActiveDocument
        if document is None:
            self._publish_pending(None)
            return
        try:
            snapshot = _snapshot(document, self.object_names,
                                 self.origin, self.axis, self.session_id, self.revision + 1)
            # Recompute notifications also occur for unchanged geometry. Compare the
            # complete payload, excluding only its transport revision; face identity
            # and exact BRep data remain unchanged.
            content = dict(snapshot)
            del content["revision"]
            if content == self.last_snapshot and self.pending_document is None:
                self.dirty = False
                self.last_error = None
                return
            self._publish_pending(document)
            self._send(snapshot)
        except (ValueError, RuntimeError, App.Base.FreeCADError) as error:
            self._publish_pending(document)
            message = str(error)
            if message != self.last_error:
                App.Console.PrintError("HornSculpt-Link: " + message + "\n")
                self.last_error = message
            return
        self.revision += 1
        self.last_snapshot = content
        self.pending_document = None
        self.document_id = str(document.Uid)
        self.dirty = False
        self.last_error = None

    def slotRecomputedDocument(self, document):
        if document == App.ActiveDocument:
            self._mark_dirty()

    def slotChangedObject(self, obj, property_name):
        if obj.Document == App.ActiveDocument and property_name in ("Shape", "Placement", "Group", "Tip"):
            self._mark_dirty()

    def slotDeletedObject(self, obj):
        if obj.Document == App.ActiveDocument:
            self._mark_dirty()

    def slotDeletedDocument(self, _document):
        self._mark_dirty()

    def stop(self):
        if self.stopped:
            return
        self.stopped = True
        App.removeDocumentObserver(self)
        self.debounce.stop()
        self.reconnect.stop()
        self.document_poll.stop()
        self.socket.abort()
        if getattr(App, "_hornsculpt_live_bridge", None) is self:
            del App._hornsculpt_live_bridge


def stop():
    """Disconnect the active link; repeated calls are harmless."""
    connection = getattr(App, "_hornsculpt_live_bridge", None)
    if connection is not None:
        connection.stop()
    if hasattr(App, "_hornsculpt_live_bridge"):
        del App._hornsculpt_live_bridge


def start(port, token, object_names=(),
          acoustic_origin_m=(0.0, 0.0, 0.0), acoustic_axis=(0.0, 0.0, 1.0)):
    """Connect the active FreeCAD document; coordinates use its global frame."""
    stop()
    if not object_names and App.GuiUp:
        import FreeCADGui as Gui
        object_names = [obj.Name for obj in Gui.Selection.getSelection() if hasattr(obj, "Shape")]
    bridge = _Bridge(port, token, object_names, acoustic_origin_m, acoustic_axis)
    App._hornsculpt_live_bridge = bridge
    return bridge
