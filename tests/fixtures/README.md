# CAD-Testdaten

`freecad-box-20x30x40mm.step` wurde mit der FreeCAD-1.1.4-API aus
`Part.makeBox(20, 30, 40)` exportiert. Die Längeneinheit der STEP-Datei ist mm;
die Maße in SI sind 0,02 × 0,03 × 0,04 m. Sie enthält einen geschlossenen
Körper mit sechs ebenen Flächen. Sie prüft
den echten asynchronen STEP-Import ohne eine laufende FreeCAD-Sitzung.

Die beiden Tetraeder-Snapshots wurden ebenfalls über die FreeCAD-API erzeugt
und enthalten jeweils einen exakten BRep-Körper in mm.

`freecad-spline-extrusion.json` enthält einen über die FreeCAD-1.1.4-Brücke
exportierten geschlossenen Körper: eine planare quadratische BSpline durch
(0,0), (10,4), (20,0) mm, drei gerade Rückkanten bis y = −10 mm und eine
Extrusion um 15 mm in Z. FreeCAD bezeichnet die gekrümmte Mantelfläche als
`SurfaceOfExtrusion`, der native Kernel als `SurfaceOfLinearExtrusion`.
Der Protokolltest prüft diese Gleichbedeutung und lehnt einen tatsächlich
abweichenden Flächentyp weiterhin ab. FreeCAD ist keine Laufzeitabhängigkeit
dieses gespeicherten Testfalls.

