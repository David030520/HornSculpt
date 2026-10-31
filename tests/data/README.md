# CAD-Regressionsdaten

`spline_boundary_cusp.brep` enthält einen gültigen geschlossenen Körper mit sechs
Flächen, Koordinaten in Metern. Er wurde aus der getrimmten Korbfläche 35 der
separaten GRAND-ORGUE-Konturanschlusskopie (Quell-BRep SHA256
`8ee771ae837ec57ab202daa590a8e838f66e421b23c81f20743edfa2f2ad238f`)
durch 1 mm Extrusion in +X und Schnitt mit einer 4-mm-Box um
`(-0.117357879186514, 0.41604817913815501, -0.24851544106888598)` erzeugt.
Die beiden ursprünglichen B-Spline-Randkurven und ihre einseitigen Tangenten
bleiben an diesem Knoten erhalten. Der konvexe Tangentenwinkel beträgt etwa
0,01338865°. Ein Tangentenwinkel allein beweist bei endlicher zulässiger
Konturabweichung keine Unmöglichkeit eines Chordnetzes; dieser gekrümmte Fall
darf daher keinen auf geraden Randrichtungen beruhenden Quickcheck auslösen.
Keine weiteren Modellteile oder externen Quelldateien werden im Test benötigt.

`trimmed_spline_boundary.brep` enthält ein geschlossenes, gültiges Prisma mit
sechs Flächen, alle Koordinaten in Metern. Die gekrümmte Grundfläche stammt aus
einer einzelnen getrimmten B-Spline-Fläche des vom Nutzer bereitgestellten
GRAND-ORGUE-Modells. Sie wurde mit Open CASCADE 7.9.3 um 1 mm in +X extrudiert;
andere Modellteile sind nicht enthalten. Ihre TopExp-Flächennummer im Prisma ist 5.

Die Randprobe `(-0.045261334628069562, 0.33703602764109225,
-0.18250201538169269)` liegt auf einer echten topologischen Randkante. Eine
unabhängige Projektion auf die ungetrimmte B-Spline-Oberfläche findet dort mit
`GeomAPI_ProjectPointOnSurf`/MIN keine Lösung; MINMAX liefert eine 222,229 mm
entfernte Lösung mit falscher Normale. Die Flächen-pcurve des von BRepExtrema
gefundenen Randträgers liefert dagegen die korrekte UV-Lage und Normale. Der
CAD-Test prüft Position und äußere Normale über den öffentlichen Geometrievertrag.
Die Datei ist ausschließlich Testeingabe, keine Produkt- oder Laufzeitabhängigkeit.

`spline_edge_parameter.brep` ist ein geschlossenes 27-Flächen-Prisma in Metern,
gebildet aus einer einzelnen weiteren Modellfläche mit 1 mm Extrusion in +X.
Die ursprüngliche Fläche hat hier die TopExp-Nummer 26. An der geprüften
Randposition ist die direkte SameParameter-Näherung der pcurve um 236,813 nm
ungenau, obwohl Kante und Fläche 100 nm Modellierungstoleranz tragen. Eine
lokale Inversion derselben pcurve reduziert den Rest auf 34,804 nm. Der Test
prüft den öffentlichen Projektionsvertrag: tatsächlicher Flächenpunkt und
zugehörige Normale, während die gemeinsame 3D-Randkurve separat erhalten bleibt.
Eine weitere Probe hat 411,084 nm Kurven-/Flächenabstand bei 401,946 nm
nativer Kantentoleranz und prüft genau diese Trennung ohne Toleranzänderung.
Der Test benötigt das Gesamtmodell nicht.

`spline_topological_junction.brep` ist ein gültiger geschlossener Körper mit neun
Flächen in Millimetern. Er wurde mit OCCT 7.9.3 durch Schnitt der getrennten
GRAND-ORGUE-Konturanschlusskopie mit einer achsenparallelen 4-mm-Box um den
betroffenen Knoten erzeugt. Andere Modellteile sind nicht enthalten; das
Gesamtmodell wird im Test nicht benötigt. Die Flächen 1, 6 und 7 teilen den
Knoten `(-68.188134021083008, 389.05774564574302, -290.75482811787301)` mm.
Seine gespeicherte Modellierungstoleranz beträgt 652,619 nm. Die eigene
Vertex-UV-Repräsentation der Fläche 6 rekonstruiert ihn mit 652,613 nm Abstand;
eine zusätzliche allgemeine Flächenprojektion liefert dagegen 691,051 nm.
Der Test verlangt die ursprüngliche gemeinsame Topologieposition, prüft ihre
lokale Abfrage und lehnt entfernte oder nicht gemeinsam inzidente Flächen ab.
Die Toleranzen und Parameter der Quelldaten werden nicht verändert.


`spline_boundary_projection.brep` ist ein gültiger geschlossener Körper mit acht
Flächen in Metern. Er stammt aus einer einzelnen Korbwand der separaten
GRAND-ORGUE-Konturanschlusskopie und wurde mit OCCT 7.9.3 um 1 mm in +X
extrudiert. Die Flächen 3 und 7 teilen eine B-Spline-Randkurve vom Grad 8 mit
93 Polen. Bei der dokumentierten Probe liefert die allgemeine Shape-Distanzsuche
nur einen 1,277110 mm entfernten Randvertex. Die direkte Suche auf genau derselben
getrimmten Kurve findet dagegen bei Parameter 0,63991455978852374 den nächsten
Punkt mit 6,529595 µm Abstand. Der Test prüft diese Lösung, beide Endpunkte und
die Unabhängigkeit von der Reihenfolge der beteiligten Flächen. Das Gesamtmodell
und externe Programme werden für den Regressionstest nicht benötigt.
