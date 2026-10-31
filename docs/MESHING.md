# Oberflächennetz

`hs_meshing` erzeugt ein isotropes Dreiecksnetz durch eine eigene
Frontal-Delaunay-Implementierung. CGAL 6.2.1 liefert exakte geometrische
Prädikate, die beschränkte 2D-Triangulierung und räumliche Suchstrukturen.
Der neutrale Vertrag enthält SI-Koordinaten, Dreiecke, Flächen-IDs,
Randbedingungen und optionale feste Punkte.

## Verfahren

Eingabe ist die bereits geschlossene physische Oberfläche einschließlich
Antriebsfläche. Der Mesher erfindet keine Abschlussflächen.

1. Ein Flächenatlas erzeugt injektive lokale Koordinaten. Native CAD-Parameter
   werden direkt übernommen: UV-stetige Komponenten bleiben zusammen, echte
   periodische Nahtverwendungen bleiben getrennt. Facettennormalen der Vorschau
   erzeugen darin keine zusätzlichen Projektionsschnitte. Orientierung und
   Einfachheit der nativen UV-Ränder werden weiterhin geprüft.
   Ohne native Parameter werden ebene Bereiche mit Bohrungen direkt projiziert,
   entwickelbare Flächen können längentreu
   aufgefaltet werden. Eine nicht direkt projizierbare Ringfläche mit genau
   zwei Randzyklen wird entlang eines kürzesten Quellkantenpfads zwischen
   diesen Rändern geöffnet. So erzeugt die topologische Diskbildung keine
   unnötige Zickzacknaht mit spitzen künstlichen Ecken. Physische Ränder,
   Flächenrollen und Quellgeometrie bleiben erhalten; beide Seiten des Schnitts
   verwenden dieselben räumlichen Netzknoten. Andere Bereiche verwenden geprüfte Projektionen mit
   begrenzter metrischer Verzerrung. Künstliche Schnitte erzeugen keine
   physischen Flächen oder zusätzlichen Krümmungssaaten. Spitze künstliche
   Schnittohren können in die benachbarte Karte wechseln, wenn beide Karten
   zusammenhängend und injektiv bleiben und die betroffenen Eckverletzungen
   strikt abnehmen. Quellfacetten, physische Tags und geschützte Knicke bleiben
   unverändert; der Normalenkegel bleibt auf 60° begrenzt.
2. Eine gemeinsame Kurvenregistrierung tastet Flächenränder nach dem skalaren
   Größenfeld ab. Nachbarn und beide Seiten periodischer Schnitte teilen
   dieselben räumlichen Knoten. Vorschauknoten sind keine obligatorischen
   inneren BEM-Knoten.
3. Frontpunkte entstehen aus der lokalen physischen Oberflächenmetrik:
   größenabhängige Off-Center-Punkte und Umkreismittelpunkte. Die beschränkte
   Delaunay-Triangulierung erhält Ränder. Zu nahe Punkte und Verletzungen
   geschützter Randsegmente werden abgewiesen. Ihre globale Koordinatenmetrik
   ist das nach physischer Facettenfläche gewichtete Gram-Mittel. Jeder neue
   Punkt muss die ausgewählte schlechte Zelle tatsächlich ersetzen.
   Native Karten verwenden dafür exakte D1 am UV-Schwerpunkt jeder
   Referenzfacette einschließlich der Kartenorientierung. Nur PL-Karten
   rekonstruieren Ableitungen aus Vorschaufacetten; dünne native
   Vorschaufacetten dürfen die reguläre CAD-Metrik nicht verfälschen.
   Native Kandidaten verwenden die tatsächlichen räumlichen Dreiecksknoten in
   der lokalen Tangentenebene. Begrenzte gedämpfte Gauss-Newton-Schritte führen
   diese physischen Ziele in den unveränderten nativen Parameterbereich zurück.
   Unterschiedliche UV-Geschwindigkeiten sind keine physischen Kantenlängen.
   Bei adaptiver Geometrietreue liefern unzureichend aufgelöste Sehnen oder
   Zellen zusätzlich native Mittelpunktkandidaten. Der Fehlerindikator misst
   die normale Abweichung; tangentiale Unterschiede zwischen UV- und räumlichem
   Mittelpunkt sind keine Krümmungsfehler.
4. Physische Diagonalwechsel und begrenzte Verschiebungen freier Innenknoten
   dürfen die dimensionslose Winkel-, Form-, Größen- und Geometrieverletzung streng
   verringern. Die Abstiegsrichtung entsteht aus dieser Energie in der
   orthonormalen physischen Tangentenmetrik; eine begrenzte Liniensuche prüft
   jeden Schritt auf der tatsächlichen Oberfläche. Gemeinsame Ränder und fixe
   Punkte bleiben unverändert. UV- und räumliche Orientierung, Randbedingungen
   bleiben verbindlich. Ein Zwischenschritt darf noch die Detailabweichung
   verletzen, wenn die gesamte Verletzungsenergie strikt sinkt; die abschließende
   Prüfung darf dies nicht. Die Optimierung endet bei erfüllten Grenzen oder
   tatsächlichem Stillstand; ein Budget von 100.000 Innenknotenversuchen pro
   Karte begrenzt nicht konvergierende Fälle. Die Liniensuche wählt den besten
   geprüften Schritt. Erst am Stillstand sucht sie feiner bis zur Auflösung
   der numerischen Ableitung; diese Suchgenauigkeit lockert keine Netzgrenze.
   Nach Diagonalwechseln und Verschiebungen können zusätzliche Punkte für
   Größe, Geometrietreue oder Formqualität nötig sein. Diese werden in der
   beschränkten Triangulierung ohne Delaunay-Annahme eingefügt. Die Verbindung
   vom Ausgangsdreieck zum Kandidaten darf keine geschützte Randkante kreuzen.
   Beide Seiten einer betroffenen Innenkante, räumliche Orientierung,
   Punktabstände und Randsegment-Schutz bleiben geprüft. Neue Randpunkte
   entstehen ausschließlich durch koordinierte Teilung des gemeinsamen Atlas.
   Jede Einfügung muss die beanstandete Zelle tatsächlich ersetzen. Eine
   innere Kantenhalbierung ersetzt beide angrenzenden Zellen; ein numerisch
   knapp neben der Kante liegender Mittelpunkt wird mit seiner tatsächlichen
   Lage eingefügt und die alte innere Diagonale anschließend entfernt.
   UV- und räumliche Orientierung aller vier resultierenden Dreiecke werden
   vorab geprüft. Nach langsamer Größenkonvergenz werden solche Teilungen
   gebündelt, bevor erneut freie Knoten verschoben werden.
   Abschließend muss jedes Dreieck die aktivierten Qualitätsgrenzen erfüllen.
   Ist ein exakter Geometrieanbieter vorhanden, werden zusätzlich die räumlichen
   Sehnenmittelpunkte und Dreiecksschwerpunkte gegen seine getrimmte Referenzfläche
   und die aktive Detailabweichung geprüft. Ein nativ ausgewerteter Punkt,
   dessen Lage strikt innerhalb des tatsächlichen Flächentrims klassifiziert ist,
   liefert eine sichere obere Distanzschranke. Liegt diese bereits innerhalb
   der Detailabweichung, ist keine globale Nächstpunktsuche nötig:
   `d_min ≤ |Sehnenpunkt − Flächenpunkt| ≤ Detailabweichung`. Bei unzureichender
   Schranke, Randpunkten oder fehlender Trim-Evidenz bleibt die genaue Projektion
   auf die begrenzte Fläche erforderlich. Eine Randsehne kann
   außerhalb des natürlichen Trägerbereichs liegen; ihre Distanz zum physischen
   Flächenrand bleibt eine gültige Messung. Die davon getrennte Innenabbildung
   verlangt weiterhin einen eindeutigen stationären Trägerpunkt. Dieselbe
   physische Distanz löst bereits in der Frontverfeinerung und Qualitätsreparatur
   lokale Verfeinerung aus; sie ist keine ausschließlich abschließende Schranke.
   Diese Stichproben
   sind unabhängig von der nativen UV-Geschwindigkeit, aber kein globaler
   Hausdorff-Nachweis.
   Reine Abstandsprüfungen verwenden dafür einen eigenen skalaren Vertrag:
   mehrere gleich nahe Punkte am Trim haben denselben Abstand und sind kein
   Fehler. Eine Positions-/Normalenabbildung darf daraus weiterhin keinen
   willkürlichen Träger auswählen; ihre Mehrdeutigkeitsprüfung bleibt erhalten.
   Auswertung von Position/D1 und aufwendige Trim-Klassifikation sind getrennte
   Abfragen; reine Metrik- und Randbindungsauswertungen klassifizieren nicht.
   OCCT verwendet konservative Begrenzungsboxen zur Beschleunigung seiner
   unveränderten Trimprüfung. Ein auf 65.536 exakte Parameterabfragen begrenzter,
   threadsicherer Cache gehört ausschließlich zum unveränderlichen CAD-Snapshot;
   weder Parameter-Rundung noch Toleranzaufweitung finden statt.
   Qualitätswerte gehören zur aktuellen Zelle und
   den Positionsrevisionen ihrer drei Knoten. Einfügen, Diagonalwechsel und
   jede versuchte oder zurückgenommene Knotenverschiebung machen betroffene
   Einträge ungültig. Bereits bewiesene identische Stichproben werden in der
   Abschlussprüfung wiederverwendet; sämtliche Qualitäts- und Topologiegrenzen
   bleiben erhalten. Lange Prüfschleifen beachten den Auftragsabbruch.
5. Erforderliche Randteilungen werden zwischen Frontdurchläufen
   zentral auf der ursprünglichen Kurve ausgeführt, für alle Nachbarn zugleich.
   Fixierte Dreiecksränder werden niemals einseitig geteilt. Überschrittene
   Arbeitsbudgets liefern einen Fehler und kein Teilnetz.
   Das Budget zählt höchstens 100.000 zusätzlich eingefügte gemeinsame
   Randknoten. Jede Runde muss neue Knoten und betroffene Karten liefern.
   Eine feste Zahl globaler Runden wäre kein belastbares Arbeitsmaß: sehr
   dünne reale CAD-Streifen können viele Teilungsstufen benötigen, obwohl
   die übrigen Karten bereits abgeschlossen sind. Solche Merkmale werden
   weder entfernt noch mit vergrößerten Toleranzen übergangen.
   Zusätzlich wird je Karte die dimensionslose Summe der verbleibenden
   Qualitäts-, Größen- und Geometrieverletzungen verfolgt. Sechs erneute
   Randverfeinerungen ohne mindestens ein Prozent Verbesserung beenden den
   Auftrag mit der gemessenen Restverletzung. Das ist eine Begrenzung
   nicht konvergierender Arbeit und kein Beweis geometrischer Unmöglichkeit.
   Auch ein außerhalb der Karte liegender Qualitätskandidat kann eine
   Randteilung erfordern: Angefragt wird nur ein vom Zellinneren aus zuerst
   gekreuztes und vom Kandidaten verletztes Randsegment. Der Kandidat selbst
   wird außerhalb der Karte niemals eingefügt.
   Nach einer Teilung werden nur betroffene Karten einschließlich aller
   Nachbarseiten neu vernetzt. Erfolgreiche unveränderte Karten bleiben
   innerhalb desselben Auftrags erhalten. Alte gemeinsame Knoten-IDs und
   Positionen, Referenzgeometrie und Größenfeld bleiben dabei unverändert;
   die abschließende Prüfung gilt weiterhin für das gesamte Netz.

Methodische Grundlage ist [Rebays Frontal-Delaunay-Verfahren](https://doi.org/10.1006/jcph.1993.1097),
das auch [Gmsh als Algorithmus 6](https://gmsh.info/doc/texinfo/) verwendet.
Gmsh 4.15.2 bleibt externe Vergleichsreferenz; automatischer Wechsel auf
MeshAdapt ist dabei deaktiviert (`Mesh.AlgorithmSwitchOnFailure=0`).
Kein Gmsh-Quellcode und keine Gmsh-Bibliothek werden übernommen. Allgemeine
Qualitäts- oder Geschwindigkeitsparität ist damit nicht nachgewiesen.

## Größenfeld und Geometrietreue

Globale Zielkante, Mindestkante, Detailabweichung, Wachstum und lokale
Zielkanten je physischer Fläche gelten für Draft und CAD gleichermaßen.
Das Größenfeld ist skalar. Unbekannte Flächen und lokale Größen außerhalb
der globalen Grenzen werden abgewiesen.

Die Mindestkante begrenzt die lokale **Zielgröße des Größenfelds**, nicht jede
tatsächliche Dreieckskante. Krümmungs-, Naht- und Randabstands-Saaten werden auf
das Intervall zwischen Mindestkante und globaler Zielkante begrenzt; eine
kleinere Detailabweichung senkt dieses Feld nicht unter die Mindestkante.
Unabhängig davon können die geometrische Konturabweichung, erhaltene feste
Randpunkte, Referenzlinsen und konforme Qualitätsverfeinerungen kürzere
Randsegmente verlangen. Insbesondere künstliche Atlasränder können dadurch
zusätzliche Knoten behalten. Die Mindestzielkante hat aktuell keinen Vorrang
vor diesen Geometriebedingungen. Eine harte Untergrenze jeder tatsächlichen
Kante wäre bei unveränderlichen kleinen Geometriedetails nicht allgemein
erfüllbar.

Schmale physische Bereiche werden anhand von Fläche und Umfang begrenzt.
Krümmungen und gekrümmte Nahtkonturen erzeugen lokale Größen. Kurze
Vorschaukanten allein erzwingen keine feine Vernetzung. Unvermeidbare physische
Randsegmente zwischen festen Punkten geben dagegen eine lokale Zielgröße vor;
das vorhandene Größenwachstum führt von dort ins Flächeninnere. Das sind
Größenvorgaben, keine zusätzlichen Winkel-/Formgrenzen. Das Feld wächst
entlang der Oberfläche; nahe gegenüberliegende Wände übernehmen nicht allein
durch räumliche Nähe die Feinheit der Gegenseite. Eine konforme
Feldunterteilung erlaubt breite Flächen mit gröberem Inneren, auch bei einer
CAD-Quelle aus langen Fächerdreiecken.

Gröber abgetastete gekrümmte Ränder können kleine Linsen gegenüber der
Quellpolygonlinie erzeugen. Erweiterungen der Referenzparametrisierung sind
nur bei explizit geprüfter Geometrie zulässig: einfache konvexe Linse,
keine Überschneidung oder verschluckte Bohrung, begrenzte Konturabweichung
und reguläre lokale Metrik. Mehrdeutige Abschnitte behalten Quellknoten.
Es gibt keine unbegrenzte Extrapolation.

Detailabweichung steuert Krümmungsgrößen und diese lokalen Prüfungen.
Sie ist keine nachgewiesene globale Hausdorff-Fehlergrenze des Netzes.
Bei reinen Polygonquellen ist die übergebene Oberfläche die Referenz.
Für akustische Genauigkeit zusätzlich den Netzkonvergenz-Auftrag verwenden.

## Exakte CAD-Quelle

FreeCAD-Link und STEP-Import liefern exakte BRep-Volumenkörper. Open CASCADE
7.9.3 prüft Einheiten, geschlossene Topologie und die explizite Flächenbindung.
Eine native Tessellation stellt die Topologie und Referenzdaten bereit; sie ist
kein fertiges BEM-Netz.

`CadGeometry::meshingSource` erzeugt Netzquelle und unveränderliche
`SurfaceParameterization` gemeinsam. Die Flächenkoordinaten werden pro
Dreiecksecke **vor** dem topologischen Zusammenführen erfasst. Quellrandkanten
tragen eine räumliche Kurven-ID, eine eigene Randverwendungs-ID und ein
gerichtetes Intervall des kanonischen Kurvenparameters. Der Anbieter besitzt
seinen vollständigen Quellsnapshot und bleibt nach Freigabe des ursprünglichen
CAD-Objekts gültig. Ein früherer Vorschauaufruf ist kein impliziter Cache.

Die zulässige Parameterdomäne gehört zum unveränderlichen Quellsnapshot.
Native pcurves und Tessellations-UVs können bei gültigen CAD-Modellen geringfügig
über die nominellen Basisgrenzen reichen. Solche Quellbereiche werden aus den
tatsächlichen Daten begrenzt übernommen; Abfragen erweitern sie niemals.
Jede dort verwendete Auswertung liefert unverändert `S(rawUV)` und dessen echtes
`D1`. Ihr räumlicher Abstand zum entsprechenden Basisrand muss sowohl unter der
zugehörigen nativen Modellierungstoleranz als auch unter der angeforderten
Detailgenauigkeit liegen. Jede pcurve-Hülle trägt ihre unveränderte Kanten-/
Flächentoleranz; echte Endpunkte tragen zusätzlich ihre eigene Vertextoleranz.
Nur Evidenz, deren UV-Hülle die Abfrage enthält, darf diese Grenze bestimmen. Es gibt weder ein
globales Toleranzmaximum noch eine Addition zum Nutzerbudget.
Positionen und Ableitungen werden nicht geklemmt. Weiter entfernte,
unbelegte oder ungenaue Auswertungen bleiben Fehler.
Der Atlas prüft die Übereinstimmung von Positionen, Dreiecksindizes,
Flächen-IDs und festen Punkten mit diesem Snapshot. Der Aufrufer darf Wall und
Driven neu zuweisen. Im Desktop und bei Studien entsteht die native Netzquelle
mit der aktuellen Detailabweichung im Hintergrundauftrag. Die archivierte
Anzeigevorschau wird dadurch nicht ersetzt; deren Revision und Normalenmaske
bleiben getrennte Zustände.

Die vorhandene Front wertet native Karten direkt als `S(u,v)` mit den echten
Ableitungen `dS/du` und `dS/dv` aus. Parameterkoordinaten haben die Einheiten des
Anbieters, Positionen und Ableitungen sind Meter beziehungsweise Meter je
Parameter. Bei umgekehrter Flächenorientierung wird die V-Kartenkoordinate
reflektiert; die nativen Ableitungen des Anbieters bleiben unverändert. Die
physische Metrik und das skalare Größenfeld bestimmen weiterhin die
Frontpunktplatzierung und alle Dreiecksgrenzen.

Physische gemeinsame Randpunkte werden einmal als `C(t)` erzeugt. Jede
Nachbarseite verwendet ihre eigene pcurve `gamma_use(t)` desselben Parameters.
Periodische Nähte besitzen zwei Randverwendungen und verschiedene UVs derselben
räumlichen Kurve. Nur echte Perioden rechtfertigen einen konstanten UV-Lift.
Künstliche Kartenschnitte werden auf derselben Fläche ausgewertet und erzeugen
keine physische Flächengruppe. Zentrale Randteilungen bestimmen den Mittelpunkt
nach räumlicher Bogenlänge und aktualisieren alle Seiten atomar. Eine separate
Nächstpunktprojektion je Nachbarseite findet im parametrischen Netzpfad nicht
statt.

Der neutrale Anbieter liefert mit `curvePartition` den vollständigen, streng
geordneten nativen Trim einschließlich aller Polynom-/Kontinuitätsspannen.
Beim CAD-Anbieter werden dafür OCCT-CN-Intervalle verwendet, also auch innere
Spline-Knoten mit höherer Stetigkeit. Jede Teilbogenlänge wird vor der adaptiven
GL5-Quadratur an diesen Spannen partitioniert, auch bei umgekehrter Laufrichtung.
Das gemeinsame numerische Fehlerbudget und die Bisektionsgrenzen bleiben
unverändert. Ein Vergleich eines ganzen Quadraturpanels mit seinen Hälften
allein könnte eine schmale Spanne in beiden Schätzungen übersehen und einen
falschen räumlichen Mittelpunkt liefern. Die Regression prüft eine gerade
Kante mit lokalisierter kubischer Parametergeschwindigkeit in beiden Richtungen.

`SameParameter` ist eine native Approximation. Der Mesher misst den Rest
`|S(gamma_use(t)) - C(t)|` an den verwendeten Parametern gegen die gewünschte
Geometriegenauigkeit. Die gespeicherte Modellierungstoleranz wird getrennt
übertragen und nicht zum Nutzerbudget addiert. Ein unverändert erhaltener
Topologievertex kann einen eigenen nativen Modellrest besitzen. Nicht endliche
Daten, inkonsistente Quellbindungen und nicht unterstützte Parameterachsen
liefern einen Fehler; Toleranzen werden nicht erhöht und es gibt keinen stillen
Wechsel auf einen anderen CAD-Mesher. Diese lokalen Prüfungen sind keine globale
Hausdorff-Garantie und beweisen keine Qualitätsfeasibility beliebiger CAD-Ecken.

Eine nachgewiesene glatte Rotationskappe erhält reguläre kartesische
Radialkoordinaten statt einer am Pol kollabierenden Winkelachse. Der aktuelle
Nachweis gilt für nicht rationale B-Spline-Generatrizen mit streng monotonem,
planarem radialem Kontrollpolygon, einem bis auf radial skalierte Rundung
kollabierten Endpunkt und verschwindender axialer/radialer Polsteigung. Damit
ist der Radius auf allen Polynomspannen eindeutig invertierbar. Position und
Ableitungen, Quell-Dreieckskoordinaten sowie pcurve-Werte und ihre Ableitungen
verwenden denselben Koordinatenwechsel; kanonische räumliche Kurven und
Flächenidentitäten ändern sich nicht. Die frühere Winkelnaht ist dadurch keine
physische Netzrandbedingung. Andere Parameterfamilien werden durch diesen
Nachweis nicht als reguläre Kappen ausgegeben.

`isInsideTrim` ist eine konservative Anbieterevidenz und standardmäßig
falsch. Der CAD-Anbieter setzt sie ausschließlich bei strikt `IN`
klassifizierten nativen Parametern, nicht bei `ON` oder nur natürlichem
UV-Domänenbesitz. Für solche Punkte ist die volle euklidische Distanz zwischen
`S(uv)` und Sehnenpunkt ein oberer Wert der Distanz zur getrimmten Fläche.
Liegt bereits dieser obere Wert innerhalb der Detailabweichung, darf die
Front-/Reparaturprüfung die globale Nächstpunktsuche auslassen. Unbekannte
Trim-Lage und größere obere Werte verwenden weiterhin die physische Suche.
Die abschließende Prüfung darf das bereits bewiesene Distanzkriterium nur für
unveränderte Knotenidentitäten und Positionsrevisionen übernehmen.

Der allgemeine `SurfaceGeometry`-Projektionsvertrag bleibt für explizite
geometrische Abfragen und Anbieter ohne native Parameterdaten verfügbar.
`projectSurface` prüft den tatsächlich gefundenen getrimmten Support einschließlich
seiner pcurve- oder Vertex-Parameter. `projectSupportingSurface` verlangt ein
eindeutiges reguläres stationäres Minimum im natürlichen Trägerbereich;
endliche B-Spline-Bereiche werden nicht erweitert. `projectBoundary` prüft die
getrimmte räumliche Kurve, ihre Extrema und Endpunkte. `projectJunction` verlangt
einen gemeinsamen nativen Vertex innerhalb seiner gespeicherten Toleranz.
Mehrdeutige oder entfernte Zweige bleiben Fehler. Diese Abfragen rekonstruieren
keine fehlenden nativen Snapshot-UVs.

Die native Trimmklassifikation verwendet pro unveränderlicher Fläche einen
`BRepTopAdaptor_FClass2d`. Ein unabhängiger konkaver Polygonprisma zeigt, dass
der strahlbasierte OCCT-Klassifizierer nach SI-Skalierung echte Innenpunkte
mit mehreren Millimetern Randabstand als außen einstufen kann. Die neue
Strategie wird auch an Außenpunkten, Löchern und periodischen Trägern geprüft;
ein `OUT` wird nicht pauschal übergangen. Bei positivem Abstand konkurriert
ein stationärer Trägerpunkt mit der gebundenen Flächen-/Randsuche. Ein ferneres
stationäres Minimum darf einen näheren Trimmrand nicht verdrängen. Der skalare
Abstand berücksichtigt beide zulässigen Zeugen; die Geometrietoleranz bleibt
unverändert. Die Regressionen erzeugen ihre Prüfkörper selbst ohne externe
Horn-Dateien.

Für einen nachgewiesenen vollen Umlauf einer planaren, einseitig zur Achse
liegenden B-Spline-Generatrix wird die Nächstpunktsuche exakt auf die
Radius-/Axialebene reduziert. Positive Gewichte und das Kontrollpolygon
zertifizieren diese Geometrieeigenschaften unabhängig von der Vorschau.
Die Kurvensuche berücksichtigt jeden nativen CN-Abschnitt, dessen Endpunkte
und stationäre Minima; die Auswahl erfolgt anhand der tatsächlichen
dreidimensionalen Distanz. Anschließend gelten dieselbe Trimklassifikation,
Orientierung und Stationaritätsprüfung wie für andere Flächen. Ein gesuchter
Punkt auf der Rotationsachse besitzt bei positivem nächstem Radius einen
ganzen Ring gleich naher Punkte und wird ausdrücklich als mehrdeutig abgewiesen.
Ein eindeutiger kollabierter Pol wird davon unterschieden. Andere Flächen
behalten ihre passende vorhandene Projektionsstrategie.

Die zugehörige 17-Pol-Regression enthält einen kurzen Sickenabschnitt, bei dem
die zweidimensionale OCCT-Extremensuche sowohl mit Grad- als auch Tree-Strategie
ein Minimum von 223,223 µm meldete, obwohl ein strikt innerhalb der Fläche
liegender Punkt bereits 92,667 µm entfernt war. Der tatsächliche Abstand von
62,712 µm wird durch die Meridianprüfung einschließlich starrer Platzierung
nachgewiesen. Die Korrektur betrifft die Projektion selbst; Genauigkeitsgrenzen
und die unabhängige finale Netzkontrolle bleiben unverändert.

Wiederverwendbare OCCT-Projektionszustände werden während jedes Aufrufs exklusiv
geliehen und auch bei Fehlern zurückgegeben; es gibt keinen globalen Cache.
Reine Polygonquellen verwenden weiterhin ihre übergebene Oberfläche als
Referenz. Vorschau und Vernetzung bleiben verschiedene Datenstände. BRep-
Artefakte im Projekt erlauben Offline-Remeshing ohne laufendes FreeCAD.

## Feste Anschlüsse und Prüfung

`SurfaceMesh.fixedVertices` erhält Quellkoordinaten exakt. Vollständig fixierte
physische Teilflächen bleiben unverändert; andere Flächen können nach innen
gröber werden. Draft verwendet dieselbe räumliche Detailabweichung über die
gesamte komponierte Innen- und Außenfläche. Ein adaptives Referenzgitter prüft
Kurvensehnen sowie Dreieckinnenpunkte, einschließlich Off-Axis-Positionen,
Leitkurven, LUT und Mund-Morphing. Jede Kontur wird anschließend unter demselben
Fehlerbudget lokal vereinfacht; unterschiedliche Knotenzahlen werden konform
verbunden und erneut gegen die komponierte Oberfläche geprüft. Eine große oder
detailreiche Kontur erzwingt so keine unnötigen Randpunkte auf kleinen Konturen.
Es enthält nur die
fünf physischen Flächengruppen; allein die Driven-Randkontur ist fixiert.
Es gibt weder einen künstlichen Wandstreifen noch ein tangentenwinkelbasiertes
Halsbudget oder einen zweiten fixierten Ring. Driven-Flächeninhalt und Umfang
werden separat gegen die Vorschau geprüft. Der geometrische Fehler ist durch
Abtastungen geprüft, keine analytische Hausdorff-Garantie. Ein beliebiger
Tangentenwinkelfehler wird damit ebenfalls nicht garantiert.

Vor und nach der Vernetzung werden Koordinaten, Indizes, Duplikate,
Degeneration, Kanteninzidenz, Orientierung, Komponenten, Vertex-Links und
Selbstüberschneidungen geprüft. Komponentenzahl, Eulercharakteristik,
physische Flächen-IDs, Randbedingungen und feste Koordinaten bleiben erhalten.
Ein zugesetzter Durchgang darf nicht allein wegen Formschluss akzeptiert werden.

Die optionalen Standardgrenzen verlangen mindestens 20° Innenwinkel und höchstens 4 für
längste Kante / kleinste Höhe. Ein gleichseitiges Dreieck hat das Verhältnis
`2 / sqrt(3)`. Die Zielkante ist ein Diskretisierungsmaß, keine harte Obergrenze
jeder Kante. Winkel- und Formgrenzen sind im gemeinsamen Netz-Panel einstellbar
und abschaltbar. Ausgeschaltet treiben diese Werte weder Frontverfeinerung noch
Endabnahme; Größenfeld und nicht entartete Elemente bleiben erforderlich.
Kein Retry ändert die gewählten Einstellungen. Projektformat 7 speichert sie;
Version 5 erhält ihre Grenzen, ältere Projekte verwenden die Standardgrenzen.
Beim Lesen gespeicherter Netze
werden gültige geschlossene Randflächen geprüft, ohne pauschale 20°/4-Sperre.

Nach dem Atlasaufbau kann ein belegter ebener Linienwinkel einen Auftrag vor
der Frontverfeinerung ablehnen. Der neutrale Geometrieanbieter weist dafür
eine exakte Ebene und zwei analytische, koplanare Linienintervalle nach.
Zusätzlich müssen Eckpunkt, beide Randabschnitte und ihre aktuellen
Registry-/UV-Verwendungen im tatsächlich benutzten Atlas erhalten sein.
Ein bloßer CAD-Tangentenwinkel oder eine grobe Vorschausehne genügt nicht.
Bei einem solchen unveränderlichen ebenen Fächer hat mindestens ein Dreieck
einen Winkel höchstens so groß wie der konvexe Eckwinkel. Mit
`Seitenverhältnis >= 1 / sin(Winkel)` kann auch die aktivierte Formgrenze
einen Widerspruch belegen. Reflexe, gekrümmte und nicht belegbare Fälle werden
dadurch nicht gesperrt; diese notwendige Bedingung garantiert kein Netz.
Bei ausgeschalteten Qualitätsgrenzen entfällt die Prüfung vollständig.
CAD-Modelltoleranzen und Konturabweichung vergrößern den Beleg nicht.

Insbesondere beweist ein kleiner gekrümmter GRAND-ORGUE-Tangentenwinkel bei
erlaubter endlicher Konturabweichung keine allgemeine Unmöglichkeit eines
endlichen Chordnetzes. Eine native Nachprüfung des lokalen CAD-Arbeitskandidaten
hat außerdem einen unabhängigen Quellfehler belegt: Auf Fläche 35 stimmen
kanonische räumliche Kurve und pcurve
trotz gesetztem `SameParameter` nicht innerhalb der angeforderten 0,25 mm
überein. Der Mesher meldet Fläche, Kurve, Randverwendung, Parameter und Rest;
er erhöht weder das Budget noch korrigiert er die CAD-Quelle still.
Vorschautopologie und `SameParameter`-Flag allein sind keine Genauigkeitsnachweise.

Flächenjobs laufen auf bis zu vier CPU-Threads. Gemeinsame Ränder ändern sich
nur zwischen vollständigen Durchläufen. Stabile IDs bestimmen Einfüge- und
Ausgabeordnung unabhängig von Speicheradressen und Threadreihenfolge.
Abbruch liefert kein Teilnetz.

## GPU und Nachweise

Atlas-Projektionen werden vor dem Zusammenführen und der Konditionierung
künstlicher Ecken auf globale Eindeutigkeit geprüft und bei Bedarf zerlegt.
Eine spätere Zerlegung würde neue spitze Rechenecken erzeugen, die sich durch
Randunterteilung nicht beseitigen lassen. Normalenoktanten werden exakt über
die Komponentenvorzeichen gewählt; eine feste Score-Toleranz darf nahezu ebene
Bereiche nicht in künstliche Inseln zerlegen. Ein spitzer künstlicher Sektor
wird als vollständiger Quellfacettenfächer zwischen benachbarten Rechencharts
derselben physischen Fläche übertragen. Das ist nur über glatte Nachbarschaft
zulässig, wenn beide verbleibenden Charts zusammenhängend und projektiv gültig
bleiben und ihre künstlichen Eckwinkel sich verbessern. Physische Ränder,
Quellknicke und Geometriepositionen ändern sich dabei nicht.
Die Draft-Regression prüft
100 × 10 mm Entry mit Exponenten 4, 8 und 12 sowie 40 × 10 mm mit Exponent 12,
jeweils Querschnittsexponent 8 und Zielkante 16 mm. Alle vier zuvor scheiternden
Fälle erfüllen die unveränderten Grenzen 20°/4 und erhalten den geschlossenen,
schnittfreien Rand einschließlich Fläche, Umfang und festen Driven-Knoten.
Ein runder Entry mit Querschnittsexponent 12 prüft zusätzlich einen solchen
Mehrfacettenfächer am Mund und erfüllt dieselben Bedingungen bei 16 mm Zielkante.
Zusätzlich werden geometrisch unveränderte Randabfragen an dünnen
Referenzdreiecken innerhalb ihrer Koordinaten-Rundungsgrenze stabil ausgewertet;
echte Abfragen außerhalb der Referenzfläche bleiben abgewiesen.

Der native CAD-Import prüft die Quelltessellation auch dann auf vollständige
gegenorientierte Randkanten, wenn OCCT den Tessellierungslauf als beendet
meldet. Dünne gekrümmte Ringe können bei grober Sehnenapproximation fehlende
Dreiecke erzeugen. Nur solche unvollständigen Tessellierungen werden automatisch
mit halbiertem Sehnenziel wiederholt, bis die Prüfung besteht oder Modellpräzision,
Geometriebudget beziehungsweise Abbruch den Lauf begrenzen. Die angeforderte
Genauigkeit wird dabei strenger eingehalten; CAD-Geometrie und finale
isotrope Qualitätsvorgaben ändern sich nicht. Die Quellpunktzahl ist weiterhin
keine Vorgabe für die gesamte Dichte des endgültigen Netzes.

Periodische native Flächen dürfen ihre äquivalenten Parameter außerhalb des
Grundintervalls auswerten. Diese Fortsetzung darf nicht als unzulässige
Spline-Extrapolation geklemmt werden. Die Regression prüft einen Zylinder mit
nach `[2π,4π]` verschobenen Nahtparametern einschließlich echter Positionen
und Ableitungen; nichtperiodische Modellabweichungen bleiben geprüft.

Eine echte CUDA-Mikromessung auf der RTX 5070 Laptop GPU verglich neunmal
dieselbe Dreiecksqualitätsauswertung einschließlich Host-/Device-Transfers.
Median: 10.000 Dreiecke CPU 0,287 ms / GPU 0,290 ms;
100.000 Dreiecke 2,860 / 1,139 ms; eine Million 29,995 / 11,586 ms.
Das misst keinen vollständigen GPU-Mesher. Bei kleinen BEM-Netzen rechtfertigt
diese Stufe derzeit keinen produktiven GPU-Pfad. Der BEM-Solver verwendet CUDA.

Tests prüfen Bohrungen, periodische Nähte, Flächenmetrik, fixe Anschlüsse,
unterschiedliche Quelltessellationsdichten, lokale Größen, Kugeln,
mehrere Draft-Formen, physische Anschlüsse, Wiederholbarkeit und Abbruch.
Native Regressionen prüfen zusätzlich räumlich identische Ebenen mit stark
veränderlicher Parametergeschwindigkeit, begrenzte optionale Kandidatenabfragen
und eine reguläre Fläche mit kollabierender Kantenprojektion in der lokalen
Tangentenebene. Ein exakter CAD-Quader von 30 × 20 × 0,0115 mm erzeugt bei
8 mm Zielkante 41.036 Dreiecke unter unveränderten 20°/4-Grenzen; Volumen,
Formschluss und Schnittfreiheit bleiben geprüft. Reale dünne Merkmale dürfen
somit zusätzliche Randteilungsstufen benötigen, ohne pauschal als unlösbar zu
gelten.

Die aktuelle Tritonia-Datei mit 49 Flächen und SHA-256
`126b20324bf7c91d24122a07bd1455818f2e85c885aa9cef1fd2c67d2390fae5`
wurde bei 8 mm Zielkante und 0,25 mm Detailabweichung vollständig geprüft.
Mit ausdrücklich ausgeschalteten Winkel-/Formgrenzen entstehen 9.542 Dreiecke
in 28,59 s; Import und Netzquelle benötigen zusätzlich 0,185 s und 0,068 s.
Alle 49 Flächen, eine geschlossene orientierte Komponente, Eulercharakteristik
−6 und Schnittfreiheit bleiben erhalten. Mindestwinkel 0,113932° und
Seitenverhältnis 557,414 kennzeichnen die schmalen Elemente dieses Diagnosefalls.
Er ist kein Nachweis der Standardqualität: Bei aktivierten 20°/4-Grenzen
beendet die Fortschrittsprüfung den Lauf nach 101,162 s mit den verbleibenden
Formverletzungen auf Fläche 12; sie akzeptiert kein Netz. Dort treffen eine
Gerade und eine native Ellipse tangential aufeinander. Zusätzlich enthält
das Modell einen etwa 11,48 µm breiten realen Streifen. Die Geometrie wird
weder entfernt noch durch vergrößerte Toleranzen ersetzt. Ein früherer
58-Flächen-Dateistand hat einen anderen Inhalt; seine Laufzeiten sind nicht
mit diesem Snapshot vergleichbar.
Aktuelle CAD- und Gesamtergebnisse stehen in STATUS.md.

Abgeschaltete Winkel-/Formgrenzen schalten die numerischen BEM-Toleranzen nicht
ab. Der in 0.1.9 beobachtete Nahquadratur-Abbruch beim kleinen Draft
(4 × 4 mm Entry, 20 mm Länge, 1 mm Wand, 8 mm Zielkante, Winkel-/Formgrenzen aus)
ist mit dem aktuellen Geometrie-/Meshingpfad nicht mehr reproduzierbar.
Bei H20°/V15° bestehen 852 Dreiecke die CPU-Nahplanung bei 100 Hz und 20 kHz;
bei H40°/V25° bestehen 1198 Dreiecke die Planung bei 100 Hz. Toleranz `1e-4`
und Unterteilungstiefe 6 bleiben unverändert. Ein echter 100-Hz-CUDA-Solve
des ersten Falls besteht mit Residuum `4,52e-6` bei Vorgabe `5e-5`.
Die Regression `bem_cuda_tests --coarse-draft` erzeugt und berechnet diesen Fall
ohne eine feste Dreieckszahl zu verlangen. Die Prüfung belegt keinen
Netzkonvergenz- oder Hochfrequenz-Genauigkeitsnachweis; eine erfolgreiche
20-kHz-Nahplanung ist keine akustische Freigabe bis 20 kHz.
PH170 besteht ohne Winkel-/Formgrenzen mit 14.332 Dreiecken, mindestens 3,263°,
Aspect 17,648 und unveränderter geschlossener Topologie (46 Flächen, Euler −12).
Dieser Netzaufbau dauerte 1,42 s; eine CAD-CUDA-Gesamtabnahme folgt daraus nicht.
Aktuelle Projekte enthalten keine bereichsspezifische Eintritts-Netzeinstellung.
V4/V5-Werte werden beim Import verworfen. Gespeicherte aktuelle Draft-Quellen
werden mit der aktuellen deterministischen Tessellation verglichen und bei
abweichender Darstellung in der Anwendung erneuert;
zugehörige alte Netze und akustische Ergebnisse bleiben als veraltet sichtbar.
