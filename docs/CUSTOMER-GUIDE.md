# HornSculpt · Kurzanleitung

Das Windows-Paket vollständig entpacken und **Start-HornSculpt.cmd** starten.
Für Draft wird kein FreeCAD benötigt. **Datei → Öffnen** lädt vorhandene
`.hsproj`-Projekte. **Speichern unter …** legt einen eigenen Projektstand an.
**Hilfe → Kurzanleitung** und **Hilfe → Wissensbasis** öffnen die integrierte
Offline-Hilfe. Quellen im Internet werden erst beim Anklicken geöffnet.

Für einen ersten Überblick **Datei → Öffnen** wählen und im entpackten Paket
`examples/HF-Draft.hsproj` öffnen. Das Beispiel enthält eine fertige Berechnung.
Mit **Speichern unter …** einen eigenen Stand anlegen, bevor Änderungen
gespeichert werden.

Das Speichern läuft im Hintergrund. Änderungen währenddessen bleiben als
ungespeichert markiert; bei Neu, Öffnen und Schließen wird gegebenenfalls
nochmals gefragt. Projektdateien sind auf 256 MiB begrenzt; größere Projekte
werden mit einer Fehlermeldung abgewiesen.

Unter **Studie → Ergebnisstände und Bericht …** lässt sich das aktuelle
Ergebnis benennen. Mehrere markierte Stände werden auf gemeinsamen Achsen
verglichen; der erste markierte Stand in Listenreihenfolge ist die Referenz.
Für den PDF-Bericht wird der gewünschte Stand ausdrücklich ausgewählt.
Benannte Ergebnisstände werden mit dem Projekt gespeichert. Nur ein aktuelles,
vollständiges Ergebnis kann als neuer Stand übernommen werden. Jeder Stand
belegt zusätzlichen Platz in der Projektdatei; die Grenze von 256 MiB gilt für
das gesamte Projekt. Projekte mit benannten Ständen verwenden Format 15 und
lassen sich in HornSculpt 0.1.52 und älteren Versionen nicht öffnen. Ohne
benannte Stände bleibt das bisherige Format erhalten.

Unter **Studie → Parameter-Sweep …** lassen sich die Studienachsen zunächst
in Ruhe einstellen. Erst **Start** fordert die GPU für die Studie an; eine
laufende Desktop-Berechnung wird dafür beendet. Die Übergabe kann mit
**Abbrechen** abgebrochen werden.

Die angebotenen Profilparameter richten sich nach dem Draft: Im symmetrischen
Modus gelten H und V jeweils für beide Seiten, im asymmetrischen Modus lassen
sich rechts, links, oben und unten getrennt wählen. Beim Öffnen einer
gespeicherten Studie gilt deren eigene Geometriebasis.

Markierte berechnete Varianten können in der Studie als **Polars H/V**
verglichen werden. Für **Polarmaps H/V** eine Variante ausdrücklich auswählen;
beide Ansichten bieten 0°- oder Maximumbezug. Fehlende Beobachterdaten bleiben
Lücken. **Im Hauptfenster analysieren** öffnet eine ausgewählte berechnete
Variante als unbenannte Kopie mit ihren gespeicherten Ergebnissen. Vorher wird
bei ungespeicherten Änderungen gefragt; die Studiendatei bleibt erhalten.

**Studie → Netzkonvergenz …** fragt zunächst die Prüffrequenzen ab. Danach
grobe und feinere Zielkante sowie mindestens zwei Stufen einstellen. Alle
Stufen verwenden dieselben Prüffrequenzen; der Vergleich belegt Konvergenz nur
für die jeweils berechneten Fälle.

## Systemvoraussetzungen

Windows 10 22H2 oder Windows 11 22H2-SV2, 23H2, 24H2 oder 25H2,
jeweils 64 Bit. Die 3D-Ansicht benötigt OpenGL 3.3 Core. Für die akustische
GPU-Berechnung ist eine NVIDIA-Grafikkarte mit CUDA Compute Capability 7.5
oder höher und ein passender NVIDIA-Treiber erforderlich; empfohlen wird
Treiberzweig R610 oder neuer. Die benötigten Qt- und CUDA-Laufzeiten sind im
Paket enthalten, ein CUDA-SDK ist nicht erforderlich. Der benötigte
GPU-Speicher hängt vom Modell und der Berechnung ab.

Geprüft wurde Windows 11 Pro 64 Bit, Build 26200, mit einer NVIDIA GeForce
RTX 5070 Laptop GPU und Treiber 610.88. Weitere Karten und Windows-Versionen
sind damit nicht einzeln geprüft. Herstellerquellen stehen im
Lieferstand.

## Erstes Draft-Horn

**Horn-Assistent …** in der Projekt-Werkzeugleiste oder im Menü **Geometrie**
öffnet den geführten Entwurf. Der Assistent fragt den freien runden
Halsdurchmesser, axiale Länge, radiale Wandstärke, geometrische H/V-Vollwinkel
und die Querschnittsgrundform ab. Die Konturvorschau und Prüfseite zeigen die
gewählten Werte; auf der vierten Seite folgen Netzvorgabe und Frequenzbereich.

Die Winkel beschreiben das OS-SE-Grundprofil. Akustische Abstrahlwinkel
ergeben sich erst aus einer Berechnung. Optional startet der Assistent nach
der Übernahme Netzaufbau und Berechnung; eine automatische Klangoptimierung
findet nicht statt.
Die Kreis-/Ellipsenform verwendet Exponent 2; Superellipsen verwenden
Exponent 4 oder 6. Der runde Hals geht glatt zwischen 0 % und 65 % der
Hornlänge in die Zielform über.

**Draft übernehmen** ersetzt die Grundgeometrie des aktiven Drafts als einen
Bearbeitungsschritt. Leitkurve, Mund-Morphing, Asymmetrie und Diffraktionsspalt
sind beim neuen Grundhorn ausgeschaltet. Baugruppenplatzierungen und
CAD-Anteile bleiben erhalten. **Abbrechen** übernimmt keine Assistentenwerte.
**Strg+Z** nimmt die Übernahme zurück, **Strg+Y** wiederholt sie. Die Datei
wird erst mit **Speichern** geschrieben. Im FreeCAD-Link zuerst auf **Draft**
wechseln, um den Assistenten zu verwenden.

## Vom Entwurf zum Ergebnis

Die Leiste über dem Modell zeigt den aktuellen Arbeitsschritt.

1. **Entwurf** öffnet die Parameter beziehungsweise die CAD-Flächenzuordnung.
   Nach einer Änderung prüft HornSculpt die Geometrie im Hintergrund.
2. **Netz erzeugen** öffnet die Netzeinstellungen und startet den gemeinsamen
   Qualitätsnetz-Auftrag. Während des Netzaufbaus ist **Netz abbrechen** in
   der Werkzeugleiste erreichbar. Fehlerdetails stehen im Status und an den
   betroffenen Einstellungen; ein altes Netz gilt nicht als aktuelles Netz.
3. **Berechnung …** zeigt die Berechnungseinstellungen und Ergebnisse unten.
   Frequenzen, Anregung und Rechenmodus prüfen, dann dort **Berechnen** starten.
   Der Rotationsmodus verwendet seine eigene geprüfte Meridianquelle.

Die Leiste unterscheidet laufende Aufträge, ungültige Geometrie,
Netz-/Berechnungsfehler und aktuelle Ergebnisse. Ergebnisse behalten ihre
Revision und werden nach fachlichen Änderungen als veraltet behandelt.
Nach einem unerwarteten Programmende bietet HornSculpt eine vorhandene
automatische Sicherung zur Wiederherstellung an. Änderungen werden dafür
ungefähr alle 60 Sekunden im Hintergrund gesichert. Der wiederhergestellte
Stand hat zunächst keinen Dateinamen und wird über **Speichern unter …**
gesichert; die ursprüngliche Projektdatei bleibt separat erhalten.

## Ansicht und Platz

Die Vorschau beginnt ohne Dreieckskanten. **Kanten** schaltet das Overlay
ein; die Wahl bleibt über Neustarts erhalten. **Einpassen** oder **F** passt
das ganze Modell ein. Das erste erzeugte Modell und ein neu übernommenes
Grundhorn werden einmal automatisch eingepasst; spätere Parameteränderungen
erhalten die Kamera.

Das Standardlayout reserviert ungefähr 35 % der Fensterhöhe für die unteren
Ergebnisdocks. Gespeicherte Benutzerlayouts bleiben erhalten.
**Ansicht → Docklayout zurücksetzen** stellt die Standardanordnung wieder her,
ohne Modellwerte oder Ergebnisse zu ändern. Auf schmalen Fenstern ordnet das
Erstlayout Projekt, Netz und Parameter als seitliche Register an; im CAD-Modus
steht dort das CAD-Dock bereit. Docks lassen sich anschließend frei vergrößern,
verschieben und ein- oder ausblenden. **Ansicht → Ergebnisse vergrößern**
(Strg+Umschalt+E) gibt den Ergebnisdiagrammen mehr Platz; erneutes Wählen
stellt die 3D-Ansicht wieder her.

## Akustische Ergebnisse

Nach einer erfolgreichen Berechnung unten die gewünschte Ansicht wählen:
Frequenzgang, H/V-Polars, Isobaren, Winkelvergleich, Phase oder Schallfeld.
**Bereich** und **Ansicht** ordnen die Ergebnisse; **Qualität** zeigt den
Ergebniszustand und erläutert die verfügbaren Angaben.
Die angezeigte Frequenz und Normierung gehören zur Ergebnisdarstellung.
Absolutpegel hängen von der tatsächlichen Anregung ab. Eine Normierung auf
0° beschreibt die Richtwirkung relativ zur Achse und ist kein Absolutpegel.

Für **CTA-2034** ein Fernfeldergebnis berechnen; die Standardrechnung legt
das benötigte H/V-Winkelraster an. Im Reiter **CTA-2034** die Anregung wählen
und bei Bedarf **CTA-CSV …** für die ausgewerteten Kurven verwenden. Die
Kurven benötigen jeweils unterschiedliche Richtungen: Für Schallleistung und
Raumprognose sind vollständige 360°-H/V-Daten erforderlich. Fehlende Werte
erscheinen als Lücken in den Diagrammen und im CTA-CSV. Der allgemeine
CSV-Export enthält die komplexen Rohdaten.

Im Ergebnisreiter **Schallfeld** H- oder V-Schnitt, Frequenz, Anregung und
**Pegel** oder **Phase** auswählen. Über **Optionen** Bereich und Raster
prüfen und anschließend **Feld berechnen** drücken. Die berechnete Ebene erscheint oben im 3D-Modell;
**3D einpassen** richtet die gemeinsame Ansicht von Modell und Feld aus. Ist
die Ebene ausgeblendet, **Ansicht → Schallfeldebene** einschalten. Die
Feldebene muss an ihren tatsächlichen Weltpunkten separat berechnet werden;
beim Öffnen eines Projekts startet diese Berechnung nicht automatisch. Eine
H/V-Fernfeldkurve ist keine räumliche Nahfeldkarte. Fehlende oder ungültige
Werte dürfen keine scheinbar vollständige Karte bilden. Die Feldfarben sind
abgestuft: Pegel in 3-dB-Schritten, Phase in 36°-Schritten und der Realteil in
zehn Stufen über den angezeigten Bereich. Die Pegelpalette entspricht der
Darstellung der Isobaren. Bei geringer Fensterhöhe öffnet **Werkzeuge öffnen**
die Arbeitsablauf- und 3D-Ansichtssteuerung in einem eigenen Fenster.
**Messpunkt** ermöglicht die Abfrage eines berechneten Feldpunkts in der
3D-Ansicht. Die Animation zeigt Phasen desselben berechneten Feldes. Der
Zwischenspeicher hält höchstens acht berechnete Ebenen; für eine nicht mehr
gespeicherte Frequenz ist **Feld berechnen** erneut erforderlich.

## CAD verwenden

Für eine STEP-Datei **Datei → STEP importieren …** wählen. Für eine laufende
FreeCAD-Verbindung oben **FreeCAD-Link** öffnen und im CAD-Dock
**Workbench installieren …** verwenden. Nach dem FreeCAD-Neustart dort die
Workbench **HornSculpt** wählen, die gewünschten geschlossenen Körper
auswählen und **Verbinden** drücken.

Flächen im Modell oder der Flächentabelle auswählen; Strg erweitert die
Auswahl. **Wall** bedeutet schallhart. **Driven** erhält eine explizite
komplexe Normalgeschwindigkeit in m/s RMS und muss mit **Zuweisen** bestätigt
werden. Positive Geschwindigkeit zeigt aus dem Festkörper in die Luft.
Nach Topologieänderungen als **prüfen** markierte oder entfallene Flächen
gezielt kontrollieren. Anschließend denselben Netz- und Berechnungsablauf
verwenden. Details stehen in [FreeCAD-LINK.md](FreeCAD-LINK.md).

## Wenn ein Schritt nicht weitergeht

| Befund | Nächster Schritt |
| --- | --- |
| Modell zu klein oder nicht sichtbar | **F** drücken; bei verstellten Docks **Ansicht → Docklayout zurücksetzen**. |
| Netzaktion gesperrt | Abschluss der Vorschau abwarten; Geometriediagnose rechts beziehungsweise CAD-Zuordnungen prüfen. |
| Letzte gültige Vorschau sichtbar | Die geänderten Parameter sind ungültig. Diagnose korrigieren oder mit **Strg+Z** zurückgehen. |
| Netz fehlgeschlagen | Fehlermeldung und Netzeinstellungen prüfen. Zielkante, Mindestkante und Detailabweichung müssen zur kleinsten tatsächlichen Geometrie passen. |
| Berechnung gesperrt | Aktuelles geprüftes Netz beziehungsweise gültige Rotationsquelle sowie mindestens eine unterstützte Anregung bereitstellen. |
| Ergebnis veraltet | Nach Geometrie-, Anregungs- oder Solveränderung den passenden Netz-/Berechnungsschritt erneut starten. |
| CAD-Verbindung unterbrochen | FreeCAD-Dokument neu berechnen und Verbindung in der Workbench prüfen; gespeicherte Daten bleiben offline lesbar. |

Die Vorschau ist kein Qualitätsnetz und ein erfolgreiches Netz kein
Konvergenznachweis. Lineares BEM weist weder Hochpegeleignung noch thermische
Belastbarkeit eines Treibers nach. Bekannte konkrete Einschränkungen
und Prüfnachweise stehen im kurzen Lieferstand und
ausführlich in STATUS.md.

## Teilgebiete und `.msh`

**Datei → Gmsh-Netz importieren …** liest ASCII MSH 2.2/4.1 mit linearen
Dreiecken und ebenen konvexen Vierecken. Dateieinheit ausdrücklich wählen;
Entity-IDs und Physical Names bleiben im Projekt erhalten. In der
Flächentabelle die tatsächlichen Antriebsflächen auf **Driven** setzen.
**Importnetz prüfen** prüft bei dieser Quelle das Originalnetz; es wird nicht
unbemerkt neu vernetzt. Eine aktive FreeCAD-Verbindung vor dem Import trennen.

Unter **Akustik → Berechnung → Teilgebiete …** benannte Innenräume anlegen
und physische Flächen sowie vorhandene Interfaceflächen getrennt zuordnen.
Eine Interfacefläche verbindet zwei Gebiete, sie ist keine Wall- oder
Driven-Fläche. Ohne vorgegebene Interfaceflächen lassen sich aus einer bereits
geschlossenen physischen Quelle ebene Öffnungsabschlüsse erzeugen.
**3D** und den **dichten GPU-Solver** verwenden; pro LEM-Port nur physische
Antriebsflächen desselben Gebiets wählen. **Interfaces** zeigt die wirklichen
Schnittflächen im Modell. Details, Normalenregeln und Grenzen stehen unter
[Akustische Teilgebiete](SUBDOMAINS.md).
