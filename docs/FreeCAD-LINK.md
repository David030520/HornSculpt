# FreeCAD-Link

FreeCAD läuft als eigener Prozess. HornSculpt empfängt die physische Randfläche
und verwendet dafür denselben Mesher und CUDA-Solver wie für Draft. FreeCAD ist
für den Draft-Modus keine Voraussetzung. Die kleine neue Python-Brücke greift
ausschließlich innerhalb der FreeCAD-API auf das CAD-Modell zu.

## Bedienung

1. In FreeCAD das Dokument erfolgreich neu berechnen. Die zu übertragenden
   geschlossenen Körper auswählen. Eine sichtbare PartDesign-Body-Auswahl ist
   einer Auswahl sämtlicher Zwischenstufen vorzuziehen.
2. In HornSculpt oben **FreeCAD-Link** wählen. Das CAD-Dock startet den lokalen
   Empfänger. **Workbench installieren …** installiert die neue HornSculpt-
   Workbench im Benutzerverzeichnis von FreeCAD. FreeCAD danach einmal neu
   starten, die Workbench **HornSculpt** wählen und dort **Verbinden** drücken.
   Bei mehreren laufenden HornSculpt-Instanzen die passende Instanz auswählen.
3. Nach dem Verbindungsaufbau erscheinen Dokument, Revision und Flächenliste.
   Flächen durch einen kurzen Linksklick im Modell oder über die Liste
   auswählen. Strg erweitert oder schaltet die Auswahl um, Umschalt wählt in
   der Liste einen Bereich. Die Auswahl wird im Modell hervorgehoben.
4. **Wall** oder **Driven** sowie bei Driven den Real- und Imaginärteil der
   Normalgeschwindigkeit einstellen und **Zuweisen** betätigen. Die Zuweisung
   gilt für sämtliche ausgewählten Flächen. Die Werte sind
   RMS in m/s; positive Geschwindigkeit zeigt aus dem Festkörper in die Luft.
   Alle neuen Flächen beginnen als Wall. Eine Simulation benötigt mindestens
   eine Driven-Fläche.
5. Die **Zielkante** in mm wählen, **Netz erzeugen** starten und danach im
   Akustikdock **Berechnen** verwenden. HornSculpt erzeugt die Vorschau aus der
   übertragenen exakten CAD-Geometrie. Die Vorschau ist kein freigegebenes
   isotropes Qualitätsnetz.

Die zusätzliche Tabellenspalte **Netz [mm]** legt eine lokale Zielkante pro
Fläche fest. Doppelklick oder F2 öffnet die Eingabe; bei Mehrfachauswahl wird
der Wert gemeinsam übernommen. **0 / global** entfernt die lokale Vorgabe.
Erlaubt sind Werte zwischen Mindestkante und globaler Zielkante im gemeinsamen
Netzpanel. Die Zuweisungen nutzen denselben Größenfeld-Vertrag, werden im
Projekt gespeichert und lassen sich mit Strg+Z/Y zurücknehmen/wiederholen.

Mit **Adaptiv** werden Bohrungen, gekrümmte Konturen und dünne Flächen lokal
feiner vernetzt. **Zielkante** begrenzt die Zielgröße in großen Bereichen;
**Min. Kante** begrenzt die lokale Verfeinerung nach unten. **Detailabweichung**
steuert das Krümmungs-Sehnenziel und **Wachstum** den Übergang zu gröberen
Bereichen. Kleinere Detailabweichungen erzeugen mehr Elemente. Diese Werte
steuern das isotrope Größenfeld; die Prüfung auf Mindestwinkel, Seitenverhältnis
und Wasserdichtigkeit bleibt unverändert. Für neue Live-Snapshots ist die
exakte BRep-Geometrie die Referenz für Flächen und gemeinsame Randkurven.

Änderungen und Neuberechnungen in FreeCAD erzeugen neue Snapshots ohne
Export-/Importschleife. Die beim Verbinden ausgewählte Objektmenge bleibt
fest. Ohne Auswahl verwendet die Brücke sichtbare Körper und vermeidet den
doppelten Export eines sichtbaren Containers samt sichtbarem Unterobjekt.
Ein Wechsel der Objektmenge lässt sich über **Trennen**, die passende
FreeCAD-Auswahl und erneutes **Verbinden** festlegen.

Die Workbench findet lokale Empfänger und ihre Sitzungsschlüssel automatisch
über kurzlebige Dateien im Benutzer-Datenverzeichnis. Nach einem
Anwendungsneustart erneut in der Workbench **Verbinden** drücken.
Bei einer unterbrochenen TCP-Verbindung zum weiterhin laufenden Empfänger
versucht die Brücke im Abstand von einer Sekunde die Wiederverbindung.

## Geometrie und Koordinaten

Übertragen werden gültige geschlossene Volumenkörper als OpenCASCADE-BRep und
ihre Flächenmetadaten. FreeCADs Anzeige-Triangulation wird vor dem Export aus
der Kopie entfernt; Knoten- und Dreieckslisten gehören nicht zum Live-Protokoll.
Offene Schalen und einzelne ungeschlossene Antriebsflächen werden
zurückgewiesen. Die Brücke erfindet keine Abschlussfläche. Das Modell muss den
vollständigen physischen Rand einschließlich Antrieb und Rückseite beschreiben.
Mehrere Körper benötigen getrennte geschlossene, geometrisch zulässige Ränder.

Der BRep-Inhalt verwendet ausdrücklich FreeCAD-mm und enthält bereits die
globale Objektplatzierung. Der native C++-Import skaliert einmal mit 0,001 auf
SI. Metadaten und Objektplatzierungen verwenden bereits Meter; lokale
Flächensignaturen werden zur Importprüfung in den globalen Rahmen transformiert.
Die normale Dreieckswindung der abgeleiteten Vorschau zeigt aus dem
Festkörper in die Luft. Die Brücke liest FreeCADs persistente Dokument- und
Objektidentitäten und verändert dafür keine Modelleigenschaften.

Die Brückenfunktion `start` akzeptiert `acoustic_origin_m` und `acoustic_axis`
im globalen Modellkoordinatensystem. Vorgaben sind `(0,0,0)` und `+Z`.
Der empfangene Achsursprung und die normierte Achse bestimmen die akustischen
Messpunkte; die Kamera verschiebt diese Fachkoordinaten nicht. Im Draft-Modus
liegt der akustische Messursprung am Hornmund.

HornSculpt erzeugt die Anzeige-Vorschau nativ mit 0,25 mm Abweichungsziel.
Sie dient der Anzeige und Flächenauswahl; sie ersetzt die exakte Quelle nicht.
Der native Import prüft geschlossene orientierte Körper und die vollständige
Flächenbindung. Fläche, Schwerpunkt, Oberflächentyp und orientierte Normale
müssen zu den übertragenen Flächensignaturen passen.
Axiale Draft-Diagramme werden für allgemeine CAD-Geometrie geleert und als
nicht verfügbar bezeichnet.

## Flächenzuordnung und Revisionen

`FaceN` dient ausschließlich als Name innerhalb eines einzelnen Quellsnapshots.
Es ist keine dauerhafte Identität. Zuordnungen werden in dieser Reihenfolge
abgeglichen: eindeutiger nativer FreeCAD-Topologiename, eindeutige geometrische
Signatur und bei unverändertem erfassten Topologiefingerabdruck eine eindeutige
strukturelle Signatur im lokalen Objektkoordinatensystem. Dieser Fingerabdruck
erfasst Flächenarten sowie Kanten- und Drahtzahlen; die strukturelle Signatur
ergänzt normierte Flächenlage und Normale. Uneindeutige Kandidaten werden nicht
automatisch zugeordnet. Die Zuordnung behält gegebenenfalls sowohl Wall/Driven
als auch die komplexe Normalgeschwindigkeit.

Nach einem nicht eindeutig auflösbaren Topologiewechsel tragen neue Flächen den
Status **prüfen**. **Zuweisen** oder **Zuordnung bestätigen** bestätigt die
ausgewählten Flächen. Verschwundene frühere Flächen erscheinen in der Anzahl
**entfallen** und benötigen **Entfallene Flächen bestätigen**. Vorher werden
deren Antriebe nicht auf andere Flächen übertragen; Netzfreigabe und Simulation
bleiben gesperrt.

Der Knopf **Entfallene Flächen bestätigen** erscheint nur, wenn frühere
Flächenzuordnungen tatsächlich entfallen sind. Fehlt er, gibt es dort nichts
zu bestätigen. Eine deaktivierte Zuordnungsbestätigung bei bereits bestätigten
Flächen ist ebenfalls ein normaler Zustand.

Solange FreeCAD eine Änderung meldet, aber keinen gültigen vollständigen Stand
geliefert hat, erscheint **Änderung ausstehend**. Das gilt auch bei einem noch
berührten oder ungültigen CAD-Modell. Vorhandene Ergebnisse bleiben sichtbar und
werden als veraltet markiert. Verbindungsabbruch und Moduswechsel entwerten
laufende Netz-/Solveraufträge; verspätete Ergebnisse dürfen keinen neueren
Stand überschreiben. Eine reine Randbedingungsänderung kann ein bereits gültiges
Netz mit aktualisierten Flächentags wiederverwenden.

## Lokales Protokoll

Der C++-Adapter lauscht ausschließlich auf IPv4-Loopback. Die Verbindung beginnt
mit Version 2, einem zufälligen Sitzungsschlüssel und einer Sitzungs-ID.
UTF-8-JSON-Nachrichten tragen einen vorangestellten 32-Bit-Längenwert in
Netzwerkbyteordnung. Eine Nachricht ist auf 64 MiB beschränkt.

Ein Snapshot enthält Dokument-ID, Sitzung, monoton steigende Revision,
SI-Metadateneinheit, Objektplatzierungen, akustische Achse und Flächensignaturen.
Je Objekt kommen ein Base64-BRep, dessen SHA-256 und die explizite Bindung seiner
lokalen Flächendurchlaufreihenfolge an die Metadaten hinzu. Diese Reihenfolge
gilt nur innerhalb desselben Artefakts. Erst der Identitätsabgleich vergibt die
dauerhaften physischen IDs; der Import verwendet niemals `FaceN` als Identität
über Revisionen hinweg. Revisionen werden als Dezimalstrings übertragen.
Unvollständige Nachrichten bleiben unsichtbar. Parser-, Größen-, Transformations-
und native BRep-Prüfungen müssen vor der atomaren Übernahme bestehen. Die
Vorschau entsteht danach im Hintergrund aus dem zugeordneten exakten Stand.
Eine währenddessen bearbeitete Randbedingung löst die erneute Zuordnung gegen
den neuen aktuellen Metadatenstand aus und wird nicht still überschrieben. Veraltete
Revisionen derselben Sitzung werden verworfen. Eine neue authentifizierte Sitzung
darf ihren Revisionszähler neu beginnen; bekannte Dokumentzuordnungen können
innerhalb der laufenden HornSculpt-Sitzung erhalten bleiben.

## Ausgeführte Prüfung und Grenzen

Die C++-Protokolltests prüfen Fragmentierung, atomare Übernahme, Größenlimit,
Authentifizierung, veraltete Revisionen, unvollständige BRep-Flächenbindung, Änderung von
Flächenreihenfolgen, Zuordnungserhalt, Topologiewechsel sowie Trennung und
Wiederverbindung. Diese Prüfungen bestanden im nativen Build am 19.10.2026.
Zusätzlich wurde die echte FreeCAD-1.1.4-API mit Qt-Ereignisschleife und TCP-Brücke
gestartet: Ein geschlossener Körper wurde von 10 auf 20 mm geändert und neu
berechnet; der neue exakte SI-Stand und die erhaltene komplexe Driven-Zuordnung
wurden geprüft. Dazu erhält `freecad_tests` die drei Argumente
FreeCAD-Python, FreeCAD-Bibliotheksordner und vollständiger Pfad zu `bridge.py`.

Zusätzlich bestand ein unsichtbarer zusammenhängender Linktest mit FreeCADCmd
1.1.4 und der echten Bridge: Ein vom öffentlichen ST260-30 abgeleitetes
geschlossenes STEP mit bündigem Idealkolben wurde live übertragen. Eine
10-mm-Platzierungsänderung erzeugte eine neue Revision; nach TCP-Abbruch
und automatischer Wiederverbindung blieb eine komplexe Driven-Zuordnung
erhalten. Das Projekt wurde nach Beenden von FreeCAD offline wieder geöffnet.
Der Test bediente keine Workbench und ist kein akustischer Solvernachweis.

Der separate Workbench-Test bestand ebenfalls mit einer echten FreeCAD-GUI und
ihren registrierten Verbinden-/Trennen-Kommandos in einem isolierten Prozess.
Zehn Neuberechnungen ohne Geometrieänderung erzeugten keine zusätzlichen
Nachrichten; die tatsächliche Änderung lieferte genau einen neuen BRep-Snapshot.
Prüfsumme, nativer BRep-Inhalt und Abmessungen wurden kontrolliert. Die laufende
Benutzersitzung blieb unverändert.

CAD-Dock und gemeinsame Netz-/Solveranbindung sind im nativen Build enthalten.
Der portable PH170-BRep-Test besteht mit eigener exakter Neuvernetzung,
zwei CUDA-Sweeps über je sieben Frequenzen und vollständigem Ergebnis-Save/Open.
Live-Übertragung, Flächenzuweisung und Solverablauf werden durch eigene
Prüfungen abgesichert. Dies ist keine allgemeine Hochfrequenz- oder
Leistungsfreigabe.

`.hsproj` ab Format 7 speichert die exakten BRep-Artefakte samt Einheiten,
Objekt- und Flächenbindung zusätzlich zu CAD-Snapshots, Vorschau, Qualitätsnetz und
akustische Ergebnisse einschließlich des ursprünglichen Rechenauftrags. Die
Datei lässt sich ohne laufendes FreeCAD öffnen. Gespeicherte Ergebnisse werden
ohne erneute Simulation dargestellt; veraltete Ergebnisse bleiben als solche
gekennzeichnet. Für eine feste, nachvollziehbare Exportmenge die gewünschten
Körper beim Verbinden auswählen.

Dateien aus Format 1–6 bleiben lesbar. Enthält ihr CAD-Stand nur ein gespeichertes
Netz, wird dadurch keine exakte CAD-Quelle erfunden. Ein neuer Live-Snapshot
liefert die BRep-Quelle. Gespeicherte Ergebnisse bleiben auch ohne diese Quelle
lesbar. Die Persistenztests bestanden ebenfalls: exakte BRep-Daten und
Zuordnungen bleiben beim Speichern/Öffnen erhalten; die Vorschau lässt sich
ohne FreeCAD aus der gespeicherten Quelle neu erzeugen. Beschädigte Daten,
Prüfsummen, Einheiten und widersprüchliche native Flächensignaturen werden
zurückgewiesen.

Bei nativen gekrümmten Flächenanschlüssen kann eine grobe Vorschaufacettierung
schneiden, obwohl die exakte BRep-Geometrie gültig ist. Der Import beginnt bei
0,25 mm Sehnenabweichung und versucht ausschließlich bei solchen Schnitten
und ansonsten vollständig bestandener Topologieprüfung höchstens drei
strengere Stufen: 0,125 / 0,0625 / 0,03125 mm. Die Schnittprüfung bleibt
verbindlich; nicht auflösbare Schnitte und andere Netzfehler werden abgewiesen.
Die exakte Quelle, Flächenkennungen, Randbedingungen und SI-Messgrenzen werden
dabei nicht geändert. Diese Vorschauprüfung ersetzt die spätere isotrope
Vernetzung und ihre vollständige BEM-Netzprüfung nicht.

`freecad_protocol` bestand am 24. Oktober 2026 einschließlich des gespeicherten
nativen Rezessfalls in 6,77 s: 24 Schnittpaare bei 0,25 mm wurden bei 0,125 mm
zu null Schnittpaaren bei 25.376 Vorschaudreiecken. BRep-Quelldaten und sieben
Flächengruppen blieben erhalten; bewusst falsche Flächenmetadaten wurden mit
präziser SI-Abweichung abgewiesen. Der echte native Import und eine
Ringänderung um 0,5 mm bestanden außerdem die Prüfung erhaltener
Kalotten-/Antriebs- und Aperturidentitäten sowie Speichern/Öffnen des Projekts.
