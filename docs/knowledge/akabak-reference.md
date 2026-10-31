# AKABAK: externe Funktions- und Bedienreferenz

Die Untersuchung betrifft die lokal installierte AKABAK-Version 3.3.2b144.
Sie ist ergänzende Referenz neben dem primär beauftragten Boundary Lab; kein
Solver, Quellcode oder Runtime-Baustein von AKABAK wird in HornSculpt
übernommen. Die nachfolgenden Aussagen stammen aus der installierten Hilfe
und statischen Ressourcen, nicht aus einem Vergleich zweier Rechenläufe.
Öffentliche Herstellerinformationen: [AKABAK](https://www.randteam.de/AKABAK/AKABAK.html)
und [VACS-Dokumentation](https://www.randteam.de/VACS/VACS-Docs.html).

## Prüfumfang und Identität

- EXE 52.755.968 Byte, SHA-256
  `605b802202249268d40a653c5793e4d5d3e4f1acbbfadf71cf6c6f8013610c11`.
- Installierte CHM-Hilfe, SHA-256
  `9399290b533672592100a0ced11c233ac1020789fee72e15834962937d5dc28d`.
  Extrahierte Kapitel im externen QA-Ordner `akabak-help-3.3.2/Chapters`:
  `Observation-Part/Form-Fields (Nodes).html`,
  `Form-Fields (Balloon).html`, `Form-Range (frame).html`,
  `Desktop/Page-Frequencies.html`, `Desktop/Form-FieldOutput.html`,
  `Appendix/Boundary Element Analysis.html` und
  `Appendix/Radiation Impedance.html`.
- REA-PE-Ressourceninventar
  `ev_a46c3a564382b2c9fda76359ea5633c337ef7e373980f77148bbed806065f1ff`:
  786 Ressourcen und 1582 Inventareinträge. 133 Stringtabellenblöcke mit
  2085 Strings wurden lokal als `akabak-resource-strings.json` dokumentiert.
  Die Strings belegen vorhandene Ressourcen, keine tatsächlich ausgeführte
  Funktion. Eine Ghidra-Suche lief in einen Timeout; keine Dekompilierung oder
  verifizierte interne Solverimplementierung wird behauptet.

## Beobachtung und Konventionen

Die Hilfe unterscheidet eine Lösungs- und eine Beobachtungsphase. Feldpunkte
lassen sich auf Polygonflächen oder räumlichen Schalen definieren. Für die
Darstellung stehen Quellenwahl, Frequenz sowie Skalargröße und Farbbereich
zur Verfügung; die Knotendichte ist ein eigenes Abtastproblem. Die Hilfe
beschreibt komplexen Feldexport und eine RMS-Mittelung bei Auswahl mehrerer
Frequenzen. Daraus folgt nicht, dass HornSculpt bereits alle diese
Bedienoptionen besitzt.

AKABAK verwendet laut seiner Hilfe Spitzenwerte und für den Pegel die
Referenz `√2 · 20 µPa`. HornSculpt verwendet komplexe RMS-Drücke mit
`20 µPa`; bei konsistenter Umrechnung liefern beide denselben dB-SPL.
Normalen zeigen laut AKABAK-Hilfe in das betrachtete akustische Gebiet.
Beim Vergleich müssen Normalrichtung, Quellenzeichen, Einheiten, Zeitphase,
Geometrie und Medienangaben deshalb explizit angepasst werden.

AKABAKs dokumentierte Fernfeldoption evaluiert in endlicher großer Entfernung,
standardmäßig dem Hundertfachen der Strahlergröße, und skaliert anschließend
Amplitude und Phase auf den angegebenen Bezug. Das ist nicht automatisch
identisch mit HornSculpts asymptotischem Fernfeld. Ein Vergleich benötigt
gleiche Beobachtungsrichtung, Abstandsnormierung und Phasenbezug.

## Kopplung und Grenzen

Die Hilfe beschreibt einen Nicht-Eindeutigkeitskompensator `NUC`, ohne hier
dessen mathematische Implementierung nachzuweisen. Insbesondere wird aus
dieser Bezeichnung kein Burton–Miller-Verfahren abgeleitet. Für die Kopplung
beschreibt die Hilfe flächengewichtete Eigen- und gegenseitige
Strahlungsimpedanzen; bei nicht einheitlichen Geschwindigkeitsformen ist die
leistungsbezogene Definition wichtig. Diese fachliche Referenz ersetzt weder
einen eigenen Porttest noch einen Vergleich identischer Daten.

## Subdomain Modeling: Befund und HornSculpt-Vertrag

Prüfstand: 31.10.2026. Untersucht wurden die installierte Hilfe, das lokale
Beispielpaket, das erneut mit REA geprüfte PE-Ressourceninventar und die
bestehenden HornSculpt-Verträge. Das ist eine Funktionsanalyse und ein
Implementierungsauftrag. Der damalige Analysebefund ist unten erhalten; der
inzwischen implementierte und geprüfte Stand ist unter
[Akustische Teilgebiete](../SUBDOMAINS.md) dokumentiert. Die REA-Analyse allein
ist kein numerischer Nachweis.

### Quellen und Evidenzgrenzen

Die CHM-Quellen sind relativ zum Wurzelverzeichnis der installierten Hilfe:

- `Chapters/Appendix/Sub-Domain Modelling.html`
- `Chapters/Appendix/Normals and Volumes.html`
- `Chapters/Appendix/Normals and Volumes - Interface Elements.html`
- `Chapters/Appendix/Meshing.html`
- `Chapters/Appendix/Infinite Baffle.html`
- `Chapters/Appendix/Boundary Element Analysis.html`
- `Chapters/BEM-Part/Form-Subdomain.html`, `Form-Interface.html` und
  `Form-Infinite-Baffle.html`
- `Chapters/Desktop/Form-NewComponent.html` und `Form-CreateSubInt.html`.

Die Hilfe hat den oben festgehaltenen SHA-256. Das installierte Readme nennt
3.3.2 b144 vom 26.06.2026. Das erneut ausgeführte REA-Inventar ergibt dieselbe
EXE-Identität und dieselbe Evidenz-ID
`ev_a46c3a564382b2c9fda76359ea5633c337ef7e373980f77148bbed806065f1ff`.
Die daraus bereits dekodierten RT_STRING-Daten enthalten insbesondere
String-ID 63475 an Datei-Offset 40535566 (verschiedene Teilgebiete erforderlich)
und ID 63486 an Offset 40536516 (fehlendes Teilgebiet). Die Ressourcen belegen
Diagnosetexte; ihre Ausführung oder konkrete Validierungslogik folgt daraus
nicht. Eine zusätzliche Ghidra-Stringabfrage wurde ohne auswertbaren Befund
abgebrochen; interne Kopplungsmatrizen werden nicht behauptet.
Der REA-Zugriff auf das OneDrive-ZIP scheiterte mit „target changed while being
read“; die Beispielbefunde stammen aus den vorhandenen entpackten Dateien.

Die offizielle [Einführung](https://randteam.de/AKABAK3/AKABAK-Introduction.html)
bestätigt Subdomain Modeling und vollständig gekoppelte LE-/BE-Teile. Die
[Studienseite](https://randteam.de/AKABAK3/AKABAK-Studies.html) nennt unter anderem
SP38 mit Innenraum, Port und Außenraum sowie ein Basshorn mit sieben
Teilgebieten. Einige dortige Studien betreffen ausdrücklich den Vorgänger
ABEC; daraus folgt kein Versions- oder Laufnachweis für aktuelle Beispiele.

In der vorhandenen AKABAK-Sitzung war ein anderes Nutzerprojekt mit einem
Solve-Protokoll sichtbar. Das Projekt wurde nicht verdrängt. Die folgenden
Dialogschritte stammen aus der Hilfe, nicht aus einer ausgeführten
Bedienprüfung des Horn-Beispiels. Es wurden keine Beispielrechnungen,
Laufzeitvergleiche oder Ergebnisgleichheiten nachgewiesen.

### Fachmodell laut Hilfe

Ein Interior-Teilraum beschreibt ein endliches berandetes akustisches Volumen,
ein Exterior-Teilraum den unbeschränkten Strahlungsraum. Ein Interface ist
eine eigene geometrische Komponente mit genau zwei Subdomain-Referenzen.
Seine Fläche ist akustisch transparent: Druck und physische
Normalgeschwindigkeit sind stetig. Zugeordnete Interface-Dreiecke gehören
beiden Teilräumen; gewöhnliche Randdreiecke gehören einem Teilraum.
Wall oder Driven auf Interface-Kindkomponenten wird laut Hilfe ignoriert.
HornSculpt soll solche widersprüchlichen Rollen stattdessen sichtbar ablehnen.

AKABAK orientiert gewöhnliche Randnormalen in den betrachteten Luftraum.
Interface-Normalen zeigen einheitlich zur ersten referenzierten Subdomain.
Ein Tausch beider Referenzen zusammen mit Normalenumkehr ist äquivalent.
Diese Konvention darf nicht still mit HornSculpts vorhandener Körpernormalen-
und Driven-Konvention vermischt werden.

Der Graph aus Teilräumen und BE-Interfaces bildet zusammenhängende akustische
Domains. Getrennte Domains können weiterhin über das LE-Netz gekoppelt sein.
Mehrere Körper im selben Außenraum sind daher etwas anderes als mehrere
Luftgebiete. Eine Schallwand ist ebenfalls kein allgemeiner Interface-Ersatz.
Eine Infinite Baffle gehört zu einem Exterior-Teilraum und verwendet laut
Hilfe eine besondere Green-Funktion. Geometrie hinter dieser Ebene kann eine
weitere Teilraumaufteilung verlangen; daraus folgt keine Schallwandpflicht
für jeden Hornmund.

### Dokumentierter Bedienablauf

| Schritt | AKABAK laut Hilfe | Folgerung für HornSculpt |
| --- | --- | --- |
| 1 | Subdomain im BEM-Baum anlegen, Interior/Exterior wählen | Benannte Luftgebiete mit stabiler ID und explizitem Typ |
| 2 | Geometrie unter dem gewählten Parent erzeugen oder über Change Subdomain zuordnen | Draft und CAD verwenden denselben Flächen-/Gebietsvertrag |
| 3 | Interface erzeugen, Subdomain 1 und 2 sowie eigene Geometrie zuordnen | Zwei geprüfte Seiten einer gemeinsamen transparenten Fläche |
| 4 | Struktur und Normalen darstellen, Gebietsabschlüsse prüfen | Auswahl und hervorgehobene Gebietsseite samt Normalpfeilen |
| 5 | Globale und lokale Meshing-Vorgaben verwenden | Gemeinsames Interface-Netz, lokale skalare Größenfelder |
| 6 | Lösen, danach Beobachtung und LE-Kopplung auswerten | Revisionsgebundene Randdaten und gezielte Wiederverwendung |

Die Hilfe beschreibt lokale Netzeinstellungen für Subdomains und Interfaces
mit globalem Rückfallwert. Ein Interface besitzt ein gemeinsames Netz; daraus
ist kein Beleg für beliebige nichtkonforme Kopplung zweier separat definierter
Netze abzuleiten. Laut Hilfe werden zuerst Teilräume und anschließend
Interfaces gelöst. Das ist eine Funktionsbeschreibung, kein Nachweis eines
bestimmten Schurkomplement-, Iterations- oder Burton–Miller-Verfahrens.
Auch die dokumentierte spätere LE-Kopplung ist keine pauschale Zusage, dass
beliebige Parameteränderungen ohne erneute BEM-Rechnung möglich sind.

### Konkrete lokale Beispiele

Alle Pfade sind relativ zum externen Beispielverzeichnis
`AKABAK-Examples`; die proprietären Projekte werden nicht ins Tool kopiert.

| Beispiel | Beobachtet | Geeignete Prüfung |
| --- | --- | --- |
| `BEM/Waveguides/Horn Simple IB.akp` | Bereits External und Interface-IB vorhanden | Basis für die zusätzliche Hornunterteilung |
| `BEM/Waveguides/Horn Simple IB (subdomains).akp` | Horn throat, Horn bell, External; Throat-Bell und Bell-External | Zusätzliche transparente Aufteilung ohne Änderung der Physik |
| `BEM/Loudspeakers/TL Subdomains (free).akp` | Transmissionline-Geometrie, Front-/Rückseitenreferenzen und LE-/BE-Feldnamen | Mehrgebiets- und Membrankopplung als späterer Folgefall |
| `Tutorials/Speaker Cabinet SP20 Sealed/(1) LE only.akp` bis `(5) LE BEM (box + IB).akp` | Benannte stufenweise Tutorialfolge | Getrennte Innen-/Außenluft mit gemeinsamer Treiberlast |
| `BEM/Scripting/Waveguide H800/H800_Driver.txt` und `H800_Interface-Exterior.txt` | Explizite Driving-Gruppe/Referenz und Interface-Geometrie | Flächengruppen und nachvollziehbare Geometriezuordnung |

Die AKP-Dateien sind binär. Lesbare Namen belegen ihre enthaltenen Begriffe,
nicht die vollständig rekonstruierte Verdrahtung oder komplexe Randlösung.
Im H800-Treibertext stehen die Driving-Gruppe in Zeile 15 und ihre Referenz
in Zeilen 22–24, im Interface-Text beginnt die Elementdefinition bei Zeile 15.

SHA-256 der drei Hauptfälle, in Tabellenreihenfolge:

1. `0ba0ac6feaf4d77065f48a3b5259925d15aa0001365ce1b1a1ab04f870d7d5d4`
2. `bc9c5832c592586d31504e8fd0267b604fe5bcb4c1262afa71476145f90dad7d`
3. `6bdec97f852f152781fd1ef5ada3468d6f8490a6f05e9d567b0d53a734579ab6`

### Abgleich mit HornSculpt

`modules/core/include/hornsculpt/mesh.hpp` beschreibt ein SurfaceMesh mit
Wall/Driven, aber keine Gebietsseite. Der Vertrag in
`modules/bem/galerkin/exterior_domain.hpp` verlangt ausdrücklich denselben
unbeschränkten Luftraum für alle Körper und weist verschachtelte Hüllen und
geschlossene Luftkammern ab. `modules/bem/include/hornsculpt/bem.hpp` bindet
einen Solverauftrag an ein Medium und ein Netz. Die Baugruppenkomposition ist
keine akustische Gebietszerlegung.

`modules/lem/include/hornsculpt/lem.hpp` besitzt gewichtete Ports und eine volle
gegenseitige Impedanzmatrix. `modules/bem/lem_coupling.hpp` bindet sie an den
bestehenden einzelnen BEM-Raum; der dortige Interfaceport ist kein allgemeines
BE↔BE-Interface. `apps/desktop/project_file.hpp` und die zugehörigen
Simulations-/LEM-JSON-Verträge benötigen für Mehrgebiete zusätzliche Identität
und Revisionsdaten. Die bisherige Außenraumvalidierung bleibt verbindlich,
bis ein echter Ersatz samt Prüfnachweisen vorliegt.

### Eigener Vertrag und erste Abnahme

Die folgenden Punkte sind HornSculpt-Entscheidungen, keine rekonstruierten
AKABAK-Interna:

1. Typisierte Gebiets- und Interface-IDs; eindeutiger Interior-/Exterior-Typ,
   Medium und orientierte Randseiten. Für den ersten Vergleich ein homogenes
   Luftmedium verwenden; allgemeine Medienübergänge separat herleiten und prüfen.
2. Ein Interface referenziert zwei verschiedene Gebiete und genau eine
   gemeinsame Geometrie/Trace-Diskretisierung mit zwei orientierten Seiten.
   Bei nach außen aus jedem Luftgebiet gerichteten Normalen gilt
   `p_A = p_B` und `v_n,A + v_n,B = 0`. Die Abbildung bestehender Körpernormalen
   auf diese Seiten wird explizit gespeichert bzw. abgeleitet und getestet.
   Eine skalare uniforme LEM-Portgeschwindigkeit ersetzt die räumliche
   Interface-Lösung nicht.
3. Jede endliche Gebietshülle einschließlich Interfaces muss geschlossen,
   orientiert und schnittfrei sein. Beim Außenraum gilt zusätzlich die
   Strahlungsbedingung im Unendlichen. Eine interne Interfacefläche ist ein
   bewusster Gebietsabschluss, keine künstliche Reparatur einer Modelllücke.
   Die Vereinigung aller Flächen ist nicht als gewöhnliche einzelne
   geschlossene SurfaceMesh-Hülle zu prüfen: gemeinsame Interfaces müssen
   in den getrennten Gebietsrandansichten validiert werden.
4. Zuerst konforme gemeinsame Interfaces; keine stillen Nächstpunkt-
   Zuordnungen. Nichtkonforme Kopplung erst mit eigener geprüfter Projektion
   und Leistungserhaltung. Keine Übernahme der in AKABAK-Hilfe diskutierten
   Tolerierung kleiner Netzlücken.
5. Innenraum-/Außenraumoperatoren und Kopplungsbedingungen gemeinsam herleiten.
   Burton–Miller bleibt der produktive Außenraumpfad. Den bisherigen Operator
   nicht ohne Herleitung auf geschlossene Innenräume anwenden: physische
   Innenraumresonanzen von fiktiven Außenraumeigenfrequenzen unterscheiden.
6. Erste nutzbare Strecke: einfaches Horn mit einer zusätzlichen transparenten
   Schnittfläche teilen und gegen dasselbe ungeteilte physische Modell prüfen.
   Danach Kammer–Port–Außenraum und beidseitig LEM-gekoppelte Membran prüfen.
   Netz-/Quadraturkonvergenz, komplexen Druck, Impedanz, lokale Druck-/Flussreste,
   Energiebilanz und Orientierungstausch prüfen. Toleranzen vorher begründen;
   alle Teilgebiete zusammen höchstens 10.000 Dreieckelemente je Prüffall.
7. LEM-Ports an Gebiet und Randseite binden; volle gegenseitige Last und
   bekannte Anregungen erhalten. Unabhängige Domains dürfen über das Netzwerk
   gekoppelt sein. Keine doppelte Flächenanregung oder doppelte Luftlast.
8. Modellbaum, Flächenauswahl, lokale Netzvorgaben, Undo/Redo, Speichern/Öffnen,
   CAD-Topologieabgleich und Ergebnisse gemeinsam umstellen. Gebiets-/
   Interfaceänderungen invalidieren passende Solverstände; fehlende
   Zuordnungen verlangen sichtbare Korrektur. Feldpunkte gehören explizit
   einem Gebiet, Fernfeld ausschließlich dem unbeschränkten Außenraum.

Als erste betroffene Prüfstellen dienen `tests/exterior_domain_tests.cpp`,
`tests/bem_lem_tests.cpp` und `tests/project_file_tests.cpp`, ergänzt um echte
Mehrgebiets-, GUI- und CUDA-Prüfungen. Bestehende Ein-Gebietsprojekte müssen
reproduzierbar bleiben. Zusätzliche Testdateien und Produktionsmodule erst
mit tatsächlich verwendeter Implementierung anlegen.

### Meshimport: dokumentierter AKABAK-Ablauf und eigener Gmsh-Vertrag

**Beobachtet – mitgelieferte Hilfe.** `AKABAK.chm`, Abschnitte
`Appendix/Mesh-File-GMSH.html`, `General-Part/Form-MeshFile.html` und
`BEM-Part/Form-Elements (Mesh-File).html`: Gmsh erzeugt ein 2D-Oberflächennetz;
AKABAK dokumentiert den Export als MSH 2.2 ASCII. Im Bereich General wird ein
Mesh-File-Objekt geöffnet, seine Ausdehnung geprüft und nötigenfalls nach SI
skaliert. Tags einschließlich Physical Names helfen, Teilmengen für BEM-Elemente
auszuwählen und Subdomains oder Interfaces zuzuordnen; Flächennormalen können
umgekehrt werden. Die Hilfe unterscheidet das importierte Flächennetz vom
anschließend in AKABAK verfeinerten BEM-Netz. Die Eingangskoordinaten werden
ohne Skalierung als SI aufgefasst. Die Hilfe nennt im 3D-Modus neben Dreiecken
auch Vierecke.

**Beobachtet – REA.** Das installierte `AKABAK.exe` ist ein x86-64-PE mit
SHA-256 `605b802202249268d40a653c5793e4d5d3e4f1acbbfadf71cf6c6f8013610c11`.
Die statische PE-Ressourceninventur lieferte Evidence-ID
`ev_a46c3a564382b2c9fda76359ea5633c337ef7e373980f77148bbed806065f1ff`.
Sie belegt weder Parseralgorithmus noch Laufzeitverhalten; es wurde keine
Ghidra-Funktionsanalyse oder Ausführung des Meshimports vorgenommen. Die
untersuchte Hilfedatei hat SHA-256
`9399290b533672592100a0ced11c233AC1020789FEE72E15834962937D5DC28D`.

**Eigener Vertrag.** Der neue Importer orientiert sich an der
[offiziellen MSH-4.1-Spezifikation](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format)
und der [MSH-2-Spezifikation](https://gmsh.info/doc/texinfo/gmsh.html#MSH-file-format-version-2-_0028Legacy_0029).
Unterstützt werden ausschließlich ASCII-Dateien 2.2 und 4.1 mit linearen
Oberflächendreiecken und planaren streng konvexen Vierecken. Punkt- und Linienelemente werden geprüft und gezählt, jedoch nicht als BEM-Flächen verwendet. Andere Elemente und nicht implementierte Abschnitte werden
mit Fehler abgewiesen. 4.1 ist eine eigene Erweiterung und keine behauptete
AKABAK-Kompatibilität. Knoten- und Element-IDs, Entity-Tags, Physical-Tags
und Namen bleiben erhalten. Der Umrechnungsfaktor nach Metern wird explizit
angegeben. Dreiecksorientierungen bleiben unverändert; offene, nicht
mannigfaltige, widersprüchlich orientierte und entartete Netzstellen werden
diagnostiziert, nicht repariert. Eine Diagnose ersetzt die vollständige
Netz- und Selbstschnittprüfung vor der BEM-Berechnung nicht.
