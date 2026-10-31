# Akustische Teilgebiete

HornSculpt kann ein akustisches Modell in endliche Innenräume und einen
unbeschränkten Außenraum mit gemeinsamen, transparenten Schnittflächen
aufteilen. Jedes Dreieck einer Schnittfläche wird einmal gespeichert und
gehört zu genau zwei verschiedenen Gebieten. Die beiden Gebietsseiten verwenden
entgegengesetzte Normalen. Die Grenzfläche ist weder eine Materialwand noch
eine vorgeschriebene Anregung: Der Druck ist stetig und die nach außen
orientierten Normalableitungen sind gegensinnig.

Ein gemeinsames konformes Dreiecksnetz beschreibt beide Seiten. Physische
Flächen gehören einem Gebiet oder sind ausdrücklich als beidseitige physische
Flächen zugeordnet. Deren Drücke bleiben auf beiden Seiten unabhängig.
Der gemeinsame Partitionsprüfer kontrolliert
Gebietszuordnung, Schnittflächen, Orientierung, geschlossene Gebietsränder und
geometrische Überschneidungen. Eine Schnittfläche muss ausdrücklich vorgegeben
sein: Der Abschluss einer Modelllücke wird nicht stillschweigend als Interface
interpretiert. Für offene, ausdrücklich beauftragte planare Randkonturen kann
der Partitioner eine konturtreue Schnittfläche vernetzen; importierte MSH-Netze
können ihre Schnittflächen ebenfalls explizit ausweisen.

## Bedienung und Gmsh-Import

**Datei → Gmsh-Netz importieren …** öffnet eine `.msh`-Datei. Die Dateieinheit
Meter oder Millimeter muss ausdrücklich gewählt werden; intern werden Meter
gespeichert. Der Import unterstützt ASCII MSH 2.2 und 4.1 mit linearen
Dreiecken, ebenen streng konvexen Vierecken sowie geprüften Hilfspunkten und
Linien. Vierecke werden orientierungstreu entlang einer Diagonalen geteilt;
die Ursprungs-Element-ID bleibt an beiden Dreiecken erhalten. Hilfselemente
werden gezählt und nicht als BEM-Flächen verwendet. Binary-, Volumen-,
höhergradige, periodische und partitionierte Eingaben werden ausdrücklich
abgewiesen. Eine solche Datei vor dem Import als unterstütztes Oberflächennetz
exportieren; es erfolgt keine stille Geometrie- oder Formatänderung.

Für eine eindeutige Flächenauswahl müssen alle Dreiecke derselben geometrischen
Entity dieselben Physical-Tags besitzen. Andernfalls den Export in getrennte
Entities aufteilen; der Import meldet den Konflikt. Der Projektvertrag begrenzt
den Import auf eine Million Knoten und zwei Millionen Einträge je Abschnitt
sowie zwei Millionen resultierende Dreiecke. Das Einlesen läuft im Hintergrund;
ein Einheitenwechsel verwirft den noch laufenden Import.

Das importierte Oberflächennetz ist eine polygonale Offline-Quelle im
bestehenden Import-/FreeCAD-Link-Arbeitsbereich. Es wird weder als exakte
BRep-Geometrie noch als laufende FreeCAD-Verbindung ausgegeben. Ursprüngliche
Knoten- und Element-Tags, geometrische Entity-IDs, Physical-Tags und Namen,
Skalierung und Quellnetz werden mit dem Projekt gespeichert. Die Entity-IDs
liefern die stabilen physischen Flächengruppen. Zunächst sind alle Gruppen
Wall; Antriebsflächen werden in der Flächentabelle bewusst auf Driven gesetzt.
Eine bestehende FreeCAD-Verbindung wird nicht automatisch beendet; vor dem
Wechsel auf eine Gmsh-Quelle ist sie zu trennen.

Unter **Akustik → Berechnung → Teilgebiete …** oder **BEM-Teilgebiete …**
in der Flächenansicht wird der Modellbaum bearbeitet. Innenräume, Außenraum
und Interfaces besitzen eigene Namen und IDs. Jedes Interface nennt zwei
verschiedene Gebiete. Benannte Elementgruppen unter diesen Einträgen fassen
mehrere Flächen zusammen. Die Mehrfachauswahl aus der Flächentabelle wird
vorselektiert; übernommen wird sie erst durch die ausdrückliche Zuordnung.
Gruppennamen ändern keine Randbedingung. Gruppen bleiben disjunkte Teilmengen
der tatsächlichen Flächenzuordnung.

Beim Öffnen eines Projekts mit expliziten Interfaces sind deren Drahtnetz
und Normalen sichtbar. **Interface-Normalen** zeigt cyanfarbene Pfeile `+n`
in Subdomain 1 und orangefarbene Pfeile `−n` in Subdomain 2. Die Legende nennt
die zugehörigen Gebietsnamen. Die Richtungen folgen der tatsächlichen
Dreiecksorientierung, auch nach einer bewussten Flächenumkehr; bei großen
Netzen werden die Pfeile je Flächengruppe ausgedünnt.

**Modellebenen** zeigt die konfigurierte unendliche Außenraum-Schallwand
und die beiden Symmetrieebenen als farbige Gitter. Ihre endliche gezeichnete
Ausdehnung ist symbolisch. Ein Schallwandpfeil zeigt in den abstrahlenden
Außenraum, die Symmetriepfeile in die behaltenen Halbseiten. Diese Angaben
werden bereits aus der aktuellen Modellkonfiguration angezeigt, ohne einen
Solve zu starten. Das Ein-/Ausblenden verändert weder Netz und Kamera noch
Revision oder Ergebnisgültigkeit.

**Ganzes Modell** spiegelt das Viertelmodell an beiden Symmetrieebenen und
zeigt das vollständige Horn. Bei aktiver XY-Symmetrie ist diese Ansicht
standardmäßig eingeschaltet; zum ursprünglichen Viertelnetz den Haken
entfernen. **Einpassen** passt die Kamera an die jeweils sichtbare Ausdehnung
an. Interfaces, Normalen, Farben und Flächenauswahl erscheinen auch in den
Spiegelbildern. Eine dort gewählte Fläche bleibt dieselbe Ursprungsfläche.
Die Symmetrieansicht erweitert weder das gespeicherte Netz noch den
Solverauftrag. Berechnete Schallfeldebenen bleiben an ihren tatsächlichen
Weltpositionen und werden nicht als unbelegte Feldwerte gespiegelt.

Abbrechen übernimmt keine Eingaben. Vorhandene IDs bleiben beim Bearbeiten
erhalten. Numerisch gleiche Face-IDs eines anderen CAD-Dokuments sind
keine Identität: Nach einem unklaren Flächenwechsel bleibt die Zuordnung
sichtbar zur Prüfung markiert und blockiert die neue Rechnung.

- **Geschlossene physische Quelle ohne Interfaceflächen:** Die gewählten
  Innenflächen definieren einen Luftraum; seine planaren Öffnungsränder erhalten
  gemeinsame, konturtreu vernetzte transparente Schnittflächen. Die physische
  Quelle muss vorher bereits geschlossen und schnittfrei sein.
- **Schon im Netz enthaltene Interfaceflächen:** Ihre IDs getrennt zuordnen.
  Das Interface verbindet seine explizit gewählten Gebiete. Die Dreiecksnormale
  zeigt in das als Gebiet 1 gewählte Gebiet;
  Normalenanzeige und bewusste Flächenumkehr verwenden. Alle Gebietsränder müssen
  einschließlich der angegebenen Interfaces vollständig geschlossen sein.
  Verbleibende Lücken werden nicht automatisch ergänzt.

**Importnetz prüfen** prüft bei Gmsh die tatsächliche importierte Diskretisierung
einschließlich der eingestellten Winkel-/Formgrenzen. Es verändert ihre Knoten
und Dreiecke nicht. Zielkantenwerte ersetzen daher kein entsprechend feines
Gmsh-Netz. Netzstatistik und Diagnose zeigen die tatsächlichen Kanten; eine
globale Dreifachkante an einem Interface ist von ungültigen Gebietsrändern zu
unterscheiden. Jeder Gebietsrand wird separat geprüft. Zu grobe Netze oder
widersprüchliche Vorgaben bleiben Fehler und erhalten keinen Erfolgsstatus.

Für die Rechnung **3D** und den **dichten GPU-Solver** verwenden. Nach einer
erfolgreichen Rechnung zeigt der Modellbaum die gelösten Gebiete und
Interfaces. **Interfaces** blendet die wirklichen Schnittflächen als
cyanfarbenes Drahtnetz ein; importierte Interfaceflächen bleiben über die
Flächentabelle zugänglich. Dies verändert weder Anregung noch Revision.
Berechnete Feldpunkte auf einer transparenten Fläche erhalten deren stetigen
Druck; die Fläche ist keine Materialwand. Physische LEM-Anschlüsse können
nach dem unten beschriebenen Vertrag an verschiedene Gebiete gekoppelt werden.

## Numerischer Vertrag

Die Darstellungen werden durch den OpenGL-Test `acoustic_guides_render`
geprüft: Schallwand und beide Symmetrieebenen, Normalenumkehr, Sichtbarkeit,
Entfernen der Hilfen sowie unveränderte Kamera, Flächenauswahl und Quellnetze.
Der native GUI-Test prüft die Anzeigen vor und nach einer CUDA-Rechnung,
das Wiederöffnen und den Wechsel zu einem Projekt ohne diese Hilfen.

### Viertelsymmetrie und Außenraum-Schallwand

**XY-Symmetrie** verwendet ein Viertelnetz mit Ebenen x=0 und y=0. Die
behaltenen Halbseiten werden eindeutig aus dem Netz bestimmt und angezeigt;
ein Netz auf beiden Seiten einer Ebene wird abgewiesen. Die metrische
Ebenentoleranz ist ausdrücklich einstellbar. Nur interne Arbeitskoordinaten
innerhalb dieses Budgets werden an einer Naht projiziert; die importierte
Quelle bleibt erhalten. Es werden keine künstlichen Abschlussdreiecke erzeugt.

Die **Schallwand des Außenraums** ist ein eigenes Objekt mit Ursprung,
Normale in den abstrahlenden Halbraum und Ebenentoleranz. Sie erlaubt räumliche
Interfaces mit offenem Rand auf der Ebene und ist vom bisherigen Vertrag einer
planaren Rayleigh-Apertur zu unterscheiden. Der Außenraum verwendet zusätzlich
den schallharten Spiegeloperator. Innenräume erhalten vier XY-Bildbeiträge,
der Außenraum mit Schallwand bis zu acht. Singular- und Nahpaare an Ebenennaht
und Kreuzungsachse berücksichtigen gemeinsame Bildknoten. Druckfreiheitsgrade
bleiben im Sektor; ein vervierfachtes Solvernetz wird nicht aufgebaut.

Der Desktoppfad unterstützt diese Kombination für 3D, dichten CUDA-Solver und
Fernfeld mit symmetrischer vorgegebener Anregung. Endliche Feldpunkte und
LEM-Kopplung werden für diesen speziellen Symmetrie-/Schallwandpfad abgewiesen.
Die übrigen Einraum- und Mehrgebietsverträge behalten ihre jeweils dokumentierten
Funktionen. Projektformat 18 speichert Gebiete, Interfaces, Gruppen, Ebenen und
Ergebnissnapshot; ältere Formate werden beim Lesen migriert.

Explizit beidseitige starre Driven-Flächen bewegen sich körperlich gemeinsam:
die in die Luft gerichtete Geschwindigkeit der zweiten Seite erhält das
Gegenzeichen. Ihre Last ist die Druckdifferenz, ihre Fläche wird physisch einmal
gezählt. Im Viertelmodell beziehen sich Gesamtfluss und Volumenimpedanz auf
das vollständige symmetrische Modell. Eine LEM-Druckdifferenzkopplung solcher
beidseitigen Ports ist nicht freigegeben.

Die Zeitkonvention ist `exp(+i ωt)`, die auslaufende Grundlösung lautet
`G(r) = exp(−ikr)/(4πr)`. Die Normale `n` zeigt jeweils aus dem betrachteten
Luftgebiet heraus. Mit Druck `p`, Normalableitung `q = ∂p/∂n`, Wellenzahl `k`
und den Randoperatoren `V`, `K`, `K′` und `W` lauten die verwendeten
Gebietsgleichungen:

```text
C(p,q) = (½ + K)p − Vq = 0
N(p,q) = Wp + (K′ − ½)q = 0
C − (i/k)N = 0
```

Die letzte Zeile ist die kombinierte Galerkin-Gleichung nach
[Burton und Miller, 1971](https://doi.org/10.1098/rspa.1971.0097). Der
[Einraumvertrag](BURTON-MILLER.md) dokumentiert dieselbe Außenraumphysik mit
entgegengesetzter Normalenrichtung; seine Vorzeichen sind entsprechend
umzurechnen. Die Operatoren und Singular- und Nahquadraturen werden unabhängig
von AKABAK im Projekt implementiert. Der Druck wird auf jedem Gebietsrand stetig linear (`P1`)
diskretisiert. Auf einem gemeinsamen Interface werden die zugehörigen
Druckknoten stark identifiziert; der räumlich lineare Interfacefluss besitzt
eigene gemeinsame `P1`-Freiheitsgrade. Für die beiden nach außen gerichteten
Gebietsnormale erhält er entgegengesetztes Vorzeichen. In der Matrix ist die
unbekannte Flussgröße als `q/k` skaliert. Bekannte physische Wand- und
Driven-Flüsse bleiben dreiecksweise konstant (`DP0`). Für eine physische
Geschwindigkeit `v` positiv in die Luft gilt mit der Luftdichte `ρ`
`q = +i ρωv` bezüglich der aus der Luft gerichteten Normalen.

Reguläre Panelpaare werden frequenz- und abstandsabhängig qualifiziert; nicht ausreichend qualifizierte Paare verwenden adaptive Integration. Die paarweise Integration verwendet reguläre, singuläre und adaptive
Nahquadratur. Die CUDA-Assembly verarbeitet die Quadraturpunkte in begrenzten
Paketen; das komplexe dichte Gleichungssystem wird mit CUDA-LU gelöst. Das
angegebene lineare Residuum wird anschließend aus der unskalierten Gleichung
bestimmt. Es ist kein Nachweis des Diskretisierungs- oder Quadraturfehlers.
Aus dem gelösten `P1`-Druck und dem Fluss werden gebietsbezogene Cauchy-Daten
gespeichert. Endliche Feldpunkte verwenden die Grenzdaten ihres Luftgebietes
und adaptive Quadratur. Ein Punkt unmittelbar auf einer transparenten
Schnittfläche erhält deren interpolierten `P1`-Druck; der Abstandsschutz an
physischen Wänden und Interface-Rändern bleibt erhalten. Fernfeldwerte werden
aus dem asymptotischen Integral des Außenraums berechnet.

## LEM-Kopplung

Physische `drivenSurface`-Ports können in verschiedenen BEM-Gebieten liegen;
ein einzelner Port muss vollständig einem Gebiet zugeordnet sein. Seine
Flächen und gegebenenfalls Geschwindigkeitsgewichte bestimmen eine
leistungskonjugierte Druckprojektion. Für jeden Port wird eine räumliche
Einheitsflussanregung gelöst. Alle Eigen- und gegenseitigen Druckantworten
bilden die komplexe Lastmatrix des elektrischen, mechanischen und akustischen
LEM-Netzwerks. Damit kann eine Membran mit Front- und Rückseite an zwei
verschiedene BEM-Luftgebiete gekoppelt werden.

Nach dem LEM-Solve berechnet HornSculpt die resultierenden physischen
Geschwindigkeiten und löst das BEM-System mit diesen Anregungen erneut. Das
veröffentlichte Residuum und Feld gehören zu dieser tatsächlichen Endlösung.
Ein zusätzlicher Druckschluss vergleicht ihre modal gewichteten Drücke mit den
in das Netzwerk zurückgeführten Portdrücken. Unabhängige vorgegebene
Anregungen gleichzeitig mit der LEM-Kopplung werden derzeit abgewiesen.
Virtuelle BEM↔LEM-`interfacePort`-Flächen sind im Mehrgebietspfad noch nicht
unterstützt; sie sind von transparenten BEM↔BEM-Interfaces zu unterscheiden.

## Daten und Grenzen

`BemSubdomainModel` speichert Gebiets- und Interface-IDs, Namen, Gebietsarten,
Dreieckszuordnungen und orientierte Gebietsseiten. `normalSign` multipliziert
die Dreiecksnormale und ergibt die aus dem jeweiligen Luftgebiet gerichtete
Normale. Der Solver erhält das Netz und den Modell-Snapshot beim Konstruktor;
der Request muss denselben Teilgebietsstand enthalten. Ergebnisse enthalten
pro Gebiet explizite Knoten- und Dreiecksindizes mit Druck und Normalableitung
sowie pro Interface die Knotenwerte des linearen Flusses. Veröffentlichte
Gebietstraces sind unveränderlich und werden beim Kopieren eines Ergebnisses
gemeinsam genutzt. Der ältere Einraum-`boundary()`-Vertrag bleibt getrennt.

Der aktuelle Pfad verwendet ein homogenes Medium, dreidimensionale konforme
Interfaces und den dichten CUDA-Solver. Nichtkonforme Netze, mehrere Interfaces
mit gemeinsamen Knoten, allgemeine Medienübergänge, FMM, der rotationssymmetrische
`m=0`-Modus und die unendliche Schallwand sind für diesen Mehrgebietspfad nicht
freigegeben. Der Netz- und Rechenaufwand wächst entsprechend dem dichten
System; eine Leistungszusage folgt daraus nicht. AKABAK dient ausschließlich
als externe Funktions- und Bedienreferenz und ist keine Produktabhängigkeit.

## Bisheriger Nachweis

Für die neue Modellsteuerung bestanden native Dialog-/Undo-Prüfungen und der
Windows-Ablauf eines Viertelmodells mit räumlichem Interface und Außenraum-
Schallwand bis zur CUDA-Rechnung, Speicherung und aktuellen Wiederöffnung.
Der Archivtest prüft alle 147 H/V-Beobachterplätze einschließlich gültiger
und ungültiger Halbraumrichtungen. Maskierte Richtungen werden in Diagrammen
und Export ausgelassen. Die vollständigen früheren Projektdateifälle blieben
nach der Formatmigration lesbar.

Ein offener Außenrand direkt auf der Schallwand wurde mit 64 Sektordreiecken
gegen ein Vollmodell mit 256 Dreiecken bei 180, 550 und 1.200 Hz geprüft:
komplexer Druck, Fernfeld und Impedanz stimmen innerhalb 0,1 % überein.
Zusätzliche Prüfungen betreffen getrennte Drücke an beidseitigen Flächen,
Sektorrandtopologie und unzulässige Masken. Die Übereinstimmung mit dem
Vollmodell ersetzt keinen unabhängigen absoluten Netzkonvergenznachweis.

Die folgenden Werte stammen aus tatsächlich ausgeführten Prüfungen des
jeweiligen Entwicklungsstands auf einer NVIDIA RTX 5070 Laptop GPU:

- Bei einer analytischen Kugel sank der komplexe Feldfehler über drei
  Netzstufen von etwa 38,8 % auf 13,4 % und 3,72 %.
- Für einen geteilten Cup-/Hornfall mit offenem Interface-Rand betrug die
  Abweichung zum ungeteilten Fall etwa 1,17 %; das lineare Residuum lag bei
  etwa `4,1·10⁻¹⁵`.
- Der LEM-Test mit Ports in verschiedenen Gebieten und beidseitiger
  Membranrückwirkung bestand.

Zusätzlich bestanden der Dreigebietstest mit einem Feldpunkt genau auf einem
transparenten Interface, die aktive Leistungsbilanz der Kugel bei ka ungefähr
1, π und 2π (1.024 Dreiecke insgesamt), der Vergleich von Singularordnungen
4/6 sowie engerer regulärer/Nahquadraturtoleranzen 10⁻⁴/10⁻⁶. Der unabhängige
CPU-Fernpaartest vergleicht P1×P1-Momente der Grundlösung und ihrer drei
Gradientenkomponenten mit einer Grad-14-Produktquadratur. Die nachträgliche
Feldprüfung verwendet benachbarte Quadraturgrade und dieselbe an das
Residuumziel gebundene Grenze wie der Einraumpfad; bei Nichterreichen entsteht
kein erfolgreicher Ergebnisstand.

| Prüfgröße | Gemessener relativer Unterschied |
| --- | --- |
| Horn, reguläre/Nahquadraturtoleranz 10⁻⁴/10⁻⁶ | höchstens 1,07·10⁻⁸ |
| Horn, Singularordnung 4/6, komplexer Druck | höchstens 0,059 % |
| Kugel, Feld geteilt/ungeteilt über drei Frequenzen | höchstens 0,164 % |
| Kugel, aktive Eingangs-/Fernfeldleistung | höchstens 0,083 % |

Die Kugelfrequenzen sind etwa 545,9, 1.715 und 3.430 Hz, der physische
Kugelradius 0,1 m und der transparente Interface-Radius 0,15 m.

Der native Gmsh-Quellenablauf bis CUDA-Rechnung, Speichern und Wiederöffnen,
separater LEM-Archivpfad, ungültige Interface-Traces, Dialoge, Undo/Redo und
der gesamte Projektarchivtest bestanden. Der OpenGL-Test prüft sichtbares
Interface-Drahtnetz, Umschalten, unveränderte physische Flächenauswahl und
Aufräumen unabhängiger Kontexte. Bestehende Einraum-CUDA-, FMM-, LEM- und
Feldtests bestanden ebenfalls.

Diese Ergebnisse belegen die genannten Prüffälle, keine allgemeine
Konvergenz- oder Leistungsfreigabe und keine numerische AKABAK-Parität.
Horn- und BEM-Prüffälle bleiben auf insgesamt höchstens 10.000
Dreieckelemente je Fall begrenzt.
