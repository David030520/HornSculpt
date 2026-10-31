# BEM und gekoppelte Modelle

## Randwertproblem und Diskretisierung

Für ein homogenes lineares Medium löst die frequenzabhängige BEM das
Helmholtz-Problem über Randgrößen. Das spart ein großes Volumennetz im
Außenraum, führt aber zu dichten Fernwechselwirkungen. Eine offene
Abstrahlungsrechnung braucht eine passende auslaufende Lösung; ein
Innenraumproblem und eine abstrahlende Außenfläche sind nicht austauschbar.
COMSOL dokumentiert diese Abwägung ausdrücklich und verwendet dort Costabels
symmetrische Kopplung mit ACA und GMRES, nicht HornSculpts Burton–Miller-Pfad
[S15]. Gleiche Diagramme bedeuten deshalb nicht gleiche Diskretisierung.

Bei Galerkin werden Test- und Ansatzfunktionen über Flächenpaare integriert.
Der Burton–Miller-Ansatz kombiniert Randgleichung und Normalableitung, um
nichtphysikalische Mehrdeutigkeiten an inneren Eigenfrequenzen zu vermeiden
[S13, S14]. Die Kombination beseitigt nicht automatisch Quadratur-,
Diskretisierungs- oder Modellfehler. Vorzeichen von Normalen, Green-Funktion,
Zeitkonvention und Kopplungsparameter gehören zum mathematischen Vertrag;
Werte aus unterschiedlichen Konventionen erst nach Umrechnung vergleichen.

Oden und Mitautoren leiten für ihr gekoppeltes Streuproblem eine schwache
Galerkin-Form her: die hypersinguläre Normalableitung wird in Terme mit
Oberflächenableitungen und schwächer singulärem Kernel überführt [S36,
S. 29–32]. Die dortigen Vorzeichen gehören zu ihrer Green-Funktion und
Normalenrichtung; sie dürfen nicht ungeprüft als Rezept für eine andere
Implementierung eingesetzt werden. Praktische Operatorprüfungen müssen die
einzelnen Blöcke sowie die zusammengesetzte Lösung einschließen.

## Singular- und Nahintegration

Gleiche, kanten- und eckenbenachbarte Elemente verlangen eine geeignete
Behandlung der singulären Operatoren; nahe, aber getrennte Flächen können
ebenfalls sehr steile Integranden erzeugen. Eine reguläre Quadraturordnung,
die für entfernte Flächen genügt, ist kein Nachweis für enge Spalte. Die
hypersinguläre Form muss mathematisch konsistent regularisiert werden [S13,
S14]. Prüfen: identische, benachbarte und getrennte Paare gegen unabhängige
Referenz; Abstand und Quadraturordnung systematisch variieren. Es gibt keine
allgemein gültige einzelne Nahschwelle.

## Beschleunigung und Symmetrie

FMM und ACA approximieren geeignete Fernwechselwirkungen, behalten aber
Nahwechselwirkungen gesondert. Näherungstoleranz, Iterationen, Speicher,
Setup und Gesamtlaufzeit sind gemeinsam gegen einen ausreichend genauen
direkten Rechenweg zu messen. COMSOLs ACA-Strategie ist eine andere Methode
als FMM [S15]; daraus folgt kein gemessener Geschwindigkeitsvorteil eines
bestimmten HornSculpt-Modus.

Ein FMM-Fehler und der Fehler der Quellenapproximation sind unabhängig
veränderbare Beiträge. Langston und Mitautoren zeigen dies für einen
adaptiven Volumenintegral-Löser und behandeln Nahwechselwirkungen mit
separater Quadratur [S37, S. 1–3]. Das Modell dort ist keine Validierung einer
hochfrequenten akustischen BEM und nennt keine HornSculpt-Laufzeit. Für
HornSculpt daher die FMM-Toleranz bei festem Netz gegen dichte Ergebnisse
testen, anschließend Netz und Nahquadratur unabhängig verfeinern.

Eine Spiegelebene ist nur zulässig, wenn Geometrie und sämtliche Quellen die
entsprechende Symmetrie oder Antisymmetrie erfüllen. Das Bildquellenmodell
hat außerdem geometrische Gültigkeitsbedingungen; COMSOL weist darauf hin,
dass eine gewöhnliche unendliche Symmetrieebene nicht beliebig eine
Durchgangsöffnung zwischen zwei Gebieten ersetzt [S15]. Für eine physische
Schallwand mit Apertur muss das tatsächliche gekoppelte Problem formuliert
werden.

## LEM und Ports

Akustische Impedanz ist Druck durch Volumenstrom in Pa·s/m³; die
flächenbezogene Größe Druck durch Geschwindigkeit hat eine andere Einheit.
Mit Bezugsfläche S wird die Volumenimpedanz durch ρc/S normiert.
Die Kopplung muss Druck und Volumenstrom mit konsistenter Orientierung und
Leistungsbilanz austauschen. Ein virtueller Interfaceport ist weder
Materialwand noch zusätzliche unabhängige Quelle. Gegenseitige Lasten
mehrerer Ports müssen in der gekoppelten Lösung enthalten sein. COMSOL zeigt
den Unterschied zwischen vorgeschriebener Bewegung und rückgekoppeltem
elektromechanischem Modell [S21, S22]; seine FEM-Beispiele sind fachliche
Orientierung, kein Nachweis für HornSculpts konkrete BEM-Implementierung.

Eine Schallfeldkarte stellt endliche Feldpunkte dar; das asymptotische
Fernfeld ist davon zu unterscheiden. Drucknullstellen besitzen keine
verlässliche Phase. Randnahe Punkte brauchen kontrollierte Nahintegration;
maskierte oder nicht konvergierte Punkte dürfen nicht als gültig erscheinen.

## Was eine Konvergenzprüfung tatsächlich verändert

Ein Galerkin-Eintrag koppelt zwei Flächenintegrale. Bei zusammenfallenden
Flächen liegt die Singularität innerhalb des Integrationsgebiets; bei Kanten-
oder Eckkontakt liegt sie an dessen Rand. Bei disjunkten Flächen ist der Kernel
zwar endlich, kann aber für einen kleinen Spalt sehr stark variieren. Sauter
und Schwab unterscheiden geometrische Fallklassen und betrachten den Einfluss
des Matrixintegrationsfehlers auf die vollständig diskrete Konvergenz [S29,
S. 12–15]. Ihre analytischen Voraussetzungen sind eng; daraus folgt keine
Konvergenzgarantie für beliebige Helmholtz-Spalte.

COMSOL 6.4 beschreibt für beinahe berührende 3D-Elemente eine zusätzliche
Koordinatentransformation und eine einstellbare höhere Quadraturordnung,
aktiviert über einen relativen Abstand [S15, „Quadrature“]. Diese
Herstelleroption belegt die Problemklasse, nicht die Gleichwertigkeit einer
HornSculpt-Einstellung. Für einen reproduzierbaren Vergleich Elementgröße,
absoluten Spalt, Verhältnis Spalt/Elementgröße, verwendete Quadratur und
Referenzfehler angeben. Erst bei stabilen Quadraturergebnissen die Netzfolge
als Diskretisierungskonvergenz auswerten.

## Feld, Rand und Stabilisierung auseinanderhalten

COMSOL beschreibt eine Interpolation nahe dem Rand, weil die direkte
Kerndarstellung im Volumen dort numerisch empfindlich sein kann [S15,
„Postprocessing Settings“]. Das ist eine Auswertungsentscheidung;
Randlösung und Volumengrenzwert sind getrennt zu testen. Bei HornSculpt sind
Gültigkeitsmaske, Randabstand und Ergebnisrevision sichtbar zu halten; eine
farbig interpolierte Fläche macht nicht-konvergierte Punkte nicht gültig.

COMSOLs Stabilisierung setzt einen weiteren Parameter und zusätzliche
Freiheitsgrade ein; die Anleitung unterscheidet große von kleinen Modellen und
nennt Einschränkungen für innere Randbedingungen [S15, „Stabilization“].
Diese Stabilisierung ist nicht gleichbedeutend mit Burton–Miller, FMM oder
einer beliebigen iterativen Vorkonditionierung. Vergleichbar sind erst
erreichte Feldfehler und Gesamtkosten unter derselben Aufgabe.

## Modalität und Schnittstellen

Kolbreks modale BEM/MMM-Kopplung bildet für jede Geschwindigkeitsmode eine
Druckantwort und projiziert sie auf alle berücksichtigten Druckmoden.
Dadurch entsteht eine Matrix mit gegenseitigen Lastanteilen [S07, S. 213–216].
Die Arbeit zeigt in einem Hybridbeispiel, dass grobe räumliche Abtastung an der
Schnittfläche höhere Moden verfälschen kann, obwohl die dominierende
Halsimpedanz vergleichsweise gut stimmt [S07, S. 218–220]. Bei HornSculpts
skalaren Ports ist das kein Nachweis einer allgemeinen modalen
BEM↔BEM-Kopplung. Es ist ein Grund, die angenommene Geschwindigkeitsform und
die resultierende Druckprojektion gezielt zu prüfen.


## COMSOL 6.4 als Methoden- und Modellierungsvergleich

Die öffentliche Dokumentation [S15] bezeichnet ihr Verfahren als direkte
BEM mit symmetrischer Costabel-Kopplung. Sie nennt adaptive cross approximation
(ACA), teilweise Matrixassemblierung und GMRES; eine stabilisierte Variante
bringt zusätzliche Freiheitsgrade mit. HornSculpts produktiver Ansatz ist
Burton–Miller-Galerkin. Diese Namen bezeichnen unterschiedliche
Formulierungs- und Beschleunigungsentscheidungen. Insbesondere ist ACA nicht
automatisch FMM, und COMSOL-Stabilisierung kein Beleg für die dortige
Verwendung von Burton–Miller. Ein Solververgleich muss deshalb dasselbe
physikalische Randwertproblem, dieselben Normierungen und nachgewiesene
Fehlerziele verwenden; gleiche Frequenz und gleiche Knotenzahl reichen nicht.

Für die Randnähe beschreibt COMSOL die Interpolation zwischen der gelösten
Randgröße und der Feld-Kernauswertung, standardmäßig über 0,5 lokale
Elementgrößen [S15, „Postprocessing Settings“]. Das ist nicht die gleiche
Größe wie eine explizite Ausschlussdistanz oder eine erhöhte Nahquadratur.
Bei Vergleichen komplexen Druck und Phase an denselben Weltpunkten prüfen,
die Wandabstände und Interpolationsoptionen dokumentieren und nahe Punkte
zusätzlich mit verfeinerter Quadratur untersuchen. Auf dem Rand liegende,
ausgeschlossene und nichtkonvergierte Punkte nicht als Luftwerte mitteln.
Die dünnspaltspezifische 3D-Transformation und der relative Grenzabstand sind
in [S15, „Quadrature“] eigenständige Maßnahmen bei beinahe benachbarten
Elementpaaren. Für konkrete Spaltfälle absoluten Abstand, lokale Kantenlänge,
Quadratur und Netzfolge separat variieren.

Für symmetrische bzw. antisymmetrische unendliche Ebenen fordert COMSOL,
dass Quellen und BEM-Ränder auf derselben Seite liegen; eine Öffnung mit
Quellen auf der gegenüberliegenden Seite ist damit nicht pauschal eine solche
Symmetrieebene [S15]. Vor Anwendung sind Geometrie **und Anregung** zu
spiegeln und anschließend Feld und Leistung gegen ein Vollmodell zu prüfen.

COMSOL koppelt BEM und FEM bidirektional über den Gesamtdruck und verwendet
FEM unter anderem für innere oder poröse Bereiche [S39]. Die
Akustik-Struktur-Kopplung überträgt strukturelle Normalbeschleunigung auf das
Fluid und dessen Drucklast zurück auf die Struktur [S40]. Bei einer innen
liegenden dünnen Schale mit Fluid beidseitig nennt die Dokumentation für BEM
eine Einschränkung und empfiehlt eine umgebende FEM-Domäne mit BEM–FEM-Rand.
Das Lüftungsgehäuse-Tutorial [S41] veranschaulicht den kombinierten Aufbau,
nicht aber eine bestätigte Kopplung in HornSculpt. HornSculpts vorhandene
LEM-Anschlüsse sind ein eigener reduzierter Vertrag; eine elastische Membran,
thermoviskose Grenzschicht oder allgemeine FEM–BEM-Kopplung darf daraus nicht
als bereits umgesetzt abgeleitet werden.

Praktischer Prüfweg: zunächst schallharte Geometrie und eindeutig orientierte
vorgegebene Geschwindigkeit vergleichen; dann bei gleichem Netzkonvergenzziel
komplexe Rand- und Feldwerte, akustische Leistung und Impedanz betrachten.
Für gekoppelte Referenzfälle zusätzlich Schnittstellen-Normalen, beidseitige
Leistungsbilanz, Frequenz, Materialparameter und verwendete Strukturdynamik
festhalten. Erst danach sind Rechenzeit und Speicher unter ausgewiesener
Hardware und Genauigkeit sinnvoll vergleichbar.
