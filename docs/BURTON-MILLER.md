# Burton–Miller: mathematischer Vertrag

Dieser Vertrag beschreibt den produktiven P1-/DP0-Galerkin-Pfad. Die
Raumtests und die analytischen CUDA-Tests wurden mit dieser Formulierung
ausgeführt; der frühere Kollokationspfad wurde aus dem Solver entfernt.

## Konventionen und Vorzeichen

Die geschlossene Randfläche Γ hat die Normale n aus dem Festkörper in das
unbeschränkte Luftgebiet. Wir verwenden RMS-Phasoren mit `exp(+iωt)`, positive
Frequenz, `k=ω/c>0`, `G(x,y)=exp(-ik|x-y|)/(4π|x-y|)` und
`q=∂ₙp=-iρωvₙ`. Die Euler-Gleichung ist `iρωv=-∇p`.
Die Potentiale sind `S[q]=∫Gq` und `D[p]=∫∂ₙᵧG p`.
Mit dieser Normalenrichtung gilt im Außenraum `p=D[p]-S[q]`.

Wir bezeichnen mit D auch den gemittelten Randoperator des Doppelschichtpotentials,
mit K′ den gemittelten Randoperator von `∂ₙₓS` und mit `H=-∂ₙₓD` den
hypersingulären Operator. Die Definition des Minuszeichens von H ist ausdrücklich
Teil des Vertrags. Die Standarddefinitionen der vier Operatoren stehen im
[Bempp-Handbuch](https://bempp.com/handbook/api/boundary_operators.html).

Die äußeren Spuren ergeben

```text
γ⁺D = ½I+D                 γ⁺S = S
∂ₙ⁺D = -H                 ∂ₙ⁺S = -½I+K′
```

Einsetzen in dieselbe Außenraumdarstellung liefert zwei Gleichungen:

```text
(½I-D)p = -Sq
Hp = -(½I+K′)q
```

Ihre Kombination mit `α=-i/k` liefert den verbindlichen diskreten Vertrag:

```text
A p = B q
A = ½ M_pp - D + α H
B = -S - α(½ M_pq + K′)
```

Die Identität ist in einer schwachen Galerkin-Form eine Massenmatrix, keine
Einheitsmatrix. Einzeln aus anderen Programmen übernommene Vorzeichen sind
unzulässig: Zeitkonvention, Normalenrichtung, Sprungrelationen und H-Definition
müssen gemeinsam stimmen. Der Ursprung der kombinierten Formulierung ist
[Burton und Miller, 1971](https://doi.org/10.1098/rspa.1971.0097).
Eine Galerkin-P1-/DP0-Darstellung mit anderer Normalen-/Zeitkonvention wird in
[Brunner et al., Abschnitt 3](https://perso.ensta.fr/~mbonnet/brunner_et_al_10.pdf)
beschrieben; deren Matrixvorzeichen werden hier deshalb nicht blind übernommen.

Für `d=x-y`, `r=|d|` sind die implementierungsnahen Ableitungen
`∂ₙᵧG=(1+ikr)G(d·nᵧ)/r²` und
`∂ₙₓG=-(1+ikr)G(d·nₓ)/r²`.
Bei einem identischen ebenen Test-/Quellpanel sind beide Normalableitungen exakt
null; das H-Integral ist trotzdem nicht null.

## Asymptotisches Fernfeld

`BemRequest.farField` wählt die analytische Grenze des Außenraumpotentials für
`x = o + r d`, `r → ∞`. `o` ist der explizite Phasenursprung. Beobachter bleiben
Weltpunkte und bestimmen ausschließlich `d = (observer-o)/|observer-o|`.
Ein Punkt gleich dem Ursprung ist ungültig. Seine Entfernung und eine etwaige
Lage innerhalb des Körpers spielen für diese Richtungsangabe keine Rolle.

Mit den obigen Potentialen ergibt die Entwicklung
`|o+r d-y| = r - d·(y-o) + O(1/r)` den ausgewerteten Vertrag:

```text
F_o(d) = 1/(4π) ∫Γ exp(+ik d·(y-o)) [ik(d·n_y) p(y) - q(y)] dS_y
p_far(d; R) = exp(-ikR)/R · F_o(d)
```

`R = referenceDistance > 0` ist standardmäßig 1 m. Das Ergebnis enthält also
die Ausbreitungsphase `exp(-ikR)` und den Pegelfaktor `1/R`. Es ist auf diesen
Abstand normierter **Fernfelddruck**, keine Berechnung des endlichen Feldes
an einem tatsächlich 1 m entfernten Punkt. Insbesondere wird die Fernfeldgrenze
nicht durch willkürlich große Beobachterkoordinaten angenähert. Die allgemeinen
Potentialdefinitionen und getrennten Fernfeldoperatoren sind auch im
[Bempp-Handbuch](https://bempp.com/handbook/api/potential_operators.html) und der
[Fernfeld-API](https://bempp-cl.readthedocs.io/en/latest/docs/bempp_cl/api/operators/far_field/helmholtz/index.html)
beschrieben; hier gelten ausdrücklich unsere positive Zeitkonvention und
obigen Vorzeichen.

Der Faktor `exp(-ikR)/R` wird getrennt vom Aperturintegral ausgewertet.
Dadurch wird die kleine Aperturphase nicht schon beim Subtrahieren eines sehr
großen Bezugsabstands aus den Quellkoordinaten ausgelöscht.

Der eigene CUDA-Kernel integriert diese glatte Darstellung mit der vorhandenen
Dreiecksregel achten Grades, stetigem P1-Druck und DP0-Fluss. Phase,
Geometrie, Summation und die gelösten Druckkoeffizienten verwenden FP64.
Netz und Phasenursprung werden gemeinsam lokal
zentriert, bevor ihre Differenzen gebildet werden. Eine gemeinsame Translation
ändert dadurch den Vertrag nicht. Druck auf der Oberfläche, Strahlungsimpedanz
und Burton–Miller-Gleichung sind von diesem Beobachtungsmodus unabhängig.
Ohne `farField` bleibt die bisherige endliche Beobachtung einschließlich
adaptiver Nahfeldintegration aktiv. Der Beobachtercache unterscheidet beide
Modi sowie Phasenursprung und Bezugsabstand.

`bem_cuda_tests --asymptotic` prüft diesen Pfad auf der echten RTX 5070 Laptop
GPU mit analytischer pulsierender Kugel und translatorischer Dipolanregung
`v_n=n_z`. Bei 100, 1000 und 1715 Hz sank der maximale komplexe Kugelfehler
von 2,493 % bei 512 Dreiecken auf 0,625 % bei 2048 Dreiecken. Der Dipol hatte
0,686 % Fehler; Gegenphase auf ±Z und die Quer-Nullstelle werden separat geprüft.
Endliche Kugelfelder, auf 1 m und dieselbe Phase zurückgerechnet, nähern sich
dem Fernfeld bei 1/10/100/10000 m mit Fehlern von rund
`2,90·10⁻⁴ / 2,90·10⁻⁵ / 2,90·10⁻⁶ / 5,4·10⁻⁸`.
Zusätzlich geprüft sind Bezugsabstand und `1/R`-Pegel, Phasenursprung,
Richtungsnormalisierung, gemeinsame Translation um `10⁹ m` auf identisch
gerundeter Geometrie, Cache-Moduswechsel, komplexe Mehrfachanregung sowie
unveränderte Oberflächendrücke und Impedanz. Diese Referenzfälle ersetzen
keinen Netzkonvergenznachweis für beliebige Hörner.

## Schwache H-Form und Räume

Druck wird mit stetigen linearen Knotenfunktionen φ dargestellt; der
Normalfluss mit unstetigen Konstanten pro Dreieck. Testfunktionen sind P1.
Für geschlossene Flächen gilt die regularisierte schwache Form

```text
H_ij = ∫Γ∫Γ G(x,y) [curlΓφ_i(x)·curlΓφ_j(y)
                    - k²(nₓ·nᵧ) φ_i(x)φ_j(y)] dS_y dS_x
curlΓφ = n × ∇Γφ
```

Damit werden schwach singuläre Kerne integriert, keine punktweise ausgewerteten
hypersingulären Selbstterme. Die globale Stetigkeit von P1 ist wesentlich für
die Regularisierung; unabhängige konstante Druckwerte je Panel würden einen
anderen Raum mit zusätzlichen Sprungtermen ergeben.
[Schwach regularisierte Helmholtz-Galerkin-Form, Gleichung 18](https://pmc.ncbi.nlm.nih.gov/articles/PMC11614954/)

Auf einem ebenen Dreieck mit Knoten a,b,c, Fläche A und Normale n sind
`∇λ₀=n×(c-b)/(2A)`, `∇λ₁=n×(a-c)/(2A)` und
`∇λ₂=n×(b-a)/(2A)`. Gradienten und Curls sind panelweise konstant.
Die exakten lokalen Massen sind `M_pp=A/12*[2 1 1;1 2 1;1 1 2]`
und `M_pq=A/3*[1;1;1]`. Druckfreiheitsgrad i bleibt Netzknoten i;
Flussfreiheitsgrad t bleibt Dreieck t. Physische Flächentags werden erhalten.

`galerkin/space.cpp` berechnet diese frequenzunabhängigen Größen einmal sowie
CSR-Knoteninzidenzen, konsistente Massendiagonale und Zeilensummen.
Die Tests prüfen eine unabhängige affine Referenz, Partition der Eins,
Massensummen, CSR-Rückverweise, Translation, beliebige eigentliche Rotation,
Skalierung, umgekehrte Wicklung und fehlerhafte Eingaben. Geschlossene,
orientierte und schnittfreie Geometrie wird zusätzlich im gemeinsamen
Netzvertrag geprüft.

## Integration und GPU-Abnahme

Reguläre, nahezu berührende, eckennachbarliche, kantennachbarliche und identische
Panelpaare benötigen unterscheidbare Integrationsregeln. Die Referenzdreiecke
haben Fläche ½; eine gemeinsame Referenz-Paarregel integriert die Konstante zu
¼. Der physische Multiplikator ist daher `4*A_test*A_trial`. Permutierte
Nachbarschaftsregeln müssen ihre baryzentrischen Koordinaten wieder auf die
ursprüngliche lokale Knotenreihenfolge abbilden.

Der Produktionspfad baut A und mehrere rechte Seiten unmittelbar auf der GPU
auf. Je Frequenz verwendet er eine pivotierte komplexe LU-Faktorisierung für
alle Anregungen und gegebenenfalls die iterative Nachkorrektur. Geometrie und
Integrationsregeln werden über den Sweep wiederverwendet. Die Matrix wird
mit den reziproken P1-Massenzeilensummen äquilibriert; das ausgewiesene relative
L2-Residuum bezieht sich auf die ursprüngliche, unskalierte Gleichung.
LU-Faktoren und Korrekturrechte-Seiten bleiben FP32. Lösung, erneute Auswertung
von `b-Ap`, Rückskalierung mit der P1-Masse und Residualnorm verwenden FP64.
So begrenzt nicht die FP32-Summation der Residualberechnung die Nachkorrektur.
Das Residuum prüft die tatsächlich assemblierte FP32-Matrix; Quadratur- und
Diskretisierungsfehler werden dadurch nicht zu Rundungsfehlern erklärt.
Singuläre Transformationen und relative Panelgeometrie werden in FP64
berechnet; ihre Kerne und Momente verwenden FP32. Bei adaptiven nahen Integralen
bleiben Plan und baryzentrische Transformationen in FP64. Für verlangte
Nahquadraturtoleranzen ab 10⁻⁵ werden Kerne und Momente in FP32 ausgewertet;
darunter wird derselbe Plan mit FP64-Kernen und -Momenten integriert.
Der Vergleich beider tatsächlich benutzten Pfade auf der Kugel ergab für
komplexen Druck und Impedanz eine maximale relative Abweichung von 2,70·10⁻⁷
bei einem unveränderten Abnahmelimit von 10⁻⁴.
Reguläre Paare verwenden eine FP64-Differenz der Panelursprünge
mit vorberechneten lokalen FP32-Offsets und FP32-Kernen und -Momenten.
Gesammelt wird in komplexen FP32-Matrizen. Der interne lokale Ursprung wird
einmal pro Solversitzung abgezogen; Beobachter erhalten denselben Versatz.
Das reguläre Feld verwendet dieselben vorbereiteten 3-/6-Punkt-FP32-Regeln
wie der gewählte reguläre Grad 2/4. Beobachterabstand, `k*r` und die
trigonometrische Phase bleiben dabei in FP64; große Abstände dürfen keine
zusätzlichen Phasenfehler durch ein gerundetes Float-Wellenzahlprodukt erzeugen.
Eine eigene CUDA-Regression prüft die radiale Propagation zwischen 1000 m und
10⁶ m bei 1715 Hz auf derselben Kugel. Der vorherige Float-Phasenpfad scheiterte
mit einem komplexen Verhältnisfehler von 2,62; der korrigierte Lauf liegt unter
3·10⁻⁷ bei unverändertem Grenzwert 10⁻³ für das facettierte Modell.
Bei einem Abstand zum Panelmittelpunkt von
höchstens zwei maximalen Panelkanten bleiben eine adaptive 16-Punkt-FP64-Regel
und ihre Fehlerkontrolle aktiv.
Ein Beobachter 1 mm außerhalb der Kugel ist in den CUDA-Tests enthalten.
Das Feld verwendet die linearen Druckfunktionen; die öffentliche
`surfacePressure`-Ausgabe ist deren exakter Mittelwert pro Dreieck.

Am 19.10.2026 bestanden die eigenen CUDA-Tests auf einer RTX 5070 Laptop GPU
mit regulärer Ordnung 4, singulärer Ordnung 4 und adaptiver Nahintegration
(`nearDistanceRatio=1`). Die Kugelnetze mit 128, 512 und 2048 Dreiecken zeigen
abnehmenden komplexen Druckfehler. Bei 2048 Dreiecken beträgt der größte
analytische Fehler über 100, 500, 1000 und 1715 Hz und drei Beobachter 0,673 %;
der größte ursprüngliche relative Gleichungsrest beträgt 2,05·10⁻⁵. 1715 Hz
entspricht bei Radius 0,1 m und c=343 m/s `ka=π`.
Der vorhandene unabhängige BEAT-CPU-Datensatz auf demselben 512er Netz stimmt
im komplexen Beobachterdruck mit relativen L2-Fehlern unter 10⁻⁶ überein.
Komplexe Mehrfachanregung, Impedanznormalisierung, Abbruch sowie Ablehnung von
inneren und auf dem Rand liegenden Beobachtern wurden ebenfalls ausgeführt.
Ein zusätzlicher 20-/50-Hz-Lauf auf dem 2048er Netz bestand mit unveränderter
Resttoleranz 5·10⁻⁵. Bei 20 Hz wurde eine iterative Nachkorrektur benötigt;
die ursprünglichen Reste lagen unter 1,7·10⁻⁵ beziehungsweise 3,9·10⁻⁵.
Translationsprüfungen vergleichen dieselbe bereits gerundete Geometrie lokal
und bei gemeinsamen Versätzen von 10⁹ beziehungsweise 10¹² m. Sie prüfen
Druck und komplexe Impedanz, ohne Genauigkeit verlorener Eingabebits zu behaupten.
Diese Läufe belegen noch keine Laufzeitparität auf dem Produktionshorn.

Der CPU-Interaktionsplan verwendet räumliche Kandidaten, begrenzte parallele
Auswertung und exakt sortierte Nahpaar-Zeilen. Ein ungeordnetes Nahpaar wird
einmal geplant; die Rückrichtung tauscht die baryzentrischen Patchkoordinaten.
Nur tatsächlich nahe disjunkte Paare werden vom regulären GPU-Pfad ausgenommen.
Falsche Kandidaten erhalten weiterhin die reguläre Regel. Solversitzungen halten
Quadratur, validierte Beobachter und dimensionierte Matrix-/LU-Arbeitspuffer
über wiederholte Aufrufe vor; Matrixwerte werden je Frequenz neu berechnet.
Die nahe Tensorregel bereitet je Patch ihre beiden Punktlisten samt
P1-Koordinaten einmal in FP64 im gemeinsamen GPU-Speicher vor. Reguläre
Nahpatches verwenden je 16 Punkte; bei kritischen Patches am räumlichen
Unterteilungslimit sind es je 64 Punkte. CPU-Fehlerschätzer und GPU-Auswertung
verwenden denselben gespeicherten Grad je Patch.

Für einen Laufzeitvergleich sind identische Nahintegrationsoptionen nötig.
`nearDistanceRatio=0` schaltet ausschließlich die adaptive Integration
disjunkter Panelpaare ab; berührende und identische Paare behalten ihre
singulären Regeln. Die genaue Nahintegration wird separat bewertet.

Die Abnahme benötigt echte CUDA-Rechnungen: komplexe Kugellösung einschließlich
Phase und Impedanz, Netz-/Quadraturkonvergenz, Eigenfrequenzumgebung, unabhängige
Mehrfachanregungen und Parität mit derselben externen Boundary-Lab-/BEAT-Version
auf identischen Netzen. Weder der Methodenname noch ein kleiner linearer Rest
ersetzt diese Nachweise. Die unabhängige Implementierung enthält keinen
übernommenen BEAT- oder Bempp-Produktivcode.

## Unabhängige Prüfung der Paarquadratur

Die CPU-Paarquadratur wird unabhängig vom linearen Solver geprüft. Identische,
kantenbenachbarte und eckenbenachbarte Dreiecke verwenden jeweils sechs, fünf
und zwei Sektoren einer Sauter–Schwab-/Duffy-Transformation. Die Abbildungen
wurden aus den Integrationsgebieten und den veröffentlichten Transformationen
in [Abschnitt 4 der Sauter–Schwab-Darstellung](https://link.springer.com/article/10.1007/s10915-023-02318-y)
hergeleitet. Die Tests prüfen die Jacobi-Determinanten numerisch, alle neun
P1-Paarmomente, die Rückpermutation bei entgegengesetzter Kantenorientierung,
Reziprozität, Skalierung und Konvergenz. Für das rechtwinklige Einheitsdreieck
ist außerdem das unabhängige analytische Integral
`∫T∫T 1/|x-y| = (2+√2)*asinh(1)/3` enthalten. Sein absoluter Fehler fällt
von `6,96·10⁻³` bei Gaußordnung 2 auf `6,57·10⁻¹⁴` bei Ordnung 16.

Disjunkte Nahpaare verwenden positive Dreiecksregeln der Grade 6 und 8 zur
Fehlerschätzung. Ein gemeinsames Fehlerbudget für das gesamte Panelpaar steuert
die Unterteilung des jeweils größten Fehlerbeitrags. Ein geometrischer Abstandstest
darf eine bereits erreichte Fehlertoleranz nicht überstimmen. Der Plan enthält
affine Teilpanelpaare und verwendet normalerweise Grad 8 für die
Auswertung. Diese Regeln beziehen sich auf
[Dunavant, 1985](https://doi.org/10.1002/nme.1620210612).
Am räumlichen Unterteilungslimit erhöht ein noch kritisches Teilpanelpaar
gezielt die Integrationsordnung: positive Tensor-Gauß-Regeln der Grade 12/14
mit 49/64 Punkten liefern Schätzer und Auswertung. Das ist eine zusätzliche
Genauigkeitsanpassung innerhalb derselben Teilflächen. Fehlertoleranz und
räumliche Tiefe bleiben gleich. Nicht erreichte Toleranz oder ein zu großer
Plan führen weiterhin zu einem Fehler; die Grenze von 50.000 Teilpanelpaaren
und die Abbruchprüfung bleiben erhalten.
Ein reproduzierter Fall mit nahe beieinanderliegenden Ecken erreicht gegenüber
unabhängigen Tensor-Gauß-Regeln der Ordnungen 16/20 einen relativen Fehler von
1,72·10⁻⁸ bei angeforderten 10⁻⁴. Toleranz und Unterteilungslimit wurden dafür
nicht erhöht. Nicht aufgelöste sehr enge Parallelflächen werden weiterhin
abgelehnt; die Fehlermeldung nennt Panels, Flächen, Abstand und Wellenzahl.

Ein weiterer zuvor abgebrochener Parallelflächenfall mit Kante 1 m und
Abstand 0,1 m besteht bei unverändert 10⁻⁴/Tiefe 6 mit 1.684 Patches,
davon 194 höherer Ordnung. Unabhängige gleichmäßige Unterteilung mit
Tensor-Gauß-Ordnungen 12/16 bestätigt S, D, K′ und H bis höchstens
1,41·10⁻⁷ relativ. Die echte CUDA-Regression `bem_cuda_thin_body` berechnet
eine geschlossene Box 100 × 100 × 10 mm bei 20/1000 Hz mit zwei komplexen
Anregungen. Der bisherige Nahpaarabbruch entfällt; Residuen sind höchstens
2,38·10⁻⁶. Gegen einen strengeren Nahplan (10⁻⁶/Tiefe 8) unterscheiden sich
die komplexen Drücke höchstens um 8,14·10⁻⁸. Impedanz, Phase und Symmetrie
werden ebenfalls geprüft. Dieser Vergleich auf demselben groben Netz ist
ein Quadratur-/Lösungstest, kein physikalischer Netzkonvergenznachweis.

Am 19.10.2026 wurden zusätzlich die einzelnen komplexen S-, D-, K′- und
H-Blöcke bei `k=0,7 m⁻¹` mit einer isolierten externen CPU-Referenz verglichen:
BEAT Engine 0.4.0, Commit
`96cb77e5eb5b8c15c6bbf50b73c4737028ab6cea`, Boundary Lab
`3a721c31c1d03211beb14fb729116c0b3505c645`, Julia 1.12.7.
Beide Rechnungen verwendeten reguläre Ordnung 4, singuläre Ordnung 8 und
`nearDistanceRatio=0`, damit disjunkte Paare dieselbe feste Integrationsregel
hatten. Die vier Fälle sind ein identisches rechtwinkliges Dreieck, zwei
Dreiecke mit entgegengesetzt orientierter gemeinsamer Kante, zwei mit genau
einem gemeinsamen Knoten und zwei räumlich getrennte Dreiecke.

Als Vergleichstoleranz gilt `10⁻⁶` relativer Frobeniusfehler je
nicht verschwindendem Operatorblock. Der größte gemessene Fehler beträgt
`2,31·10⁻⁷` für K′ des Kantenpaars; die anderen Topologien liegen unter
`1,02·10⁻¹⁴`. D und K′ des identischen ebenen Panels sind in beiden Rechnungen
exakt null. Dieser Vergleich prüft Quadratur, Ableitungszeichen und schwache
H-Form auf der CPU; er ist keine zusätzliche GPU-Abnahme.

Die externen Rohdaten liegen außerhalb des Produkts unter
`work/reference-env/` im ursprünglichen Arbeitsverzeichnis. SHA-256:

| Datei | SHA-256 |
|---|---|
| `pair-operators-beat.json` | `119fc659927c4cd4ab7df0511f5f1a53a2e727ba3f670b4e6315a201b40a2cc5` |
| `pair-operators-hornsculpt.json` | `9ab0b00f9c731f4590a0de6b686fe558df69d5b9c4437da2f684cc180b58cd8a` |
| `pair-operator-comparison.json` | `6b10f6df2cd0e908f90591483ff23ac06976bca83abc196929ea61027679ffe8` |

Eine weitere CPU-Probe bildet einen hornähnlichen Abstand ab: zwei versetzte,
entgegengesetzt orientierte gleichseitige Dreiecke mit 8 mm Kantenlänge und
6 mm Abstand. Mit der Produkttoleranz `10⁻⁴` und maximaler Tiefe 6 erzeugt
der bei 16 kHz vorbereitete Plan 16 Teilpanelpaare bei Tiefe 2. Für 100, 500,
1000, 2000, 4000, 8000 und 16000 Hz beträgt der größte relative Blockfehler
von S, D, K′ oder H `2,27·10⁻⁹` gegenüber unabhängiger Tensor-Gauß-Quadratur
der Ordnung 20. Ordnung 16 und 20 stimmen dabei bis `1,81·10⁻¹⁴` überein.
Die separate Integration eines nahen Beobachterfelds mit komplexem linearem
Druck und konstantem Fluss weicht höchstens um `2,57·10⁻⁵` ab.
Der Datensatz `near-horn-8mm-6mm.json` hat SHA-256
`af7aafbcd5427523178e10d3facded457a12c0c12be177c47f963a4b1a59f664`.
Diese lokale Geometrieprobe ersetzt weder die Fehlerprüfung des gesamten
Produktionshorns noch einen Laufzeitvergleich.

## Unterstützter Außenraum

Nach der Prüfung auf geschlossene, orientierbare und schnittfreie Netze prüft
`validateExteriorComponents` jede zusammenhängende Oberfläche separat. Ihr
vorzeichenbehaftetes Volumen wird mit eigenem lokalem Ursprung und kompensierter
Summation berechnet und muss positiv sein. Ein positives Gesamtvolumen allein
kann einen kleineren, falsch gewickelten Körper verdecken und reicht deshalb
nicht aus. Die Prüfung verändert weder Flächen noch Wicklungen.

Mehrere getrennte, nach außen orientierte Festkörper im selben unbeschränkten
Luftgebiet sind zulässig. Verschachtelte Oberflächen werden über ihren Raumwinkel
erkannt und explizit abgelehnt. Eingeschlossene Luftkavitäten und die Kopplung
getrennter Innen- und Außenluftgebiete sind damit nicht unterstützt. Ein offenes
Horninnere mit Verbindung durch den Mund zum Außenraum ist davon nicht betroffen.
Die CPU-Tests decken getrennte Körper, invertierte Teilkörper trotz positivem
Gesamtvolumen, verschachtelte gleichgerichtete Oberflächen, überlappende
Begrenzungsboxen ohne Verschachtelung und CAD-Platzierungen von `10¹² m` ab.
Paarabstände werden vor der Interpolation aus Ursprungsdifferenz und lokalen
Kanten gebildet. Die zusätzliche Quadraturregression vergleicht identische,
bereits im Speicher gerundete Koordinaten vor und nach lokaler Zentrierung;
S-, D-, K′- und H-Momente für identische, kantenbenachbarte und nahe disjunkte
Paare stimmen bei gemeinsamen Platzierungen von `10⁹ m` und `10¹² m` bis
`2·10⁻¹²` relativ überein.

Diese CPU-Verträge belegen keine uneingeschränkte Platzierungsreichweite des
gesamten Solvers: Die gemeinsame Netzprüfung bestimmt ihre Degenerationsschwelle
noch aus der globalen Begrenzungsbox. Extrem weit voneinander entfernte kleine
Körper können deshalb bereits vor der Außenraumprüfung als numerisch degeneriert
abgelehnt werden. Die Tests der lokalen Komponentenvolumina und Paarquadratur
umgehen diese Einschränkung nicht im Produktionspfad.
