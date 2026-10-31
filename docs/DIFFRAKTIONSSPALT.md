# Diffraktionsspalt für PA-Hörner

Stand: 29. Oktober 2026. Optionale Draft-Halsbaugruppe mit zwei unabhängig
einstellbaren Radien, gespeichert mit Projektformat 14. Projekte vor Format 13
öffnen mit deaktiviertem Spalt; Format 13 behält Superellipse und axiale Glättung.

## Empfehlung und fachlicher Zweck

Eine optionale Halsbaugruppe in Draft verbindet den tatsächlichen
Treiberanschluss über einen Formadapter mit einem schmalen Spalt und der
anschließenden OS-SE-Aufweitung. Der Spalt ist eine offene Stelle des
zusammenhängenden Luftwegs. Seine Wände, Tiefe und Übergänge werden als reale
Geometrie modelliert. Am Spalt wird keine zusätzliche Schallquelle angesetzt.

Ein vertikaler, also in Y hoher und in X schmaler Spalt ist ein sinnvoller
Ausgangspunkt für breite horizontale Abstrahlung. Das ist ein Entwurfsansatz,
keine Zusage für einen bestimmten Winkel oder Frequenzbereich. Die anschließende
Aufweitung, Mundgröße, Phasenverteilung und der Treiber bleiben entscheidend.
Die klassischen CD-Hörner und ihre unterschiedlichen Kompromisse sind in
[JBL Technical Note 1/4](https://www.cieri.net/Documenti/JBL/Technical%20Notes/JBL%20Technical%20Note%20-%20Vol.1%2C%20No.4%20-%20Constant%20Directivity%20Horns.pdf)
dokumentiert; der Link enthält das originale Herstellerdokument als Archivkopie.

Ein Spalt ist nicht grundsätzlich für jedes PA-Horn erforderlich. JBL beschreibt
bei traditionellen Spaltübergängen Nachteile für Impedanz, Frequenzgang und
Verzerrungen und stellt ihnen kontinuierliche Übergänge gegenüber. Daher gehört
zum Entwurf immer ein Vergleich mit dem unveränderten OS-SE-Horn und einer
geglätteten Spaltvariante. Das ist keine pauschale Wertung aller Spalthörner.
[JBL Technical Note 1/31](https://jblpro.com/en-US/site_elements/tech-note-progressive-transitiontm-pt-waveguides)

ATH besitzt mit `Slot.Length` bereits ein anfängliches gerades Profilsegment,
dessen Länge als Ausdruck über dem Azimut angegeben werden kann. Ein solcher
Profilversatz und ein expliziter Rund-zu-Spalt-Adapter sind unterschiedliche
Operationen. Eine einzelne Eingabe „Slot.Length“ würde Breite, Höhe, Fläche und
Übergang des hier vorgeschlagenen PA-Spalts nicht ausreichend beschreiben.
[ATH 4.8.2 User Guide, Abschnitt 4.1.1, S. 18](https://at-horns.eu/release/Ath-4.8.2-UserGuide.pdf)

## Geometrie und Bedienung

Die folgende Gestaltung ist eine eigene Parameterisierung. Die Quellen
begründen den fachlichen Kontext, nicht diese konkrete Parameterisierung.

Der Treiberanschluss bleibt bei `z = 0`; die Ausbreitung erfolgt nach `+Z`.
Intern gelten Meter und Radiant. Im Reiter **Spalt** aktiviert
**Spalt vor OS-SE aktivieren** die Halsbaugruppe.
Die OS-SE-Parameter bleiben im Profilabschnitt.

| Eingabe | Bedeutung |
| --- | --- |
| Anschluss B/H [mm], Form | Vorhandener physischer Treiberanschluss; Kreis bei B = H und Exponent 2. |
| Spaltbreite X [mm] | Volle lichte Ausdehnung in X, kein Halbradius. |
| Spalthöhe Y [mm] | Volle lichte Ausdehnung in Y. Breite und Höhe können für die zweite Orientierung getauscht werden. |
| Spaltkontur | Superellipse mit Formexponent oder Rechteck mit echtem Eckradius. |
| Eckradius [mm] | Kreisradius der vier Querschnittsecken, unabhängig vom Ausgang. 0 ergibt ein scharfes Rechteck; maximal die halbe kleinere lichte Spaltabmessung. |
| Adapterlänge [mm] | Länge vom Treiberanschluss bis zum Beginn des Spalts. |
| Spaltlänge [mm] | Länge des geraden Abschnitts mit konstantem Spaltquerschnitt; 0 bedeutet direkte Anschlusskante. |
| Spaltausgang | Axiale Glättung oder Kreisbogen mit echtem Radius im radialen Längsschnitt. |
| Glättungslänge [mm] | Im Modus Axiale Glättung: positive axiale Übergangslänge; 0 erhält die ursprüngliche Profiltangente. |
| Ausgangsradius [mm] | Im Kreisbogenmodus: echter Kreisradius zum Horn hin, unabhängig vom Eckradius. 0 erhält den ungerundeten Anschluss. |

Die Gesamtlänge bleibt die vorhandene Länge vom Treiberanschluss bis zum Mund.
Der OS-SE-Abschnitt erhält die verbleibende Länge. Änderungen am Adapter
verändern deshalb nachvollziehbar die verfügbare Aufweitung; sie verschieben
weder den Treiberanschluss noch heimlich den Mund. Die Anzeige nennt zusätzlich
die abgeleitete Spaltposition, verbleibende Aufweitungslänge und Spaltfläche.

Unterstützt wird ein einzelner mittiger Spalt, ausgerichtet nach X/Y.
Mehrfachspalte, Stege, Phase-Plugs und beliebige Verdrehungen sind eigene spätere
Geometrieaufgaben. Sie werden nicht durch inaktive Schaltflächen angekündigt.
Ein Spalt darf bei gleichzeitig größerer Höhe dieselbe oder eine größere Fläche
als der Anschluss besitzen; eine Flächenverengung wird nicht vorausgesetzt.

## Anschluss an die bestehende Architektur

`EntryParameters` beschreibt den Eintritt der gesamten Draft-Geometrie.
Dieser Vertrag bleibt der physische Treiberanschluss. Die Spaltkontur erhält
eigene Daten; einfach `entry.width` zu verkleinern würde den Treiberanschluss
und seine Driven-Fläche fälschlich verändern.

Die Komposition in `modules/geometry` verbindet drei axiale Segmente:
Formadapter, optionalen geraden Spalt und den anschließenden Profilabschnitt.
OS-SE bleibt eine `LongitudinalProfile`-Implementierung. Sie erhält für den
letzten Abschnitt die Spalt-Halbachsen und dessen lokale axiale Koordinate.
Der Adapter interpoliert die polaren Radien mit
`h(x) = 10x³ − 15x⁴ + 6x⁵`, `x = z / Adapterlänge`.
Er bleibt zwischen den positiven Endradien und besitzt an beiden Enden
verschwindende erste und zweite axiale Ableitungen. Damit treffen sich Adapter
und gerader Spalt C2-stetig. Die endliche Superellipse ist keine scharfkantige
Rechteckkontur; ihr Exponent ist mindestens 2.

Bei einer geometrischen Kante stimmen die Positionen beider Segmente überein;
abweichende einseitige Tangenten sind ausdrücklich Teil des Modells. Ein glatter
Anschluss verwendet über die Länge `t` die lokale Profilkoordinate
`g(u) = t·(6x³ − 8x⁴ + 3x⁵)`, `x = u/t`, `u = z − z_Spaltausgang`.
Ab `u = t` gilt wieder `g(u) = u`. Es gelten
`g(0) = g′(0) = g″(0) = 0` und `g(t) = t`, `g′(t) = 1`, `g″(t) = 0`.
Die Koordinate wächst monoton; der Anschluss ist C2-stetig, ohne die
verbleibende Profillänge oder Mundposition zu ändern.
Eine Kante darf nicht durch Vorschau-Normalenglättung als verrundet erscheinen.
Die axiale Glättung ist ausdrücklich kein Kreisradius.

Der alternative Rechteckquerschnitt verwendet vier Geraden und tangential
anschließende Kreisbögen. Für volle lichte Maße B und H und Eckradius r gilt
`A = B·H − (4−π)·r²`. Die folgende Aufweitung skaliert diese Kontur mit den
H/V-Profilradien; dabei können aus den Kreisecken elliptische Bögen werden.
Der eingestellte Eckradius bezeichnet die konstante lichte Spaltkontur.

Der Ausgangsradius R erzeugt in jedem radialen Längsschnitt einen echten
Kreisbogen. Seine Tangenten treffen den geraden Spalt und das unveränderte
nachfolgende Profil. Mit Profilradius f(u), Spalt-Polarradius a und lokaler
Koordinate u ab nominalem Spaltausgang erfüllt die Anschlussstelle t
`f(t)−a = R·(1−1/sqrt(1+f′(t)²))`; die axiale Kreismitte liegt bei
`zc = t−R·f′(t)/sqrt(1+f′(t)²)`. Der Kreis beginnt bei zc im geraden Spalt.
Die Gleichung wird aus Kreisgeometrie und Tangentengleichheit hergeleitet;
es ist kein allgemeiner räumlicher Rollkugelradius.

Der Radius wird unverändert übernommen. Fehlende Tangentiallösungen,
Profilsprünge oder ein Rückschnitt bis in den Adapter werden abgelehnt.
Ein ohnehin tangentialer Profilbeginn allein ist keine Radiusverrundung.
Die Anschlussstellen können je Umfangswinkel verschieden liegen; der
Formübergang muss hinter dem letzten Anschluss beginnen.

Leitkurvenposition und prozentuale Formübergänge beziehen sich weiterhin auf die
gesamte Hornlänge. Für die erste Stufe muss eine aktive Leitkurve hinter dem
Spaltübergang liegen. Die Profilinversion verwendet dort die lokale Position
`z - z_Spaltausgang`; der Treffer wird anschließend an der vollständig
komponierten Geometrie geprüft. Formübergang und Mund-Morphing dürfen den
definierten Anschluss oder Spalt nicht nachträglich verändern. Überlappende
Bereiche werden mit einer konkreten Eingabediagnose abgelehnt.
**Formbeginn hinter Spalt setzen** setzt nur den Beginn des glatten
Prozentverlaufs. Das Ende muss dahinter liegen. Die LUT bleibt unverändert;
sie muss bis zum Anschlussende Formanteil 0 halten. Im Kreisbogenmodus wird
der tatsächliche Anschlussbereich aus der gültigen Geometrie übernommen.
Solange keine gültige Geometrie vorliegt, nennt die Fehlerdiagnose die Ursache;
der Einstellknopf setzt keine geschätzte Position.

Analyse, Vorschau, STEP-Export und isotropes Netz verwenden dieselbe endgültige
Kontur. `A(z)`, `V(z)` und Wandwinkel werden daraus ausgewertet. An einer
bewussten Kante sind einseitige Wandwinkel nötig; ein undefinierter Ableitungswert
darf nicht als geglätteter Messpunkt erscheinen. Parameter, Undo/Redo,
Projektdatei, Baugruppeninstanzen und Sweep-Varianten verwenden denselben
unveränderlichen Parametersatz. Die Analyse nimmt die Segmentgrenzen explizit
auf; bei einem Knick stehen vor und nach einer Ableitungslücke die beiden
einseitigen Grenzwerte. STEP exportiert abschnittsweise gefittete Flächen mit
gemeinsamen Ringnähten und vorgeschriebenen Randtangenten. Der Export enthält
keine Deckelfläche im inneren Spalt. Die bestehende adaptive STEP-Prüfung bleibt
ein Stichprobennachweis, kein globaler Hausdorff-Beweis.
Eckradius 0 bleibt für Vorschau und Netz zulässig. Der STEP-Export lehnt diese
scharfe Rechteckkontur ausdrücklich ab: Deren bei der Aufweitung wandernde
Kantennähte werden vom aktuellen Flächenexport nicht unterstützt. Positive
Eckradien verwenden den regulären toleranzgeprüften Export.

## Netz und Akustik

Physische Abschnittskanten erhalten gemeinsame Randkurven und beständige
Flächenidentitäten im neutralen Geometrievertrag. Alle Spaltwände sind `Wall`.
Die lokalen Draft-IDs 1–5 behalten ihre Rollen. Formadapter (6), Spaltwand (7)
und Halsaußenwand (8) sind separat auswählbar, einfärbbar und mit lokalen
Zielkanten belegbar. Bei gerader Spaltlänge 0 entfällt Fläche 7.
Beim Kreisradius endet die Spaltgruppe am frühesten Kreisbogenbeginn über
alle Umfangswinkel; die anschließende Rundung gehört zur Hornwand. Damit
liegt keine nominelle Flächennaht mitten im glatten Übergang. Der nominelle
Spaltausgang bleibt eine Profilkoordinate, keine zusätzliche Netzkante.
Die einzige vorgegebene Anregung bleibt am realen Treiberanschluss; bei einem
gekoppelten LEM-Treiber gelten dessen bestehende Portregeln. Die innere offene
Spaltquerschnittsfläche erhält weder einen Deckel noch eine neue Driven-Fläche.
Die geschlossene Materialoberfläche mit vorhandenem Anregungsabschluss bleibt
über den bisherigen BEM-Pfad prüfbar.

Kleine Abstände und hohe Krümmung müssen über allgemeine Geometriegenauigkeit,
physische Randkurven und das skalare Größenfeld aufgelöst werden. Dafür gibt es
keine Hals-Sondervernetzung und keinen künstlichen Kragen. Unvereinbare
Mindestkante oder Formabweichung werden verständlich gemeldet. Eng gegenüberliegende
Flächen benötigen auch dann geeignete Nahintegration, wenn das lineare Residuum
klein ist. Netz- und Quadraturverfeinerung werden getrennt geprüft.

Ein rechteckiger Spalt ist nicht rotationssymmetrisch; der m=0-Modus ist dafür
unzulässig. Zulässige Spiegelsymmetrien folgen aus Geometrie und Anregung zusammen.
Einzelhörner und gekoppelte Instanzen nutzen weiterhin dieselben Solverobjekte.

Der lineare BEM-Pfad kann komplexen Druck, Phase, Impedanz und H/V-Richtwirkung
vergleichen. Er liefert keine belastbare Aussage zu Turbulenz, Strömungsablösung,
Kompression, harmonischen Verzerrungen oder maximalem PA-Pegel. Auch ein
Querschnittsflächenverhältnis ist kein Belastbarkeitsnachweis. Für solche Aussagen
sind zusätzliche physikalische Modelle oder Messungen erforderlich. Ein idealer
Kolben bildet die reale Hochton-Phasenverteilung eines Kompressionstreibers nur
näherungsweise ab.

## Prüfung und Grenzen

Die eigenen Tests prüfen Geometrie, Analyse, Projekt-Roundtrip, Undo/Redo,
Sweeps, Vorschau und STEP-Rückimport. Der Mesh-Test verwendet 6 × 18 mm lichte
Spaltmaße bei 60 mm Gesamtlänge, 12 mm Adapter, 3 mm geradem Spalt und 0 bzw.
4 mm Glättung. Er prüft geschlossene orientierte 2-Mannigfaltigkeit,
Schnittfreiheit, Quellfläche und Flächentags bei aktiven Qualitätsgrenzen.

Der CUDA-Test vergleicht ohne Spalt, mit Knick und mit Glättung bei gleichem
12-mm-Treiberanschluss, Medium und Anregungsbezug. Er vergleicht komplexen
Druck auf Achse und H/V bei 45° sowie die Eingangsimpedanz bei 1 und 4 kHz.
Netzverfeinerung und strengere Nahquadratur werden getrennt vom linearen
Residuum beurteilt. Das ist ein begrenzter numerischer Regressionstest, keine
Abnahme eines bestimmten PA-Horns oder des gesamten Hochtonbereichs.

Für einen konkreten Entwurf bleiben H/V-Isobaren, On-/Off-Axis-Verläufe,
komplexe Impedanz und Schallfeld gemeinsam zu bewerten. Die sechs bisherigen Spaltwerte und beide Radien
sind im vorhandenen Sweep-Pfad verfügbar; eine breitere H-Kurve allein
entscheidet nicht über die Eignung. Baugruppen bleiben gemäß dem bestehenden
Sweep-Vertrag von Parameterstudien ausgeschlossen.

Geometrische Pflichtfälle sind deaktivierter Spalt (identische bestehende
Geometrie), gleiche Anschluss-/Spaltkontur, Spaltlänge null, vertauschte X/Y-Maße,
erhaltener Leitkurventreffer und widersprüchliche Übergangsbereiche. Zu lange
Adapter, nichtpositive lichte Maße, nichtendliche Radien und Selbstschnitte
werden abgelehnt. Der tatsächlich ausgeführte Abnahmestand steht in
STATUS.md.
