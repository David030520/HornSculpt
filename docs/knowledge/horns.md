# Hornentwicklung: Form, Abstrahlung und Grenzen

## Profil und Querschnitt

Das Längsprofil beschreibt die Öffnung vom Hals zum Mund; eine
Querschnittskurve beschreibt den Umfang an einer axialen Position. Beide sind
unterschiedliche Freiheitsgrade. OS-SE verbindet eine profilbildende
hyperbolische Komponente mit einer glatten Termination; die veröffentlichten
Parameter erlauben u. a. den konischen Grenzfall. Der geometrische
Öffnungswinkel ist kein garantierter akustischer Abstrahlwinkel [S01, S02].
Die gewählte Termination, Schallwand und Treiberanbindung beeinflussen die
tatsächliche Lösung.

Eine Leitkurve ist ein geschlossener räumlicher Zwang an einer **bestimmten
axialen Position**: die Profile müssen diese Kurve treffen. Sie ist weder bloß
eine Mundkontur noch das nachträgliche Mund-Morphing. Das ATH-Handbuch trennt
explizite Profilvorgabe, Leitkurve und Morphing [S02]. In HornSculpt müssen
Kurvenform, physische Größe und Rotation sowie axialer Formübergang getrennt
und das Treffen der Kurve nach allen Transformationen geprüft werden.
Unlösbare Kombinationen dürfen nicht unsichtbar korrigiert werden.

## Wirkungsgrad, Mündung und Richtwirkung

Keele untersucht den Einfluss der Mundgröße auf Reflexionen bei exponentiellen
Hörnern und beschreibt die Schwierigkeit konstanter Richtwirkung mit rein
exponentiellen Profilen [S03, S04]. In seinem eindimensionalen Exponentialhorn-
Modell ergibt die Optimierung über `f_c` bis `10 f_c` für die angenommene
Halbraumlast `k_c a_M = 0,935`, mit `k_c = 2πf_c/c` und Mundradius `a_M`
[S03, PDF-S. 13 und 20]. Andere Abstrahlmodelle und Frequenzbereiche liefern
andere Werte; das Modell sagt selbst keine Richtwirkung voraus. Die untersuchten
konisch-exponentiellen Hörner verwenden einen vollen −6-dB-Winkel als
Coverage-Definition [S04, PDF-S. 4 und 20].
Daraus folgt kein universeller optimaler Mund für alle Profile und
Einbausituationen. Die freie Mündung, endliche Schallwand und nahe Wände sind
verschiedene Randwertprobleme. Kolbrek zeigt bei reflektierenden Begrenzungen
die Bedeutung gegenseitiger Strahlungsimpedanzen; deren Wirkung hängt vom
Abstand und vom Durchlassbereich des Horns ab [S10].

Polars und Isobaren ergeben sich aus dem komplexen berechneten Feld. Für einen
Winkelvergleich Bezug, Frequenz und Schwellenwert nennen; die erste
Schwellenkreuzung der Hauptkeule ist nicht identisch mit jedem möglichen
Coverage-Maß. Mehrere Quellen können interferieren: ihre komplexen Drücke
müssen unter gleicher Zeitkonvention addiert werden, bevor daraus Pegel
gebildet werden. Geometrische Schnitte allein liefern kein räumliches
Abstrahlverhalten [S09].

## Beugung und höhere Moden

Kanten, Mundumgebung und abrupte Übergänge können reflektierte und gebeugte
Anteile erzeugen; deren Abstand und Phase zur Direktwelle verändern
Frequenzgang und Richtwirkung. Eine starke Glättung der Ergebnisdarstellung
kann diese Effekte verdecken. Vergleichsrechnungen sollten Mundradius,
Schallwandabmessung und Umgebung getrennt variieren und Netzkonvergenz an den
betroffenen Kanten prüfen. Ein ideales unendliches Baffle ist kein endliches
Gehäuse [S10, S16, S17].

Kolbrek erläutert, dass die Schlussfolgerung zur optimalen Mündungsgröße im
Keele-Modell auf ebener Wellenfront und Kolbenlast in unendlicher Schallwand
beruht; bei seiner Betrachtung gekrümmter Wellenfronten verändert sich das
Ergebnis [S34, S. 9–10]. Die Regel `kr≈1` ist daher eine Orientierung für
die untersuchte Klasse, kein allgemeiner Entwurfsbeweis. Das ist eine Einordnung durch Kolbrek; Keele benennt die Grenzen seiner
angenommenen Mündungslast selbst [S03, PDF-S. 13].
JBL vergleicht in seiner PT-Notiz unterschiedliche Horn-Treiber-Paare
[S27, S. 6–9]; der Vergleich ist deshalb keine isolierte Messung allein
der Hornform.

Ein eindimensionales Hornmodell beschreibt eine nützliche
Niederfrequenznäherung, enthält aber die räumliche Mundverteilung nicht
allgemein. Höhere Quermoden können oberhalb ihrer jeweiligen Grenzfrequenz
ausbreitungsfähig werden; auch darunter können abklingende Moden an Übergängen
wichtig sein. Mode Matching koppelt die Moden benachbarter geometrisch
einfacher Abschnitte und verlangt Konvergenz sowohl der Segmentierung als auch
der Modenzahl [S08, S09, S10]. Eine 2D-Rotationsrechnung mit m=0 setzt
rotationssymmetrische Geometrie **und** Anregung voraus; sie repräsentiert keine
allgemeinen azimutalen Moden.

## Treiber und Tiefton

Keele verknüpft in einem vereinfachten Tiefton-Ersatzmodell Treiberparameter,
Halsfläche sowie Vorder- und Rückkammer [S05, PDF-S. 5–7, 11–13]. Die
angenommene Verlustvereinfachung und Hornlast begrenzen die Gleichungen;
maximale nominale Effizienz ist nicht gleich maximale Bandbreite. Eine
Hubgrenze allein beschreibt weder thermische Grenzen noch nichtlineare
Kompression.

Eine feste Normalgeschwindigkeit untersucht die akustische Last, beantwortet
aber nicht automatisch die Frage nach einem spannungsgetriebenen Treiber.
Das elektrische und mechanische Ersatznetz muss die Druckrückwirkung
aufnehmen; Kolbrek betrachtet endliche Hornlänge, Last und
Parameterempfindlichkeit [S11], COMSOL stellt feste Bewegung und gekoppelte
Treibermodelle gegenüber [S21]. Daten eines ähnlichen Treibers bleiben
gekennzeichnete Annahmen. Ein kleines Rückvolumen ist bei genügend hoher
Frequenz nicht zwangsläufig durch eine einzelne Nachgiebigkeit beschrieben.

Kolbreks Tieftonarbeit [S11, S. 4–9] zeigt am vereinfachten Modell, dass eine
Kompensation der Reaktanz an nur einer Frequenz im übrigen Durchlassbereich
ungünstig sein kann. Sie untersucht außerdem endliche Hornlasten und
Parameterempfindlichkeit; in ihrem konkreten Beispiel reagieren Ergebnisse
unter anderem auf Kraftfaktor, Gleichstromwiderstand und bewegte Masse.
Der erwärmte Schwingspulenwiderstand ist deshalb ein relevanter Modellparameter.
Diese Fallstudie ersetzt weder eine Messung der konkreten Einheit noch eine
nichtlineare thermische Simulation.
