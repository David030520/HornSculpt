# Schallfeld und Phase

Das Register **Phase** und die Schallfeldebene im oberen **3D-Modell** zeigen
ausschließlich komplexe Ergebnisdaten eines bestimmten BEM-Snapshots.
Das Register **Schallfeld** enthält dessen Berechnungs- und Farbsteuerung.
Darstellungsparameter werden im
Projekt gespeichert. Ihre Änderung verändert weder Anregung noch Solverrevision.
Veraltete Ergebnisse tragen die Kennzeichnung ihrer ursprünglichen Revision.

## Schallfeld

1. Eine Simulation berechnen. Neu berechnete Ergebnisse enthalten zusätzlich
   den vollständigen kontinuierlichen P1-Druck an den Netzpunkten und die
   stückweise konstante DP0-Normalableitung auf den Dreiecken. Ältere Dateien
   ohne diese Randwerte bleiben lesbar; für ihre Feldansicht ist eine neue
   Simulation erforderlich.
2. Im Register **Schallfeld** Frequenz, Anregung und H- oder V-Schnitt auswählen.
   Quer- und Axialbereich beziehen sich auf den gespeicherten Beobachtungsrahmen
   des Ergebnisses. Die Kamera und spätere Achsänderungen verschieben diesen
   Ergebnisrahmen nicht. Frequenz und Anregung werden als Darstellungswahl
   gespeichert; fehlt die gespeicherte Auswahl im Ergebnis, wird eine verfügbare
   gewählt. Nach dem Öffnen wird die Ebene erst auf Anforderung neu berechnet.
3. Raster und Randabstand einstellen und **Feld berechnen** wählen. Der
   abbrechbare Hintergrundjob löst keine neue BEM-Randwertaufgabe. Er berechnet
   das Helmholtz-Integral an den gewählten endlichen Beobachtungspunkten.
4. Zwischen **Pegel**, **Phase** und **Re(p·eⁱᵠ)** wechseln. Der Phasenregler φ
   dreht die gespeicherten komplexen RMS-Druckwerte für den Realteil. Die
   Spitzenauslenkung des harmonischen Drucks beträgt √2 mal den RMS-Realteil.
5. Die berechnete Ebene erscheint mit dem CAD im oberen 3D-Fenster. Drehen,
   Verschieben und Zoomen verwenden dieselbe Kamera und dieselben tatsächlichen
   Weltkoordinaten. Beim ersten Feld eines geöffneten Projekts wird eine
   schräge gemeinsame Ansicht eingepasst. Frequenz- und Farbänderungen erhalten
   danach die manuelle Kamera. **3D einpassen** stellt die Feldansicht erneut
   ein; **Ansicht → Schallfeldebene** blendet sie ein oder aus.

Der Pegel verwendet 20 µPa als Bezug. Phase hat eine zyklische Farbskala von
−180° bis +180°; Realteil eine vorzeichenbehaftete Skala in Pa RMS. Die
Darstellung besteht aus echten 3D-Dreiecken zwischen den berechneten
Weltpunkten. Der Fragmentshader interpoliert zuerst den komplexen Druck und
wertet daraus die Farbe aus. Er mischt keine gewickelten Phasenwinkel oder
Phasenfarben über den Sprung bei ±180°. Diese Darstellungsinterpolation
erzeugt keine zusätzlichen akustischen Stützstellen.
Vor der Darstellung werden auch maskierte Weltpunkte auf Übereinstimmung mit
Schnittrahmen, Quer-/Axialbereich und Rasterindizes geprüft. Widersprüchliche
Koordinaten werden abgewiesen; die Anzeige verschiebt keine Proben.

Die Feldebene ist opak und unbeleuchtet, damit ihre Farben der eingeblendeten
Skala entsprechen. CAD und Feld verwenden denselben Tiefenpuffer: Vorderes
CAD verdeckt dahinterliegende Feldflächen; vordere Feldflächen verdecken
hinteres CAD. Die Ebene wird nicht als transparente Überlagerung durch feste
Geometrie hindurch gezeichnet. Frequenz, H/V-Schnitt, Revision und eine
gegebenenfalls veraltete Revision stehen im 3D-Fenster.

Festkörper, ausgeschlossene Randnähe und nicht konvergierte Quadraturpunkte
werden ausgespart. Dreiecke mit einer ungültigen Rasterstützstelle werden
nicht gezeichnet; darunter liegende Geometrie bleibt sichtbar.
Auf der unendlichen Schallwandebene werden keine singulären Grenzwerte erzeugt.
Außerhalb der Luftkavität ist hinter der Schallwand kein abstrahlender Halbraum.
Die Anzeige nennt die Anzahl der nicht konvergierten Punkte. Ein kleinerer
Randabstand kann mehr Nahquadratur benötigen; eine ausgesparte Stelle ist kein
berechneter Schalldruck von null.

Bei einem berechneten komplexen Druck von exakt null ist die Phase undefiniert
und bleibt ebenfalls ausgespart. Derselbe gültige Nullwert erscheint im Pegel
an der unteren Skala und im Realteil weiß.

Mindestens etwa sechs Rasterpunkte je Wellenlänge sind eine Orientierung für
die Abbildung des laufenden Feldes, kein Konvergenznachweis. Enge Spalte und
hohe lokale Gradienten benötigen feinere Schnitte. Ein dichteres Darstellungs-
raster verbessert nicht die Genauigkeit eines unzureichenden BEM-Netzes.

**CSV …** exportiert die tatsächlichen komplexen Rasterwerte, Weltkoordinaten
und eine Gültigkeitsmaske. **Ansicht → 3D-Bild exportieren …** exportiert den
aktuellen 3D-Framebuffer einschließlich CAD, Feld und Skala als PNG.
Das komplexe Gitter wird im Arbeitsspeicher gehalten; Projektdateien speichern
die vollständigen BEM-Randwerte und die gewählten Anzeigeparameter. Ein
geöffneter Feldschnitt lässt sich daraus erneut auswerten.

Für einen ausdrücklich angeforderten Feldstart beim Öffnen dient die Option
`--evaluate-field`. Beispielsweise öffnet
`HornSculpt.exe Modell.hsproj --view field --evaluate-field` die Feldansicht
und berechnet den Schnitt mit den gespeicherten Anzeigeparametern. Ohne diese
Option wird kein Feldjob automatisch gestartet. Zusammen mit
`--screenshot Feld.png` wartet die Aufnahme auf das abgeschlossene Feld; ein
fehlendes Ergebnis, Abbruch oder Berechnungsfehler führt zu einem Fehlerstatus
anstelle eines Screenshots des leeren Bereitschaftszustands.

## Integral und Konventionen

Es gelten RMS-Phasoren mit `exp(+iωt)` und der auslaufenden Greenschen Funktion
`G = exp(−ikr)/(4πr)`. `q = ∂p/∂ν` gehört stets zur äußeren Normale des
geschlossenen Netzes. Die gespeicherten Cauchy-Daten verwenden die volle
Netzreihenfolge auch bei einer rotationssymmetrisch reduzierten BEM-Lösung.

- Außenraum eines Festkörpers: `p(x) = D p − S q`.
- Inneres der geschlossenen Luftkavität: `p(x) = S q − D p`.
- Abstrahlender Halbraum der unendlichen Schallwand:
  `p(x) = −2 S_A q_A`, ausschließlich über die tatsächliche Apertur.

Die Feldnachauswertung verwendet eigene adaptive Dreiecksquadratur mit
Grad-4-/Grad-8-Vergleich, begrenzter räumlicher Unterteilung und zusätzlicher
Wellenlängenauflösung. Bei überschrittener Tiefe wird der Punkt ausdrücklich
ausgespart. Die Luft-/Festkörperzuordnung verwendet den Raumwinkel des
geschlossenen Netzes. Der dargestellte Feldschnitt ist eine endliche
Nah-/Fernraumauswertung, kein aus Polardaten hochgerechnetes Bild.

## Phase

Das Register **Phase** verwendet den komplexen Druck des ausgewählten Winkels
für H und V. Zwischen berechneten Winkelstützpunkten wird komplex linear
interpoliert. Die Option **Bezug auf Winkel** zieht die Phase des einstellbaren
Bezugswinkels derselben Ebene ab. Fehlende Werte und exakte Drucknullstellen
bleiben Lücken. Sprünge der gewickelten Phase werden ebenfalls unterbrochen.

**Entwickeln** wählt zwischen benachbarten Frequenzstützstellen jeweils den
nächsten 360°-Zweig. Ein grobes Frequenzraster kann vollständige Umläufe
verbergen; diese sind aus den vorhandenen Daten nicht rekonstruierbar.
**Laufzeit** addiert `360°·f·t` zur dargestellten Phase. Die Kompensation ist
eine Anzeigeoperation; sie verändert die gespeicherten Drücke nicht.

## Referenz und Prüfungen

Als Bedienreferenz wurde die lokal ausgelieferte AKABAK-Hilfe gelesen:
**Observation / Fields (Nodes)** beschreibt regelmäßig verteilte
Beobachtungspunkte auf Ebenen sowie eine frequenzabhängige Druckkontur.
**Export Fields** beschreibt Koordinaten und Real-/Imaginärwerte.
Die eigenen Ansichten und Integralverfahren wurden unabhängig implementiert;
AKABAK ist keine Laufzeit- oder Testabhängigkeit.

Die REA-Prüfung des ausgelieferten AKABAK-Programms identifizierte das PE-Artefakt
mit SHA-256 `605b802202249268d40a653c5793e4d5d3e4f1acbbfadf71cf6c6f8013610c11`
und dem Evidenzbeleg
`ev_a46c3a564382b2c9fda76359ea5633c337ef7e373980f77148bbed806065f1ff`.
Sie inventarisierte statische Ressourcen. Daraus wird keine Aussage über
ausgeführte AKABAK-Berechnungen oder dessen interne Numerik abgeleitet.

Die gezielten Tests prüfen den Außenraum gegen einen analytischen Monopol,
den Aperturhalbraum gegen das exakte Achsenintegral einer Kreisscheibe und
die Kavitätsorientierung gegen eine ebene Welle. Weitere Prüfungen betreffen
Geometriemasken, Quadraturgrenzen, Abbruch, komplexe Phasenbezüge, Entwicklung,
asynchrone Veröffentlichung und den Projekt-Roundtrip der Randwerte.
Den tatsächlich ausgeführten Prüfstand nennt die Abschlussmeldung.

Geprüfter Teilstand vor der Umstellung auf die 3D-Ebene: Release-Build der CPU-Feldbibliothek, Feld-/Phasenwidgets
und Projektpersistenz ohne neue Warnungen; CTest `field_evaluation` (0,05 s),
`field_phase` (0,17 s) und `project_file` (22,73 s) bestanden. Zusätzlich wurden
die drei Feldgrößen und die Phasenkurven mit dem nativen Windows-Renderer
visuell geprüft. Das Feldbild dieser Widgetprüfung verwendet ausdrücklich
analytische Cauchy-Testdaten und ist kein Coax-Simulationsergebnis.
Die CPU-Nachauswertung wurde bei der gemeinsamen Integration zusätzlich gegen
die tatsächlich berechneten CUDA-Aperturwerte an endlichen Beobachtern geprüft.
Die 3D-Ebene besitzt gezielte Rendererprüfungen für Tiefenverdeckung,
Masken, Phasenumschlag, Realteil und Kameraerhalt.

Am 30. Oktober 2026 wurde der produktive Release-`MeshViewport` in einem
separaten unsichtbaren QA-Prozess mit echten gespeicherten BEM-Randwerten
geprüft: am HF-Draft sowie am CAD-Netz des ST260-30 mit separat modelliertem
idealem Halskolben. Eine 1-kHz-Feldebene wurde jeweils aus den gespeicherten
Spuren berechnet und mit Geometrie, Pegel- und Phasenlegende in einen
1000 × 700 Pixel großen OpenGL-Framebuffer gerendert und zurückgelesen.
Win32-Prüfungen bestätigten während der Läufe, dass kein prozesseigenes
Fenster sichtbar war und das Vordergrundfenster unverändert blieb. Die Bilder
wurden visuell geprüft; ein Pixelvergleich ist kein physikalischer Nachweis
der Farbskala.

Eine separate warme Render-/Kamera-/Readback-Messung auf Intel Arc 140T
verwendete je Fall zehn Vorlauf- und 40 Messbilder. Die 95. Perzentile der
vier Draft-/CAD- und Pegel-/Phasenfälle lagen zwischen 6,93 und 8,07 ms.
Enthalten sind Kamera-/Einpassungsarbeit, Rendern und GPU-zu-CPU-Readback,
nicht Solver, Feldberechnung oder Dateispeicherung. Rohdaten, Artefakthashes
und Versuchsgrenzen liegen im externen QA-Nachweis `hidden-gl/README.md`.
Dies ist keine Messung sichtbarer FPS oder Desktop-Eingabelatenz und keine
Renderabnahme des installierten Desktopprozesses oder der NVIDIA-Grafik.
