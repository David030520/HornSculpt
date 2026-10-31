# Simulations- und Validierungspraxis

## Vor der Rechnung

1. Quellenmodell, akustische Achse, Umgebung, Medium, Frequenzband und
   Zeitsignatur festlegen. Schalldruck als RMS oder Spitze und Bezug 20 µPa
   ausdrücklich nennen; 1 m Fernfeldnormierung ist kein endlicher 1-m-Punkt.
2. Geometrie auf physische Anschlüsse, geschlossene orientierte Topologie,
   Selbstüberschneidungen und Flächenzuordnung prüfen. Ein äußerer Luftspalt
   darf nicht durch ein erfundenes Wandstück geschlossen werden.
3. Wellenlänge λ = c/f und kleinste relevante geometrische Details betrachten.
   COMSOL nennt λ/N mit N ungefähr 5–10 als Startregel, abhängig von
   Diskretisierung [S18]; das ist keine Genauigkeitsgarantie für HornSculpt.
4. Treibergeschwindigkeit, komplexe Phase und Orientierung sowie elektrische
   Anregung und Rückwirkung unterscheiden. Unbekannte Daten als Annahme
   kennzeichnen.

## Nach der Rechnung

| Frage | Prüfverfahren |
|---|---|
| Algebra gelöst? | Relatives Residuum und Iterations-/Abbruchstatus prüfen. |
| Integrale gelöst? | Singularordnung und Nahquadratur getrennt variieren, enge Abstände einbeziehen. |
| Netz ausreichend? | Mindestens zwei, besser drei Netze mit gleicher Geometrie/Anregung vergleichen; komplexen Druck, Impedanz und Winkelschnitte betrachten. |
| Physik richtig? | Analytische Kugel und Kolben mit bekannter Normierung; Resonanzdurchlauf, Symmetrie und Superposition [S16, S17]. |
| Modell passend? | Treiber, Gehäuse, Öffnungen und Umgebung gegen reale Einbausituation prüfen; Messungen mit Abstand und Fensterung dokumentieren. |
| Referenz vergleichbar? | Identisches Netz, Medium, Frequenz, Anregung, Achsen, Phasor- und Impedanzkonventionen verwenden. |
| Leistung vergleichbar? | Gleiche Hardware, Genauigkeit und Warmzustand; Aufbau, Lösung, Nachverarbeitung, Gesamtzeit und Spitzenspeicher getrennt messen. |

Das lineare Residuum ist kein Ersatz für die übrigen Zeilen. Bei tiefen
Pegeln sind relative Fehler und Phase numerisch empfindlich; zusätzlich
absolute komplexe Fehler auswerten. Keine erfolgreiche GPU-Berechnung
behaupten, wenn nur CPU oder ein Gerätecheck lief.

## HornSculpt: vorhandener Nachweis und offene Fälle

Der dokumentierte Vergleich mit gepinntem Boundary Lab/BEAT umfasst einen
512-Dreiecks-Kugelfall und Operatorblöcke; ein feinerer eigener analytischer
Kugelfall prüft Konvergenz. Referenz-CUDA-Piloten für Kugel und ein Horn sind
dokumentiert. Die Einzelheiten und Einschränkungen stehen in
Paritätsstatus und
Solver-Leistungsmessung. Dies ist keine
vollständige Paritäts- oder Leistungsfreigabe für beliebige Modelle.

Für die weitere Abnahme: Horn und mindestens ein asymmetrisches CAD-Modell mit
identischen Netzen und mehreren Anregungen vergleichen; komplexe Fern- und
Nahfelder, Flächendrücke und normierte Impedanzen auswerten; Resonanzen und
Nullstellen gezielt abdecken; Netz und Quadratur unabhängig verfeinern.
Coax-Spalt und BEM↔LEM-Port zusätzlich auf Druck-/Volumenstrombilanz und
Sensitivität gegenüber angenommener Kammergeometrie prüfen. Messungen sind
notwendig, um Treiberparameter und die reale Konstruktion zu validieren.

Die CTA-Auswertung mittelt definierte H/V-Richtungen und verwendet eine
70-Richtungs-Näherung für die Schallleistung, keine vollständige
Kugelabtastung. Mehr dazu in [CTA-2034](../CTA-2034.md). Raumprognose und
akustischer Winkel sind keine Garantien für einen konkreten Raum.
