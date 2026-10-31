# Diagrammskalierung

Das Kontextmenü jedes Linien-, Polar- und Isobarendiagramms enthält
**Skalierung …** und **PNG exportieren …**. Die Eingaben verwenden die
angezeigte Achseneinheit, bei Frequenzen Hz und bei Winkeln Grad.

- Linien- und Isobarenachsen unterstützen automatische oder feste Grenzen
  sowie eine eigene Schrittweite. Ein leerer Schritt bedeutet automatisch.
  Frequenzachsen bleiben logarithmisch; feste Frequenzgrenzen müssen positiv sein.
- Polars besitzen eine dB-Skala und einen Winkelschritt. Die dB-Vorgabe ist
  −30 bis +3 dB mit 3-dB-Schritten. Der äußere Ring zeigt das Maximum.
- Isobaren verwenden dieselben dB-Grenzen und Farbstufen. Außerhalb der
  vorhandenen Frequenz-/Winkeldaten bleibt die Fläche leer; es wird nicht
  über den Messbereich extrapoliert.
- Fehlende Daten (`NaN`) unterbrechen Kurven. Bei fester dB-Skala wird `−∞ dB`
  bei verschwindendem Druck an der unteren Skala dargestellt; +∞ an der oberen.

Die gemeinsame akustische dB-Skala gehört zu den Simulationseinstellungen im
Projekt. Einzelne Achseneinstellungen und Polar-Winkelschritte werden über
`QSettings` am stabilen `objectName` des Diagramms gespeichert. Neue Daten
setzen diese Einstellungen nicht zurück. Die Anwendung setzt ihre Vorgaben
vor dem `objectName`; spätere explizite API-Aufrufe überschreiben die Skala.

`PlotAxisScale` enthält optional gemeinsame Min-/Max-Grenzen und einen optionalen
positiven Schritt. `PlotDbScale` verlangt endliche, aufsteigende Grenzen und
höchstens 1000 Teilstriche einschließlich Enden. Die Eingabe akzeptiert auch
wissenschaftliche Schreibweise und rundet nicht auf eine feste Dezimalzahl.
Bei sehr dichten Achsen werden Beschriftungen ausgedünnt; die gespeicherten
Grenzen und Schritte bleiben erhalten. PNG-Ausgaben übernehmen die sichtbare Skala.
