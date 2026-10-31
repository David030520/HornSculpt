# Wissensbasis Hornentwicklung und Simulation

Stand 29.10.2026. Diese Übersicht richtet sich an Entwickler und Anwender.
Sie ist eine Auswahl, kein vollständiges Literaturverzeichnis. Die einzelnen
Quellen, ihre Zugänglichkeit und der tatsächlich geprüfte Umfang stehen im
[Quellenverzeichnis](sources.md). Referenzen [S01] usw. beziehen sich darauf.

- [Hornprofile, Leitkurven, Richtwirkung und Moden](horns.md)
- [BEM, Quadratur, Beschleunigung und Kopplung](bem.md)
- [Akustische Teilgebiete, Interfaces und Gmsh-Import](../SUBDOMAINS.md)
- [Simulation, Vergleich und Verifikation](validation.md)
- Entscheidungsweg und durchgerechnete Beispiele
- [AKABAK als externe Bedien- und Funktionsreferenz](akabak-reference.md)
- [Quellen und Lesestatus](sources.md)

**Grundsatz:** Ein glattes geometrisches Profil ist keine Vorhersage der
Abstrahlung. Ein kleines lineares Residuum beweist weder ein ausreichend
feines Netz noch die richtige physikalische Modellierung. Ergebnisse gelten
für den dokumentierten Geometrie-, Randbedingungs- und Solverstand.

Boundary Lab von JWSound/BEAT ist die beauftragte externe Vergleichsreferenz,
isoliert und auf die in parity-reference-status.md
genannten Commits gepinnt. Analytische Fälle bleiben eine unabhängige zweite
Referenz. Weder ein bestandener Kugelfall noch eine installierte
Referenzumgebung belegen vollständige Parität oder Geschwindigkeit.

Hinweise auf kommerzielle Produkte beschreiben öffentliche Dokumentation oder
einen externen Vergleich; sie sind keine Laufzeitabhängigkeit von HornSculpt.
Urheberrechtlich geschützte Artikel werden nicht in Projekt oder Paket kopiert.
