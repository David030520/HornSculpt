# HornSculpt

Native C++20-/Qt-6-Anwendung fuer Hornentwurf, isotropes Oberflaechenmeshing,
CUDA-BEM mit Burton-Miller-Galerkin und einen optionalen lokalen FreeCAD-Link.
Die Arbeitsbereiche Draft und FreeCAD-Link verwenden gemeinsame Geometrie-,
Netz-, Simulations- und Ergebnisvertraege.

Dieses Repository enthaelt ausschliesslich das Tool mit Builddateien, Tests,
Lizenzhinweisen und ausgewaehlter Bedienungs-/Methodendokumentation.
Kundenmodelle, konkrete Horn-/Coax-Projekte, Ergebnisse, Laufzeitpakete,
persoenliche Daten und die private Projekthistorie werden nicht veroeffentlicht.
Eine Regression mit einem privaten CAD-Modellausschnitt bleibt samt Fixture
im privaten Projekt; die allgemeinen Tool-Tests sind enthalten.

## Bauen

Windows x64, Visual Studio 2022, CMake ab 3.25, Qt 6.8.3, CUDA Toolkit 13.3,
CGAL 6.2.1, Boost 1.88.0 und Open CASCADE 7.9.3. Abhaengigkeiten explizit
vorbereiten; der Build laedt keine Pakete herunter.

```powershell
.\build.ps1 -QtRoot <Qt-Prefix> -CgalRoot <CGAL-Prefix> `
  -BoostRoot <Boost-Prefix> -CudaRoot <CUDA-Prefix> `
  -OcctRoot <OCCT-SDK-Prefix> -OcctRuntimeDirectory <OCCT-Laufzeitordner> `
  -BuildDirectory <Buildordner-ausserhalb-des-Repositories> -Test
```

Danach den vorhandenen Build inkrementell aktualisieren. Versionen werden
ueber Git-Commits und inkrementelle Pushes gepflegt; keine vollstaendigen
Tool-Kopien als Historie anlegen. Eine Installation mit `cmake --install
<Buildordner> --config Release --prefix <Paketordner>` ist fuer konkrete
Auslieferungen vorgesehen. Start: `Start-HornSculpt.cmd` im installierten Paket.

Bedienung: [Kurzanleitung](docs/CUSTOMER-GUIDE.md),
[FreeCAD-Link](docs/FreeCAD-LINK.md), [Wissensbasis](docs/knowledge/README.md).
Lizenzhinweise fuer mitgelieferte Drittanbieterbestandteile stehen unter `licenses/`.
