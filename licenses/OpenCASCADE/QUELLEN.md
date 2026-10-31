# CAD-Abhängigkeiten

Der native STEP-/BRep-Pfad verwendet Open CASCADE **7.9.3**. Der eigene
Frontal-Delaunay-Mesher bleibt unabhängig von dessen Vorschautessellation.
FreeCAD ist weiterhin ein optionaler externer Prozess.

Der Build erwartet einen vorbereiteten SDK-Ordner und einen separaten Ordner
mit dessen Laufzeit-DLLs. `build.ps1` erhält beide Pfade explizit; CMake lädt
keine Pakete herunter. Header und Importbibliotheken werden nur im Build
benötigt. Das portable Paket enthält die benötigten dynamischen Laufzeiten.

Geprüfte Windows-x64-Artefakte aus den offiziellen OCCT-Releases:

| Artefakt | SHA-256 |
|---|---|
| [OCCT 7.9.3 VC14 x64 SDK](https://github.com/Open-Cascade-SAS/OCCT/releases/download/V7_9_3/opencascade-7.9.3-vc14-64.zip) | `2A4F5E187992427674869E89477882D1ACE4D67695D4B27173830D93F8CF11BE` |
| [Passende Drittanbieter-Laufzeiten](https://github.com/Open-Cascade-SAS/OCCT/releases/download/V7_9_3/3rdparty-vc14-64.zip) | `F25BE51AF6160A2DB74D19A4ADC0C08728517822167B01459F8E775927A2A2F6` |

Die STEP-Bibliothek importiert im vorkompilierten SDK auch XDE-/Darstellungs-
Bibliotheken. Deren vollständige Loader-Abhängigkeiten gehören zum Paket,
auch wenn HornSculpt seine eigene Qt-/OpenGL-Ansicht verwendet: TBB, jemalloc,
FreeImage, FreeType, OpenVR und die mitgelieferten FFmpeg-Laufzeiten. Letztere
melden über ihre tatsächliche Bibliotheks-API LGPL 2.1 oder später.
Debug-CRT-Laufzeiten werden nicht verteilt. Die explizite Laufzeitliste steht
in `modules/geometry/cad/CMakeLists.txt` und wird mit reduzierter PATH-Variable
im portablen Paket geprüft.

OCCT steht unter LGPL 2.1 mit eigener zusätzlicher Ausnahme. Das Paket enthält
die jeweiligen Lizenztexte und Quellenverweise unter `licenses/OpenCASCADE`.
Gmsh ist keine Laufzeit- oder Build-Abhängigkeit. Der CAD-Kern wird nicht aus
einer bestehenden HornSculpt-Installation übernommen.
