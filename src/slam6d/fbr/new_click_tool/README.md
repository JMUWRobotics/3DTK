# SoftwarePraktikum
Alle Materialen zur Organisatorischen Führung des Softwarepraktikum SS2026 Klick-Tool

## Struktur

- `MVP/`: C++ Prototyp des Klick-Tools
- `MVP/assets/`: Beispielbilder
- `MVP/coordinates/`: Beispielausgaben mit gespeicherten Koordinaten
- `docs/pflichtenheft/`: Pflichtenheft als LaTeX-Datei und PDF
- `docs/templates/`: Vorlagen für Projektdokumentation



Grafiktool: ImGUI -> enable. PhD-Student that can (possibly help):Fabian Arzberger

1.	Load one Image (png, jpg) and show/return coordinates of clicked pixel.
	  Up to 0.5 subpixel-accuracy (maybe 1/8 if we double the pixels with each zoom-step?)
	  Option to create png from 3DTK (scan to panorama – function) and load as png automatically

2. 	Load two images like in 1. on top of each other or next to each other. Maybe compare picture size to screen size to choose layout automatically.
	  Click, display and return corresponding coordinates.
	  Thin line between points, to see corresponding points.

3. 	Integrate. Commandline-parameters + description, see project icpFixpoint as example
