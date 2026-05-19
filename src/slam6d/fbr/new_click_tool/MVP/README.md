# ImGui Klick-Tool

Ein C++ Tool mit Dear ImGui zur visuellen Auswahl und Speicherung von Bildkoordinaten.

---

## Voraussetzungen (Requirements)

Damit das Projekt gebaut werden kann, wird ein C++ Compiler, **CMake** und **Git**. Wähle die passenden Befehle für dein Betriebssystem:

### macOS
Nutze den Paketmanager [Homebrew](https://brew.sh/):
1. C++ Compiler (Apple Command Line Tools) installieren:
   `xcode-select --install`
2. Werkzeuge via Brew installieren:
   `brew install cmake pkg-config git`

### Linux: Ubuntu / Debian
`sudo apt update` /
`sudo apt install build-essential cmake git curl zip unzip tar pkg-config pkg-config libxinerama-dev libxcursor-dev xorg-dev libglu1-mesa-dev`

### Linux: Arch
`sudo pacman -Syu` /
`sudo pacman -S base-devel cmake git curl zip unzip tar pkg-config`

---

## Installation und Bauen

**1. Repository klonen & Ordner öffnen**
`git clone <repo-link>` /
`cd 3DTK/src/slam6d/fbr/new_click_tool/MVP`

**2. Projekt bauen**
`make`

**3. Tool mit Beispielbild starten**
`./build/tool assets/K2_pic.jpg`

Alternativ kann ein anderes Bild angegeben werden:
`./build/tool assets/mein_Bild.jpg`

Weitere Makefile-Regeln:
- `make clean`: Build-Artefakte innerhalb des Build-Ordners entfernen
- `make fclean`: kompletten Build-Ordner entfernen
- `make re`: Build-Ordner entfernen und Projekt neu bauen
