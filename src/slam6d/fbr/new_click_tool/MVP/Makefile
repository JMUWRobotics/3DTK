# Variablen
BUILD_DIR = build
EXEC = new_click_tool

# Default Target
all: configure
	@echo "Kompiliere das Projekt..."
	cmake --build $(BUILD_DIR) -j 8

configure:
	@echo "Konfiguriere CMake (Lokale 3rdparty-Abhaengigkeiten)..."
	cmake -B $(BUILD_DIR) -S .

clean:
	@echo "Loesche Build-Ordner..."
	rm -rf $(BUILD_DIR)

# fclean (Full clean) löscht zusätzlich den generierten Koordinaten-Ordner
fclean: clean
	@echo "Loesche generierte Koordinaten..."
	rm -rf Koordinaten

re: fclean all

# Run-Befehl 
run: all
	@echo "Starte $(EXEC)..."
	./$(BUILD_DIR)/$(EXEC) $(IMAGE)

.PHONY: all configure clean fclean re run
