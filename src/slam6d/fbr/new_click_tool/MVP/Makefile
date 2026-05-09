VCPKG_ROOT ?= $(HOME)/vcpkg

all: configure
	cmake --build build

configure:
	cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE=$(VCPKG_ROOT)/scripts/buildsystems/vcpkg.cmake

clean:
	rm -rf build

re: fclean all

.PHONY: all configure fclean re
