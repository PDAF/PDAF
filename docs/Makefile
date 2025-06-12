BUILD_DIR = ./_build
SRC_DIR_LIB = ../src
SRC_DIR_INT = ../interface

all: docs src-browser-lib src-browser-interface

docs:
	jupyter-book build -W -n --keep-going .

src-browser-lib:
	ford -d $(SRC_DIR_LIB) -o $(BUILD_DIR)/html/src FORD_options.md

src-browser-interface:
	ford -d $(SRC_DIR_INT) -o $(BUILD_DIR)/html/interface FORD_options_int.md

.PHONY: clean

clean:
	jupyter-book clean .
