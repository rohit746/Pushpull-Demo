CXX ?= c++
PKG_CONFIG ?= pkg-config

CXXFLAGS += -std=c++17 -O2 -Wall -Wextra -pedantic $(shell $(PKG_CONFIG) --cflags sdl2)
LDLIBS += $(shell $(PKG_CONFIG) --libs sdl2)

.PHONY: all clean run

all: pushpull_pp

pushpull_pp: src/main.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDLIBS)

run: pushpull_pp
	./pushpull_pp

clean:
	rm -f pushpull_pp
