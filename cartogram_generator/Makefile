#CXXFLAGS += -I include -std=c++14 -Wall -Wextra -Wshadow -Werror -O3 -fPIC
CXXFLAGS += -I include -std=c++14 -Wall -Wextra -O3 -fPIC

MASON ?= mason
VARIANT = variant 1.1.4
GEOMETRY = geometry 1.1.0
RAPIDJSON = rapidjson 1.1.0

DEPS = `$(MASON) cflags $(VARIANT)` `$(MASON) cflags $(GEOMETRY)`
RAPIDJSON_DEP = `$(MASON) cflags $(RAPIDJSON)`

default: cartogram

mason_packages/headers/geometry:
	$(MASON) install $(VARIANT)
	$(MASON) install $(GEOMETRY)
	$(MASON) install $(RAPIDJSON)

build:
	mkdir -p build

CFLAGS += -fvisibility=hidden

build/geojson.o: src/mapbox/geojson.cpp include/mapbox/geojson.hpp include/mapbox/geojson_impl.hpp include/mapbox/geojson_value_impl.hpp build mason_packages/headers/geometry Makefile
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(DEPS) $(RAPIDJSON_DEP) -c $< -o $@

build/libgeojson.a: build/geojson.o
	$(AR) -rcs $@ $<

cartogram: main.cpp build/libgeojson.a
	$(CXX) $(CFLAGS) $(CXXFLAGS) $(DEPS) $(RAPIDJSON_DEP) $< -Lbuild -lgeojson -lgmp -o $@
