cartogram: main.cpp cgal_typedef.h constants.h fill_with_density.cpp\
	fill_with_density.h geo_div.cpp geo_div.h map_state.cpp map_state.h\
	read_csv.cpp read_csv.h read_geojson.cpp read_geojson.h rescale_map.cpp\
	rescale_map.h write_eps.cpp write_eps.h
	g++ -O3 -o cartogram main.cpp fill_with_density.cpp geo_div.cpp\
	  map_state.cpp read_geojson.cpp read_csv.cpp rescale_map.cpp write_eps.cpp\
		-lboost_program_options -lgmp -std=c++17 -pthread
