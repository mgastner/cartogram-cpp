cartogram: main.cpp constants.h geo_div.cpp geo_div.h map_state.cpp\
	map_state.h print_to_eps.cpp print_to_eps.h read_csv.cpp read_csv.h\
	read_geojson.cpp read_geojson.h rescale_map.cpp rescale_map.h
	g++ -O3 -o cartogram main.cpp geo_div.cpp map_state.cpp print_to_eps.cpp\
	  read_geojson.cpp read_csv.cpp rescale_map.cpp -lboost_program_options\
		-lgmp
