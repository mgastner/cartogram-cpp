// TO DO: positional matching of argument flags

#include "constants.h"
#include "map_state.h"
#include "blur_density.h"
#include "fill_with_density.h"
#include "flatten_density.h"
#include "densify.h"
#include "project.h"
#include "read_csv.h"
#include "read_geojson.h"
#include "rescale_map.h"
#include "write_eps.h"
#include "check_topology.h"
#include "write_to_json.h"
#include "geo_div.h"
#include "densification_points.h"
#include <boost/program_options.hpp>
#include <iostream>
#include <cmath>

// Functions that are called if the corresponding command-line options are
// present
void on_geometry(const std::string geometry_file_name)
{
  std::cerr << "Using geometry from file " << geometry_file_name << std::endl;
  return;
}

void on_visual_variable_file(const std::string geometry_file_name)
{
  std::cerr << "Using visual variables from file "
            << geometry_file_name
            << std::endl;
  return;
}

int main(const int argc, const char *argv[])
{

  using namespace boost::program_options;
  std::string geo_file_name;

  // Default number of grid cells along longer Cartesian coordinate axis.
  int long_grid_side_length = default_long_grid_side_length;

  // World maps need special projections. By default, we assume that the
  // input map is not a world map.
  bool world;

  // Other boolean values that are needed to parse the command line arguments
  bool input_polygons_to_eps,
       density_to_eps;

  // Parse command-line options. See
  // https://theboostcpplibraries.com/boost.program_options
  variables_map vm;
  try {
    options_description desc{"Options"};
    desc.add_options()(
      "help,h", "Help screen"
      )(
      "geometry,g",
      value<std::string>(&geo_file_name)->required()->notifier(on_geometry),
      "GeoJSON file"
      )(
      "visual_variable_file,v",
      value<std::string>()->notifier(on_visual_variable_file),
      "CSV file with ID, area, and (optionally) colour"
      )(
      "id,i",
      value<std::string>(),
      "Column name for IDs of geographic divisions (default: 1st CSV column)"
      )(
      "area,a",
      value<std::string>(),
      "Column name for target areas (default: 2nd CSV column)"
      )(
      "color,c",
      value<std::string>(),
      "Column name for colors (default: 3rd CSV column if it exists)"
      )(
      "long_grid_side_length,l",
      value<int>(&long_grid_side_length),
      "Number of grid cells along longer Cartesian coordinate axis"
      )(
      "world,w",
      value<bool>(&world)->default_value(false)->implicit_value(false),
      "Boolean: is input a world map in longitude-latitude format?"
      )(
      "input_polygons_to_eps",
      value<bool>(&input_polygons_to_eps)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: make EPS image input_polygons.eps?"
      )(
      "density_to_eps",
      value<bool>(&density_to_eps)
      ->default_value(false)
      ->implicit_value(true),
      "Boolean: make EPS images input_*.eps?"
      );
    store(parse_command_line(argc, argv, desc), vm);
    if (vm.count("help") || vm.empty()) {
      std::cerr << desc << '\n';
      return EXIT_SUCCESS;
    } else {
      notify(vm);  // Triggers notifier functions such as on_geometry()
    }
  } catch (const error &ex) {
    std::cerr << "ERROR: " << ex.what() << std::endl;
    return EXIT_FAILURE;
  }
  MapState map_state(vm["visual_variable_file"].as<std::string>(),
                     world,
                     density_to_eps);

  // Read visual variables (e.g. area, color) from CSV
  read_csv(vm, &map_state);

  // Read geometry
  try {
    read_geojson(geo_file_name, &map_state);
  } catch (const std::system_error& e) {
    std::cerr << "ERROR: "
              << e.what()
              << " ("
              << e.code()
              << ")"
              << std::endl;
    return EXIT_FAILURE;
  }

  try {
    holes_inside_polygons(&map_state);
  } catch (const std::system_error& e) {
    std::cerr << "ERROR: " << e.what() << " (" << e.code() << ")" << std::endl;
    return EXIT_FAILURE;
  }

  // Rescale map to fit into a rectangular box [0, lx] * [0, ly].
  rescale_map(long_grid_side_length, &map_state);
  if (input_polygons_to_eps) {
    std::cout << "Writing input_polygons.eps" << std::endl;
    write_map_to_eps("input_polygons.eps", &map_state);
  }

  // File name without extension
  std::string map_name = geo_file_name.substr(0, geo_file_name.size() - 8);

  // Map after rescaling, wihtout projecting
  std::string precart_name = map_name + "_pre-cartogram.geojson";
  json precart_json = cgal_to_json(&map_state);
  write_to_json(precart_json, geo_file_name, precart_name);

  // Round all points in cartogram
  round_points(&map_state);

  // Map after rounding points
  std::string postround_name = map_name + "_pre-cartogram_post-rounding.geojson";
  json cart_json_1 = cgal_to_json(&map_state);
  write_to_json(cart_json_1, geo_file_name, postround_name);

  std::map<std::string, double> initial_area_errors = map_state.area_err();

  std::map<std::string, std::vector<double>> debug_area_error;

  // Start map integration
  while (map_state.n_finished_integrations() < 20 && // max_integrations
         map_state.max_area_err() > max_permitted_area_error) {

    std::cout << std::endl
              << "Integration number "
              << map_state.n_finished_integrations()
              << std::endl;

    // Rounding Points
    round_points(&map_state);

    // Run set_zero_target_area if running for first time
    if (map_state.n_finished_integrations() == 0) {
      map_state.set_zero_target_area();
    }

    // Filling density
    fill_with_density(&map_state, map_name);

    // Blurring map
    if (map_state.n_finished_integrations() == 0) {
      blur_density(5.0, &map_state);
    } else {
      double blur_width =
      (std::pow(2.0, 3 - int(map_state.n_finished_integrations())) > 1e-1) ? std::pow(2.0, 3 - int(map_state.n_finished_integrations())) : 0.0;
      //std::cout << "Blur condition: " << std::pow(2.0, 3 - (int)(map_state.n_finished_integrations())) << "with" << map_state.n_finished_integrations() << "\n";
      std::cout << "Blur width: " << blur_width << "\n";
      blur_density(blur_width, &map_state);
    }

    // Flattening density
    flatten_density(&map_state);

    //point_search(&map_state, 197.485, 197.486, 257.620, 257.621);

    // Densify
    // Densify map
    map_state.set_geo_divs(densify(map_state.geo_divs()));

    // Choosing diaganols that are inside graticule cells
    choose_diag_4(&map_state);

    // Projecting with Triangulation
    project_with_triangulation(&map_state);

    // Rounding map again
    round_points(&map_state);

    // Printing cartogram
    std::string cartogram_file_name =
      map_name +
      "_cartogram_" +
      std::to_string(map_state.n_finished_integrations()) +
      ".geojson";
    json cart_json = cgal_to_json(&map_state);
    write_to_json(cart_json, geo_file_name, cartogram_file_name);

    // Update debug_area_error

    std::map<std::string, double> current_area_errors = map_state.area_err();

    for (auto gd : map_state.geo_divs()){
      debug_area_error[gd.id()].push_back(current_area_errors[gd.id()]);
    }

    map_state.inc_integration();

  }

  std::cout << std::endl
            << "Final cartogram: cartogram_"
            << map_state.n_finished_integrations()
            << ".geojson"
            << std::endl << std::endl;

  // Export population debug map

  std::map<std::string, std::vector<double>> &population_debug =
    *map_state.debug_population();

  FILE *debug_pop = fopen("population_debug.csv", "w");

  if (debug_pop == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }

  std::cout << "writing debug file...\n";

  fprintf(debug_pop, "Region,Population");
  for (unsigned int i = 0; i < map_state.n_finished_integrations(); i++){
    fprintf(debug_pop, ",Integration %d", i);
  }

  for (auto gd : map_state.geo_divs()){
    fprintf(debug_pop, "\n%s", gd.id().c_str());
    fprintf(debug_pop, ",%f", map_state.target_areas_at(gd.id()));
    for (unsigned int i = 0; i < map_state.n_finished_integrations(); i++){
      fprintf(debug_pop, ",%f", population_debug[gd.id()][i]);
    }
  }

  fclose(debug_pop);



  // std::map<std::string, std::vector<double>> &debug_area_error =
  //   *map_state.debug_area_error();

  FILE *debug_area_error_file = fopen("area_error_debug.csv", "w");

  if (debug_area_error_file == NULL)
  {
    printf("Error opening file!\n");
    exit(1);
  }

  std::cout << "writing debug file...\n";

  fprintf(debug_area_error_file, "Region,Initial Area Errors");
  for (unsigned int i = 0; i < map_state.n_finished_integrations(); i++){
    fprintf(debug_area_error_file, ",Integration %d", i);
  }

  for (auto gd : map_state.geo_divs()){
    fprintf(debug_area_error_file, "\n%s", gd.id().c_str());
    fprintf(debug_area_error_file, ",%f", initial_area_errors[gd.id()]);
    for (unsigned int i = 0; i < map_state.n_finished_integrations(); i++){
      fprintf(debug_area_error_file, ",%f", debug_area_error[gd.id()][i]);
    }
  }

  fclose(debug_area_error_file);


/*
  fill_with_density(&map_state);
    if (map_state.n_finished_integrations() == 0) {
      blur_density(5.0, &map_state);
    } else{
      blur_density(0.0, &map_state);
    }
    flatten_density(&map_state);
    //project(&map_state);
    //map_state.set_geo_divs(densify(map_state.geo_divs()));
    // project_graticule_centroids(&map_state);
    // project_with_triangulation(&map_state);
    choose_diag_2(&map_state);
    project_with_triangulation(&map_state);
    map_state.inc_integration();
*/
/*
  json cart_json = cgal_to_json(&map_state);
  write_to_json(cart_json, geo_file_name, "cartogram.geojson");
*/
  return EXIT_SUCCESS;
}
