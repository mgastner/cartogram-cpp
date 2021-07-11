#include "cgal_typedef.h"
#include "write_to_json.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "csv.hpp"
#include <iostream>
#include <fstream>

nlohmann::json cgal_to_json(CartogramInfo *cart_info)
{
  nlohmann::json container;
  for (auto &inset_state : *cart_info->ref_to_inset_states()) {
    for (auto gd : inset_state.geo_divs()) {
      nlohmann::json gd_container;
      for (auto pwh : gd.polygons_with_holes()) {

        // Get exterior ring of Polygon_with_holes
        Polygon ext_ring = pwh.outer_boundary();
        nlohmann::json polygon_container;
        nlohmann::json er_container;
        for (unsigned int i = 0; i < ext_ring.size(); i++) {

          // Get exterior ring coordinates
          double arr[2];
          arr[0] = ext_ring[i][0];
          arr[1] = ext_ring[i][1];
          er_container.push_back(arr);
        }
        polygon_container.push_back(er_container);

        // Get holes of Polygon_with_holes
        for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); hci++) {
          Polygon hole = *hci;
          nlohmann::json hole_container;
          for (unsigned int i = 0; i < hole.size(); i++) {

            // Get hole coordinates
            double arr[2];
            arr[0] = hole[i][0];
            arr[1] = hole[i][1];
            hole_container.push_back(arr);
          }
          polygon_container.push_back(hole_container);
        }
        gd_container.push_back(polygon_container);
      }
      container.push_back(gd_container);
    }
  }
  return container;
}

int cartogram_id_from_csv(const boost::program_options::variables_map vm,
                                  nlohmann::json geo_div_properties) {
  // Get name of CSV file from vm
  std::string csv_name;
  if (vm.count("visual_variable_file")) {
    csv_name = vm["visual_variable_file"].as<std::string>();
  } else {
    std::cerr << "ERROR: No CSV file given!" << std::endl;
    _Exit(15);
  }

  // Do we need to run try/catch checks like we did in read_csv()?

  // Opening CSV Reader
  csv::CSVReader reader(csv_name);

  // Finding index of column with IDs
  std::string id_header;
  int id_col = 0;
  if (vm.count("id")) {
    id_header = vm["id"].as<std::string>();
    id_col = reader.index_of(id_header);
  } else {
    id_header = reader.get_col_names()[0];
  }

  // Iterate through CSV rows
  int row_num = 1; // To be used as cartogram_id
  for (auto row : reader) {
    std::string csv_id = row[id_col].get();
    std::string json_id = geo_div_properties[id_header];

    // Return cartogram_id
    if (csv_id == json_id) {
      return row_num;
    }

    row_num++;
  }
}

void write_to_json(nlohmann::json container,
                   std::string old_geo_fn,
                   std::string new_geo_fn,
                   std::ostream &new_geo_stream,
                   bool output_to_stdout,
                   const boost::program_options::variables_map vm)
{
  std::ifstream i(old_geo_fn);
  nlohmann::json old_j;
  i >> old_j;
  nlohmann::json newJ;

  // Loop over multipolygons in the container
  for (int i = 0; i < (int) container.size(); i++) {
    newJ["features"][i]["properties"] = old_j["features"][i]["properties"];
    newJ["features"][i]["id"] = old_j["features"][i]["id"];
    newJ["features"][i]["type"] = "Feature";
    newJ["features"][i]["geometry"]["type"] = "MultiPolygon";

    // loop over Polygon_with_holes in the multipolygon
    for (int j = 0; j < (int) container[i].size(); j++) {

      // Loop over exterior ring and holes in the Polygon_with_holes
      for (int k = 0; k < (int) container[i][j].size(); k++) {
        newJ["features"][i]["geometry"]["coordinates"][j][k] =
          container[i][j][k];
      }
    }

    // Get cartogram_id for current GeoDiv based on CSV row numbers
    if (!newJ["features"][i]["properties"]["cartogram_id"]) {
      newJ["features"][i]["properties"]["cartogram_id"] =
        cartogram_id_from_csv(vm, newJ["features"][i]["properties"]);
    }
  }
  newJ.push_back({"type", old_j["type"]});
  if (output_to_stdout) {
    new_geo_stream << newJ << std::endl;
  } else {
    std::ofstream o(new_geo_fn);
    o << newJ << std::endl;
  }
  return;
}
