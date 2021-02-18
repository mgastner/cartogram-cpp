#include "area_error.h"
#include <algorithm>

double max_area_err(MapState *map_state){

  // Formula for relative area error:
  // area_on_cartogram / target_area - 1

  double sum_target_area = 0.0;
  double sum_cart_area = 0.0;

  for (auto gd : map_state->geo_divs()){
    sum_target_area += map_state->target_areas_at(gd.id());
    sum_cart_area += gd.area();
  }

  double mae = 0.0;

  for (auto gd : map_state->geo_divs()){
    double obj_area =
      map_state->target_areas_at(gd.id()) * sum_cart_area / sum_target_area;
    double relative_area_error = gd.area() / obj_area - 1;
    if (relative_area_error < 0){
      mae = std::max(mae, -relative_area_error);
    }
    else{
      mae = std::max(mae, relative_area_error);
    }
  }

  std::cout << "max. area err: " << mae << "\n";

  return mae;
}
