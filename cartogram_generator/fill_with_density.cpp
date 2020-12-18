#include "map_state.h"
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Bbox_2.h>

// Struct to store intersection data
struct intersection {
	double x; // x co-ordinate of intersection
	double target_density; // geo-div's target_density
	bool direction; // whether intersection enters or exits

	// from https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
	// Overloading operator for comparing this data type
	bool operator < (const intersection& rhs) const {
		if (x < rhs.x) {
			return true;
		}
		else if (x == rhs.x && direction < rhs.direction) {
			return true;
		}
		return x < rhs.x;
	}

};

void fill_with_density(MapState* map_state) {
	double total_current_area = 0.0;
	for (auto gd: map_state->geo_divs()) {
		total_current_area += gd.area();
	}
	double total_target_area = 0.0;
	for (auto gd: map_state->geo_divs()) {
		total_target_area += map_state->target_areas_at(gd.id());
	}

	double mean_density = total_target_area / total_current_area;
	FTReal2d & rho_init = *map_state->ref_to_rho_init();

	// Assigning all densities to 0 initially
	for (unsigned int i = 0; i < map_state->lx(); i++) {
		for (unsigned int j = 0; j < map_state->ly(); j++) {
			rho_init(i, j) = 0;
		}
	}

	// Resolution to set density, must be power of 2
	double res = 16;

	// An array (map_intersections) to store vectors of intersections
	std::vector<intersection> map_intersections[(int) (map_state->ly() * res)];

	// Iterating through geodivs in map_state
	for (auto gd: map_state->geo_divs()) {
		std::cout << "Working on gd with ID: " << gd.id() << std::endl;

		// Associative area only called once to find out target density
		double target_density = map_state->target_areas_at(gd.id()) / gd.area();
		target_density /= res;

		// Iterating through polygons with holes in map_state
		for (int j = 0; j < gd.n_polygons_with_holes(); j++) {
			std::cout << "Polygon " << j << " in GeoDiv" << std::endl;
			Polygon_with_holes pwh = gd.polygons_with_holes()[j];
			CGAL::Bbox_2 bb = pwh.bbox();

			// Cylcling through y co-ordinates in bounding box of pwh
			for (unsigned int k = floor(bb.ymin()) - 1; k <= ceil(bb.ymax()) + 1;
			     k++) {

				// Cycling through each of res number of lines in one cell
				for (double line_y = k + (1.0/res)/2; line_y < k + 1;
				     line_y += (1.0/res)) {
					Polygon ext_ring = pwh.outer_boundary();
					double prev_point[2];
					prev_point[0] = ext_ring[0][0];
					prev_point[1] = ext_ring[0][1];
					// Temporary vector of intersections for that particular line
					std::vector <intersection> intersections;

					// The following algorithm works by iterating through "res" lines in
					// each cell. For each line, it iterates through every edge in a
					// polygon and stores any intersections. Finally, once all
					// intersections have been stored, we iterate between intersections,
					// and add the appropriate densities.

					// Running algorithm on exterior ring
					for (unsigned int l = 1; l <= ext_ring.size(); l++) {
						double curr_point[2];
						curr_point[0] = ext_ring[(l)%ext_ring.size()][0];
						curr_point[1] = ext_ring[(l)%ext_ring.size()][1];
						// Checking if intersection present
						if ((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
						    (curr_point[1] >= line_y && prev_point[1] <= line_y)) {

							if (curr_point[1] == prev_point[1]) {
								continue; // Pre-condition to ignore grazing incident
							}
							// Adding an epsilon in case of intersection with point
							else if (curr_point[1] == line_y) {
								curr_point[1] += 0.00001 * (1/res);
							}
							else if (prev_point[1] == line_y) {
								prev_point[1] += 0.00001 * (1/res);
							}
							// Creating an intersection and storing it in a vector
							intersection temp;
							temp.x = (curr_point[0] * (prev_point[1] - line_y) +
							          prev_point[0] * (line_y - curr_point[1])) /
							         (prev_point[1] - curr_point[1]);
							temp.target_density = target_density;
							temp.direction = false; // temporary value
							intersections.push_back(temp);
						}
						prev_point[0] = curr_point[0];
						prev_point[1] = curr_point[1];
					}

					// Running algorithm on each hole
					for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
						Polygon hole = *hci;
						prev_point[0] = hole[0][0];
						prev_point[1] = hole[0][1];

						for (unsigned int l = 1; l <= hole.size(); l++) {
							double curr_point[2];
							curr_point[0] = hole[(l)%hole.size()][0];
							curr_point[1] = hole[(l)%hole.size()][1];
							if ((curr_point[1] <= line_y && prev_point[1] >= line_y) ||
							    (curr_point[1] >= line_y && prev_point[1] <= line_y)) {

								if (curr_point[1] == prev_point[1]) {
									continue;
								}
								else if (curr_point[1] == line_y) {
									std::cout << "curr_point[1] == line_y" << std::endl;
									curr_point[1] += 0.00001 * (1/res);
								}
								else if (prev_point[1] == line_y) {
									std::cout << "prev_point[1] == line_y" << std::endl;
									prev_point[1] += 0.00001 * (1/res);
								}
								intersection temp;
								temp.x = (curr_point[0] * (prev_point[1] - line_y) +
								          prev_point[0] * (line_y - curr_point[1])) /
								         (prev_point[1] - curr_point[1]);
								temp.target_density = target_density;
								temp.direction = false; // temporary value
								intersections.push_back(temp);
							}
							prev_point[0] = curr_point[0];
							prev_point[1] = curr_point[1];
						}
					}

					// Checking if odd number of intersections
					if (intersections.size() % 2 != 0) {
						std::cout << "Incorrect Topology" << std::endl;
						std::cout << "Number of intersections: " << intersections.size();
						std::cout << std::endl;
						std::cout << "Y-coordinate: " << line_y << std::endl;
						std::cout << "Intersection points: " << std::endl;
						for (unsigned int l = 0; l < intersections.size(); l++) {
							std::cout << intersections[l].x << std::endl;
						}
					}
					else {
						std::sort(intersections.begin(), intersections.end());
						// Adding sorted vector of intersections to array map_intersections
						for (unsigned int l = 0; l < intersections.size(); l++) {
							intersections[l].direction = (l%2 == 0);
							map_intersections[(int) ((line_y - (1/res)/2) * res)].
							push_back(intersections[l]);
						}
					}

				}
			}
		}
	}

	// Cycling through y co-ordinates in map_state
	for (unsigned int k = 0; k < map_state->ly(); k++) {
		// Cycling through each of res number of lines in one cell
		for (double line_y = k + (1.0/res)/2; line_y < k + 1; line_y += (1.0/res)) {

			std::vector<intersection> intersections =
				map_intersections[(int) ((line_y - (1/res)/2) * res)];

			// Sorting vector in ascending order of intersection
			sort(intersections.begin(), intersections.end());

			// Filling lines with no intersections with mean_density
			if (intersections.size() == 0) {
				for (unsigned int l = 0; l < map_state->lx(); l++) {
					rho_init(l, k) += mean_density/res;
				}
			}
			else {
				// Filling from first coordinate upto first GeoDiv
				for (unsigned int l = 1; l <= ceil(intersections[0].x); l++) {
					if (l == ceil(intersections[0].x)) {
						rho_init(l - 1, k) += ( (mean_density/res) *
						                        (intersections[0].x - floor(intersections[0].x)) );
					}
					else {
						rho_init(l - 1, k) += mean_density/res;
					}
				}
				// Filling any empty spaces between GeoDivs
				for (unsigned int l = 1; l < intersections.size() - 1; l += 2) {
					double left_x = intersections[l].x;
					double right_x = intersections[l + 1].x;
					for (unsigned int m = ceil(left_x); m <= ceil(right_x); m++) {
						if (left_x == right_x) {
							// Pre-condition, to ensure different intersecting points
							continue;
						}

						if (m == ceil(left_x)) {
							rho_init(m - 1, k) += ( (mean_density/res) *
							                        (ceil(left_x) - left_x) );
						}
						else if (m == ceil(right_x)) {
							rho_init(m - 1, k) += ( (mean_density/res) *
							                        (right_x - floor(right_x)) );
						}
						else {
							rho_init(m - 1, k) += (mean_density/res);
						}
					}
				}
				// Filling from last GeoDiv upto last coordinate
				for (unsigned int l = ceil(intersections.back().x); l <= map_state->lx(); l++) {
					if (l == ceil(intersections.back().x)) {
						rho_init(l - 1, k) += ( (mean_density/res) *
						                        (ceil(intersections.back().x) - intersections.back().x) );
					}
					else {
						rho_init(l - 1, k) += mean_density/res;
					}

				}
			}

			// Filling GeoDivs by iterating through intersections
			for (unsigned int l = 0; l < intersections.size(); l += 2) {
				double left_x = intersections[l].x;
				double right_x = intersections[l + 1].x;
				// Checking for Intersection of Polygons/Holes/Geodivs
				if (intersections[l].direction == intersections[l + 1].direction) {
					std::cout << "Invalid Geometry!" << std::endl;
					std::cout << "Intersection of Polygons/Holes/Geodivs" << std::endl;
					std::cout << "Y-coordinate: " << k << std::endl;
					std::cout << "Left X-coordinate: " << left_x << std::endl;
					std::cout << "Right X-coordinate: " << right_x << std::endl;
					continue; // to highlight where intersection is present
				}
				// Filling each cell between intersections
				for (unsigned int m = ceil(left_x); m <= ceil(right_x); m++) {
					if (m == ceil(left_x)) {
						rho_init(m - 1, k) += ( intersections[l].target_density *
						                        (ceil(left_x) - left_x) );
					}
					else if (m == ceil(right_x)) {
						rho_init(m - 1, k) += ( intersections[l].target_density *
						                        (right_x - floor(right_x)) );
					}
					else {
						rho_init(m - 1, k) += intersections[l].target_density;
					}
				}
			}

		}
	}

	map_state->execute_fwd_plan();
	return;

}
