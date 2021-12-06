#ifndef DENSIFY_H_
#define DENSIFY_H_

#include <list>
#include "cgal_typedef.h"
#include "geo_div.h"

std::vector<GeoDiv> densified_geo_divs(std::vector<GeoDiv>,
                                       const unsigned int,
                                       const unsigned int);

#endif
