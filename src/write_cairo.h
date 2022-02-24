#ifndef WRITE_CAIRO_H
#define WRITE_CAIRO_H

void write_map_to_ps(const std::string, const bool, const bool, InsetState*);
void write_graticule_heatmap_to_ps(const std::string, const bool, InsetState*);
void write_density_to_ps(const std::string, const double*, InsetState*);


#endif
