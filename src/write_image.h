#ifndef WRITE_IMAGE_H
#define WRITE_IMAGE_H

void write_map_image(const std::string, const bool, const bool, InsetState *);
void write_graticule_heatmap_image(const std::string, const bool, InsetState *);
void write_density_bar_image(const std::string);
void write_density_image(const std::string, const double *, InsetState *, const bool = false);

#endif
