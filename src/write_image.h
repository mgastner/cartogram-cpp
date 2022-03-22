#ifndef WRITE_IMAGE_H
#define WRITE_IMAGE_H

void write_map_image(const std::string, const bool, const bool, const bool, InsetState *);
void write_graticule_heatmap_image(const std::string, const bool, const bool, InsetState *);
void write_density_bar_image(const std::string, const bool);
void write_density_image(const std::string, const double *,
                         const bool, const bool, InsetState *);

#endif
