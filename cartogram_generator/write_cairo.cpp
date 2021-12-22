#include <cairo/cairo.h>
#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"

void write_cairo_map_to_eps(std::string eps_name, bool plot_graticule, InsetState* inset_state) {
    cairo_surface_t *surface;
    cairo_t *cr;

    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 390, 60);
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 0, 0, 0);
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 40.0);

    cairo_move_to(cr, 10.0, 50.0);
    cairo_show_text(cr, "Drawing using Cairo");

    cairo_surface_write_to_png(surface, "image.png");

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}
