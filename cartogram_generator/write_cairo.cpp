#include <cairo/cairo.h>
#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"

void write_cairo_polygons_to_eps(const char *filename,
                           bool fill_polygons,
                           bool colors,
                           bool plot_graticule,
                           InsetState *inset_state)
{
    cairo_surface_t *surface;
    cairo_t *cr;

    surface = cairo_ps_surface_create(filename, 800, 800);
    cr = cairo_create(surface);

    cairo_set_source_rgb(cr, 0, 0, 0);
    unsigned int line_width = std::min(inset_state->lx(), inset_state->ly());
    cairo_set_line_width(cr, 0.001 * line_width);

    for (auto gd : inset_state->geo_divs()) {
        for (auto pwh : gd.polygons_with_holes()) {
            Polygon ext_ring = pwh.outer_boundary();
            // Move to starting coordinates
            cairo_move_to(cr, ext_ring[0][0], ext_ring[0][1]);
            // Plot each point in exterior ring
            unsigned int n = ext_ring.size();
            for (unsigned int i = 1; i <= n; ++i) {
                cairo_line_to(cr, ext_ring[i % n][0], ext_ring[i % n][1]);
            }

            // Plot holes
            for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
                Polygon hole = *hci;
                cairo_move_to(cr, hole[0][0], hole[0][1]);
                unsigned int hole_size = hole.size();
                for (unsigned int i = 1; i <= hole_size; ++i) {
                    cairo_line_to(cr, hole[i % hole_size][0], hole[i % hole_size][1]);
                }
            }
        }
    }
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

void write_cairo_map_to_eps(std::string eps_name, bool plot_graticule, InsetState* inset_state) {
    // solid: clockwise, holes: anti-clockwise
    const char* fname = eps_name.c_str();
    cairo_surface_t *surface;
    cairo_t *cr;

    // surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 390, 60);
    // cr = cairo_create(surface);

    // cairo_set_source_rgb(cr, 0, 0, 0);
    // cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
    //                        CAIRO_FONT_WEIGHT_NORMAL);
    // cairo_set_font_size(cr, 40.0);

    // cairo_move_to(cr, 10.0, 50.0);
    // cairo_show_text(cr, "Drawing using Cairo");

    // cairo_surface_write_to_png(surface, "image.png");

    // cairo_destroy(cr);
    // cairo_surface_destroy(surface);

    Check whether the has all GeoDivs colored
    bool has_colors = (inset_state->colors_size() == inset_state->n_geo_divs());
    write_cairo_polygons_to_eps(fname,
                        true,
                        has_colors,
                        plot_graticule,
                        inset_state);
}
