#include <cairo/cairo.h>
#include <cairo/cairo-ps.h>
#include <cairo/cairo-pdf.h>
#include "constants.h"
#include "cartogram_info.h"
#include "inset_state.h"
#include "iostream"

void write_polygon_to_cairo_surface(cairo_t *cr,
                            int height,
                            bool fill_polygons,
                            bool colors,
                            bool plot_graticule,
                            InsetState *inset_state)
{
    unsigned int line_width = std::min(inset_state->lx(), inset_state->ly());
    cairo_set_line_width(cr, 0.001 * line_width);
    // draw the shapes
    for (auto gd : inset_state->geo_divs()) {
        for (auto pwh : gd.polygons_with_holes()) {
            Polygon ext_ring = pwh.outer_boundary();
            // Move to starting coordinates
            cairo_move_to(cr, ext_ring[0][0], height - ext_ring[0][1]);
            // Plot each point in exterior ring
            unsigned int n = ext_ring.size();
            for (unsigned int i = 1; i < n; ++i) {
                cairo_line_to(cr, ext_ring[i][0], height - ext_ring[i][1]);
            }

            // Plot holes
            for (auto hci = pwh.holes_begin(); hci != pwh.holes_end(); ++hci) {
                Polygon hole = *hci;
                cairo_move_to(cr, hole[0][0], height - hole[0][1]);
                unsigned int hole_size = hole.size();
                for (unsigned int i = 1; i <= hole_size; ++i) {
                    cairo_line_to(cr, hole[i % hole_size][0], height - hole[i % hole_size][1]);
                }
            }

            if (colors || fill_polygons) {
                if (inset_state->is_input_target_area_missing(gd.id())) {
                    // Fill path with dark-grey
                    cairo_set_source_rgb(cr, 0.9375, 0.9375, 0.9375);
                } else if (colors) {
                    // Get color
                    Color col = inset_state->colors_at(gd.id());
                    double red = col.r / 255.0;
                    double green = col.g / 255.0;
                    double blue = col.b / 255.0;

                    // Fill path
                    cairo_set_source_rgb(cr, red, green, blue);
                } else if (fill_polygons) {
                    // Fill path with default color
                    cairo_set_source_rgb(cr, 0.96, 0.92, 0.70);
                }
            }
            cairo_fill(cr);
        }
    }

    //add the labels
    for (auto gd : inset_state->geo_divs()) {
        // go to a specific coordinate to place the label
        cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 14.0);
        Point label_coordinate = gd.point_on_surface();
        cairo_move_to(cr, label_coordinate.x(), height - label_coordinate.y());
        //get the label
        std::string label = inset_state->labels_at(gd.id());
        const char* label_char = label.c_str();
        cairo_show_text(cr, label_char);
    }

    //plot the graticule
    if (plot_graticule) {
        boost::multi_array<XYPoint, 2> &cum_proj =
        *inset_state->ref_to_cum_proj();
        unsigned int graticule_line_spacing = 7;

        // Set line width of graticule lines
        cairo_set_line_width(cr, 0.0005 * line_width);
        cairo_set_source_rgb(cr, 0, 0, 0);

        // Vertical graticule lines
        for (unsigned int i = 0; i <= inset_state->lx(); i += graticule_line_spacing) {
            cairo_move_to(cr, cum_proj[i][0].x, cum_proj[i][0].y);
            for (unsigned int j = 1; j < inset_state->ly(); ++j) {
                cairo_line_to(cr, cum_proj[i][j].x, cum_proj[i][j].y);
            }
            cairo_stroke(cr);
        }

        // Horizontal graticule lines
        for (unsigned int j = 0; j <= inset_state->ly(); j += graticule_line_spacing) {
            cairo_move_to(cr, cum_proj[0][j].x, cum_proj[0][j].y);
            for (unsigned int i = 1; i < inset_state->lx(); ++i) {
                cairo_line_to(cr, cum_proj[i][j].x, cum_proj[i][j].y);
            }
            cairo_stroke(cr);
        }
    }
}
void write_cairo_polygons_to_png(std::string fname,
                            bool fill_polygons,
                            bool colors,
                            bool plot_graticule,
                            InsetState *inset_state)
{
    const char* filename = fname.c_str();
    unsigned int width = inset_state->lx();
    unsigned int height = inset_state->ly();

    cairo_surface_t *surface;
    cairo_t *cr;

    surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, width, height);
    cr = cairo_create(surface);

    write_polygon_to_cairo_surface(cr, height, fill_polygons, colors, plot_graticule, inset_state);

    cairo_surface_write_to_png(surface, filename);
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
}

void write_cairo_polygons_to_ps(std::string fname,
                           bool fill_polygons,
                           bool colors,
                           bool plot_graticule,
                           InsetState *inset_state)
{
    const char* filename = fname.c_str();
    cairo_surface_t *surface;
    cairo_t *cr;
    unsigned int width = inset_state->lx();
    unsigned int height = inset_state->ly();

    surface = cairo_ps_surface_create(filename, width, height);
    cr = cairo_create(surface);

    //add comments
    std::string title = "%%Title: " + fname;
    cairo_ps_surface_dsc_comment(surface, title.c_str());
    cairo_ps_surface_dsc_comment(surface, "%%Creator: Michael T. Gastner et al.");
    cairo_ps_surface_dsc_comment(surface, "%%For: Humanity");
    cairo_ps_surface_dsc_comment(surface, "%%Copyright: License CC BY-NC-ND 2.0");
    cairo_ps_surface_dsc_comment(surface, "%%Magnification: 1.0000");

    write_polygon_to_cairo_surface(cr, height, fill_polygons, colors, plot_graticule, inset_state);

    cairo_show_page(cr);

    cairo_surface_destroy(surface);
    cairo_destroy(cr);
}

// This function outputs both png and ps files
void write_cairo_map(std::string filename, bool plot_graticule, InsetState* inset_state) {
    // solid: clockwise, holes: anti-clockwise
    std::string png_filename = filename + ".png";
    std::string ps_filename = filename + ".ps";

    //Check whether the has all GeoDivs colored
    bool has_colors = (inset_state->colors_size() == inset_state->n_geo_divs());
    write_cairo_polygons_to_png(png_filename,
                        true,
                        has_colors,
                        plot_graticule,
                        inset_state);
   write_cairo_polygons_to_ps(ps_filename,
                       true,
                       has_colors,
                       plot_graticule,
                       inset_state);
}
