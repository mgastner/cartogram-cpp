#ifndef RESCALE_MAP_H_
#define RESCALE_MAP_H_

void rescale_map(unsigned int, InsetState*, bool);
void normalize_inset_area(InsetState *inset_state, double);
void unscale_map(InsetState*);
void shift_inset_to_target_position(InsetState *inset_state,
                                    std::string pos, 
                                    std::map <std::string, CGAL::Bbox_2>);

#endif
