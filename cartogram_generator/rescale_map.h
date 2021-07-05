#ifndef RESCALE_MAP_H_
#define RESCALE_MAP_H_

void rescale_map(unsigned int, InsetState*, bool);
void normalize_inset_area(InsetState *inset_state, double);
void unscale_map(InsetState*);
void shift_insets_to_target_position(CartogramInfo*);

#endif
