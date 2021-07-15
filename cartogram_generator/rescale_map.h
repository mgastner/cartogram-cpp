#ifndef RESCALE_MAP_H_
#define RESCALE_MAP_H_

void rescale_map(unsigned int, InsetState*, bool);
void normalize_inset_area(InsetState *inset_state, double, bool = false);
void shift_insets_to_target_position(CartogramInfo*);

#endif
