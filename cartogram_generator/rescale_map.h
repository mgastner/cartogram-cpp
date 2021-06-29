#ifndef RESCALE_MAP_H_
#define RESCALE_MAP_H_

void rescale_map(unsigned int, InsetState*, bool);
void rescale_to_each_other(InsetState *inset_state);
void unscale_map(InsetState*);
void rescale_to_position(InsetState *inset_state, std::string pos);

#endif
