#pragma once

#include <vector>
#include "pam.h"

struct read_info;
struct write_info;

double remap_to_palette(
	const f_pixel* input, size_t width, size_t height,
	write_info* output_image,
	std::vector<colormap_item>& map, double min_opaque_val
	);
	
void remap_to_palette_floyd(
	const f_pixel* input, size_t width, size_t height,
	write_info* output_image,
	const std::vector<colormap_item>& map, double min_opaque_val,
	const double* edge_map
	);

