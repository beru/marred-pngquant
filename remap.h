#pragma once

#include <vector>
#include "pam.h"

struct read_info;
struct write_info;

int best_color_index(
	const std::vector<colormap_item>& map,
	f_pixel px,
	double* dist_out
	);

void viter_init(
	const std::vector<colormap_item>& map,
	f_pixel* average_color, double* average_color_count, f_pixel* base_color, double* base_color_count
	);

void viter_update_color(
	f_pixel acolor, double value, std::vector<colormap_item>& map, int match,
	f_pixel* average_color, double* average_color_count,
	const f_pixel* base_color, const double* base_color_count
	);

void viter_do_interation(
	const std::vector<hist_item>& hist,
	std::vector<colormap_item>& map,
	double min_opaque_val
	);

void viter_finalize(
	std::vector<colormap_item>& map,
	const f_pixel* average_color, const double* average_color_count
	);

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

