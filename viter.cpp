#include "pam.h"
#include "viter.h"
#include "nearest.h"
#include <stdlib.h>
#include <string.h>

/*
 * Voronoi iteration: new palette color is computed from weighted average of colors that map to that palette entry.
 */
void viter_init(
	const colormap* map,
	viter_state* average_color
	)
{
	memset(average_color, 0, sizeof(average_color[0])*map->colors);
}

void viter_update_color(
	const f_pixel acolor,
	const double value,
	const colormap* map,
	int match,
	viter_state state[]
	)
{
	assert(match < map->colors);
    state[match].a += acolor.alpha * value;
    state[match].r += acolor.r * value;
    state[match].g += acolor.g * value;
    state[match].b += acolor.b * value;
    state[match].total += value;
}

void viter_finalize(
	colormap* map,
	const viter_state average_color[]
	)
{
	for (size_t i=0; i<map->colors; i++) {
        double a=0, r=0, g=0, b=0, total=0;
		colormap_item& pal = map->palette[i];
		const viter_state& s = average_color[i];
		a += s.a;
		r += s.r;
		g += s.g;
		b += s.b;
		total += s.total;
		if (total) {
			pal.acolor = f_pixel(a,r,g,b) / total;
		}
		
		pal.popularity = total;
		
	}
}

double viter_do_iteration(
	histogram* hist,
	colormap* const map,
	double min_opaque_val,
	viter_callback callback
	)
{
	std::vector<viter_state> average_color(map->colors);
//	viter_init(map, &average_color[0]);
	nearest_map* const n = nearest_init(map);
	hist_item* const achv = hist->achv;
	const int hist_size = hist->size;

	double total_diff=0;
	for (size_t j=0; j<hist->size; ++j) {
		double diff;
		unsigned int match = nearest_search(n, achv[j].acolor, min_opaque_val, &diff);
		total_diff += diff * achv[j].perceptual_weight;

		viter_update_color(achv[j].acolor, achv[j].perceptual_weight, map, match, &average_color[0]);

		if (callback) callback(&achv[j], diff);
	}
	nearest_free(n);
	viter_finalize(map, &average_color[0]);

	return total_diff / hist->total_perceptual_weight;

}

