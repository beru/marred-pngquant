#include "remap.h"

#include <stddef.h>
#include "png.h"	/* libpng header; includes zlib.h */
#include "rwpng.h"

int best_color_index(
	f_pixel px,
	const std::vector<colormap_item>& map,
	double min_opaque_val, double* dist_out
	)
{
	int ind = 0;
	double dist = colordifference(px, map[0].acolor);

	for (size_t i=1; i<map.size(); ++i) {
		double newdist = colordifference(px, map[i].acolor);
		if (newdist < dist) {
			ind = i;
			dist = newdist;
		}
	}

	if (dist_out) *dist_out = dist;
	return ind;
}

/*
 * Voronoi iteration: new palette color is computed from weighted average of colors that map to that palette entry.
 */
void viter_init(
	const std::vector<colormap_item>& map,
	f_pixel* average_color, double* average_color_count,
	f_pixel* base_color, double* base_color_count
	)
{
	for (size_t i=0; i<map.size(); i++) {
		average_color_count[i] = 0;
		average_color[i] = f_pixel(0,0,0,0);
	}

	// Rather than only using separate mapping and averaging steps
	// new palette colors are computed at the same time as mapping is done
	// but to avoid first few matches moving the entry too much
	// some base color and weight is added
	if (base_color) {
		for (size_t i=0; i<map.size(); i++) {
			const colormap_item& nmi = map[i];
			double value = 1.0 + nmi.popularity/2.0;
			base_color_count[i] = value;
			base_color[i] = nmi.acolor * value;
		}
	}
}

void viter_update_color(
	f_pixel acolor, double value, std::vector<colormap_item>& map, int match,
	f_pixel* average_color, double* average_color_count,
	const f_pixel* base_color, const double* base_color_count
	)
{
	f_pixel& ac = average_color[match];
	double& acc = average_color_count[match];
	const f_pixel& bc = base_color[match];
	const double& bcc = base_color_count[match];
	ac += acolor * value;
	acc += value;

	if (base_color) {
		map[match].acolor = (ac + bc) / (acc + bcc);
	}
}

void viter_finalize(
	std::vector<colormap_item>& map,
	const f_pixel* average_color, const double* average_color_count
	)
{
	for (size_t i=0; i<map.size(); i++) {
		colormap_item& pal = map[i];
		const double acc = average_color_count[i];
		if (acc) {
			pal.acolor = average_color[i] / acc;
		}
		pal.popularity = acc;
	}
}

void viter_do_interation(
	const std::vector<hist_item>& hist,
	std::vector<colormap_item>& map,
	double min_opaque_val
	)
{
	std::vector<f_pixel> average_color(map.size());
	std::vector<double> average_color_count(map.size());

	viter_init(map, &average_color[0], &average_color_count[0], NULL,NULL);
	for (size_t j=0; j<hist.size(); j++) {
		const hist_item& hi = hist[j];
		int match = best_color_index(hi.acolor, map, min_opaque_val, NULL);
		viter_update_color(hi.acolor, hi.perceptual_weight, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
	}

	viter_finalize(map, &average_color[0], &average_color_count[0]);
}

double remap_to_palette(
	const read_info* input_image, write_info* output_image,
	std::vector<colormap_item>& map, double min_opaque_val
	)
{
	const rgb_pixel** input_pixels = (const rgb_pixel**) input_image->row_pointers;
	unsigned char** row_pointers = output_image->row_pointers;
	int rows = input_image->height;
	int cols = input_image->width;
	double gamma = input_image->gamma;

	int remapped_pixels=0;
	double remapping_error=0;
	
	int transparent_ind = best_color_index(f_pixel(0,0,0,0), map, min_opaque_val, NULL);

	std::vector<f_pixel> average_color(map.size());
	std::vector<double> average_color_count(map.size());
	viter_init(map, &average_color[0], &average_color_count[0], NULL, NULL);

	for (int row=0; row<rows; ++row) {
		const rgb_pixel* input_row = input_pixels[row];
		unsigned char* output_row = row_pointers[row];
		for (int col=0; col<cols; ++col) {

			f_pixel px = to_f(gamma, input_row[col]);
			int match;

			if (px.a < 1.0/256.0) {
				match = transparent_ind;
			}else {
				double diff;
				match = best_color_index(px, map,min_opaque_val, &diff);

				remapped_pixels++;
				remapping_error += diff;
			}

			output_row[col] = match;

			viter_update_color(px, 1.0, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
		}
	}

	viter_finalize(map, &average_color[0], &average_color_count[0]);

	return remapping_error / MAX(1,remapped_pixels);
}

template <typename T>
T limitValue(T val, T min, T max)
{
	if (val < min) return min;
	if (max < val) return max;
	return val;
}

/**
  Uses edge/noise map to apply dithering only to flat areas. Dithering on edges creates jagged lines, and noisy areas are "naturally" dithered.
 */
void remap_to_palette_floyd(
	const read_info* input_image, write_info* output_image,
	const std::vector<colormap_item>& map, double min_opaque_val,
	const double* edge_map
	)
{
	const rgb_pixel** input_pixels = (const rgb_pixel**) input_image->row_pointers;
	unsigned char** row_pointers = output_image->row_pointers;
	int rows = input_image->height;
	int cols = input_image->width;
	double gamma = input_image->gamma;
	
	int ind = 0;
	int transparent_ind = best_color_index(f_pixel(0,0,0,0), map, min_opaque_val, NULL);

	f_pixel* thiserr = NULL;
	f_pixel* nexterr = NULL;
	double sr=0, sg=0, sb=0, sa=0;
	int fs_direction = 1;

	/* Initialize Floyd-Steinberg error vectors. */
	thiserr = (f_pixel*) malloc((cols + 2) * sizeof(*thiserr));
	nexterr = (f_pixel*) malloc((cols + 2) * sizeof(*thiserr));
	srand(12345); /* deterministic dithering is better for comparing results */

	for (int col=0; col<cols+2; ++col) {
		const double rand_max = RAND_MAX;
		thiserr[col].r = ((double)rand() - rand_max/2.0)/rand_max/255.0;
		thiserr[col].g = ((double)rand() - rand_max/2.0)/rand_max/255.0;
		thiserr[col].b = ((double)rand() - rand_max/2.0)/rand_max/255.0;
		thiserr[col].a = ((double)rand() - rand_max/2.0)/rand_max/255.0;
	}

	for (int row=0; row<rows; ++row) {
		memset(nexterr, 0, (cols + 2) * sizeof(*nexterr));

		int col = (fs_direction) ? 0 : (cols - 1);

		const rgb_pixel* input_row = input_pixels[row];
		unsigned char* output_row = row_pointers[row];
		do {
			f_pixel px = to_f(gamma, input_row[col]);
            double dither_level = edge_map ? edge_map[row*cols + col] : 0.9;
			
			/* Use Floyd-Steinberg errors to adjust actual color. */
			sr = px.r + thiserr[col + 1].r * dither_level;
			sg = px.g + thiserr[col + 1].g * dither_level;
			sb = px.b + thiserr[col + 1].b * dither_level;
			sa = px.a + thiserr[col + 1].a * dither_level;

			sr = limitValue(sr, 0.0, 1.0);
			sg = limitValue(sg, 0.0, 1.0);
			sb = limitValue(sb, 0.0, 1.0);
			sa = limitValue(sa, 0.0, 1.0);

			if (sa < 1.0/256.0) {
				ind = transparent_ind;
			}else {
				ind = best_color_index(f_pixel(sa,sr,sg,sb), map, min_opaque_val, 0);
			}

			output_row[col] = ind;

			f_pixel xp = map[ind].acolor;
			f_pixel err;
			err.r = (sr - xp.r);
			err.g = (sg - xp.g);
			err.b = (sb - xp.b);
			err.a = (sa - xp.a);
			
			// If dithering error is crazy high, don't propagate it that much
			// This prevents crazy geen pixels popping out of the blue (or red or black! ;)
			if (err.r*err.r + err.g*err.g + err.b*err.b + err.a*err.a > 8.0/256.0) {
				dither_level *= 0.5;
            }
			double colorimp = (3.0 + map[ind].acolor.a)/4.0 * dither_level;
			err.r *= colorimp;
			err.g *= colorimp;
			err.b *= colorimp;
			err.a *= dither_level;
			
			/* Propagate Floyd-Steinberg error terms. */
			if (fs_direction) {
				thiserr[col + 2] += err * 7.0f / 16.0;
				nexterr[col	   ] += err * 3.0f / 16.0;
				nexterr[col + 1] += err * 5.0f / 16.0;
				nexterr[col + 2] += err        / 16.0;
			}else {
				thiserr[col	   ] += err * 7.0f / 16.0;
				nexterr[col	   ] += err	       / 16.0;
				nexterr[col + 1] += err * 5.0f / 16.0;
				nexterr[col + 2] += err * 3.0f / 16.0;
			}

			if (fs_direction) {
				++col;
				if (col >= cols) break;
			}else {
				--col;
				if (col < 0) break;
			}
		}while (1);
		std::swap(thiserr, nexterr);
		fs_direction = !fs_direction;
	}
}

