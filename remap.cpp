#include "remap.h"

#include <stddef.h>
#include "png.h"	/* libpng header; includes zlib.h */
#include "rwpng.h"
#include "dxor.h"

int best_color_index(
	const std::vector<colormap_item>& map,
	f_pixel px,
	double* dist_out
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
		int match = best_color_index(map, hi.acolor, NULL);
		viter_update_color(hi.acolor, hi.perceptual_weight, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
	}

	viter_finalize(map, &average_color[0], &average_color_count[0]);
}

double remap_to_palette(
	const f_pixel* input, size_t width, size_t height,
	write_info* output_image,
	std::vector<colormap_item>& map, double min_opaque_val
	)
{
	unsigned char** row_pointers = output_image->row_pointers;
	
	int remapped_pixels = 0;
	double remapping_error = 0;
	
	int transparent_ind = best_color_index(map, f_pixel(0,0,0,0), NULL);

	std::vector<f_pixel> average_color(map.size());
	std::vector<double> average_color_count(map.size());
	viter_init(map, &average_color[0], &average_color_count[0], NULL, NULL);
	
	const f_pixel* pLine = input;
	for (size_t y=0; y<height; ++y) {
		unsigned char* output_row = row_pointers[y];
		for (size_t x=0; x<width; ++x) {
			f_pixel px = pLine[x];
			int match;
			if (px.alpha < 1.0/256.0) {
				match = transparent_ind;
			}else {
				double diff;
				match = best_color_index(map, px, &diff);
				++remapped_pixels;
				remapping_error += diff;
			}
			output_row[x] = match;
			viter_update_color(px, 1.0, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
		}
		pLine += width;
	}
	viter_finalize(map, &average_color[0], &average_color_count[0]);
	return remapping_error / max(1,remapped_pixels);
}

/**
  Uses edge/noise map to apply dithering only to flat areas. Dithering on edges creates jagged lines, and noisy areas are "naturally" dithered.
 */
void remap_to_palette_floyd(
	const f_pixel* input, size_t width, size_t height,
	write_info* output_image,
	const std::vector<colormap_item>& map, double min_opaque_val,
	const double* edge_map
	)
{
	unsigned char** row_pointers = output_image->row_pointers;
	
	int ind = 0;
	int transparent_ind = best_color_index(map, f_pixel(0,0,0,0), NULL);

	f_pixel* thiserr = NULL;
	f_pixel* nexterr = NULL;
	int fs_direction = 1;

	/* Initialize Floyd-Steinberg error vectors. */
	thiserr = (f_pixel*) malloc((width + 2) * sizeof(*thiserr));
	nexterr = (f_pixel*) malloc((width + 2) * sizeof(*thiserr));
	sdxor156(12345); /* deterministic dithering is better for comparing results */

	const double INVFACTOR = 1.0 / 255.0;
	for (size_t x=0; x<width+2; ++x) {
		thiserr[x].r = (dxor156() - 0.5) * INVFACTOR;
		thiserr[x].g = (dxor156() - 0.5) * INVFACTOR;
		thiserr[x].b = (dxor156() - 0.5) * INVFACTOR;
		thiserr[x].alpha = (dxor156() - 0.5) * INVFACTOR;
	}

	const f_pixel* pInputLine = input;
	for (size_t y=0; y<height; ++y) {
		memset(nexterr, 0, (width + 2) * sizeof(*nexterr));

		int x = (fs_direction) ? 0 : (width - 1);
		unsigned char* output_row = row_pointers[y];
		const double* pEdge = &edge_map[y*width];
		do {
			f_pixel px = pInputLine[x];
			double dither_level = pEdge[x];
			
			/* Use Floyd-Steinberg errors to adjust actual color. */
			f_pixel tmp = px + thiserr[x + 1] * dither_level;
			tmp.r = limitValue(tmp.r, 0.0, 1.0);
			tmp.g = limitValue(tmp.g, 0.0, 1.0);
			tmp.b = limitValue(tmp.b, 0.0, 1.0);
			tmp.alpha = limitValue(tmp.alpha, 0.0, 1.0);

			if (tmp.alpha < 1.0/256.0) {
				ind = transparent_ind;
			}else {
				ind = best_color_index(map, tmp, 0);
			}

			output_row[x] = ind;

			f_pixel err = tmp - map[ind].acolor;
			
			// If dithering error is crazy high, don't propagate it that much
			// This prevents crazy geen pixels popping out of the blue (or red or black! ;)
			if (err.r*err.r + err.g*err.g + err.b*err.b + err.alpha*err.alpha > 8.0/256.0) {
				dither_level *= 0.5;
			}
			double colorimp = (3.0 + map[ind].acolor.alpha)/4.0 * dither_level;
			err.r *= colorimp;
			err.g *= colorimp;
			err.b *= colorimp;
			err.alpha *= dither_level;
			
#if 1
			// changed kernel after reading the paper : Reinstating FloydSteinberg: Improved Metrics for Quality Assessment of Error Diffusion Algorithms (Sam Hocevar, Gary Niger)
			/* Propagate Floyd-Steinberg error terms. */
			if (fs_direction) {
				thiserr[x + 2] += err * 7.0 / 16.0;
				nexterr[x + 0] += err * 4.0 / 16.0;
				nexterr[x + 1] += err * 5.0 / 16.0;
				nexterr[x + 2] += err * 0.0 / 16.0;
			}else {
				thiserr[x + 0] += err * 7.0 / 16.0;
				nexterr[x + 0] += err * 0.0 / 16.0;
				nexterr[x + 1] += err * 5.0 / 16.0;
				nexterr[x + 2] += err * 4.0 / 16.0;
			}
#else
			/* Propagate Floyd-Steinberg error terms. */
			if (fs_direction) {
				thiserr[x + 2] += err * 7.0 / 16.0;
				nexterr[x	   ] += err * 3.0 / 16.0;
				nexterr[x + 1] += err * 5.0 / 16.0;
				nexterr[x + 2] += err		   / 16.0;
			}else {
				thiserr[x	   ] += err * 7.0 / 16.0;
				nexterr[x	   ] += err		   / 16.0;
				nexterr[x + 1] += err * 5.0 / 16.0;
				nexterr[x + 2] += err * 3.0 / 16.0;
			}
#endif

			if (fs_direction) {
				++x;
				if (x >= width) break;
			}else {
				--x;
				if (x < 0) break;
			}
		}while (1);
		std::swap(thiserr, nexterr);
		fs_direction = !fs_direction;

		pInputLine += width;
	}
}

