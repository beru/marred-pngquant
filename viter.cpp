#include "pam.h"
#include "viter.h"

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

double viter_do_interation(
	const std::vector<hist_item>& hist,
	std::vector<colormap_item>& map,
	double min_opaque_val
	)
{
	std::vector<f_pixel> average_color(map.size());
	std::vector<double> average_color_count(map.size());

	viter_init(map, &average_color[0], &average_color_count[0], NULL,NULL);
	
    double total_diff=0, total_weight=0;
	for (size_t j=0; j<hist.size(); j++) {
		const hist_item& hi = hist[j];
		double diff;
		int match = best_color_index(map, hi.acolor, &diff);
        total_diff += diff * hi.perceptual_weight;
        total_weight += hi.perceptual_weight;
		
		viter_update_color(hi.acolor, hi.perceptual_weight, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
	}
	
	viter_finalize(map, &average_color[0], &average_color_count[0]);
	
	return total_diff / total_weight;
}

