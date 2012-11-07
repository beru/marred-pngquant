
#include "pam.h"
#include "nearest.h"
#include "mempool.h"
#include <stdlib.h>

struct color_entry {
	f_pixel color;
	double radius;
	unsigned int index;
};

struct sorttmp {
	double radius;
	unsigned int index;
};

struct head {
	f_pixel center;
	double radius;
	unsigned int num_candidates;
	color_entry* candidates;
};

struct nearest_map {
    head* heads;
    mempool* mempool;
    unsigned int num_heads;
};

static
int compareradius(const void* ap, const void* bp)
{
	double a = ((const sorttmp*)ap)->radius;
	double b = ((const sorttmp*)bp)->radius;
	return a > b ? 1 : (a < b ? -1 : 0);
}

head build_head(
	f_pixel px,
	const colormap* map,
	int num_candidates,
	mempool* m,
	unsigned int skip_index[],
	unsigned int *skipped
	)
{
	std::vector<sorttmp> colors(map->colors);
	unsigned int colorsused=0;

	for (unsigned int i=0; i<map->colors; i++) {
		if (skip_index[i]) continue;
		colors[colorsused].index = i;
		colors[colorsused].radius = colordifference(px, map->palette[i].acolor);
		colorsused++;
	}

	qsort(&colors[0], colorsused, sizeof(colors[0]), compareradius);
	assert(colorsused < 2 || colors[0].radius <= colors[1].radius);

	num_candidates = min<int>(colorsused, num_candidates);

    head h;
    h.candidates = (color_entry*) mempool_new(m, num_candidates * sizeof(h.candidates[0]));
    h.center = px;
	h.num_candidates = num_candidates;
	for (unsigned int i=0; i<num_candidates; ++i) {
		color_entry entry;
		entry.color = map->palette[colors[i].index].acolor;
		entry.index = colors[i].index;
		entry.radius = colors[i].radius;
		h.candidates[i] = entry;
	}
	h.radius = colors[num_candidates-1].radius/4.0; // /2 squared

	for (unsigned int i=0; i<num_candidates; ++i) {
		assert(colors[i].radius <= h.radius*4.0);
        // divide again as that's matching certain subset within radius-limited subset
        // - 1/256 is a tolerance for miscalculation (seems like colordifference isn't exact)
        if (colors[i].radius < h.radius/4.0 - 1.0/256.0) {
			skip_index[colors[i].index]=1;
            (*skipped)++;
		}
	}
	return h;
}

static
colormap* get_subset_palette(const colormap* map)
{
    // it may happen that it gets palette without subset palette or the subset is too large
    int subset_size = (map->colors+3)/4;

    if (map->subset_palette && map->subset_palette->colors <= subset_size) {
        return map->subset_palette;
    }

    const colormap *source = map->subset_palette ? map->subset_palette : map;
    colormap *subset_palette = pam_colormap(subset_size);

    for(unsigned int i=0; i < subset_size; i++) {
        subset_palette->palette[i] = source->palette[i];
    }

    return subset_palette;
}


nearest_map* nearest_init(const colormap* map)
{
    mempool* m = NULL;
    nearest_map* centroids = (nearest_map*) mempool_new(m, sizeof(*centroids));
    centroids->mempool = m;

    unsigned int skipped=0;
	std::vector<unsigned int> skip_index(map->colors);

    colormap* subset_palette = get_subset_palette(map);
    const int selected_heads = subset_palette->colors;
    centroids->heads = (head*) mempool_new(centroids->mempool, sizeof(centroids->heads[0])*(selected_heads+1)); // +1 is fallback head

    unsigned int h=0;
    for(; h < selected_heads; h++)
    {
        unsigned int num_candiadtes = 1+(map->colors - skipped)/((1+selected_heads-h)/2);

        centroids->heads[h] = build_head(subset_palette->palette[h].acolor, map, num_candiadtes, centroids->mempool, &skip_index[0], &skipped);
        if (centroids->heads[h].num_candidates == 0) {
            break;
        }
    }

    centroids->heads[h].radius = MAX_DIFF;
	f_pixel px; px.alpha = px.r = px.g = px.b = 0.0;
    centroids->heads[h].center = px;
    centroids->heads[h].num_candidates = 0;
    centroids->heads[h].candidates = (color_entry*) mempool_new(centroids->mempool, (map->colors - skipped) * sizeof(centroids->heads[h].candidates[0]));
    for (unsigned int i=0; i<map->colors; i++) {
        if (skip_index[i]) continue;
		color_entry entry;
		entry.color = map->palette[i].acolor;
		entry.index = i;
		entry.radius = 999;
        centroids->heads[h].candidates[centroids->heads[h].num_candidates++] = entry;
    }
    centroids->num_heads = ++h;

    // get_subset_palette could have created a copy
    if (subset_palette != map->subset_palette) {
        pam_freecolormap(subset_palette);
    }

    return centroids;
}

unsigned int nearest_search(const struct nearest_map *centroids, const f_pixel px, const double min_opaque_val, double* diff)
{
    const int iebug = px.alpha > min_opaque_val;

    const struct head *const heads = centroids->heads;
    for (unsigned int i=0; i<centroids->num_heads; i++) {
        double headdist = colordifference(px, heads[i].center);

        if (headdist <= heads[i].radius) {
            assert(heads[i].num_candidates);
            unsigned int ind=heads[i].candidates[0].index;
            double dist = colordifference(px, heads[i].candidates[0].color);

            /* penalty for making holes in IE */
            if (iebug && heads[i].candidates[0].color.alpha < 1.0) {
                dist += 1.0/1024.0;
            }

            for (unsigned int j=1; j<heads[i].num_candidates; j++) {
                double newdist = colordifference(px, heads[i].candidates[j].color);

                /* penalty for making holes in IE */
                if (iebug && heads[i].candidates[j].color.alpha < 1.0) {
                    newdist += 1.0/1024.0;
                }

                if (newdist < dist) {

                    dist = newdist;
                    ind = heads[i].candidates[j].index;
                }
            }
            if (diff) *diff = dist;
            return ind;
        }
    }
    assert(0);
    return 0;
}

void nearest_free(struct nearest_map *centroids)
{
    mempool_free(centroids->mempool);
}
