
#include "pam.h"
#include "nearest.h"

struct color_entry {
	f_pixel color;
	double radius;
	int index;
};

struct sorttmp {
	double radius;
	int index;
};

struct head {
	f_pixel center;
	double radius;
	int num_candidates;
	color_entry candidates[256];
};

static int compareradius(const void* ap, const void* bp)
{
	double a = ((sorttmp*)ap)->radius;
	double b = ((sorttmp*)bp)->radius;
	return a > b ? 1 : (a < b ? -1 : 0);
}

head heads[256];
int num_heads;
int skipped = 0;

struct head build_head(f_pixel px, const colormap* map, int num_candidates, int skip_index[])
{
	sorttmp colors[map->colors];
	int colorsused=0;

	for (int i=0; i<map->colors; i++) {
		if (skip_index[i]) continue;
		colors[colorsused].index = i;
		colors[colorsused].radius = colordifference(px, map->palette[i].acolor);
		colorsused++;
	}

	qsort(&colors, colorsused, sizeof(colors[0]), compareradius);
	assert(colorsused < 2 || colors[0].radius <= colors[1].radius);

	head h;
	h.center = px;
	num_candidates = MIN(colorsused, num_candidates);
	h.num_candidates = num_candidates;
	for (int i=0; i<num_candidates; i++) {
		h.candidates[i] = (color_entry) {
			.color = map->palette[colors[i].index].acolor,
			.index = colors[i].index,
			.radius = colors[i].radius,
		};
	}
	h.radius = colors[num_candidates-1].radius/4.0f; // /2 squared

	for (int i=0; i<num_candidates; i++) {
		assert(colors[i].radius <= h.radius*4.0f);
		if (colors[i].radius < h.radius/4.f) { // divide again as that's matching certain subset within radius-limited subset
			skip_index[colors[i].index]=1;
			skipped++;
		}
	}
	return h;
}


void nearest_init(const colormap* map)
{
	int skip_index[map->colors];
	for (int j=0; j<map->colors;j++) {
		skip_index[j] = 0;
	}
	
	int max_heads = map->subset_palette->colors;
	int h = 0;
	for (; h<max_heads; h++) {
		int num_candiadtes = 1+(map->colors-skipped)/((1+max_heads-h)/2);
		int idx = best_color_index(map->subset_palette->palette[h].acolor, map, 1.0, NULL);
		heads[h] = build_head(map->palette[idx].acolor, map, num_candiadtes, skip_index);
		if (heads[h].num_candidates == 0) {
			break;
		}
	}
	heads[h].radius = 9999999;
	heads[h].center = (f_pixel){0,0,0,0};
	heads[h].num_candidates = 0;
	for (int i=0; i<map->colors; i++) {
		if (skip_index[i]) continue;
		heads[h].candidates[heads[h].num_candidates++] = (color_entry) {
			.color = map->palette[i].acolor,
			.index = i,
			.radius = 999,
		};
	}
	num_heads = ++h;
}

int nearest_search(f_pixel px)
{
	for (int i=0; i<num_heads; i++) {
		double headdist = colordifference(px, heads[i].center);

		if (headdist <= heads[i].radius) {
			assert(heads[i].num_candidates);
			int ind=heads[i].candidates[0].index;
			double dist = colordifference(px, heads[i].candidates[0].color);
			for (int j=1; j<heads[i].num_candidates; j++) {
				double newdist = colordifference(px, heads[i].candidates[j].color);
				if (newdist < dist) {
					dist = newdist;
					ind = heads[i].candidates[j].index;
				}
			}
			return ind;
		}
	}
	assert(0);
	return 0;
}

