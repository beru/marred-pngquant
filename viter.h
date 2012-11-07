#pragma once


struct viter_state {
    double a, r, g, b, total;
};

typedef void (*viter_callback)(hist_item* item, double diff);

void viter_init(
	const colormap* map,
	viter_state state[]
);

void viter_update_color(
	const f_pixel acolor,
	const double value,
	const colormap* map,
	int match,
	viter_state state[]
);

void viter_finalize(
	colormap* map,
	const viter_state state[]
);

double viter_do_iteration(
	histogram* hist,
	colormap* const map,
	const double min_opaque_val,
	viter_callback callback
);


