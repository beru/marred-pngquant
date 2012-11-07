#pragma once

colormap*
mediancut(
	histogram* hist,
	const double min_opaque_val,
	unsigned int newcolors,
	const double target_mse
);