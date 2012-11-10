#pragma once

colormap*
mediancut(
	std::vector<hist_item>& hist,
	const double min_opaque_val,
	uint newcolors,
	const double target_mse
);

