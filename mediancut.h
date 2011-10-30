#pragma once

#include <vector>

std::vector<colormap_item> mediancut(std::vector<hist_item>& hist, double min_opaque_val, int reqcolors);
