./stator_gen profile | nc_contour_pocket --tool_r 0.5 | (./stator_gen drill && cat) | (turtlecam ../toolpaths/stator_outside.lua && cat) | (turtlecam ../toolpaths/spacer.lua && cat) | nc_backplot
