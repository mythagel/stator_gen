./stator_gen profile | nc_svgpath -f 200 | nc_contour_pocket --tool_r 0.5 --cut_z -1.5 -d 0.5 | (./stator_gen drill && cat) | (turtlecam ../toolpaths/stator_outside.lua && cat) | (turtlecam ../toolpaths/spacer.lua && cat) | nc_backplot
