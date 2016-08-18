local cam = require("cam")

local od = 90;
local id = od - 4;

f = 3000;

local center = {}
center.x = 0;
center.y = 0;

-- TODO need to plunge to depth + ensure moves are on clearance plane

cam.polygon(id/2, 128, nil, center, f);
cam.polygon(od/2, 128, nil, center, f);
