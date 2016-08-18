local cam = require("cam")

local od = 90;

f = 50;

local center = {}
center.x = od/2;
center.y = od/2;

-- TODO need to plunge to depth + ensure moves are on clearance plane

move_to(center.x, center.y, nil);
cam.polygon(od/2, 128, nil, center, f);
