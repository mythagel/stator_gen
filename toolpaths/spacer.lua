local cam = require("cam")

local od = 50;
local id = od - 4;
local tool_r = 1/2

f = 300;

local center = {}
center.x = od/2;
center.y = od/2;

-- TODO need to plunge to depth + ensure moves are on clearance plane

move_to(center.x, center.y, nil);
cam.polygon((id/2)-tool_r, 128, nil, center, f);
move_to(center.x, center.y, nil);
cam.polygon((od/2)+tool_r, 128, nil, center, f);
