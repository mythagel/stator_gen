
//$fn = 64;

module spacer() {
    od = 90;
    id = od - 4;
    h = 1;
    difference() {
        cylinder(r=od/2, h=h);
        translate([0,0,-0.5]) cylinder(r=id/2, h=h+1);
    }
}

module diaphragm() {
    od = 90;
    % cylinder(r=od/2, h=0.003);
}

module stator() {
    od = 90;
    id = od - 10;
    h = 2;
    
    difference() {
        cylinder(r=od/2, h=h);
        for(x = [1 : 15 : 90]) {
            for(y = [1 : 15 : 90]) {
                translate([x-45, y-45, -0.5]) cylinder(r=10/2, h = h+1);
            }
        }
    }
    
    difference() {
        cylinder(r=od/2, h=h);
        translate([0,0,-0.5]) cylinder(r=id/2, h=h+1);
    }
}

module stator_assy() {
    translate([0,0,2.1+1.1+1.1]) stator();
    translate([0,0,2.1+1.1]) spacer();
    % translate([0,0,2.1+1.001]) diaphragm();
    translate([0,0,2.1]) spacer();
    stator();
}

module body() {
    od = 100;
    stator_id = 90;
    id = stator_id - 10;
    h = 10;
    difference() {
        // main body
        intersection() {
            translate([0,0,h/2]) scale([1,1,0.5]) sphere(r=od/2);
            cylinder(r=od/2, h=h);
        }
        
        // through hole
        translate([0,0,-0.5]) cylinder(r=id/2, h=h+1);
        
        // stator shelf
        translate([0,0,h-7.1]) cylinder(r=stator_id/2, h=6+2);
        
        // cover shelf
        translate([0,0,h-1]) cylinder(r=(stator_id+6)/2, h=1+1);
    }
}

translate([0,0,10-7]) stator_assy();
body();