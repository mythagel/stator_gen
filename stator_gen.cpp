// Nicholas Gill (c) 2016
#include <random>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <string>


#define JC_VORONOI_IMPLEMENTATION
#include "jc_voronoi.h"

inline std::string r6(double v) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << v;
    auto s = ss.str();
    
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if(s.back() == '.') s.pop_back();
    return s;
}


struct point {
    double x;
    double y;
};
struct polar_point {
    double r;
    double t;
};

double distance(const point& p0, const point& p1) {
    auto x = p1.x - p0.x;
    auto y = p1.x - p0.y;
    return std::abs(std::sqrt(x*x + y*y));
}

static const double PI = 3.14159265359;

point lerp(const point& p0, const point& p1, float t) {
    return { (1-t)*p0.x + t*p1.x, (1-t)*p0.y + t*p1.y };
}

double f = 50;
point pos;
point bezier_cp;
void bezier_curve_to(bool abs, float x1, float y1, float x, float y) {
    point p[3];
    p[0] = pos;
    p[1] = {abs ? x1 : pos.x + x1, abs ? y1 : pos.y + y1};
    p[2] = {abs ? x : pos.x + x, abs ? y : pos.y + y};
    pos = p[2];
    bezier_cp = p[1];

    for(float t = 0.0; t < 1.0; t += 0.03) {
        auto ab = lerp(p[0], p[1], t);
        auto bc = lerp(p[1], p[2], t);
        auto point = lerp(ab, bc, t);

        std::cout << "G01 X" << r6(point.x) << " Y" << r6(point.y) << " F" << r6(f) << "\n";
    }
    std::cout << "G01 X" << r6(pos.x) << " Y" << r6(pos.y) << " F" << r6(f) << "\n";
}

int main() {
    std::vector<polar_point> points;
    double r = 60;  // mm

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> theta(0, 2*PI);
    std::uniform_real_distribution<> radius(0, r);
    for (int n = 0; n < 100; ++n) {
        auto t = theta(gen);
        auto r = radius(gen);

        points.push_back({r, t});
    }


    if (false) {
        std::sort(begin(points), end(points), [](const polar_point& p0, const polar_point& p1) -> bool {
            return std::tie(p0.t, p0.r) < std::tie(p1.t, p1.r);
        });

        for (auto& p : points) {
            auto x = p.r * std::cos(p.t);
            auto y = p.r * std::sin(p.t);
            std::cout << "G83 X" << std::fixed << x << " Y" << y << " Z-1 R1 Q0.5 F50" << '\n';
        }
        std::cout << '\n';
    } else {
        std::vector<jcv_point> jcv_points;

        if (false) {
            for (auto& p : points) {
                auto x = p.r * std::cos(p.t) + 70;
                auto y = p.r * std::sin(p.t) + 70;
                jcv_points.push_back({x, y});
            }
        } else {

            for (unsigned i = 0; i < 18; ++i) {
                int pointoffset = 0;
                jcv_points.emplace_back();
                jcv_points.back().x = (float)(pointoffset + rand() % (70-2*pointoffset));
                jcv_points.back().y = (float)(pointoffset + rand() % (70-2*pointoffset));
            }
        }

        jcv_diagram diagram;
        std::memset(&diagram, 0, sizeof diagram);
        jcv_diagram_generate(jcv_points.size(), jcv_points.data(), 80, 80, &diagram);

        auto sites = jcv_diagram_get_sites(&diagram);
        for (unsigned i = 0; i < diagram.numsites; ++i) {
            auto site = &sites[i];
            for (auto e = site->edges; e; e = e->next) {
                auto p = point{site->p.x, site->p.y};
                auto m0 = lerp(p, {e->pos[0].x, e->pos[0].y}, 0.92);
                auto m1 = lerp(p, {e->pos[1].x, e->pos[1].y}, 0.92);
                e->pos[0].x = m0.x;
                e->pos[0].y = m0.y;
                e->pos[1].x = m1.x;
                e->pos[1].y = m1.y;
            }
        }
        for (unsigned i = 0; i < diagram.numsites; ++i) {
            auto site = &sites[i];

            auto first = lerp({site->edges->pos[0].x, site->edges->pos[0].y}, {site->edges->pos[1].x, site->edges->pos[1].y}, 0.5);
            std::cout << "G0" << " X" << std::fixed << first.x << " Y" << first.y << "\n";
            pos = first;

            for (auto e = site->edges; e; e = e->next) {
                auto m0 = lerp({e->pos[0].x, e->pos[0].y}, {e->pos[1].x, e->pos[1].y}, 0.5);
                auto next = e->next ? e->next : site->edges;
                auto m1 = lerp({next->pos[0].x, next->pos[0].y}, {next->pos[1].x, next->pos[1].y}, 0.5);
                bezier_curve_to(true, e->pos[1].x, e->pos[1].y, m1.x, m1.y);
            }
        }

    }
}