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
#define JCV_REAL_TYPE double
#define JCV_FABS fabs
#define JCV_ATAN2 atan2
#include "jc_voronoi.h"

#include "clipper.hpp"
namespace cl = ClipperLib;


bool output_profile = false;
bool output_drill = true;

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
    auto y = p1.y - p0.y;
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
    double r = 90/2;  // mm

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
                jcv_points.back().x = (float)(pointoffset + rand() % ((int)(r*2)-2*pointoffset));
                jcv_points.back().y = (float)(pointoffset + rand() % ((int)(r*2)-2*pointoffset));
            }
        }

        jcv_diagram diagram;
        std::memset(&diagram, 0, sizeof diagram);
        jcv_diagram_generate(jcv_points.size(), jcv_points.data(), r*2, r*2, &diagram);

        auto sites = jcv_diagram_get_sites(&diagram);

        bool reduce = false;
        if (reduce) {
            for (int i = 0; i < diagram.numsites; ++i) {
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
        }

        bool bezier = false;
        if (bezier) {
            for (int i = 0; i < diagram.numsites; ++i) {
                auto site = &sites[i];

                auto first = lerp({site->edges->pos[0].x, site->edges->pos[0].y}, {site->edges->pos[1].x, site->edges->pos[1].y}, 0.5);
                std::cout << "G0" << " X" << std::fixed << first.x << " Y" << first.y << "\n";
                pos = first;

                for (auto e = site->edges; e; e = e->next) {
                    auto next = e->next ? e->next : site->edges;
                    auto m1 = lerp({next->pos[0].x, next->pos[0].y}, {next->pos[1].x, next->pos[1].y}, 0.5);
                    bezier_curve_to(true, e->pos[1].x, e->pos[1].y, m1.x, m1.y);
                }
            }
        } else {

            double scale = 10e12;
            auto scale_point = [&](const jcv_point& p) -> cl::IntPoint {
                return cl::IntPoint(p.x * scale, p.y * scale);
            };
            auto unscale = [&](const cl::IntPoint& p) -> point {
                return {static_cast<double>(p.X) / scale, static_cast<double>(p.Y) / scale};
            };

            cl::Path circle_clip;
            for (double t = 0; t < 2*PI; t += 0.03) {
                auto x = r * std::cos(t);
                auto y = r * std::sin(t);
                circle_clip.push_back(scale_point({x + r, y + r}));
            }

            double offset = -1.5;

            if (output_profile) {
                for (int i = 0; i < diagram.numsites; ++i) {
                    auto site = &sites[i];

                    /* create clipper path for each site, offset inwards, then output */

                    cl::Paths paths;
                    paths.emplace_back();
                    auto first = site->edges->pos[0];
                    paths.back().push_back(scale_point(first));

                    for (auto e = site->edges; e; e = e->next) {
                        paths.back().push_back(scale_point(e->pos[1]));
                    }

                    cl::Clipper clipper;
                    clipper.AddPaths(paths, cl::ptSubject, true);
                    clipper.AddPath(circle_clip, cl::ptClip, true);
                    paths.clear();
                    clipper.Execute(cl::ctIntersection, paths);

                    cl::ClipperOffset co;
                    co.AddPaths(paths, cl::jtRound, cl::etClosedPolygon);
                    co.ArcTolerance = 0.1 * scale;
                    cl::Paths solution;
                    co.Execute(solution, offset * scale);
                    cl::CleanPolygons(solution);

                    for(auto& path : solution) {

                        std::cout << std::fixed << "G00 X" << static_cast<double>(path.begin()->X)/scale << " Y" << static_cast<double>(path.begin()->Y)/scale << "\n";
                        for(auto& p : path) {
                            std::cout << std::fixed << "G01 X" << static_cast<double>(p.X)/scale << " Y" << static_cast<double>(p.Y)/scale << " F50\n";
                        }
                        std::cout << std::fixed << "G01 X" << static_cast<double>(path.begin()->X)/scale << " Y" << static_cast<double>(path.begin()->Y)/scale << "\n";
                        std::cout << "\n";
                    }
                }
            }
            if (output_drill) {
                for (int i = 0; i < diagram.numsites; ++i) {
                    auto site = &sites[i];

                    double drill_diameter = 1;

                    cl::Paths clip;
                    clip.emplace_back();
                    auto first = site->edges->pos[0];
                    clip.back().push_back(scale_point(first));
                    for (auto e = site->edges; e; e = e->next) {
                        clip.back().push_back(scale_point(e->pos[1]));
                    }

                    cl::Clipper clipper;
                    clipper.AddPaths(clip, cl::ptSubject, true);
                    clipper.AddPath(circle_clip, cl::ptClip, true);
                    clip.clear();
                    clipper.Execute(cl::ctIntersection, clip);

                    if (clip.empty())
                        continue;

                    cl::ClipperOffset co;
                    co.AddPaths(clip, cl::jtRound, cl::etClosedPolygon);
                    co.ArcTolerance = 0.1 * scale;
                    clip.clear();
                    co.Execute(clip, offset * scale);
                    cl::CleanPolygons(clip);

                    // random depth - number of inscribed circles for maximum voronoi cell dimensions
                    double drill_offset = drill_diameter;
                    for(unsigned depth = 0; depth < 100; ++depth) {

                        cl::PolyTree pt;
                        {
                            cl::Paths paths;
                            paths.emplace_back();

                            double dt = 0;
                            unsigned n_holes = (2*PI * (drill_offset/2)) / drill_diameter;
                            auto theta = 2*PI / n_holes;
                            if (depth % 2 == 0) {
                                dt = theta/2;
                            }
                            auto t = dt;
                            for (unsigned i = 0; i < n_holes; ++i, t += theta) {
                                auto x = (drill_offset/2) * std::cos(t);
                                auto y = (drill_offset/2) * std::sin(t);
                                paths.back().push_back(scale_point({x + site->p.x, y + site->p.y}));
                            }

                            cl::Clipper clipper;
                            clipper.AddPaths(paths, cl::ptSubject, false);
                            clipper.AddPaths(clip, cl::ptClip, true);
                            clipper.Execute(cl::ctIntersection, pt);
                        }

                        for(auto node = pt.GetFirst(); node; node = node->GetNext()) {
                            auto& path = node->Contour;
                            for(auto& point : path) {
                                auto p = unscale(point);
                                std::cout << "G83 X" << std::fixed << p.x << " Y" << p.y << " Z-1 R1 Q0.5 F50" << '\n';
                            }
                            std::cout << "\n";
                        }

                        drill_offset += drill_diameter;
                    }
                }
            }
        }

    }
}
