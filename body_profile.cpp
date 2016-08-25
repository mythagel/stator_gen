#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>

#include "clipper.hpp"
namespace cl = ClipperLib;

double PI = 3.14159265359;

struct point {
    double x;
    double y;
};

inline std::string r6(double v) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(6) << v;
    auto s = ss.str();
    
    s.erase(s.find_last_not_of('0') + 1, std::string::npos);
    if(s.back() == '.') s.pop_back();
    return s;
}

int main() {
    double scale = 10e12;
    auto scale_point = [&](double x, double y) -> cl::IntPoint {
        return cl::IntPoint(x * scale, y * scale);
    };
    auto unscale = [&](const cl::IntPoint& p) -> point {
        return {static_cast<double>(p.X) / scale, static_cast<double>(p.Y) / scale};
    };

    unsigned n_segments = 8;
    double delta_theta = (2*PI) / n_segments;

    cl::Paths paths;

    double r0 = (2.5*25.4)/2;
    double r1 = 100/2;
    double theta = 0;
    for (unsigned segment = 0; segment < n_segments; ++segment, theta += delta_theta) {
        if (segment == 0)
            continue;

        paths.emplace_back();
        auto& path = paths.back();

        double t0 = theta - delta_theta/2;
        double t1 = theta + delta_theta/2;

        double x0 = r0 * std::cos(t0);
        double y0 = r0 * std::sin(t0);
        path.push_back(scale_point(x0, y0));

        for (double t = t0; t < t1; t += 0.03) {
            double x1 = r1 * std::cos(t);
            double y1 = r1 * std::sin(t);
            path.push_back(scale_point(x1, y1));
        }

        double x1 = r1 * std::cos(t1);
        double y1 = r1 * std::sin(t1);
        path.push_back(scale_point(x1, y1));

        double x2 = r0 * std::cos(t1);
        double y2 = r0 * std::sin(t1);
        path.push_back(scale_point(x2, y2));

        for (double t = t1; t > t0; t -= 0.03) {
            double x3 = r0 * std::cos(t);
            double y3 = r0 * std::sin(t);
            path.push_back(scale_point(x3, y3));
        }
    }

    double offset = -4.0;
    cl::ClipperOffset co;
    co.AddPaths(paths, cl::jtRound, cl::etClosedPolygon);
    co.ArcTolerance = 0.1 * scale;

    cl::Paths solution;
    co.Execute(solution, offset * scale);
    cl::CleanPolygons(solution);

    for(auto& path : solution) {

        auto first = unscale(*path.begin());
        std::cout << "M " << r6(first.x) << " " << r6(first.y) << " ";
        for(auto& point : path) {
            auto p = unscale(point);
            std::cout << "L " << r6(p.x) << " " << r6(p.y) << " ";
        }
        std::cout << "L " << r6(first.x) << " " << r6(first.y) << " ";
    }
}
