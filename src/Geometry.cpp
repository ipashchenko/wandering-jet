#include "Geometry.h"
#include "utils.h"
#include "MyExceptions.h"
#include <math.h>
#include <iostream>
#include <list>

using Eigen::Vector3d;


std::list<double>
intersection(Vector3d R0, Vector3d Rd, double A, double B, double C, double D, double E, double F, double G, double H,
             double I, double J) {
    std::list<double> result;

    double x0 = R0[0];
    double y0 = R0[1];
    double z0 = R0[2];

    double xd = Rd[0];
    double yd = Rd[1];
    double zd = Rd[2];

    double Aq = A*xd*xd + B*yd*yd + C*zd*zd + D*xd*yd + E*xd*zd + F*yd*zd;

    double Bq = 2*A*x0*xd + 2*B*y0*yd + 2*C*z0*zd + D*(x0*yd + y0*xd) + E*x0*zd + F*(y0*zd + yd*z0) + G*xd + H*yd +
            I*zd;
    double Cq = A*x0*x0 + B*y0*y0 + C*z0*z0 + D*x0*y0 + E*x0*z0 + F*y0*z0 + G*x0 + H*y0 + I*z0 + J;
    double Dscr = Bq*Bq - 4*Aq*Cq;

    if (Aq == 0) {
        result = std::list<double>{-Cq/Bq};
    }

    else if (Dscr < 0) {
        result = std::list<double>{};
    }

    else {
        double t0 = (-Bq - sqrt(Dscr))/(2*Aq);

        if (t0 > 0) {
            result = std::list<double>{t0};
        }
        else {
            double t1 = (-Bq + sqrt(Dscr))/(2*Aq);
            result = std::list<double>{t0, t1};
        }
    }

    return result;
}

std::pair<Vector3d, Vector3d> Geometry::full_infinite_path(Ray &ray) const {
    Vector3d diff = ray.origin() - origin();
    Vector3d ref_point = diff - diff.dot(ray.direction()) * ray.direction();
    Vector3d point_in = ref_point - big_scale() * ray.direction();
    Vector3d point_out = ref_point + big_scale() * ray.direction();
    return std::pair<Vector3d,Vector3d>{point_in, point_out};
}

std::pair<Vector3d, Vector3d> Geometry::half_infinite_path(Ray &ray, const Vector3d &point) const {
    Vector3d check_point = point-ray.direction();
    if (is_within(check_point)) {
        return std::pair<Vector3d,Vector3d>{point-big_scale()*ray.direction(), point};
    }
    else {
        return std::pair<Vector3d,Vector3d>{point, point+big_scale()*ray.direction()};
    }
}
