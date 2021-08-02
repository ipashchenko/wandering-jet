#include <iostream>
//#include "Intersection.h"
#include "utils.h"


Vector3d Intersection::direction() const {
    return direction_;
}

void Intersection::set_direction(const Vector3d& direction) {
    direction_ = direction;
}

void Intersection::set_borders(const Vector3d& point_in, const Vector3d& point_out) {
    borders_ = std::make_pair(point_in, point_out);
}

const std::pair<Vector3d, Vector3d>& Intersection::get_path() const {
    return borders_;
}

// Ctor for finite intersections
Intersection::Intersection(const Ray& ray, const Vector3d& new_point_in, const Vector3d& new_point_out) {
    init(ray, new_point_in, new_point_out);
    is_finite_ = true;
}

Intersection::Intersection(const Ray& ray, const std::pair<Vector3d, Vector3d>& borders) {
    set_direction(ray.direction());
    set_borders(borders.first, borders.second);
    is_finite_ = true;
}

// Ctor for full infinite intersections
Intersection::Intersection(const Ray& ray, const Geometry &geo) : Intersection() {
    Vector3d diff = ray.origin() - geo.origin();
    Vector3d ref_point = diff - diff.dot(ray.direction()) * ray.direction();
    Vector3d point_in = ref_point - geo.big_scale() * ray.direction();
    Vector3d point_out = ref_point + geo.big_scale() * ray.direction();
    init(ray, point_in, point_out);
    is_finite_ = false;
}

// Ctor for half infinite intersections
Intersection::Intersection(const Ray& ray, const Vector3d& point, const Geometry &geometry) : Intersection() {
    // Check if ``ray`` was inside ``geometry`` before coming to ``point``.
    Vector3d check_point = point-ray.direction();
    if (geometry.is_within(check_point)) {
        init(ray, point-geometry.big_scale()*ray.direction(), point);
    } else {
        init(ray, point, point+geometry.big_scale()*ray.direction());
    }
    is_finite_ = false;
}

bool Intersection::is_finite() const {
    return is_finite_;
}

// Cause C++ can't call ctors inside ctor (internally created object will be lost once exiting its scope)
void Intersection::init(const Ray& ray, const Vector3d& point_in, const Vector3d& point_out) {
    set_direction(ray.direction());
    set_borders(point_in, point_out);
}

void Intersection::set_point_in(const Vector3d& point) {
    borders_.first = point;
}

void Intersection::set_point_out(const Vector3d& point) {
    borders_.second = point;
}

double Intersection::length_pc() const {
    return (borders_.first - borders_.second).norm()/pc;
}
