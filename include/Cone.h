#ifndef BK_TRANSFER_CONE_H
#define BK_TRANSFER_CONE_H

#include "Geometry.h"


class Cone: public Geometry {
    public:
        Cone(const Vector3d &origin, const Vector3d &direction, const double &angle, const double &scale);
        const Vector3d& origin() const override ;
        const Vector3d& direction() const ;
        const double angle() const;
        std::list<Intersection> hit(Ray &ray) const override;
        bool is_within(Vector3d& point) const override;
        const double big_scale() const override ;
        double radius_at_given_distance(const Vector3d& point) const override;

    private:
        Vector3d origin_;
        Vector3d direction_;
        double angle_;
        double cos_angle_;
        double big_scale_;
};


#endif //BK_TRANSFER_CONE_H
