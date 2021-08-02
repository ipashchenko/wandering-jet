#ifndef MHD_TRANSFER_RAY_H
#define MHD_TRANSFER_RAY_H

#include <Eigen/Eigen>
#include "Geometry.h"

using Eigen::Vector3d;


class Ray {
    friend class Cone;
    public:
        Ray(Vector3d &origin, Vector3d &direction);
        Vector3d point(double t);
        const Vector3d& origin() const;
        const Vector3d& direction() const;

    protected:
        Vector3d origin_;
        // Direction from the observer
        Vector3d direction_;
};


#endif //MHD_TRANSFER_RAY_H
