#include "utils.h"
#include "NField.h"
#include "MyExceptions.h"


NField::NField(bool in_plasma_frame, ParticlesDistribution* particles, Geometry* geometry_out, Geometry* geometry_in,
			   VField* vfield) :
    in_plasma_frame_(in_plasma_frame) {
    particles_ = particles;
    geometry_in_ = geometry_in;
    geometry_out_ = geometry_out;
	vfield_ = vfield;
}

double NField::nf(const Vector3d &point) const {
    double x, y, r_point, r_border, result;
    x = point[0];
    y = point[1];
    r_point = sqrt(x*x + y*y);
	
	
	double factor = 1.0;
	
	if(geometry_out_) {
		r_border = geometry_out_->radius_at_given_distance(point);
		if (r_point > r_border) {
			factor = exp(-pow(r_point - r_border, 2.0)/l_eps_N/l_eps_N);
		}
	}
//	if(geometry_in_) {
//		r_border = geometry_in_->radius_at_given_distance(point);
//		if (r_point < r_border) {
//			factor = exp(-pow(r_point - r_border, 2.0)/l_eps_N/l_eps_N);
//		}
//	}
	return _nf(point)*factor;
	
//    if(geometry_out_) {
////        std::cout << "Zeroing N' because of the OUTER geometry!" << "\n";
//        // Find radius of outer border at given point
//        r_border = geometry_out_->radius_at_given_distance(point);
//        if (r_point > r_border) {
//            return 0.0;
//        }
//    }
//    if(geometry_in_) {
////        std::cout << "Zeroing N' because of the INNER geometry!" << "\n";
//        // Find radius of inner border at given point
//        r_border = geometry_in_->radius_at_given_distance(point);
//        if (r_point < r_border) {
//            return 0.0;
//        }
//    }
//    return _nf(point);
}

double NField::nf_plasma_frame(const Vector3d &point, double &gamma) const {
    double n = nf(point);
    if (in_plasma_frame_) {
        return n;
    } else {
        return n/gamma;
    }
}


BKNField::BKNField(double n_0, double n_n, ParticlesDistribution* particles, bool in_plasma_frame,
                   Geometry* geometry_out, Geometry* geometry_in,  VField* vfield) :
        NField(in_plasma_frame, particles, geometry_out, geometry_in, vfield),
        n_0_(n_0),
        n_n_(n_n),
        r_mean_(0.0),
        r_width_(0.0),
        background_fraction_(0.0),
        is_profile_set_(false) {}


double BKNField::_nf(const Vector3d &point) const {
    double r = point.norm();
    double raw_density = n_0_ * pow(r / pc, -n_n_);
    if(is_profile_set_){
        double x = point[0];
        double y = point[1];
        double R_cur = hypot(x, y);
        double R_out = geometry_out_->radius_at_given_distance(point);
//        if(R_cur/R_out < 0.25){
//            std::cout << R_cur/R_out << "\n";
//        }
        return (generalized1_gaussian1d(R_cur / R_out, r_mean_, r_width_, 2) + background_fraction_) * raw_density;
    }else{
        return raw_density;
    }
}


void BKNField::set_heating_profile(double r_mean, double r_width, double background_fraction) {
    if (r_mean > 1 || r_mean < 0.0 || r_width < 0.0 || background_fraction < 0.0 || background_fraction_ > 1.0){
        throw BadHeatingProfileParameters();
    }
    r_mean_ = r_mean;
    r_width_ = r_width;
    background_fraction_ = background_fraction;
    is_profile_set_ = true;
}


EquipartitionBKNfield::EquipartitionBKNfield(ParticlesDistribution *particles,
											 std::vector<VectorBField*> vbfields,
											 Geometry *geometry_out, Geometry *geometry_in, VField *vfield,
											 double fraction):
	NField(true, particles, geometry_out, geometry_in, vfield),
	fraction_(fraction) {
	vbfields_ = vbfields;
}

double EquipartitionBKNfield::_nf(const Vector3d &point) const {
	auto v = vfield_->vf(point);
	double b = 0.0;
	for(auto bfield: vbfields_) {
		Vector3d vb = bfield->bf_plasma_frame(point, v);
		b += vb.norm();
	}
	double raw_density = fraction_*particles_->get_equipartition_bsq_coefficient()*b*b;
	return raw_density;
}
