#ifndef JETPOL_NFIELD_H
#define JETPOL_NFIELD_H

#include <Eigen/Eigen>
#include "Geometry.h"
#include "ParticlesDistribution.h"
#include "BField.h"
#include "VField.h"

using Eigen::Vector3d;

const double l_eps_N = 0.0001*pc;


class NField {
        friend class Jet;
    public:
        virtual double _nf(const Vector3d &point) const = 0;
        double nf(const Vector3d &point) const;
        double nf_plasma_frame(const Vector3d &point, double &gamma) const;


    protected:
        explicit NField(bool in_plasma_frame,
                        ParticlesDistribution* particles,
                        Geometry* geometry_in = nullptr,
                        Geometry* geometry_out = nullptr,
						VField* vfield = nullptr);
        bool in_plasma_frame_;
        ParticlesDistribution* particles_;

        // Inner border
        Geometry* geometry_in_;
        // Outer border
        Geometry* geometry_out_;
		VField* vfield_;
};


class BKNField: public NField {
    public:
        BKNField(double n_0, double n_n, ParticlesDistribution* particles, bool in_plasma_frame,
                 Geometry* geometry_out, Geometry* geometry_in = nullptr, VField* vfield = nullptr);
        double _nf(const Vector3d &point) const override;
        void set_heating_profile(double r_mean, double r_width, double background_fraction = 0.0);
    private:
        double n_0_;
        double n_n_;
        double r_mean_;
        double r_width_;
        double background_fraction_;
        bool is_profile_set_;
};


// Here we need V-field to transfer B-field to plasma frame via B' = B_lab/Gamma. Here we assume dominating transverse
// component.
class EquipartitionBKNfield : public NField {
	public:
		EquipartitionBKNfield(ParticlesDistribution* particles, std::vector<VectorBField*> vbfields,
							  Geometry* geometry_out, Geometry* geometry_in = nullptr, VField* vfield = nullptr,
							  double fraction = 1.0);
		double _nf(const Vector3d &point) const override;
	private:
		std::vector<VectorBField*> vbfields_;
		double fraction_;
};


//class FFNField: public PowerLawNField {
//    public:
//        FFNField(double M_BH, double A,  double nu, double s, double gamma_min);
//        double _nf(const Vector3d& point) const override;
//    private:
//        double M_BH_;
//        // Parameters of the stream function
//        double A_;
//        double nu_;
//        // Reference position
//        double z1_;
//        // Width of the Gaussian ring - 5*r_g
//        double delta_;
//        // Radius where n peaks on z = +/-z_1 - free parameters [0, 100*r_g]
//        double Rp_;
//        // Number density of nonthermal electrons at (R, z) = (R-p, +/-z_1)
//        double n0_;
//};
//
//

#endif //BK_TRANSFER_NFIELD_H
