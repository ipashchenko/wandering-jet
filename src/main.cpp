#include <string>
#include <iostream>
#include <fstream>
#include "ImagePlane.h"
#include "Observation.h"
#include "BField.h"
#include "VField.h"
#include "Cone.h"
#include "Cylinder.h"
#include "Parabaloid.h"
#include "Jet.h"
#include "utils.h"
#include <cmath>
#include "NField.h"
#include "Pixel.h"
#include <ctime>
#include <chrono>

using Eigen::Vector3d;
using Eigen::Matrix3Xd;
using std::vector;
using std::pair;
using namespace boost::numeric::odeint;
namespace ph = std::placeholders;

typedef std::chrono::high_resolution_clock Clock;


void run_on_analytic(double los_angle, double cone_half_angle, std::string epoch) {
	auto t1 = Clock::now();
	std::clock_t start;
	start = std::clock();

	// FIXME: z = 0 case leads to NaN intersections
    // M87
//    double redshift = 0.00436;
	// Tested
//	double redshift = 0.1;
	// 3C 454.3
	double redshift = 0.59;
//    double los_angle = 17.0*M_PI/180.0;

    // FIXME: tested
//    double los_angle = 0.75*M_PI/180.0;

	// Observed frequencies in GHz
    std::vector<double> nu_observed_ghz{15.4};
    std::vector<double> total_fluxes;
    // Frequencies in the BH frame in Hz
    std::vector<double> nu_bh;
    for(auto nu_obs_ghz : nu_observed_ghz) {
        nu_bh.push_back(nu_obs_ghz*1E+09*(1+redshift));
    }

    // Setting geometry ================================================================================================
    Vector3d origin = {0., 0., 0.};
    Vector3d direction = {0., 0., 1.};
    double big_scale = 1000*pc;
//	double cone_half_angle = 0.5*M_PI/180.0;
    Cone geometry(origin, direction, cone_half_angle, big_scale);

    // Setting components of B-fields ==================================================================================
  	// HelicalConicalBField jetbfield(1.00, 1, 89.*M_PI/180., true, 0.0, &geometry);
	  // ToroidalBField jetbfield(1.0, 1, true, 0.0, &geometry);
    ReversedPinchConicalBField jetbfield(0.5, 1.0, &geometry);
    std::vector<VectorBField*> vbfields;
    vbfields.push_back(&jetbfield);

    // Setting components of N-fields ==================================================================================
    // Exponent of the power-law particle Lorenz factor distribution
    double s = 2.5;
    double ds = 0.0;
    double gamma_min = 10.0;
    PowerLaw particles(s, gamma_min, "pairs");
    // Value at r = 1 pc
	double K_1 = 100;
    // Exponent of the decrease
	double n = 1.5;


    // Setting V-field =================================================================================================
    VField* vfield;
    bool central_vfield = true;
	double Gamma = 10.0;
    if (central_vfield) {
        vfield = new ConstCentralVField(Gamma, &geometry, 0.0);
    } else {
        vfield = new ConstFlatVField(Gamma, &geometry, 0.0);
    }
	
    // BKNField power_law_nfield_spine(K_1, n, &particles, true, &geometry, nullptr, vfield);
	// Equipartition particles density
	EquipartitionBKNfield power_law_nfield_spine(&particles, vbfields, &geometry, nullptr, vfield);
	std::vector<NField*> nfields;
	nfields.push_back(&power_law_nfield_spine);
	
	
    Jet bkjet(&geometry, vfield, vbfields, nfields);

    // FIXME: Put inside frequency loop for dep. on frequency
    // Setting parameters of pixels and image ==========================================================================
	int number_of_pixels_along = 500;
	int number_of_pixels_across = 300;
    // Non-uniform pixel from ``pixel_size_mas_start`` (near BH) to ``pixel_size_mas_stop`` (image edges)
	double pixel_size_mas_start = 0.01;
	double pixel_size_mas_stop = 0.1;
    auto image_size = std::make_pair(number_of_pixels_across, number_of_pixels_along);
    auto pc_in_mas = mas_to_pc(redshift);
    std::cout << "pc_in_mas " << pc_in_mas << std::endl;
    // Log10 of pixel size in cm
    auto lg_pixel_size_start = log10(pixel_size_mas_start*pc_in_mas*pc);
    auto lg_pixel_size_stop = log10(pixel_size_mas_stop*pc_in_mas*pc);

	std::cout << "along x across : " << number_of_pixels_along << ", " << number_of_pixels_across << "\n";
    std::cout << "Setting pixel size (pc) from " << pow(10.0, lg_pixel_size_start)/pc << " to " << pow(10.0, lg_pixel_size_stop)/pc << std::endl;
    for(auto jet_side : {true, false}) {

        // Ignore CJ
        if(jet_side == false) {
            continue;
        }

        ImagePlane imagePlane(image_size, lg_pixel_size_start, lg_pixel_size_stop, los_angle, jet_side);
        // Array of pixel sizes in cm
        auto pixel_sizes = imagePlane.getPixelSizes();
        // Array of pixel solid angles in rad*rad
        std::vector<std::vector<double>> pixel_solid_angles;
        pixel_solid_angles.resize(pixel_sizes.size());

        for(unsigned long i=0; i < pixel_sizes.size(); i++) {
            pixel_solid_angles[i].resize(pixel_sizes[0].size());
            for(unsigned long j=0; j < pixel_sizes[0].size(); j++) {
                // Divide by ``pc_in_mas*pc`` to bring ``cm`` to ``mas`` at source redshift
                pixel_solid_angles[i][j] = (pixel_sizes[i][j]/(pc_in_mas*pc))*(pixel_sizes[i][j]/(pc_in_mas*pc))*mas_to_rad*mas_to_rad;
            }
        }

        // Array of scale factors. Divide resulting image on this to obtain flux density in Jy. Accounts for cosmological
        // scaling of intensity
        std::vector<std::vector<double>> scales;
        scales.resize(pixel_sizes.size());
        for(unsigned long i=0; i < pixel_sizes.size(); i++) {
            scales[i].resize(pixel_sizes[0].size());
            for(unsigned long j=0; j < pixel_sizes[0].size(); j++) {
                scales[i][j] = 1E-23*(1.+redshift)*(1.+redshift)*(1.+redshift)/pixel_solid_angles[i][j];
            }
        }

        Observation observation(&bkjet, &imagePlane);

        // FIXME: Put out of frequency loop - these do not depend on frequency
        // Setting transfer-specific parameters ========================================================================
        double tau_max = 10;
        double dt_max_pc = 0.01;
        double dt_max = pc*dt_max_pc;
        double tau_min_log10 = -10.0;
        double tau_min = pow(10.,tau_min_log10);
        int n_ = 100;
        double relerr = 1e-10;

        // Solve for all Stokes parameters ("full") or only full intensity ("I")?
//        string polarization = "I";
        string polarization = "full";

        for(int i_nu=0; i_nu < nu_observed_ghz.size(); i_nu++) {
            if(jet_side) {
                std::cout << "Running transfer for frequency " << nu_observed_ghz[i_nu] << " GHz for approaching jet" << std::endl;
            } else {
                std::cout << "Running transfer for frequency " << nu_observed_ghz[i_nu] << " GHz for counter-jet" << std::endl;
            }
            observation.run(n_, tau_max, dt_max, tau_min, nu_bh[i_nu], polarization, relerr);
            string value = "tau";
            auto image_tau = observation.getImage(value);

            value = "I";
            auto image_i = observation.getImage(value);
            for (unsigned long int i = 0; i < image_i.size(); ++i) {
                for (unsigned long int j = 0; j < image_i[i].size(); ++j) {
                    image_i[i][j] = image_i[i][j]/scales[i][j];
                }
            }

            value = "l";
            auto image_l = observation.getImage(value);

            std::fstream fs;
            // Remove trailing zeros: https://stackoverflow.com/a/46424921
            std::ostringstream oss;
            oss << std::setprecision(8) << std::noshowpoint << nu_observed_ghz[i_nu];
            std::string freq_name = oss.str();

            std::string file_tau, file_tau_fr, file_i, file_q, file_u, file_v, file_l;
            if(jet_side) {
                file_tau = "jet_image_tau_" + freq_name + ".txt";
                file_tau_fr = "jet_image_taufr_" + freq_name + ".txt";
                file_i = "jet_image_i_" + freq_name + "_" + epoch + ".txt";
                file_q = "jet_image_q_" + freq_name + "_" + epoch + ".txt";
                file_u = "jet_image_u_" + freq_name + "_" + epoch + ".txt";
                file_v = "jet_image_v_" + freq_name + "_" + epoch + ".txt";
                file_l = "jet_image_l_" + freq_name + "_" + epoch + ".txt";
            } else {
                file_tau = "cjet_image_tau_" + freq_name + ".txt";
                file_tau_fr = "cjet_image_taufr_" + freq_name + ".txt";
                file_i = "cjet_image_i_" + freq_name + ".txt";
                file_q = "cjet_image_q_" + freq_name + ".txt";
                file_u = "cjet_image_u_" + freq_name + ".txt";
                file_v = "cjet_image_v_" + freq_name + ".txt";
                file_l = "cjet_image_l_" + freq_name + ".txt";
            }

            // Remove old file
            std::remove(file_i.c_str());
            std::remove(file_q.c_str());
            std::remove(file_u.c_str());
            std::remove(file_v.c_str());
            std::remove(file_l.c_str());
            std::remove(file_tau.c_str());
            std::remove(file_tau_fr.c_str());

            fs.open(file_tau, std::ios::out | std::ios::app);
            if (fs.is_open()) {
                write_2dvector(fs, image_tau);
                fs.close();
            }

            fs.open(file_i, std::ios::out | std::ios::app);
            if (fs.is_open()) {
                write_2dvector(fs, image_i);
                fs.close();
            }

            fs.open(file_l, std::ios::out | std::ios::app);
            if (fs.is_open()) {
                write_2dvector(fs, image_l, pc);
                fs.close();
            }

            if (polarization == "full") {
                value = "Q";
                auto image_q = observation.getImage(value);
                for (unsigned long int i = 0; i < image_q.size(); ++i) {
                    for (unsigned long int j = 0; j < image_q[i].size(); ++j) {
                        image_q[i][j] = image_q[i][j]/scales[i][j];
                    }
                }

                value = "U";
                auto image_u = observation.getImage(value);
                for (unsigned long int i = 0; i < image_u.size(); ++i) {
                    for (unsigned long int j = 0; j < image_u[i].size(); ++j) {
                        image_u[i][j] = image_u[i][j]/scales[i][j];
                    }
                }

                value = "V";
                auto image_v = observation.getImage(value);
                for (unsigned long int i = 0; i < image_v.size(); ++i) {
                    for (unsigned long int j = 0; j < image_v[i].size(); ++j) {
                        image_v[i][j] = image_v[i][j]/scales[i][j];
                    }
                }

                value = "tau_fr";
                auto image_tau_fr = observation.getImage(value);

                fs.open(file_tau_fr, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_tau_fr);
                    fs.close();
                }

                fs.open(file_q, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_q);
                    fs.close();
                }

                fs.open(file_u, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_u);
                    fs.close();
                }

                fs.open(file_v, std::ios::out | std::ios::app);
                if (fs.is_open()) {
                    write_2dvector(fs, image_v);
                    fs.close();
                }
            }
        }
    }

	std::cout << "CPU Time: "
						<< (std::clock() - start) / (double) (CLOCKS_PER_SEC)
						<< " s" << std::endl;
	auto t2 = Clock::now();
	std::cout << "User time: "
						<< std::chrono::duration_cast<std::chrono::seconds>(
								t2 - t1).count()
						<< " s" << std::endl;
}

int main(int argc, char *argv[]) {
	// These values were used to obtained double humped polarization structure
//	double los_angle = 0.75*M_PI/180.0;
	if(argc != 4){
		std::cout << argc << "\n";
		std::cout << "Supply LOS-angle (rad), half-opening angle (rad) and epoch\n";
		return 1;
	}

	double los_angle = atof(argv[1]);
	double cone_half_angle = atof(argv[2]);
	std::string epoch = argv[3];
    run_on_analytic(los_angle, cone_half_angle, epoch);
    return 0;
}
