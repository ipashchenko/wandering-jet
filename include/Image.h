#ifndef MHD_TRANSFER_IMAGE_H
#define MHD_TRANSFER_IMAGE_H

#include <array>
#include <string>
#include <Eigen/Eigen>
#include <memory>
#include "Pixel.h"

using std::array;
using std::vector;
using Eigen::Vector3d;


// ``pixel_scale`` should be number of cm in one pixel. If one pixel is 0.1 mas
// and ``mas_to_pc(z)`` parsecs in one mas for given redshift ``z`` then it
// should be ``0.1*mas_to_pc*pc_to_cm`` where pc_to_cm ~ 3*10^18cm.
class Image {
    public:
        Image(pair<unsigned long int, unsigned long int> image_size, double lg_pixel_size_start_cm,
            double lg_pixel_size_stop_cm, bool jet_side);
        const unsigned long int num_of_pixels_;
        const pair<unsigned long int,unsigned long int> image_size;
        vector<vector<double>> getImage(const string& value);
        vector<vector<double>> getPixelSizes();
        Vector3d getCoordinate(unsigned long int i, unsigned long int j);
        vector<Pixel>& getPixels();

    private:
        bool jet_side_;
        pair<unsigned long int,unsigned long int> image_size_;
        // Log10 of pixel sizes in cm
        double pixel_size_start_;
        double pixel_size_stop_;
        // n_across, n_along array with pixel sizes [in cm!]
        vector<vector<double>> pixel_sizes_;
        // 2 arrays of the coordinates of pixel centers (n_across, n_along) [in cm!]
        vector<vector<double>> pixel_center_coordinates_along_;
        vector<vector<double>> pixel_center_coordinates_across_;
        vector<Pixel> pixels_;
};


#endif //MHD_TRANSFER_IMAGE_H
