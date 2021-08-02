#ifndef MHD_TRANSFER_PIXEL_H
#define MHD_TRANSFER_PIXEL_H

#include <array>
#include <string>
#include <map>
#include <Eigen/Eigen>

using std::array;
using std::map;
using std::string;
using std::pair;
using Eigen::Vector3d;


class Pixel {
    public:
        Pixel(double size_, Vector3d &coordinate_, pair<unsigned long int,unsigned long int> &ij_);
        Vector3d getCoordinate();
        double getValue(const string& value);
        void setValue(const string& value, double newvalue);
        pair<unsigned long int,unsigned  long int> getij() {
            return ij_;
        }

    private:
        double size_;
        Vector3d coordinate_;
        pair<unsigned long int,unsigned  long int> ij_;
        map<string,double> values_;
};


#endif //MHD_TRANSFER_PIXEL_H
