#include <iostream>
#include "Pixel.h"

using std::pair;


Pixel::Pixel(double size, Vector3d &coordinate, pair<unsigned long int,unsigned long int> &ij):
    size_(size),
    coordinate_(coordinate),
    ij_(ij) {};

Vector3d Pixel::getCoordinate() {
    return coordinate_;
}

double Pixel::getValue(const string& value) {
    return values_[value];
}

void Pixel::setValue(const string& value, double newvalue) {
    values_[value] = newvalue;
}
