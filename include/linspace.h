#ifndef JETPOL_LINSPACE_H
#define JETPOL_LINSPACE_H

#include <iostream>
#include <vector>

template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspaced;

    auto start = static_cast<double>(start_in);
    auto end = static_cast<double>(end_in);
    auto num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}


#endif //JETPOL_LINSPACE_H
