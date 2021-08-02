#ifndef MHD_TRANSFER_INCLUDE_LOGSPACE_H_
#define MHD_TRANSFER_INCLUDE_LOGSPACE_H_


template<typename T>
class Logspace {
    private:
        T curValue, base;

    public:
        Logspace(T first, T base) : curValue(first), base(base) {}

        T operator()() {
            T retval = curValue;
            curValue *= base;
            return retval;
        }
};

// https://stackoverflow.com/a/21429452
std::vector<double> pyLogspace(double start, double stop, int num = 50, double base = 10) {
    double realStart = pow(base, start);
    double realBase = pow(base, (stop-start)/(num-1));

    std::vector<double> retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, Logspace<double>(realStart,realBase));
    return retval;
}

#endif //MHD_TRANSFER_INCLUDE_LOGSPACE_H_
