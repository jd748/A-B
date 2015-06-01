#ifndef DP_HH
#define DP_HH
#include <vector>

class DP
{
public:

    double at(int, int, int);
    int size();
    void set(int, int, int, double);

    const int n;
    const double delta;
    const double M;

    DP(int, double, double);

private:
    std::vector<double> arr;
    int index(int, int, int);
    int num_i;
    int num_j;
};

#endif
