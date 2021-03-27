#ifndef _MX_RAND_H
#define _MX_RAND_H

#include <cstdlib>
#include <ctime>
#include <cmath>

namespace mx
{

class RNG
{

public:
    RNG() { std::srand( std::time(NULL) ); }
    double rand( double rnd_min, double rnd_max )
    {
        return (double)std::rand() / (double)RAND_MAX * (rnd_max - rnd_min) + rnd_min;
    }
    double rand_1000() { return rand( -1000, 1000 ); }
    double rand_10() { return rand( -10, 10 ); }
    double rand_1() { return rand( -1, 1 ); }
    double rand_normal() { return 0.3989422804 * std::exp( -0.5 * std::pow( rand(-1,1), 2 ) ); }
};

}

#endif
