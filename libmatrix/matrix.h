#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

namespace mx
{

class Matrix
{
    int _n_row;
    int _n_col;

public:
    Matrix()
    {
        _n_row = 0;
        _n_col = 0;
    }

    void print( std::ostream& os=std::cout );
};

}

#endif
