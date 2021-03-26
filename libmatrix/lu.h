#ifndef _MX_LU_H
#define _MX_LU_H

#include "matrix.h"

namespace mx
{

class LinearSolver
{
    Matrix lu;

public:
    void set_matrix( Matrix mat ) { lu = mat; }
    void lu_decomp();
    Matrix solve_lower_triangular( const Matrix& mat, const Matrix& b_vec );
    Matrix solve_upper_triangular( const Matrix& mat, const Matrix& b_vec );
    Matrix get_lower();
    Matrix get_upper();
};

}

#endif
