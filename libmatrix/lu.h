#ifndef _MX_LU_H
#define _MX_LU_H

#include "matrix.h"

namespace mx
{

enum LinearSolverStatus{
    EMPTY,
    MAT_SET,
    LU_SUCCESS,
    CHOLE_SUCCESS
};

class LinearSolver
{
    Matrix _mat;
    LinearSolverStatus status;
    Matrix solve_lower_triangular( const Matrix& b_vec );
    Matrix solve_upper_triangular( const Matrix& b_vec );

public:
    LinearSolver() { status = EMPTY; }
    LinearSolver( const Matrix& mat ) { _mat = mat; status = MAT_SET; }
    void set_matrix( const Matrix& mat ) { _mat = mat; status = MAT_SET; }
    int lu_decomp();
    Matrix get_lower();
    Matrix get_upper();
    LinearSolverStatus get_status() { return status; }
    Matrix solve_vec( const Matrix& b );
};

}

#endif
