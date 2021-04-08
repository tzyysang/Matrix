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

enum LinearSolverMode{
    NONE,
    PARTIAL_LU,
    COMPLETE_LU,
    CHOLE
};

class LinearSolver
{
    Matrix _mat;
    LinearSolverStatus status;
    LinearSolverMode mode;
    Matrix solve_lower_triangular( const Matrix& b_vec );
    Matrix solve_upper_triangular( const Matrix& b_vec );
    std::vector<int> perm;
    std::vector<int> q_perm;
    double abs_threshold;

public:
    LinearSolver();
    LinearSolver( const Matrix& mat );
    void set_matrix( const Matrix& mat );
    int lu_decomp();
    int lu_decomp_partial();
    int chole_decomp();
    Matrix get_lower();
    Matrix get_upper();
    Matrix get_chole();
    LinearSolverStatus get_status() { return status; }
    Matrix solve_vec( const Matrix& b );
    int find_max( int j );
    std::tuple<int,int> find_max_complete( int idx );
    Matrix permute_vec( const Matrix& b );
    Matrix permute_vec_q( const Matrix& b );
    Matrix permute();
    Matrix permute( const Matrix& mat );
    int rank();
    Matrix matrix_lu() { return _mat; }
};

}

#endif
