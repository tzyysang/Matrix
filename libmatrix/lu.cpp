#include "lu.h"

namespace mx
{

void LinearSolver::lu_decomp()
{
    auto [row, col] = lu.size();
    assert( row>0 && col>0 );
    assert( row==col );

    std::cout << "row=" << row << " col=" << col << std::endl;
}

Matrix LinearSolver::get_lower()
{
    auto [row, col] = lu.size();
    assert( row>0 && col>0 );
    assert( row==col );

    Matrix res = Eye(row);
    for( int i=1; i<row; i++ )
        for( int j=0; j<i; j++ )
            res(i,j) = lu(i,j);
    return res;
}

Matrix LinearSolver::get_upper()
{
    auto [row, col] = lu.size();
    assert( row>0 && col>0 );
    assert( row==col );

    Matrix res = Eye(row);
    for( int i=0; i<row; i++ )
        for( int j=i; j<col; j++ )
            res(i,j) = lu(i,j);
    return res;
}

Matrix LinearSolver::solve_lower_triangular( const Matrix& mat, const Matrix& b_vec )
{
    auto [row, col] = lu.size();
    assert( row>0 && col>0 );
    assert( row==col );
    assert( row==b_vec.n_row() );

    Matrix x( row, 1 );
    for( int i=0; i<row; i++ )
    {
        x(i,0) += b_vec(i,0);
        for( int j=0; j<i; j++ )
        {
            x(i,0) -= mat(i,j) * x(j,0);
        }
        x(i,0) /= mat(i,i);
    }
    return x;
}

Matrix LinearSolver::solve_upper_triangular( const Matrix& mat, const Matrix& b_vec )
{
    auto [row, col] = lu.size();
    assert( row>0 && col>0 );
    assert( row==col );
    assert( row==b_vec.n_row() );

    Matrix x( row, 1 );
    for( int i=row-1; i>=0; i-- )
    {
        x(i,0) += b_vec(i,0);
        for( int j=i+1; j<col; j++ )
        {
            x(i,0) -= mat(i,j) * x(j,0);
        }
        x(i,0) /= mat(i,i);
    }
    return x;
}

}
