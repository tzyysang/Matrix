#include "lu.h"

namespace mx
{

int LinearSolver::lu_decomp()
{
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    for( int k=0; k<row-1; k++ )
    {
        double pivot = _mat(k,k);
        if( pivot==0.0 ) return -1;
        for( int i=k+1; i<row; i++ )
            _mat(i,k) = _mat(i,k) / pivot;

        for( int i=k+1; i<row; i++ )
            for( int j=k+1; j<col; j++ )
                _mat(i,j) = _mat(i,j) - _mat(i,k)*_mat(k,j);
    }

    status = LU_SUCCESS;
    return 0;
}

Matrix LinearSolver::get_lower()
{
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    Matrix res = Eye(row);
    for( int i=1; i<row; i++ )
        for( int j=0; j<i; j++ )
            res(i,j) = _mat(i,j);
    return res;
}

Matrix LinearSolver::get_upper()
{
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    Matrix res = Zeros(row);
    for( int i=0; i<row; i++ )
        for( int j=i; j<col; j++ )
            res(i,j) = _mat(i,j);
    return res;
}

Matrix LinearSolver::solve_lower_triangular( const Matrix& mat, const Matrix& b_vec )
{
    auto [row, col] = _mat.size();
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
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );
    assert( row==b_vec.n_row() );

    Matrix x( row, 1 );
    for( int i=row-1; i>=0; i-- )
    {
        x(i,0) += b_vec(i,0);
        for( int j=col-1; j>i; j-- )
        {
            x(i,0) -= mat(i,j) * x(j,0);
        }
        x(i,0) /= mat(i,i);
    }
    return x;
}

Matrix LinearSolver::solve_vec( const Matrix& b )
{
    assert( b.n_row()==_mat.n_row() );
    return solve_upper_triangular( get_upper(), solve_lower_triangular( get_lower(), b ) );
}

}
