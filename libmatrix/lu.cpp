#include "lu.h"

namespace mx
{

void LinearSolver::set_matrix( const Matrix& mat )
{
    auto [row, col] = mat.size();

    if( row<=0 || col<=0 ) return;
    if( row!=col ) return;

    _mat = mat;
    perm.resize( row );
    for( int i=0; i<row; i++ )
        perm[i] = i;
    status = MAT_SET;
}

int LinearSolver::find_max( int j )
{
    /// find the max_abs entris in mat[ j:end, j ]
    int max_idx = j;
    double max_val = std::abs( _mat(j,j) );
    for( int i=j+1; i<_mat.n_row(); i++ )
    {
        double val = std::abs( _mat(i,j) );
        if( val > max_val )
        {
            max_val = val;
            max_idx = i;
        }
    }
    return max_idx;
}

int LinearSolver::lu_decomp()
{
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    for( int k=0; k<row-1; k++ )
    {
        int m = find_max( k );
        _mat.swap_row( k, m );
        perm[k] = m;

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

Matrix LinearSolver::solve_lower_triangular( const Matrix& b_vec )
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
            x(i,0) -= _mat(i,j) * x(j,0);
        }
    }
    return x;
}

Matrix LinearSolver::solve_upper_triangular( const Matrix& b_vec )
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
            x(i,0) -= _mat(i,j) * x(j,0);
        }
        x(i,0) /= _mat(i,i);
    }
    return x;
}

Matrix LinearSolver::permute_vec( const Matrix& b )
{
    Matrix pb = b;
    for( int i=0; i<pb.n_row()-1; i++ )
    {
        std::swap( pb(i), pb( perm[i] ) );
    }
    return pb;
}

Matrix LinearSolver::permute()
{
    /// return the permutation matrix
    assert( status == LU_SUCCESS );
    int n = perm.size();
    std::vector<int> p_vec( n );
    for( int i=0; i<n; i++ )
        p_vec[i] = i;
    for( int i=0; i<n; i++ )
        std::swap( p_vec[i], p_vec[ perm[i] ] );

    Matrix p_mat = Zeros(n);
    for( int i=0; i<n; i++ )
        p_mat(i,p_vec[i]) = 1.0;

    return p_mat;
}

Matrix LinearSolver::solve_vec( const Matrix& b )
{
    assert( b.n_row()==_mat.n_row() );
    return solve_upper_triangular( solve_lower_triangular( permute_vec( b ) ) );
}

}
