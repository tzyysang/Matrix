#include "lu.h"

namespace mx
{

LinearSolver::LinearSolver()
:   status(EMPTY),
    mode(NONE),
    abs_threshold(1e-16),
    _rank(-1)
{
}

LinearSolver::LinearSolver( const Matrix& mat )
:   status(EMPTY),
    mode(NONE),
    abs_threshold(1e-16),
    _rank(-1)
{
    set_matrix(mat);
}

void LinearSolver::set_matrix( const Matrix& mat )
{
    auto [row, col] = mat.size();

    if( row<=0 || col<=0 ) return;
    if( row!=col ) return;

    _mat = mat;
    perm.resize( row );
    for( int i=0; i<row; i++ )
        perm[i] = i;

    q_perm.resize( col );
    for( int i=0; i<col; i++ )
        q_perm[i] = i;

    _rank = -1;
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

std::tuple<int,int> LinearSolver::find_max_complete( int idx )
{
    /// find the max_abs entris in mat[ idx:end, idx:end ]
    int max_i = idx;
    int max_j = idx;
    double max_val = std::abs( _mat(idx,idx) );
    int n = _mat.n_row();
    for( int i=idx; i<n; i++ )
    {
        for( int j=idx; j<n; j++ )
        {
            double val = std::abs( _mat(i,j) );
            if( val > max_val )
            {
                max_val = val;
                max_i = i;
                max_j = j;
            }
        }
    }
    return std::make_tuple( max_i, max_j );
}

int LinearSolver::lu_decomp_partial()
{
    /// LU decompostition with partial pivoting
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
    mode = PARTIAL_LU;
    return 0;
}

int LinearSolver::lu_decomp()
{
    /// LU decomposition with complete pivoting
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    for( int k=0; k<row-1; k++ )
    {
        auto [ m, n ] = find_max_complete( k );
        _mat.swap_row( k, m );
        _mat.swap_col( k, n );
        perm[k] = m;
        q_perm[k] = n;

        double pivot = _mat(k,k);
        if( pivot==0.0 ) return -1;
        for( int i=k+1; i<row; i++ )
            _mat(i,k) = _mat(i,k) / pivot;

        for( int i=k+1; i<row; i++ )
            for( int j=k+1; j<col; j++ )
                _mat(i,j) = _mat(i,j) - _mat(i,k)*_mat(k,j);
    }

    status = LU_SUCCESS;
    mode = COMPLETE_LU;
    return 0;
}

int LinearSolver::find_max_pivot( int j )
{
    /// find max diagonal element in _mat[ j:end, j:end ]
    int n = _mat.n_row();
    int max_idx = j;
    double max_val = _mat(j,j);
    for( j=j+1; j<n; j++ )
    {
        if( _mat(j,j)>max_val )
        {
            max_val = _mat(j,j);
            max_idx = j;
        }
    }
    return max_idx;
}

int LinearSolver::chole_decomp_pivoting()
{
    /// Cholesky decomposition with pivoting
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    Matrix res = Zeros( row );

    for( int k=0; k<row; k++ )
    {
        int q = find_max_pivot( k );

        if( _mat(q,q)<=0.0 )
        {
            _rank = k;
            return -1;
        }

        res.swap_col(k, q);
        _mat.swap_col(k, q);
        _mat.swap_row(k, q);
        perm[k] = q;

        res(k,k) = std::sqrt( _mat(k,k) );
        double r = 1.0/res(k,k);
        for( int j=k+1; j<col; j++ )
        {
            res(k,j) = r*_mat(k,j);
        }
        for( int i=k+1; i<row; i++ )
            for( int j=k+1; j<col; j++ )
                _mat(i,j) = _mat(i,j) - res(k,i)*res(k,j);

    }

    for( int i=0; i<row; i++ )
        for( int j=0; j<=i; j++)
            _mat(i,j) = res(j,i);

    status = CHOLE_SUCCESS;
    return 0;
}

int LinearSolver::chole_decomp()
{
    /// Cholesky decomposition
    auto [row, col] = _mat.size();
    assert( row>0 && col>0 );
    assert( row==col );

    for( int i=0; i<row; i++ )
    {
        if( _mat(i,i)<=0.0 ) return -1;

        for( int j=0; j<i+1; j++ )
        {
            double sum = 0.0;
            for( int k=0; k<j; k++ )
                sum += _mat(i,k) * _mat(j,k);

            if( i==j )
                _mat(i,j) = sqrt( _mat(i,i) - sum );
            else
                _mat(i,j) = 1.0 / _mat(j,j) * ( _mat(i,j) - sum );
        }
    }

    status = CHOLE_SUCCESS;
    mode = CHOLE;
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

Matrix LinearSolver::get_chole()
{
    assert( status == CHOLE_SUCCESS );
    int row = _mat.n_row();
    Matrix res = Zeros(row);
    for( int i=0; i<row; i++ )
        for( int j=0; j<=i; j++ )
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
    int r = rank();
    for( int i=0; i<r; i++ )
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
    /// swap for (Pb)
    Matrix pb = b;
    for( int i=0; i<pb.n_row()-1; i++ )
    {
        std::swap( pb(i), pb( perm[i] ) );
    }
    return pb;
}

Matrix LinearSolver::permute_vec_q( const Matrix& b )
{
    /// swap for (Qb), Q is the column permutation matrix
    Matrix pb = b;
    for( int i=pb.n_row()-1; i>=0; i-- )
    {
        std::swap( pb(i), pb( q_perm[i] ) );
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

Matrix LinearSolver::permute( const Matrix& mat )
{
    /// return P * A * Q
    assert( status == LU_SUCCESS );
    int n = perm.size();

    Matrix p_mat = mat;
    for( int i=0; i<n; i++ )
        p_mat.swap_row(i,perm[i]);
    for( int i=0; i<n; i++ )
        p_mat.swap_col(i,q_perm[i]);

    return p_mat;
}

Matrix LinearSolver::permute_chole( const Matrix& mat )
{
    /// return the P^t * A * P
    assert( status == CHOLE_SUCCESS );
    int n = perm.size();

    Matrix p_mat = mat;
    for( int i=0; i<n; i++ )
        p_mat.swap_row(i,perm[i]);
    for( int i=0; i<n; i++ )
        p_mat.swap_col(i,perm[i]);

    return p_mat;
}

int LinearSolver::rank()
{
    int size = _mat.n_row();
    if( _rank!=-1 ) return _rank;
    if( mode!=COMPLETE_LU ) return size;

    double threshold = _mat(0,0) * abs_threshold * size;
    for( int i=0; i<size; i++ )
        if( std::abs(_mat(i,i)) < threshold ) return _rank = i;
    return _rank = size;
}

Matrix LinearSolver::solve_vec_chole( const Matrix& b )
{
    auto [row, col] = _mat.size();
    Matrix y( row, 1 );
    Matrix x( row, 1 );
    int r = rank();

    Matrix P_mat = Eye(row);
    for( int i=0; i<row; i++ )
    {
        P_mat.swap_col( perm[i], i );
    }

    Matrix pb = P_mat.transpose()*b;

    /// solve L
    for( int i=0; i<r; i++ )
    {
        y(i,0) += pb(i,0);
        for( int j=0; j<i; j++ )
        {
            y(i,0) -= _mat(i,j) * y(j,0);
        }
        y(i,0) /= _mat(i,i);
    }

    /// solve L^*
    for( int i=row-1; i>=0; i-- )
    {
        x(i,0) += y(i,0);
        for( int j=row-1; j>i; j-- )
        {
            x(i,0) -= _mat(j,i) * x(j,0);
        }
        x(i,0) /= _mat(i,i);
    }

    return P_mat * x;
}

Matrix LinearSolver::solve_vec( const Matrix& b )
{
    assert( b.n_row()==_mat.n_row() );
    if( status==LU_SUCCESS )
        return permute_vec_q( solve_upper_triangular( solve_lower_triangular( permute_vec( b ) ) ) );
    if( status==CHOLE_SUCCESS )
        return solve_vec_chole( b );

    assert( status==LU_SUCCESS || status==CHOLE_SUCCESS );
    return Matrix();
}

}
