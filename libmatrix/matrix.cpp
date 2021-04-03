#include "libmatrix/matrix.h"

namespace mx
{

void Matrix::print( std::ostream& os ) const
{
    for( int i=0; i<_n_row; i++ )
    {
        for( int j=0; j<_n_col; j++ )
        {
            os << std::setprecision(6) << std::setw(12) << _mat[ index(i,j) ];
        }
        os << std::endl;
    }
}

Matrix::Matrix()
{
    _n_row = 0;
    _n_col = 0;
}

Matrix::Matrix( MatrixInitilizer mx_init )
{
    switch( mx_init.matrix_type )
    {
        case MatInit::ZEROS:
            resize( mx_init.size, mx_init.size, 0.0 );
            break;
        case MatInit::EYE:
            resize( mx_init.size, mx_init.size, 0.0 );
            for( int i=0; i<mx_init.size; i++ )
                (*this)(i,i) = 1.0;
            break;
        case MatInit::RAND:
            init_mat_random( mx_init.size, []{ return _mat_rng.rand_1000(); }  );
            break;
        case MatInit::RAND_SYM:
            init_mat_rand_sym( mx_init.size, []{ return _mat_rng.rand_normal(); } );
            break;
        case MatInit::RAND_LOWTRI:
            init_mat_rand_low_tri( mx_init.size, []{ return _mat_rng.rand_10(); } );
            break;
        case MatInit::RAND_SPD:
            init_mat_rand_spd( mx_init.size, []{ return _mat_rng.rand_1(); } );
            break;
        default:
            assert( false && "Bad MatrixInitilizer type" );
    }
}

Matrix::Matrix( int row, int col, double val )
{
    resize( row, col, val );
}

Matrix::Matrix( std::tuple<int,int> s, double val )
{
    resize( std::get<0>(s), std::get<1>(s), val );
}

Matrix::Matrix( std::initializer_list<double> list )
{
    _n_row = list.size();
    _n_col = 1;
    _mat.reserve( _n_row*_n_col );
    _mat.insert( _mat.end(), list.begin(), list.end() );
}

Matrix::Matrix( std::initializer_list< std::initializer_list<double> > lists )
{
    _n_row = lists.size();
    _n_col = 0;
    for( auto list : lists )
        if( _n_col < (int)list.size() ) _n_col = (int)list.size();
    _mat.reserve( _n_row*_n_col );

    for( auto list : lists )
    {
        _mat.insert( _mat.end(), list.begin(), list.end() );
        for( int i=list.size(); i<_n_col; i++ ) _mat.push_back(0.0);
    }
}

Matrix::Matrix( std::vector< std::vector<double> > vecs )
{
    _n_row = vecs.size();
    _n_col = 0;
    for( auto vec : vecs )
        if( _n_col < (int)vec.size() ) _n_col = (int)vec.size();

    _mat.reserve( _n_row*_n_col );
    for( auto vec : vecs )
    {
        _mat.insert( _mat.end(), vec.begin(), vec.end() );
        for( int i=vec.size(); i<_n_col; i++ ) _mat.push_back(0.0);
    }
}

Matrix::Matrix( const char* file_name )
{
    read_from_file( file_name );
}

void Matrix::read_from_file( const char* file_name )
{
    std::ifstream ifs( file_name );
    if( !ifs.is_open() )
    {
        std::cerr << "Cannot open file: " << file_name << std::endl;
        assert( false && "Failed to open file" );
    }

    int row, col;
    std::string type;
    ifs >> row >> col >> type;
    assert( row>0 && col>0 );
    if( type!="DENSE" ) std::cerr << "Error: only DENSE matrix is supported now!" << std::endl;
    assert( type=="DENSE" );

    resize( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            ifs >> (*this)(i,j);
    ifs.close();
}
void Matrix::write_to_file( const char* file_name, int precision )
{
    assert( _n_row>0 && _n_col>0 );
    assert( precision>0 );

    std::ofstream ofs( file_name );
    if( !ofs.is_open() )
    {
        std::cerr << "Cannot open file: " << file_name << std::endl;
        assert( false && "Failed to open file" );
    }

    ofs << _n_row << " " << _n_col << " DENSE" << std::endl;

    for( int i=0; i<_n_row; i++ )
    {
        for( int j=0; j<_n_col; j++ )
        {
            ofs << std::setprecision(precision) << (*this)(i,j) << " ";
        }
        ofs << std::endl;
    }
    ofs.close();
}

double& Matrix::operator()( int row, int col )
{
    return _mat[ index(row,col) ];
}
double& Matrix::operator()( int idx )
{
    assert( _n_col==1 && "mat is not a vector" );
    assert( idx>=0 && idx<_n_row );
    return _mat[ idx ];
}

double Matrix::operator()( int row, int col ) const
{
    return _mat[ index(row,col) ];
}
double Matrix::operator()( int idx ) const
{
    assert( _n_col==1 && "mat is not a vector" );
    assert( idx>=0 && idx<_n_row );
    return _mat[ idx ];
}

void Matrix::resize( int row, int col, double val )
{
    _mat.clear();
    _n_row = row;
    _n_col = col;
    _mat.resize( row*col, val );
}

int Matrix::size( int dim ) const
{
    assert( dim==0 || dim==1 );
    return (dim==0) ? _n_row : _n_col;
}

Matrix Matrix::transpose() const
{
    assert( _n_row>0 && _n_col>0 );
    Matrix res(_n_col,_n_row);
    for( int i=0; i<_n_row; i++ )
        for( int j=0; j<_n_col; j++ )
            res(j,i) = (*this)(i,j);
    return res;
}

template<typename F>
void Matrix::init_mat_random( int n, F&& rand )
{
    resize( n, n );
    for( int i=0; i<n; i++ )
        for( int j=0; j<n; j++ )
            _mat[ index(i,j) ] = rand();
}

template<typename F>
void Matrix::init_mat_rand_sym( int n, F&& rand )
{
    init_mat_rand_low_tri( n, rand );

    for( int i=0; i<n-1; i++ )
        for( int j=i+1; j<n; j++ )
            _mat[ index(i,j) ] = _mat[ index(j,i) ];
}

template<typename F>
void Matrix::init_mat_rand_low_tri( int n, F&& rand )
{
    resize( n, n );
    for( int i=0; i<n; i++ )
        for( int j=0; j<=i; j++ )
            _mat[ index(i,j) ] = rand();
}

template<typename F>
void Matrix::init_mat_rand_spd( int n, F&& rand )
{
    init_mat_rand_low_tri( n, rand );
    for( int i=0; i<n; i++ )
        for( int j=0; j<i; j++ )
            for( int k=0; k<=j; k++ )
                (*this)(j,i) += (*this)(j,k) * (*this)(i,k);
    for( int i=0; i<n; i++ )
    {
        double sum = 0.0;
        for( int k=0; k<=i; k++ )
            sum += (*this)(i,k)*(*this)(i,k);
        (*this)(i,i) = sum;
    }

    for( int i=0; i<n; i++ )
        for( int j=0; j<i; j++ )
            (*this)(i,j) = (*this)(j,i);
}

double Matrix::norm( int p )
{
    /// lp-norm, p<=0 for infinite form (max of abs of entries)
    if( p<=0 ) return norm_inf();
    if( p==1 ) return norm_1();

    double res = 0.0;
    for( int i=0; i<_n_row; i++ )
        for( int j=0; j<_n_col; j++ )
            res += std::pow( (*this)(i,j), p );

    if( p==2 ) return std::sqrt( res );

    return std::pow( res, 1.0/(double)p );
}

double Matrix::norm_1()
{
    /// l1-norm is the sum of abs of all entries
    double res = 0.0;
    for( int i=0; i<_n_row; i++ )
        for( int j=0; j<_n_col; j++ )
            res += std::abs( (*this)(i,j) );
    return res;
}

double Matrix::norm_inf()
{
    /// return max of abs of entries
    return std::abs( *(std::max_element( _mat.begin(), _mat.end(), [](double x, double y){ return std::abs(x)<std::abs(y); } )) );
}

Matrix Matrix::submatrix( int r_beg, int r_end, int c_beg, int c_end )
{
    /// return submatrix as new Matrix
    if( r_beg<0 ) r_beg = _n_row + r_beg;
    if( r_end<0 ) r_end = _n_row + r_end;
    if( c_beg<0 ) c_beg = _n_col + c_beg;
    if( c_end<0 ) c_end = _n_col + c_end;
    assert( r_beg>=0 && r_beg<=r_end );
    assert( c_beg>=0 && c_beg<=c_end );
    int row = r_end-r_beg+1, col = c_end-c_beg+1;
    Matrix sub( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            sub(i,j) = (*this)( r_beg+i, c_beg+j );
    return sub;
}

void Matrix::swap_row( int i, int j )
{
    /// swap two rows
    assert( i>=0 && i<_n_row );
    assert( j>=0 && j<_n_row );
    if( i==j ) return;
    for( int k=0; k<_n_col; k++ )
    {
        std::swap( (*this)(i,k), (*this)(j,k) );
    }
}

void Matrix::swap_col( int i, int j )
{
    /// swap two columns
    assert( i>=0 && i<_n_col );
    assert( j>=0 && j<_n_col );
    if( i==j ) return;
    for( int k=0; k<_n_row; k++ )
    {
        std::swap( (*this)(k,i), (*this)(k,j) );
    }
}

}
