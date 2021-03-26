#include "libmatrix/matrix.h"

namespace mx
{

std::ostream& operator<<( std::ostream& os, const Matrix& mat )
{
    mat.print(os);
    return os;
}

void Matrix::print( std::ostream& os ) const
{
    os ;
    for( int i=0; i<_n_row; i++ )
    {
        for( int j=0; j<_n_col; j++ )
        {
            os << std::setprecision(6) << std::setw(8) << _mat[ index(i,j) ];
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
        case MatrixInitilizer::MxInitType::ZEROS:
            resize( mx_init.size, mx_init.size, 0.0 );
            break;
        case MatrixInitilizer::MxInitType::EYE:
            resize( mx_init.size, mx_init.size, 0.0 );
            for( int i=0; i<mx_init.size; i++ )
                (*this)(i,i) = 1.0;
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

}
