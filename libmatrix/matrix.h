#ifndef _MX_MATRIX_H
#define _MX_MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

namespace mx
{

struct MatrixInitilizer
{
    enum MxInitType{
        ZEROS,
        EYE
    } matrix_type;
    int size;
};

struct Zeros : public MatrixInitilizer
{
    Zeros( int s )
    {
        matrix_type = ZEROS;
        size = s;
    }
};

struct Eye : public MatrixInitilizer
{
    Eye( int s )
    {
        matrix_type = EYE;
        size = s;
    }
};

class Matrix
{
    int _n_row;
    int _n_col;
    std::vector< double > _mat;

    int index(int row, int col) const
    {
        assert( row>=0 && row<_n_row && col>=0 && col<_n_col );
        return row*_n_col + col;
    }

public:
    /* in basic.cpp */
    Matrix();
    Matrix( MatrixInitilizer mx_init );
    Matrix( int row, int col, double val=0.0 );
    Matrix( std::tuple<int,int>, double val=0.0 );
    Matrix( std::initializer_list<double> list );
    Matrix( std::initializer_list< std::initializer_list<double> > lists );
    Matrix( std::vector< std::vector<double> > vecs );
    double& operator()( int row, int col );
    double& operator()( int idx );
    double operator()( int row, int col ) const;
    double operator()( int idx ) const;
    void print( std::ostream& os=std::cout ) const;
    void resize( int row, int col, double val=0.0 );
    void reserve( int row, int col ){ _mat.reserve( row*col ); }
    std::tuple<int, int> size() const { return {_n_row, _n_col}; }
    int size( int dim ) const;
    int n_row() const { return _n_row; }
    int n_col() const { return _n_col; }
};

    /* in operation.cpp */
std::ostream& operator<<( std::ostream& os, const Matrix& mat );
Matrix operator+( const Matrix& mat1, const Matrix& mat2 );
Matrix operator*( const Matrix& mat1, const Matrix& mat2 );
Matrix operator-( const Matrix& mat1 );
Matrix operator-( const Matrix& mat1, const Matrix& mat2 );

}

#endif
