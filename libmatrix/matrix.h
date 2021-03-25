#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

namespace mx
{

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
    Matrix( int row, int col, double val=0.0 );
    Matrix( std::initializer_list<double> list );
    Matrix( std::initializer_list< std::initializer_list<double> > lists );
    Matrix( std::vector< std::vector<double> > vecs );
    double& operator()( int row, int col );
    double& operator()( int idx );
    double operator()( int row, int col ) const;
    double operator()( int idx ) const;
    void print( std::ostream& os=std::cout ) const;
    void reserve( int row, int col ){ _mat.reserve( row*col ); }
    void resize( int row, int col, double val=0.0 );
};


std::ostream& operator<<( std::ostream& os, const Matrix& mat );
}

#endif
