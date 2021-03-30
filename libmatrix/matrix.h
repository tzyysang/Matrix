#ifndef _MX_MATRIX_H
#define _MX_MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <functional>
#include <algorithm>
#include <fstream>

#include "matrix_init.h"
#include "matrix_range.h"
#include "rand.h"

namespace mx
{

class Matrix
{
    int _n_row;
    int _n_col;
    std::vector< double > _mat;
    inline static RNG _mat_rng;

    /* in basic.cpp */
public:
    Matrix();
    Matrix( MatrixInitilizer mx_init );
    Matrix( int row, int col, double val=0.0 );
    Matrix( std::tuple<int,int>, double val=0.0 );
    Matrix( std::initializer_list<double> list );
    Matrix( std::initializer_list< std::initializer_list<double> > lists );
    Matrix( std::vector< std::vector<double> > vecs );
    Matrix( const char* file_name );
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
    Matrix transpose() const;
    double norm( int p=2 );
    double norm_1();
    double norm_inf();
    Matrix submatrix( int r_beg, int r_end, int c_beg, int c_end );
    void read_from_file( const char* file_name );
    void write_to_file( const char* file_name, int precision=16 );
private:
    template<typename F>
    void init_mat_random( int n, F&& rand );
    template<typename F>
    void init_mat_rand_sym( int n, F&& rand );
    template<typename F>
    void init_mat_rand_low_tri( int n, F&& rand );
    template<typename F>
    void init_mat_rand_spd( int n, F&& rand );

    /* inline functions */
private:
    int index(int row, int col) const
    {
        assert( row>=0 && row<_n_row && col>=0 && col<_n_col );
        return row*_n_col + col;
    }

};

    /* in operation.cpp */
std::ostream& operator<<( std::ostream& os, const Matrix& mat );
Matrix operator+( const Matrix& mat1, const Matrix& mat2 );
Matrix operator*( const Matrix& mat1, const Matrix& mat2 );
Matrix operator-( const Matrix& mat1 );
Matrix operator-( const Matrix& mat1, const Matrix& mat2 );
Matrix operator*( double scalar, const Matrix& mat );
Matrix operator*( const Matrix& mat, double scalar );
Matrix operator/( const Matrix& mat, double scalar );

}

#endif
