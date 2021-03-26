#include "libmatrix/matrix.h"

namespace mx
{

Matrix operator+( const Matrix& mat1, const Matrix& mat2 )
{
    assert( mat1.size()==mat2.size() );
    auto [row, col] = mat1.size();
    Matrix res( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            res(i,j) = mat1(i,j) + mat2(i,j);

    return res;
}

Matrix operator-( const Matrix& mat1 )
{
    auto [row, col] = mat1.size();
    Matrix res( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            res(i,j) = -mat1(i,j);
    return res;
}

Matrix operator-( const Matrix& mat1, const Matrix& mat2 )
{
    return mat1 + (-mat2);
}

Matrix operator*( const Matrix& mat1, const Matrix& mat2 )
{
    assert( mat1.n_col()==mat2.n_row() );
    int row = mat1.n_row(), col = mat2.n_col(), len = mat1.n_col();
    Matrix res( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            for( int k=0; k<len; k++ )
                res(i,j) += mat1(i,k) * mat2(k,j);
    return res;
}

std::ostream& operator<<( std::ostream& os, const Matrix& mat )
{
    mat.print(os);
    return os;
}

}
