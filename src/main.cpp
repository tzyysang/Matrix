
#include "matrix.h"
#include "lu.h"
#include "Eigen/Dense"

Eigen::MatrixXd mx_to_eigen( mx::Matrix mat )
{
    auto [row, col] = mat.size();
    Eigen::MatrixXd eig_mat( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            eig_mat(i,j) = mat(i,j);
    return eig_mat;
}

int main( int argc, char* argv[] )
{
    std::cout << "hello world" << std::endl;

    mx::Matrix mat1 = mx::Eye(3);
    mx::Matrix mat2 = { {1, 4, 7},
                        {2, 5, 8},
                        {3, 6, -10} };
    std::cout << "mat1 = \n" << mat1 << std::endl;
    std::cout << "mat2 = \n" << mat2 << std::endl;
    std::cout << "2*mat2 - mat1/0.5 = \n" << 2*mat2 - mat1/0.5 << std::endl;

    std::cout << "norm1 = " << mat2.norm(1) << std::endl;
    std::cout << "norm2 = " << mat2.norm() << std::endl;
    std::cout << "norm-1 = " << mat2.norm(-1) << std::endl;

    mx::Matrix mat3( mx::Rand(7) );
    std::cout << "Rand matrix = \n" << mat3 << std::endl;


    mx::Matrix mat6 = mx::RandSPD(10);
    Eigen::MatrixXd emat = mx_to_eigen( mat6 );
    Eigen::EigenSolver<Eigen::MatrixXd> eig_solver( emat );
    std::cout << "The eigenvalues of RandSPD matrix are:\n" << eig_solver.eigenvalues().real() << std::endl;

    std::cout << "RandSPD matrix = \n" << mat6 << std::endl;

    mx::Matrix mat7 = mx::RandSPD(10);
    mx::LinearSolver ls( mat7 );
    ls.lu_decomp();
    mx::Matrix low = ls.get_lower();
    mx::Matrix up = ls.get_upper();

    std::cout << "|mat - LU| = " << (mat7 - low*up).norm() << std::endl;

    return 0;
}
