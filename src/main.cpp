
#include "matrix.h"
#include "lu.h"
#include "third_party/Eigen/Dense"

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

    mx::Matrix mat6 = mx::RandSPD(10);
    Eigen::MatrixXd emat = mx_to_eigen( mat6 );
    Eigen::EigenSolver<Eigen::MatrixXd> eig_solver( emat );
    std::cout << "The eigenvalues of RandSPD matrix are:\n" << eig_solver.eigenvalues().real() << std::endl;

    std::cout << "RandSPD matrix = \n" << mat6 << std::endl;

    int n = 10;
    mx::Matrix mat3 = mx::Rand(n);
    mx::Matrix mat7 = mx::RandSPD(n);
    mx::LinearSolver ls( mat7 );
    ls.lu_decomp();
    std::cout << "|mat - LU| = " << ( ls.permute()*mat7 - ls.get_lower()*ls.get_upper()).norm() << std::endl;

    double abs_err = 0.0;
    for( int i=0; i<n; i++ )
    {
        mx::Matrix vec = mat3.submatrix( 0,-1, i, i );
        mx::Matrix sol = ls.solve_vec( vec );
        abs_err += ( mat7*sol - vec ).norm();
    }
    std::cout << "MAE = " << abs_err/n << std::endl;

    mx::Matrix mat = { { 1, 0,-1, 2},
                       { 0, 4, 2, 2},
                       {-1, 2, 6, 1},
                       { 2, 2, 1, 7} };

    char filename[] = "test.txt";
    mat7.write_to_file( filename, 17 );
    mx::Matrix mat_read( filename );
    std::cout << "mat_read - mat7 = " << (mat_read-mat7).norm() << std::endl;


    return 0;
}
