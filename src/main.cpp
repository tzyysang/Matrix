
#include "matrix.h"
#include "lu.h"
#include "third_party/Eigen/Dense"

#include <cstring>

static Eigen::MatrixXd mx_to_eigen( mx::Matrix mat )
{
    auto [row, col] = mat.size();
    Eigen::MatrixXd eig_mat( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            eig_mat(i,j) = mat(i,j);
    return eig_mat;
}

static int LU_error_bench()
{
    /// test LU error |PA-LU| in random matrices
    std::cout << "[LU_error_benchmark]" << std::endl;

    int size = 200;
    int runs = 10;
    mx::LinearSolver ls;
    double error = 0.0;

    for( int i=0; i<runs; i++ )
    {
        mx::Matrix mat = mx::RandSPD(size);
        ls.set_matrix( mat );
        ls.lu_decomp();
        error += ( ls.permute()*mat - ls.get_lower()*ls.get_upper()).norm();
    }
    error /= (double)runs;
    std::cout << error << std::endl;

    if( error>1e-12 ) return -1;
    return 0;
}

static int run_benchmarks( int argc, char* argv[] )
{
    int status = 0;
    for( int i=1; i<argc; i++ )
    {
        if( std::strcmp( argv[i], "-bench_LU_error" ) == 0 )
            status = status || LU_error_bench();
    }

    return status;
}

int main( int argc, char* argv[] )
{

    if( argc>=2 )
    {
        return run_benchmarks( argc, argv );
    }

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

    int n = 500;
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
