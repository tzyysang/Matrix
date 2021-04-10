
#include "matrix.h"
#include "lu.h"
#include "third_party/Eigen/Dense"

#include <cstring>

static bool apprx_equal( double x, double y, double err=1e-6 )
{
    if( x==y ) return true;
    return ( std::abs(x-y)/std::abs(x) < err );
}

static Eigen::MatrixXd mx_to_eigen( mx::Matrix mat )
{
    auto [row, col] = mat.size();
    Eigen::MatrixXd eig_mat( row, col );
    for( int i=0; i<row; i++ )
        for( int j=0; j<col; j++ )
            eig_mat(i,j) = mat(i,j);
    return eig_mat;
}

static int bench_LU_error()
{
    /// test LU error |PA-LU| in random matrices
    std::cout << "[LU_error benchmark]" << std::endl;

    int size = 100;
    int runs = 10;
    mx::LinearSolver ls;
    double error = 0.0;

    for( int i=0; i<runs; i++ )
    {
        mx::Matrix mat = mx::RandSPD(size);
        ls.set_matrix( mat );
        ls.lu_decomp_partial();
        error += ( ls.permute()*mat - ls.get_lower()*ls.get_upper()).norm();
    }
    error /= (double)runs;
    std::cout << error << std::endl;

    if( error>1e-12 ) return -1;
    return 0;
}


static int bench_LU_solve_partial()
{
    /// test LU by solving Ax = b problems
    std::cout << "[LU_solve_partial benchmark]" << std::endl;

    int size = 300;
    mx::Matrix mat = mx::Rand(size);
    mx::Matrix b_vecs = mx::Rand(size);

    mx::LinearSolver ls(mat);
    ls.lu_decomp_partial();
    Eigen::MatrixXd eig_mat = mx_to_eigen(mat);

    double mae = 0.0;
    for( int i=0; i<1; i++ )
    {
        mx::Matrix b = b_vecs.submatrix(0,-1,i,i);
        mx::Matrix x = ls.solve_vec( b );

        Eigen::VectorXd bb = mx_to_eigen( b );
        Eigen::VectorXd eig_x = eig_mat.partialPivLu().solve(bb);

        for( int k=0; k<size; k++ )
        {
            mae += std::abs( x(k) - eig_x(k) );
            if( !apprx_equal( x(k), eig_x(k) ) ) return -1;
        }
    }
    std::cout << "root MAE = " << mae/size << std::endl;
    return 0;
}

static int bench_LU_complete()
{
    /// test complete LU by comparing the solve LU matrix
    std::cout << "[LU_complete benchmark]" << std::endl;
    int size = 300;
    mx::Matrix mat = mx::RandSPD(size);

    mx::Matrix b_vecs = mx::Rand(size);
    mx::LinearSolver ls(mat);
    ls.lu_decomp();
    mx::Matrix mat_LU = ls.matrix_lu();
    Eigen::MatrixXd eig_mat = mx_to_eigen(mat);
    Eigen::MatrixXd eig_LU = eig_mat.fullPivLu().matrixLU();

    double error = 0.0;
    for( int i=0; i<size; i++ )
    {
        for( int j=0; j<size; j++ )
        {
            error += std::abs( mat_LU(i,j) - eig_LU(i,j) );
            if( !apprx_equal( mat_LU(i,j), eig_LU(i,j) ) ) return -1;
        }
    }
    std::cout << "LU error = " << error << std::endl;
    return 0;
}

static int bench_LU_solve_complete()
{
    /// test complete LU by solving Ax=b problem
    // m17.txt is a non-full-rank matrix
    std::cout << "[LU_solve_complete benchmark]" << std::endl;

    int size = 17;
    mx::Matrix mat =
    {{0.6114542991615833,0.4959730689634669,0.4338442394906383,-0.6059840184538489,0.7734632143309984,0.6199151432222963,0.5957941603614909,0.4232669829797772,-0.2866920379405562,-0.5311551041133061,0.4093570365294774,-0.6517997882998695,0.6678207151630747,-0.1854401330763814,0.1905327375397907,-0.1045284468809141,-0.6210839022198178},
     {0.4959730689634669,0.5033541714337005,0.6003788224846217,-0.6950894038410733,0.8363331512150799,0.7563758664665459,0.4097944887651946,0.3442704360438816,-0.1997696728176254,-0.33425516944617,0.5889461962090665,-0.655351454469583,0.4123445047203768,-0.1226456444138011,0.2241986293213251,-0.3369505585726722,-0.5881937099473099},
     {0.4338442394906383,0.6003788224846217,1.030563493499112,-0.8644579133173537,0.9270200133832005,1.123147986120328,0.2233772074011102,0.4197185366107688,-0.04878040870166504,-0.4635191480624934,1.220576542703738,-0.6581647514664183,0.1515848829494345,0.268063359201695,0.5482187059533523,-0.818459480555451,-0.8031504168438786},
     {-0.6059840184538489,-0.6950894038410733,-0.8644579133173537,1.049576850232317,-1.266457150682637,-1.09191494320318,-0.4552092988006269,-0.3514219165548715,0.2629324838022246,0.1407978264092501,-0.7460540792816779,0.9692896555980832,-0.4021727041665062,0.3222712476427083,-0.1868171631751369,0.540573690721267,0.6941313784495257},
     {0.7734632143309984,0.8363331512150799,0.9270200133832005,-1.266457150682637,1.789952569260828,1.069804407722728,0.4103746607823333,0.4435283671822944,-0.4119233036062188,0.07778951319265202,0.9083888493269152,-1.467912539519186,0.7334424541675612,-0.764186171409888,0.03496986497172794,-0.2958940661309953,-0.5789737637414147},
     {0.6199151432222963,0.7563758664665459,1.123147986120328,-1.09191494320318,1.069804407722728,2.73705345527091,1.186696601389957,-0.2661652726016464,-1.224158226112906,-0.3887568718544729,1.177029813229312,-1.60051004931481,-0.4152828281549049,0.4195920026038625,0.7392336881054459,-2.435554762318429,-0.1631527068554603},
     {0.5957941603614909,0.4097944887651946,0.2233772074011102,-0.4552092988006269,0.4103746607823333,1.186696601389957,1.858122081750125,0.07782060700270607,-0.7912084261207234,0.3497974121148039,-0.7990946987201394,-0.2121141120423674,-0.56283627313839,-0.1199261882109738,0.910446946004995,-0.6491494843079402,-0.5122164481837134},
     {0.4232669829797772,0.3442704360438816,0.4197185366107688,-0.3514219165548715,0.4435283671822944,-0.2661652726016464,0.07782060700270607,1.940748786520396,-0.1884219935218615,-0.2985682714176112,1.003788507235018,-0.3953904126059085,-0.1511736348236684,0.2396646231240736,0.244711946550298,0.7898015038840533,-1.222439842738219},
     {-0.2866920379405562,-0.1997696728176254,-0.04878040870166504,0.2629324838022246,-0.4119233036062188,-1.224158226112906,-0.7912084261207234,-0.1884219935218615,2.193606410381424,0.1208210708065102,-0.9757430845049383,1.971187648256199,0.1032751462309084,0.8792372804506035,-0.4348857112051209,1.559542408263886,-1.274971515666162},
     {-0.5311551041133061,-0.33425516944617,-0.4635191480624934,0.1407978264092501,0.07778951319265202,-0.3887568718544729,0.3497974121148039,-0.2985682714176112,0.1208210708065102,3.455078401159608,-2.03491826559076,0.5289670748369611,-2.4149424434306,-0.1986979077151061,-0.3565229516426843,1.050178749092096,0.1733424755781697},
     {0.4093570365294774,0.5889461962090665,1.220576542703738,-0.7460540792816779,0.9083888493269152,1.177029813229312,-0.7990946987201394,1.003788507235018,-0.9757430845049383,-2.03491826559076,4.886989835605583,-1.026852423304454,2.589094504884017,1.245249159451097,-0.51913634136345,-1.256151960466295,0.1538402386663449},
     {-0.6517997882998695,-0.655351454469583,-0.6581647514664183,0.9692896555980832,-1.467912539519186,-1.60051004931481,-0.2121141120423674,-0.3953904126059085,1.971187648256199,0.5289670748369611,-1.026852423304454,4.333891751994419,0.2083905479442092,1.622258296459258,-0.7395984672645999,0.8997043450540816,-0.8047303433072759},
     {0.6678207151630747,0.4123445047203768,0.1515848829494345,-0.4021727041665062,0.7334424541675612,-0.4152828281549049,-0.56283627313839,-0.1511736348236684,0.1032751462309084,-2.4149424434306,2.589094504884017,0.2083905479442092,5.158755504204164,0.1957093025923436,-1.23832000482993,0.6018768603526936,0.03004774359702991},
     {-0.1854401330763814,-0.1226456444138011,0.268063359201695,0.3222712476427083,-0.764186171409888,0.4195920026038625,-0.1199261882109738,0.2396646231240736,0.8792372804506035,-0.1986979077151061,1.245249159451097,1.622258296459258,0.1957093025923436,3.769234982998008,-0.4783638816394558,0.7251814106392145,-1.252849164036135},
     {0.1905327375397907,0.2241986293213251,0.5482187059533523,-0.1868171631751369,0.03496986497172794,0.7392336881054459,0.910446946004995,0.244711946550298,-0.4348857112051209,-0.3565229516426843,-0.51913634136345,-0.7395984672645999,-1.23832000482993,-0.4783638816394558,3.339709603990589,-1.677822063704201,-0.3294081103904408},
     {-0.1045284468809141,-0.3369505585726722,-0.818459480555451,0.540573690721267,-0.2958940661309953,-2.435554762318429,-0.6491494843079402,0.7898015038840533,1.559542408263886,1.050178749092096,-1.256151960466295,0.8997043450540816,0.6018768603526936,0.7251814106392145,-1.677822063704201,6.354246433776926,-0.2511285936062306},
     {-0.6210839022198178,-0.5881937099473099,-0.8031504168438786,0.6941313784495257,-0.5789737637414147,-0.1631527068554603,-0.5122164481837134,-1.222439842738219,-1.274971515666162,0.1733424755781697,0.1538402386663449,-0.8047303433072759,0.03004774359702991,-1.252849164036135,-0.3294081103904408,-0.2511285936062306,4.39337363046909}};

    mx::Matrix b_vecs = mx::Rand(size);
    mx::LinearSolver ls(mat);
    ls.lu_decomp();

    Eigen::MatrixXd eig_mat = mx_to_eigen(mat);
    double error;
    for( int i=0; i<size; i++ )
    {
        mx::Matrix b = b_vecs.submatrix( 0,-1, i, i );
        mx::Matrix x = ls.solve_vec(b);

        Eigen::VectorXd bb = mx_to_eigen(b);
        Eigen::VectorXd xx = eig_mat.fullPivLu().solve(bb);

        for( int j=0; j<size; j++ )
        {
            error += std::abs( xx(j) - x(j) );
            if( !apprx_equal( xx(j), x(j) ) ) return -1;
        }
    }
    std::cout << "solved error = " << error << std::endl;
    return 0;
}

static int bench_Chole_decomp()
{
    /// test cholesky decomposition, no pivoting
    int size = 10;
    mx::Matrix mat = mx::RandSPD( size );
    mx::LinearSolver ls( mat );
    int status = ls.chole_decomp();
    auto L = ls.get_chole();
    auto LL = L * L.transpose();

    Eigen::MatrixXd eig_mat = mx_to_eigen( mat );
    Eigen::LLT<Eigen::MatrixXd> eig_ll( eig_mat );
    Eigen::MatrixXd eig_L = eig_ll.matrixL();

    Eigen::EigenSolver<Eigen::MatrixXd> es( eig_mat );
    auto V = es.eigenvalues().real();
    for( int i=0; i<V.size(); i++ )
        assert( V(i)>0.0 && "Matrix is not SPD!" );

    double error = 0.0;
    for( int i=0; i<size; i++ )
    {
        for( int j=0; j<=i; j++ )
        {
            error += std::abs( eig_L(i,j) - L(i,j) );
            if( std::isnan(L(i,j)) ) return -1;
            if( !apprx_equal( eig_L(i,j), L(i,j) ) )
                return -1;
        }
    }
    std::cout << "error = " << error << std::endl;
    return 0;
}

static int bench_Chole_decomp_pivot()
{
    /// test cholesky decomposition with pivoting
    int size = 100;
    mx::Matrix mat = mx::RandSPD( size );

    Eigen::MatrixXd eig_mat = mx_to_eigen( mat );
    Eigen::EigenSolver<Eigen::MatrixXd> es( eig_mat );
    auto V = es.eigenvalues().real();
    for( int i=0; i<V.size(); i++ )
        assert( V(i)>0.0 && "Matrix is not SPD!" );

    mx::LinearSolver ls( mat );
    int status = ls.chole_decomp_pivoting();
    auto L = ls.get_chole();
    auto LL = L * L.transpose();
    auto M = ls.permute_chole( mat );

    double error = 0.0;
    for( int i=0; i<size; i++ )
    {
        for( int j=0; j<=i; j++ )
        {
            error += std::abs( M(i,j) - LL(i,j) );
            if( std::isnan( LL(i,j) ) ) return -1;
            if( !apprx_equal( M(i,j), LL(i,j) ) ) return -1;
        }
    }
    std::cout << "LL^* error = " << error << std::endl;
    if( error > 1e-15*size*size ) return -1;

    /// Solve Ax=b problems
    double mae = 0.0;
    mx::Matrix b_vecs = mx::Rand(size);
    for( int i=0; i<size; i++ )
    {
        auto b = b_vecs.submatrix(0,-1,i,i);
        auto x = ls.solve_vec( b );
        auto bb = mat * x;

        for( int j=0; j<size; j++ )
        {
            if( std::isnan( x(j) ) ) return -1;
            //if( !apprx_equal( bb(j), b(j), 1e-2 ) ) return -1;
        }
        mae += (bb-b).norm();
    }
    mae = mae/size/size;
    std::cout << "root MAE = " << mae << std::endl;
    if( mae > 0.1 ) return -1;

    return 0;
}

static int run_benchmarks( int argc, char* argv[] )
{
    int status = 0;
    for( int i=1; i<argc; i++ )
    {
        if( std::strcmp( argv[i], "-bench_LU_error" ) == 0 )
            status = status || bench_LU_error();
        else if( std::strcmp( argv[i], "-bench_LU_solve_partial" ) == 0 )
            status = status || bench_LU_solve_partial();
        else if( std::strcmp( argv[i], "-bench_LU_complete" ) == 0 )
            status = status || bench_LU_complete();
        else if( std::strcmp( argv[i], "-bench_LU_solve_complete" ) == 0 )
            status = status || bench_LU_solve_complete();
        else if( std::strcmp( argv[i], "-bench_Chole_decomp" ) == 0 )
            status = status || bench_Chole_decomp();
        else if( std::strcmp( argv[i], "-bench_Chole_decomp_pivot" ) == 0 )
            status = status || bench_Chole_decomp_pivot();
        else
        {
            std::cerr << "invalid command: " << argv[i] << std::endl;
            return -1;
        }
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

    int n = 5;
    mx::Matrix mat3 = mx::Rand(n);
    mx::Matrix mat7 = mx::RandSPD(n);
    std::cout << "RandSPD matrix = \n" << mat7 << std::endl;
    mx::LinearSolver ls( mat7 );
    ls.lu_decomp();
    std::cout << "|mat - LU| = " << ( ls.permute( mat7 ) - ls.get_lower()*ls.get_upper()).norm() << std::endl;

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
