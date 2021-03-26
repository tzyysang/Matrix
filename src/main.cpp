
#include "matrix.h"
#include "lu.h"

int main( int argc, char* argv[] )
{
    std::cout << "hello world" << std::endl;

    mx::Matrix mat1 = {1, 2, 3};
    mx::Matrix mat2 = { {1, 4, 7},
                        {2, 5, 8},
                        {3, 6, 10} };
    mx::Matrix perm = { {0, 1, 0},
                        {0, 0, 1},
                        {1, 0, 0} };
    std::cout << "mat1 = \n" << mat1 << std::endl;
    std::cout << "mat2 = \n" << mat2 << std::endl;

    mx::Matrix mat3 = mx::Eye(3);
    std::cout << "mat3 = \n" << mat3 << std::endl;

    mx::Matrix mat4 = mat2 - mx::Eye(3);
    std::cout << "mat4 = \n" << mat4 << std::endl;

    std::cout << "mat4 = \n" << perm * mat4 << std::endl;

    mx::Matrix mat5 = { { 1, 4, 7, 2, 5},
                        { 2, 1,-3, 0, 8},
                        {-1,-8, 5, 2,-4},
                        { 7, 3, 0, 4,-1},
                        { 5,-4, 2,-3,-7} };

    std::cout << "mat5 = \n" << mat5 << std::endl;
    mx::LinearSolver ls;
    ls.set_matrix( mat5 );
    mx::Matrix low = ls.get_lower();
    std::cout << "low = \n" << low << std::endl;

    mx::Matrix b_vec = { 1, 2, 3, 4, 5 };
    std::cout << "b_vec = \n" << b_vec << std::endl;
    mx::Matrix ans = ls.solve_lower_triangular( low, b_vec );
    std::cout << "ans = \n" << ans << std::endl;

    mx::Matrix bb = low * ans;
    std::cout << "bb = \n" << bb << std::endl;

    mx::Matrix up  = ls.get_upper();
    std::cout << "up  = \n" << up << std::endl;
    mx::Matrix ans_up = ls.solve_upper_triangular( up, b_vec );
    std::cout << "ans_up = \n" << ans_up << std::endl;

    mx::Matrix cc = up * ans_up;
    std::cout << "cc = \n" << cc << std::endl;

    return 0;
}
