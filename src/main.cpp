
#include "matrix.h"

int main( int argc, char* argv[] )
{
    std::cout << "hello world" << std::endl;

    mx::Matrix mat1 = {1, 2, 3};
    mx::Matrix mat2 = { {1, 4, 7},
                        {2, 5, 8},
                        {3, 6, 10}};
    std::cout << "mat1 = \n" << mat1 << std::endl;
    std::cout << "mat2 = \n" << mat2 << std::endl;

    mat1(1) = 3.33;
    std::cout << "mat1 = \n" << mat1 << std::endl;

    mat2(2,2) = 9.99;
    std::cout << "mat2 = \n" << mat2 << std::endl;

    mx::Matrix mat3(3,4);
    std::cout << "mat3 = \n" << mat3 << std::endl;
    mat3.resize( 4, 8, 1.11 );
    std::cout << "mat3 = \n" << mat3 << std::endl;
}
