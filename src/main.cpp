
#include "matrix.h"
#include "lu.h"

int main( int argc, char* argv[] )
{
    std::cout << "hello world" << std::endl;

    mx::Matrix mat1 = mx::Eye(3);
    mx::Matrix mat2 = { {1, 4, 7},
                        {2, 5, 8},
                        {3, 6, 10} };
    std::cout << "mat1 = \n" << mat1 << std::endl;
    std::cout << "mat2 = \n" << mat2 << std::endl;
    std::cout << "2*mat2 - mat1/0.5 = \n" << 2*mat2 - mat1/0.5 << std::endl;

    mx::Matrix mat3( mx::Rand(10) );
    std::cout << "rand matrix = \n" << mat3 << std::endl;

    mx::Matrix mat6 = mx::RandSPD(10);
    std::cout << "rand SPD matrix = \n" << mat6 << std::endl;

    return 0;
}
