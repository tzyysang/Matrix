
#include "matrix.h"

int main( int argc, char* argv[] )
{
    std::cout << "hello world" << std::endl;

    mx::Matrix mat1 = {1, 2, 3};
    mx::Matrix mat2 = { {1, 4, 7},
                        {2, 5, 8},
                        {3, 6, 10}};
    mx::Matrix perm = { {0, 1, 0},
                        {0, 0, 1},
                        {1, 0, 0}};
    std::cout << "mat1 = \n" << mat1 << std::endl;
    std::cout << "mat2 = \n" << mat2 << std::endl;

    mx::Matrix mat3 = mx::Eye(3);
    std::cout << "mat3 = \n" << mat3 << std::endl;

    mx::Matrix mat4 = mat2 - mx::Eye(3);
    std::cout << "mat4 = \n" << mat4 << std::endl;

    std::cout << "mat4 = \n" << perm * mat4 << std::endl;

}
