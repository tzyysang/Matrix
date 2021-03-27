
namespace mx
{

typedef struct MatrixInitilizer
{
    enum MxInitType{
        ZEROS,
        EYE,
        RAND,
        RAND_SYM,
        RAND_LOWTRI,
        RAND_SPD
    } matrix_type;
    int size;

    MatrixInitilizer( MxInitType type, int s )
    {
        matrix_type = type;
        size = s;
    }
} MatInit;

struct Zeros : public MatrixInitilizer
{
    Zeros( int s ) : MatrixInitilizer(ZEROS, s) {};
};

struct Eye : public MatrixInitilizer
{
    Eye( int s ) : MatrixInitilizer(EYE, s) {};
};

struct Rand : public MatrixInitilizer
{
    Rand( int s ) : MatrixInitilizer(RAND, s) {};
};

struct RandSPD : public MatrixInitilizer
{
    RandSPD( int s ) : MatrixInitilizer(RAND_SPD, s) {};
};

}
