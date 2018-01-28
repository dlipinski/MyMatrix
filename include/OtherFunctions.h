#ifndef OTHERFUNCTIONS_H
#define OTHERFUNCTIONS_H
#include <stdlib.h>
#include <ctime>
class OtherFunctions
{
    public:
        OtherFunctions();
        virtual ~OtherFunctions();

    protected:

    private:
};
unsigned long long int dRand(unsigned long long int  fMin, unsigned long long int  fMax)
{
    srand(time(NULL));
    unsigned long long int  f = rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

unsigned long long int myRand()
{
    srand(time(NULL));
    return (int)rand() % 100 + 1;
}


#endif // OTHERFUNCTIONS_H
