#ifndef MYFRACTION_H
#define MYFRACTION_H
#include <iostream>

using namespace std;

long int gcd(long int x,long int y);
long int get_common_denominator(long int x, long int y);

class MyFraction
{
    public:
         MyFraction();
         MyFraction(long int x,long int y);
         void set_values(long int x,long int y){numerator = x;denuminator = y;}
         long int get_numerator () const {return numerator;};
         long int get_denuminator () const {return denuminator;};
         float get_float_value(){return (float)numerator/denuminator;}
         double get_double_value(){return (double)numerator/denuminator;}
         void print_value(){
            cout << numerator << "/" << denuminator;
         }
         void print_double_value(){
            cout << get_double_value();
         }

         MyFraction shortenFraction(MyFraction myFraction){
            long int _gcd = gcd(myFraction.get_numerator(),myFraction.get_denuminator());
            myFraction.set_values(myFraction.get_numerator()/_gcd,myFraction.get_denuminator()/_gcd);
            return myFraction;
        }
         MyFraction operator+ (const MyFraction& b){
            MyFraction myFraction;
            long int common_denuminator = get_common_denominator(denuminator,b.get_denuminator());
            myFraction.numerator =
                ((common_denuminator/this->denuminator)*(this->numerator))
                +
                ((common_denuminator/b.get_denuminator())*(b.get_numerator()));
            myFraction.denuminator = common_denuminator;
            return shortenFraction(myFraction);
            //return myFraction;
         }
         MyFraction operator- (const MyFraction& b){
            MyFraction myFraction;
            long int common_denuminator = get_common_denominator(denuminator,b.get_denuminator());
            myFraction.numerator =
                ((common_denuminator/this->denuminator)*(this->numerator))
                -
                ((common_denuminator/b.get_denuminator())*(b.get_numerator()));
            myFraction.denuminator = common_denuminator;
            return shortenFraction(myFraction);
            //return myFraction;
         }
         MyFraction operator* (const long int& b){
            MyFraction myFraction;
            myFraction.numerator = this->numerator * b;
            myFraction.denuminator = this->denuminator;
            return shortenFraction(myFraction);
            //return myFraction;
         }
         MyFraction operator* (const MyFraction& b){
            MyFraction myFraction;
            if(b.get_denuminator()!=0){
                myFraction.numerator = this->numerator * b.get_numerator();
                myFraction.denuminator = this->denuminator * b.get_denuminator();
            }
            else{
                myFraction.numerator = 0;
                myFraction.denuminator = 1;
            }
            return shortenFraction(myFraction);
            //return myFraction;
         }
         MyFraction operator/ (const MyFraction& b){
            MyFraction myFraction;
            myFraction.numerator = this->numerator * b.get_denuminator();
            myFraction.denuminator = this->denuminator * b.get_numerator();
            return shortenFraction(myFraction);
            //return myFraction;
         }
         //MyFraction operator= (const MyFraction&b){return b;}
    protected:

    private:
          long int numerator,denuminator;
};

MyFraction::MyFraction(){
    numerator  = 0;
    denuminator = 1;
};

MyFraction::MyFraction(long int x,long int y){
    numerator = x;
    denuminator = y;
}



long int gcd(long int x,long int y)
{
    for (;;)
    {
        if (x == 0) return y;
        y %= x;
        if (y == 0) return x;
        x %= y;
    }
}

long int get_common_denominator(long int x, long int y){
    long int temp = gcd(x, y);
    return temp ? (x / temp * y) : 0;
}



#endif // MYFRACTION_H
