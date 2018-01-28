#ifndef MYMATRIX_H
#define MYMATRIX_H
#include "MyFraction.h"
#include <typeinfo>
using namespace std;

template <typename type>



class MyMatrix
{

    private:
        int width,height;

    public:
        type *matrix;

        MyMatrix(){width=height=0;}

        MyMatrix(int x,int y){
            width=x;height=y;
            matrix= new type[x*y];
            for(int i=0;i<x*y;i++) matrix[i]=0;
        }

        MyMatrix(int x){width=height=x;matrix= new type[x*x];for(int i=0;i<x*x;i++) matrix[i]=0;}

        int get_width() const {return width;}

        int get_height() const {return height;}

        void set_values(type a[]){
                for(int i=0;i<height*width;i++)
                    matrix[i]=a[i];
        }

        MyMatrix<type> operator+ (const MyMatrix& b){
            MyMatrix result(width,height);
              for (int i=0; i<width; i++) {
                for (int j=0; j<height; j++) {
                  result.matrix[i*width+j] = matrix[i*width+j] + b.matrix[i*width+j];
                }
            }
            return result;
        }

        MyMatrix<type> operator* (const MyMatrix& b){
                    int rows = b.get_height();
                    int cols = b.get_width();
                    MyMatrix result(rows,cols);
                    for (int i=0; i<rows; i++) {
                        for (int j=0; j<cols; j++) {
                          for (int k=0; k<rows; k++) {
                            result.matrix[i*cols+j] += matrix[i*width+k] * b.matrix[k*cols+j];
                          }
                        }
                    }
                    return result;
        }

        MyMatrix<type>& operator= (const MyMatrix& b){
            this->set_values(b.matrix);
            return * this;
        }

    //------------------------------------------------------------------------------------------------------------GAUSS
     MyMatrix<type> Gauss_full() {
            int n = height;
            for (int i = 0; i < n; i++) {

                // znajd wiersz z maksymalnym elementem
                type maxEl = matrix[i*width*i];
                int maxRow = i;
                int maxCol = i;
                for (int k = i+1; k < n; k++) {
                    if (matrix[k*width+i] > maxEl) {
                        maxEl = matrix[k*width+i];
                        maxRow = k;
                        maxCol = i;
                    }
                }
                // zamieñ maksymalny wiersz z obecnym
                for (int k = i; k < n+1; k++) {
                    type pom = matrix[maxRow*width+k];
                    matrix[maxRow*width+k] = matrix[i*width+k];
                    matrix[i*width+k] = pom;
                }
                // zamien maksymalna kolumne z obecna
                for (int k = i; k < n; k++) {
                    type pom = matrix[k*width+maxCol];
                    matrix[k*width + maxCol] = matrix[k*width+i];
                    matrix[k*width+i] = pom;
                }

                // wyprowad zera przed obecnym wierszem
                for (int k = i+1; k <= n; k++) {
                    type c = (-matrix[k*width+i]) / matrix[i*width+i];
                    for (int j = i; j <= n; j++) {
                        if (i == j) {
                            matrix[k*width+j] = 0;
                        } else {
                            matrix[k*width+j] +=  c * matrix[i*width+j];
                        }
                    }
                }
            }
            // rozwi¹¿ Ax = B za pomoc¹ powsta³ej macierzy trójk¹tnej
            MyMatrix<type> x(1,n);
            for (int i=n-1; i>=0; i--) {
                x.matrix[i] = matrix[i*width+n] / matrix[i*width+i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k*width+n] = matrix[k*width+n] - (matrix[k*width+i] * x.matrix[i]);
                }
            }
            return x;
        }
//------------------------------------------------------------------------------------------------------------GAUSS


};



#endif // MYMATRIX_H

