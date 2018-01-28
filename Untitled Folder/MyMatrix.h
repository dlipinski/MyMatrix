#ifndef MYMATRIX_H
#define MYMATRIX_H
#include "MyFraction.h"

using namespace std;

class MyMatrix
{
    private:
        int width,height;
    public:
        MyFraction *matrix;
        MyMatrix(){width=height=0;}
        MyMatrix(int x,int y){width=x;height=y;matrix= new MyFraction[x*y];}
        MyMatrix(int x){width=height=x;matrix= new MyFraction[x*x];}
        void set_size(int i,int j){width=i;height=j;matrix= new MyFraction[i*j];}
        int get_width() const {return width;}
        int get_height() const {return height;}
        void set_values(MyFraction a[]){
            for(int i=0;i<height*width;i++)
                matrix[i]=a[i];
        }
        void zero_values(){
            for(int i=0;i<width*height;i++)
                    matrix[i].set_values(0,1);
         }
        void print_values(){
            for(int i=0;i<width*height;i++){
                if((i+1)%width==0){
                    matrix[i].print_value();
                    cout << endl;
                }
                else{
                    matrix[i].print_value();
                    cout << " ";
                }
            }
        }
        void print_double_values(){
            for(int i=0;i<width*height;i++){
                if((i+1)%width!=0){
                    matrix[i].print_double_value();
                    cout << endl;
                }
                else{
                    matrix[i].print_double_value();
                    cout << " ";
                }
            }
        }
        MyFraction& at (int i1, int i2) {return matrix[i1 * this->get_width() + i2];}

        MyMatrix swap_line(int a, int b){
            MyMatrix myMatrix(width);
            for(int i=0;i<width*width;i++)
                myMatrix.matrix[i].set_values(matrix[i].get_numerator(),matrix[i].get_denuminator());
            MyFraction line_a[width];
            MyFraction line_b[width];
            for(int i=0;i<width;i++){
                line_a[i].set_values(matrix[a*width+i].get_numerator(),matrix[a*width+i].get_denuminator());
                line_b[i].set_values(matrix[b*width+i].get_numerator(),matrix[b*width+i].get_denuminator());
            }
            for(int i=0;i<width;i++){
                myMatrix.at(b,i).set_values(line_a[i].get_numerator(),line_a[i].get_denuminator());
                myMatrix.at(a,i).set_values(line_b[i].get_numerator(),line_b[i].get_denuminator());
            }
            return myMatrix;
        }

        MyMatrix operator+ (const MyMatrix& b){
            MyMatrix myMatrix(this->get_width(),this->get_height());
            for(int i=0;i<this->get_height();i++){
                for(int j=0;j<this->get_width();j++){
                    myMatrix.at(i,j)
                    =
                    matrix[i * get_width() + j]
                    +
                    b.matrix[i * b.get_width() + j];
                }
            }
            return myMatrix;
        }
        MyMatrix operator- (const MyMatrix& b){
            MyMatrix myMatrix(this->get_width(),this->get_height());
            for(int i=0;i<this->get_height();i++){
                for(int j=0;j<this->get_width();j++){
                    myMatrix.at(i,j)
                    =
                    matrix[i * get_width() + j]
                    -
                    b.matrix[i * b.get_width() + j];
                }
            }
            return myMatrix;
        }
        MyMatrix operator* (const long int& b){
            MyMatrix myMatrix(this->get_width(),this->get_height());
            for(int i=0;i<this->get_height();i++){
                for(int j=0;j<this->get_width();j++){
                    myMatrix.at(i,j)
                    =
                    matrix[i * get_width() + j]
                    *
                    b;
                }
            }
            return myMatrix;
        }
        MyMatrix operator* (const MyMatrix& b){


            if(b.get_width()!=1){
                MyMatrix myMatrix(width,height);
                int k=0;
                int size = myMatrix.get_width();
                    for(int i=0;i<height;i++){
                        for(int j=0;j<width;j++){

                                    myMatrix.matrix[i * size + j]
                                    =
                                    myMatrix.matrix[i * size + j]
                                    +
                                    matrix[size*i+j] * b.matrix[size*k+j];
                                }

                    }
                return myMatrix;
            }
            else{
                MyMatrix myMatrix(1,height);
                for(int i=0;i<b.get_height();i++){
                    for(int j=0;j<width;j++){
                            myMatrix.matrix[i]
                            =
                            myMatrix.matrix[i]
                            +
                            matrix[height*i+j] * b.matrix[j];
                        }
                    }
                return myMatrix;
            }
        }




        void operator= (const MyMatrix& b){
            set_values(b.matrix);

        }



    protected:


};

#endif // MYMATRIX_H
