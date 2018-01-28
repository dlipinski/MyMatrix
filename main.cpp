#include <iostream>
#include <numeric>
#include "MyMatrix.h"
#include "OtherFunctions.h"
#include <ctime>
#include <fstream>
#include <cmath>
#include <ctime>
#include <Eigen/Dense>
#include <iomanip>
#include<math.h>
#include <vector>
#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <cstdlib>
#include <Eigen/Dense>
#include <iomanip>
using namespace std;

using namespace Eigen;


int _main()
{
    cout << "> START PROGRAM" << endl;
    MyFraction s(0,1);
    int matrix_size = 150;
    int matrix_amount = 3;
    int loops = 1;

    //create and open all needed text files
    ofstream time_file;
    time_file.open("time_file.txt");
    time_file << "A*X,(A+B+C)*X,A*(B*C)" <<endl<< endl;
    ofstream error_file;
    error_file.open("error_file.txt");
    error_file.precision(25);
    fstream myfile("data2.txt", ios_base::in);
    ofstream test_file;
    test_file.open("test_file.txt");
//------------------------------------------------------------------------------------------------------------
    //Crate array and generate values for MyFraction matrix
    cout << "> READ DATA FROM FILE TO INTO DATA ARRAY" << endl;
    MyFraction myFraction_data[matrix_amount+1][matrix_size*matrix_size];
    double double_data[matrix_amount+1][matrix_size*matrix_size];
    float float_data[matrix_amount+1][matrix_size*matrix_size];
    int a,b;
    for(int i=0;i<matrix_amount+1;i++)
        for(int j=0;j<(matrix_size)*matrix_size;j++){
                myfile>>a;
                myfile>>b;
                myFraction_data[i][j].set_values(a,b);
                double_data[i][j]=(double)a/b;
                float_data[i][j]=(float)a/b;
        }
        myfile.close();
    fstream my_gauss_file("data2.txt", ios_base::in);


//------------------------------------------------------------------------------------------------------------- MYMATRIX VALUES AND OPERATIONS
//--MYFRACTION
    cout << "> CREATE MYMATRIX MYFRACTION OBJECTS" << endl;
    MyMatrix<MyFraction> matrix_a(matrix_size,matrix_size,s);
    MyMatrix<MyFraction> matrix_b(matrix_size,matrix_size,s);
    MyMatrix<MyFraction> matrix_c(matrix_size,matrix_size,s);
    vector<MyFraction> vector_x(matrix_size);

    vector<MyFraction> matrix_answer_1(matrix_size);
    vector<MyFraction> matrix_answer_2(matrix_size);
    MyMatrix<MyFraction> matrix_answer_3(matrix_size,matrix_size,s);

    cout << "> PUT DATA IN MYMATRIX MYFRACTION OBJECTS" << endl;
    matrix_a=myFraction_data[0];
    matrix_b=myFraction_data[1];
    matrix_c=myFraction_data[2];
    vector_x.insert(vector_x.begin(),myFraction_data[3],myFraction_data[3]+matrix_size);;
    cout << "> MYMATRIX MYFRACTION OPERATIONS" << endl;
    //A*X
    clock_t S_myMatrix_time_1 = clock();
        for(int i=0;i<loops;i++)
            matrix_answer_1=matrix_a*vector_x;
    clock_t E_myMatrix_time_1 = clock();
    //(A+B+C)*X
    clock_t S_myMatrix_time_2 = clock();
         for(int i=0;i<loops;i++)
            matrix_answer_2 = (matrix_a+matrix_b+matrix_c) * vector_x;
    clock_t E_myMatrix_time_2 = clock();
    //A*(B*C)
    clock_t S_myMatrix_time_3 = clock();
        for(int i=0;i<loops;i++)
            matrix_answer_3=matrix_a*(matrix_b*matrix_c);
    clock_t E_myMatrix_time_3 = clock();

    double myMatrix_time_1 = E_myMatrix_time_1 - S_myMatrix_time_1;
    double myMatrix_time_2 = E_myMatrix_time_2 - S_myMatrix_time_2;
    double myMatrix_time_3 = E_myMatrix_time_3 - S_myMatrix_time_3;
    time_file <<fixed<< "MyMatrix MyFraction times:" << endl << myMatrix_time_1 << endl << myMatrix_time_2 << endl << myMatrix_time_3 << endl<< endl;
//--DOUBLE

    cout << "> CREATE MYMATRIX DOUBLE OBJECTS" << endl;
    MyMatrix<double> matrixd_a(matrix_size,matrix_size+1,1.0);
    MyMatrix<double> matrixd_b(matrix_size,matrix_size+1,1.0);
    MyMatrix<double> matrixd_c(matrix_size,matrix_size+1,1.0);
    vector<double> vectord_x(matrix_size);

    vector<double> matrixd_answer_1(matrix_size);
    vector<double> matrixd_answer_2(matrix_size);
    MyMatrix<double> matrixd_answer_3(matrix_size,matrix_size+1,1.0);

    cout << "> PUT DATA IN MYMATRIX DOUBLE OBJECTS" << endl;
    matrixd_a=double_data[0];
    matrixd_b=double_data[1];
    matrixd_c=double_data[2];
    vectord_x.insert(vectord_x.begin(),double_data[3],double_data[3]+matrix_size);;

    cout << "> MYMATRIX DOUBLE OPERATIONS" << endl;
    //A*X
    clock_t S_myMatrixd_time_1 = clock();
        for(int i=0;i<loops;i++)
            matrixd_answer_1=matrixd_a*vectord_x;
    clock_t E_myMatrixd_time_1 = clock();
    //(A+B+C)*X
    clock_t S_myMatrixd_time_2 = clock();
         for(int i=0;i<loops;i++)
            matrixd_answer_2 = (matrixd_a+matrixd_b+matrixd_c) * vectord_x;
    clock_t E_myMatrixd_time_2 = clock();
    //A*(B*C)
    clock_t S_myMatrixd_time_3 = clock();
        for(int i=0;i<loops;i++)
            matrixd_answer_3=matrixd_a*(matrixd_b*matrixd_c);
    clock_t E_myMatrixd_time_3 = clock();

    double myMatrixd_time_1 = E_myMatrixd_time_1 - S_myMatrixd_time_1;
    double myMatrixd_time_2 = E_myMatrixd_time_2 - S_myMatrixd_time_2;
    double myMatrixd_time_3 = E_myMatrixd_time_3 - S_myMatrixd_time_3;
    time_file <<fixed<< "MyMatrix double times:" << endl << myMatrixd_time_1 << endl << myMatrixd_time_2 << endl << myMatrixd_time_3 << endl<< endl;

//--FLOAT

    cout << "> CREATE MYMATRIX FLOAT OBJECTS" << endl;
    MyMatrix<float> matrixf_a(matrix_size,matrix_size+1,1.0);
    MyMatrix<float> matrixf_b(matrix_size,matrix_size+1,1.0);
    MyMatrix<float> matrixf_c(matrix_size,matrix_size+1,1.0);
    vector<float> vectorf_x(matrix_size);

    vector<float> matrixf_answer_1(matrix_size);
    vector<float> matrixf_answer_2(matrix_size);
    MyMatrix<float> matrixf_answer_3(matrix_size,matrix_size+1,1.0);

    cout << "> PUT DATA IN MYMATRIX FLOAT OBJECTS" << endl;
    matrixf_a=float_data[0];
    matrixf_b=float_data[1];
    matrixf_c=float_data[2];
    vectorf_x.insert(vectorf_x.begin(),float_data[3],float_data[3]+matrix_size);;

    cout << "> MYMATRIX FLOAT OPERATIONS" << endl;
    //A*X
    clock_t S_myMatrixf_time_1 = clock();
        for(int i=0;i<loops;i++)
            matrixf_answer_1=matrixf_a*vectorf_x;
    clock_t E_myMatrixf_time_1 = clock();
    //(A+B+C)*X
    clock_t S_myMatrixf_time_2 = clock();
         for(int i=0;i<loops;i++)
            matrixf_answer_2 = (matrixf_a+matrixf_b+matrixf_c) * vectorf_x;
    clock_t E_myMatrixf_time_2 = clock();
    //A*(B*C)
    clock_t S_myMatrixf_time_3 = clock();
        for(int i=0;i<loops;i++)
            matrixf_answer_3=matrixf_a*(matrixf_b*matrixf_c);
    clock_t E_myMatrixf_time_3 = clock();

    double myMatrixf_time_1 = E_myMatrixf_time_1 - S_myMatrixf_time_1;
    double myMatrixf_time_2 = E_myMatrixf_time_2 - S_myMatrixf_time_2;
    double myMatrixf_time_3 = E_myMatrixf_time_3 - S_myMatrixf_time_3;
    time_file <<fixed<< "MyMatrix float times:" << endl << myMatrixf_time_1 << endl << myMatrixf_time_2 << endl << myMatrixf_time_3 << endl<< endl;

//--------------------------------------------------------------------------------------------------------- EIGEN VALUES AND OPERATIONS
    cout << "> CREATE EIGEN OBJECTS" << endl;
    //declare matrixe and vextors
    MatrixXd ad(matrix_size,matrix_size);
    MatrixXd bd(matrix_size,matrix_size);
    MatrixXd cd(matrix_size,matrix_size);
    VectorXd vd(matrix_size);

    VectorXd ansd_1(matrix_size);
    VectorXd ansd_2(matrix_size);
    MatrixXd ansd_3(matrix_size,matrix_size);

    cout << "> PUT DATA IN EIGEN OBJECTS" << endl;
    //initiate
    for(int i=0;i<matrix_size;i++){
        for(int j=0;j<matrix_size;j++){
            ad(i,j) = double_data[0][i*matrix_size+j];
            bd(i,j) = double_data[1][i*matrix_size+j];
            cd(i,j) = double_data[2][i*matrix_size+j];
        }
        vd(i)=double_data[3][i];
        }
     cout << "> EIGEN DOUBLE OPERATIONS" << endl;
    //measure time for eigen double operations
    clock_t S_eigen_double_time_1 = clock();
        for(int i=0;i<loops;i++)
            ansd_1 = ad * vd;
    clock_t E_eigen_double_time_1 = clock();

    clock_t S_eigen_double_time_2 = clock();
        for(int i=0;i<loops;i++)
            ansd_2 = (ad+bd+cd) * vd;
    clock_t E_eigen_double_time_2 = clock();

    clock_t S_eigen_double_time_3 = clock();
        for(int i=0;i<loops;i++)
            ansd_3 = ad * ( bd * cd );
    clock_t E_eigen_double_time_3 = clock();

    double eigen_double_time_1 = E_eigen_double_time_1-S_eigen_double_time_1;
    double eigen_double_time_2 = E_eigen_double_time_2-S_eigen_double_time_2;
    double eigen_double_time_3 = E_eigen_double_time_3-S_eigen_double_time_3;
    time_file << "Eigen double times:" << endl << eigen_double_time_1 << endl << eigen_double_time_2 << endl << eigen_double_time_3 << endl<< endl;


//-------GAUSS OPERATIONS
cout << ">GAUSS OPERATIONS"<<endl;
VectorXd eigen_inside_v(matrix_size);
for(int i=0;i<matrix_size;i++)
    eigen_inside_v[i]=matrixd_a[i][matrix_size];//wslozenie ostatniej kolumny do wynikow dla eigena
//eigen
  clock_t s_gauss_eigen_partial_1 = clock();
      VectorXd eigen_partial_result =  ad.partialPivLu().solve(eigen_inside_v);
  clock_t e_gauss_eigen_partial_1 = clock();
//my_gauss_full
    MyMatrix<double> test1 = matrixd_a;
    MyMatrix<float> test2 = matrixf_a;
    clock_t s_gauss_double_full_1 = clock();
        vector<double> double_full_result = matrixd_a.Gauss_full();
    clock_t e_gauss_double_full_1 = clock();

    clock_t s_gauss_float_full_1 = clock();
        vector<float> float_full_result = matrixf_a.Gauss_full();
    clock_t e_gauss_float_full_1 = clock();

  clock_t s_gauss_double_partial_1 = clock();
      vector<double> double_partial_result = test1.Gauss_partial();
  clock_t e_gauss_double_partial_1 = clock();

  clock_t s_gauss_float_partial_1 = clock();
      vector<float> float_partial_result = test2.Gauss_partial();
  clock_t e_gauss_float_partial_1 = clock();

  double gauss_partial_eigen_time_1 = e_gauss_eigen_partial_1 - s_gauss_eigen_partial_1;

  double gauss_full_double_time_1 = e_gauss_double_full_1 - s_gauss_double_full_1;
  double gauss_full_float_time_1 = e_gauss_float_full_1-s_gauss_float_full_1;
  double gauss_partial_double_time_1 = e_gauss_double_partial_1 - s_gauss_double_partial_1;
  double gauss_partial_float_time_1 = e_gauss_float_partial_1 - s_gauss_float_partial_1;
 time_file <<fixed<< endl<<endl<<"Gauss times" << endl;
 time_file <<fixed<<"Eigen: "<<gauss_partial_eigen_time_1<<endl;

time_file <<fixed<< endl<<endl<<"FullGaussTimes" << endl;
 time_file <<fixed<<"Double: "<<gauss_full_double_time_1<<endl;
 time_file <<fixed<<"Float: "<<gauss_full_float_time_1<<endl;


time_file <<fixed<< endl<<endl<<"PartialGaussTimes" << endl;
 time_file <<fixed<<"Double: "<<gauss_partial_double_time_1<<endl;
 time_file <<fixed<<"Float: "<<gauss_partial_float_time_1<<endl;

//full gaus errors array
double gauss_full_double_errors[matrix_size];
double gauss_full_float_errors[matrix_size];
//partial gauss errors array
double gauss_partial_double_errors[matrix_size];
double gauss_partial_float_errors[matrix_size];
//max errors
double gauss_full_max_double_error =0;
double gauss_full_max_float_error =0;
double gauss_partial_max_double_error =0;
double gauss_partial_max_float_error =0;
//avg errors
double gauss_full_avg_double_error =0;
double gauss_full_avg_float_error =0;
double gauss_partial_avg_double_error =0;
double gauss_partial_avg_float_error =0;

for(int i;i<matrix_size;i++){
    gauss_full_double_errors[i]=abs(eigen_partial_result[i]-double_full_result[i]);
    gauss_full_float_errors[i]=abs(eigen_partial_result[i]-float_full_result[i]);

    gauss_partial_double_errors[i]=abs(eigen_partial_result[i]-double_partial_result[i]);
    gauss_partial_float_errors[i]=abs(eigen_partial_result[i]-float_partial_result[i]);

    if(gauss_full_double_errors[i]>gauss_full_max_double_error)gauss_full_max_double_error=gauss_full_double_errors[i];
    if(gauss_full_float_errors[i]>gauss_full_max_float_error)gauss_full_max_float_error=gauss_full_float_errors[i];
    if(gauss_partial_double_errors[i]>gauss_partial_max_double_error)gauss_partial_max_double_error=gauss_partial_double_errors[i];
    if(gauss_partial_float_errors[i]>gauss_partial_max_float_error)gauss_partial_max_float_error=gauss_partial_float_errors[i];

    gauss_full_avg_double_error+=gauss_full_double_errors[i];
    gauss_full_avg_float_error+=gauss_full_float_errors[i];
    gauss_partial_avg_double_error+=gauss_partial_double_errors[i];
    gauss_partial_avg_float_error+=gauss_partial_float_errors[i];

}
gauss_full_avg_double_error/=matrix_size;
gauss_full_avg_float_error/=matrix_size;

gauss_partial_avg_double_error/=matrix_size;
gauss_partial_avg_float_error/=matrix_size;



error_file << fixed << "Errors:" << endl;

error_file << fixed << "Full gauss double errors:" << endl;
error_file << fixed << "Max: " << gauss_full_max_double_error<< endl;
error_file << fixed << "Avg: " << gauss_full_avg_double_error<< endl<<endl;
error_file << fixed << "Full gauss float errors:" << endl;
error_file << fixed << "Max: " << gauss_full_max_float_error<< endl;
error_file << fixed << "Avg: " << gauss_full_avg_float_error<< endl<<endl;

error_file << fixed << "Partial gauss double errors:" << endl;
error_file << fixed << "Max: " << gauss_partial_max_double_error<< endl;
error_file << fixed << "Avg: " << gauss_partial_avg_double_error<< endl<<endl;
error_file << fixed << "Partial gauss float errors:" << endl;
error_file << fixed << "Max: " << gauss_partial_max_float_error<< endl;
error_file << fixed << "Avg: " << gauss_partial_avg_float_error<< endl<<endl;

//OPERATIONS TESTING
error_file << "A*X,(A+B+C)*X,A*(B*C)"<<endl;
    cout << "> CREATE ERRORS ARRAYS" << endl;
    double MyFraction_errors_1[matrix_size];
    double MyFraction_errors_2[matrix_size];
    double MyFraction_errors_3[matrix_size][matrix_size];

    double double_errors_1[matrix_size];
    double double_errors_2[matrix_size];
    double double_errors_3[matrix_size][matrix_size];

    double float_errors_1[matrix_size];
    double float_errors_2[matrix_size];
    double float_errors_3[matrix_size][matrix_size];

    cout << "> ERRORS OPERATIONS" << endl;
    for(int i=0;i<matrix_size;i++){
            MyFraction_errors_1[i] = abs(ansd_1[i] - matrix_answer_1[i].get_double_value_big());
            MyFraction_errors_2[i] = abs(ansd_2[i] - matrix_answer_2[i].get_double_value_big());


            double_errors_1[i] = abs(ansd_1[i] - matrixd_answer_1[i]);
            double_errors_2[i] = abs(ansd_2[i] - matrixd_answer_2[i]);

            float_errors_1[i] = abs(ansd_1[i] - matrixf_answer_1[i]);
            float_errors_2[i] = abs(ansd_2[i] - matrixf_answer_2[i]);
        for(int j=0;j<matrix_size;j++){
            MyFraction_errors_3[i][j] = abs(ansd_3(i,j) - matrix_answer_3[i][j].get_double_value_big());
            double_errors_3[i][j] = abs(ansd_3(i,j) - matrixd_answer_3[i][j]);
            float_errors_3[i][j] = abs(ansd_3(i,j) - matrixf_answer_3[i ][j ]);
        }
    }
    long double max_MyFraction_error_1=0;
    double max_MyFraction_error_2=0;
    double max_MyFraction_error_3=0;

    double avg_MyFraction_error_1=0;
    double avg_MyFraction_error_2=0;
    double avg_MyFraction_error_3=0;

    double max_double_error_1=0;
    double max_double_error_2=0;
    double max_double_error_3=0;

    double avg_double_error_1=0;
    double avg_double_error_2=0;
    double avg_double_error_3=0;

    double max_float_error_1=0;
    double max_float_error_2=0;
    double max_float_error_3=0;

    double avg_float_error_1=0;
    double avg_float_error_2=0;
    double avg_float_error_3=0;

    for(int i=0;i<matrix_size;i++){

        if(MyFraction_errors_1[i]>max_MyFraction_error_1) max_MyFraction_error_1 = MyFraction_errors_1[i];
        if(MyFraction_errors_2[i]>max_MyFraction_error_2) max_MyFraction_error_2 = MyFraction_errors_2[i];

        if(double_errors_1[i]>max_double_error_1) max_double_error_1 = double_errors_1[i];
        if(double_errors_2[i]>max_double_error_2) max_double_error_2 = double_errors_2[i];

        if(float_errors_1[i]>max_float_error_1) max_float_error_1 = float_errors_1[i];
        if(float_errors_2[i]>max_float_error_2) max_float_error_2 = float_errors_2[i];

        avg_MyFraction_error_1+=MyFraction_errors_1[i];
        avg_MyFraction_error_2+=MyFraction_errors_2[i];

        avg_double_error_1+=double_errors_1[i];
        avg_double_error_2+=double_errors_2[i];

        avg_float_error_1+=float_errors_1[i];
        avg_float_error_2+=float_errors_2[i];

        for(int j=0;j<matrix_size;j++){
            if(MyFraction_errors_3[i][j]>max_MyFraction_error_3) max_MyFraction_error_3 = MyFraction_errors_3[i][j];
            if(double_errors_3[i][j]>max_double_error_3) max_double_error_3 = double_errors_3[i][j];
            if(float_errors_3[i][j]>max_float_error_3) max_float_error_3 = float_errors_3[i][j];

            avg_MyFraction_error_3+=MyFraction_errors_3[i][j];
            avg_double_error_3+=double_errors_3[i][j];
            avg_float_error_3+=float_errors_3[i][j];
        }
    }
    avg_MyFraction_error_1/=matrix_size;
    avg_MyFraction_error_2/=matrix_size;
    avg_MyFraction_error_3/=(matrix_size*matrix_size);

    avg_double_error_1/=matrix_size;
    avg_double_error_2/=matrix_size;
    avg_double_error_3/=(matrix_size*matrix_size);

    avg_float_error_1/=matrix_size;
    avg_float_error_2/=matrix_size;
    avg_float_error_3/=(matrix_size*matrix_size);



    error_file << "Max MyFraction errors:" << endl;
    error_file << fixed << max_MyFraction_error_1 << endl << max_MyFraction_error_2 << endl << max_MyFraction_error_3 << endl<<endl;
    error_file << "Average MyFraction errors:" << endl;
    error_file << fixed << avg_MyFraction_error_1 << endl << avg_MyFraction_error_2 << endl << avg_MyFraction_error_3 << endl<<endl;
    error_file << "Max double errors:" << endl;
    error_file << fixed << max_double_error_1 << endl << max_double_error_2 << endl << max_double_error_3 << endl<<endl;
    error_file << "Average double errors:" << endl;
    error_file << fixed << avg_double_error_1 << endl << avg_double_error_2 << endl << avg_double_error_3 << endl<<endl;
    error_file << "Max float errors:" << endl;
    error_file << fixed << max_float_error_1 << endl << max_float_error_2 << endl << max_float_error_3 << endl<<endl;
    error_file << "Average float errors" << endl;
    error_file << fixed << avg_float_error_1 << endl << avg_float_error_2 << endl <<avg_float_error_3 << endl<<endl;
    error_file.close();
    test_file.close();
     time_file.close();
    cout << "> TESTS" << endl;



    return 0;
}
