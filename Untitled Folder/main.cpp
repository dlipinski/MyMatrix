#include <iostream>
#include <numeric>
#include "MyFraction.h"
#include "MyMatrix.h"
#include "OtherFunctions.h"
#include <ctime>
#include <Eigen/Dense>
#include <fstream>
#include <cmath>

using namespace std;
using namespace Eigen;


int main()
{
    cout << "> START PROGRAM" << endl;
    int matrix_size = 3;
    int matrix_amount = 3;
    int min_value = 60;
    int max_value = 50;
    int loops = 10000;

    //create and open all needed text files
    ofstream time_file;
    time_file.open("time_file.txt");
    ofstream error_file;
    error_file.open("error_file.txt");
//------------------------------------------------------------------------------------------------------------
    //Crate array and generate values for MyFraction matrix
    cout << "> CREATE AND GETERENTE INTO DATA ARRAY" << endl;
    MyFraction myFraction_data[matrix_amount+1][matrix_size*matrix_size];
    double double_data[matrix_amount+1][matrix_size*matrix_size];
    float float_data[matrix_amount+1][matrix_size*matrix_size];
    for(int i=0;i<matrix_amount+1;i++)
        for(int j=0;j<matrix_size*matrix_size;j++){
                myFraction_data[i][j].set_values(dRand(min_value,max_value),dRand(min_value,max_value));
                double_data[i][j]=myFraction_data[i][j].get_double_value();
                float_data[i][j]=myFraction_data[i][j].get_float_value();
        }
//------------------------------------------------------------------------------------------------------------- MYMATRIX VALUES AND OPERATIONS
    cout << "> CREATE MYMATRIX OBJECTS" << endl;
    MyMatrix matrix_a(matrix_size);
    MyMatrix matrix_b(matrix_size);
    MyMatrix matrix_c(matrix_size);
    MyMatrix vector_x(1,matrix_size);

    MyMatrix matrix_answer_1(1,matrix_size);
    MyMatrix matrix_answer_2(1,matrix_size);
    MyMatrix matrix_answer_3(matrix_size);

    cout << "> PUT DATA IN MYMATRIX OBJECTS" << endl;
    matrix_a.set_values(myFraction_data[0]);
    matrix_b.set_values(myFraction_data[1]);
    matrix_c.set_values(myFraction_data[2]);
    vector_x.set_values(myFraction_data[3]);

    cout << "> MYMATRIX OPERATIONS" << endl;
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
    time_file <<fixed<< "MyMatrix times:" << endl << myMatrix_time_1 << endl << myMatrix_time_2 << endl << myMatrix_time_3 << endl;
//--------------------------------------------------------------------------------------------------------- EIGEN VALUES AND OPERATIONS
    cout << "> CREATE EIGEN OBJECTS" << endl;
    //declare matrixe and vextors
    MatrixXd ad(matrix_size,matrix_size);
    MatrixXd bd(matrix_size,matrix_size);
    MatrixXd cd(matrix_size,matrix_size);
    VectorXd vd(matrix_size);

    MatrixXf af(matrix_size,matrix_size);
    MatrixXf bf(matrix_size,matrix_size);
    MatrixXf cf(matrix_size,matrix_size);
    VectorXf vf(matrix_size);

    VectorXd ansd_1(matrix_size);
    VectorXd ansd_2(matrix_size);
    MatrixXd ansd_3(matrix_size,matrix_size);

    VectorXf ansf_1(matrix_size);
    VectorXf ansf_2(matrix_size);
    MatrixXf ansf_3(matrix_size,matrix_size);
    cout << "> PUT DATA IN EIGEN OBJECTS" << endl;
    //initiate
    for(int i=0;i<matrix_size;i++){
        for(int j=0;j<matrix_size;j++){
            ad(i,j) = matrix_a.at(i,j).get_double_value();
            bd(i,j) = matrix_b.at(i,j).get_double_value();
            cd(i,j) = matrix_c.at(i,j).get_double_value();

            af(i,j) = matrix_a.at(i,j).get_float_value();
            bf(i,j) = matrix_b.at(i,j).get_float_value();
            cf(i,j) = matrix_c.at(i,j).get_float_value();
        }
        vd(i)=vector_x.at(i,0).get_double_value();
        vf(i)=vector_x.at(i,0).get_float_value();
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
    time_file << "Eigen double times:" << endl << eigen_double_time_1 << endl << eigen_double_time_2 << endl << eigen_double_time_3 << endl;
    cout << "> EIGEN FLOAT OPERATIONS" << endl;
    //measure time for eigen float values
    clock_t S_eigen_float_time_1 = clock();
        for(int i=0;i<loops;i++)
            ansf_1 = af * vf;
    clock_t E_eigen_float_time_1 = clock();

    clock_t S_eigen_float_time_2 = clock();
        for(int i=0;i<loops;i++)
            ansf_2 = (af+bf+cf) * vf;
    clock_t E_eigen_float_time_2 = clock();

    clock_t S_eigen_float_time_3 = clock();
        for(int i=0;i<loops;i++)
            ansf_3 = af * ( bf * cf );
    clock_t E_eigen_float_time_3 = clock();

    double eigen_float_time_1 = E_eigen_float_time_1-S_eigen_float_time_1;
    double eigen_float_time_2 = E_eigen_float_time_2-S_eigen_float_time_2;
    double eigen_float_time_3 = E_eigen_float_time_3-S_eigen_float_time_3;
    time_file << "Eigen float times:" << endl << eigen_float_time_1 << endl << eigen_float_time_2 << endl << eigen_float_time_3 << endl;
    time_file.close();
//------------------------------------------------------------------------------------------- CHECK VALUES
    cout << "> CREATE ERRORS ARRAYS" << endl;
    double double_errors_1[matrix_size];
    double double_errors_2[matrix_size];
    double double_errors_3[matrix_size][matrix_size];

    float float_errors_1[matrix_size];
    float float_errors_2[matrix_size];
    float float_errors_3[matrix_size][matrix_size];

    cout << "> ERRORS OPERATIONS" << endl;
    for(int i=0;i<matrix_size;i++){
            double_errors_1[i] = abs(ansd_1(i) - matrix_answer_1.at(i,0).get_double_value());
            double_errors_2[i] = abs(ansd_2(i) - matrix_answer_2.at(i,0).get_double_value());
            float_errors_1[i] = abs(ansf_1(i) - matrix_answer_1.at(i,0).get_float_value());
            float_errors_2[i] = abs(ansf_2(i) - matrix_answer_2.at(i,0).get_float_value());
        for(int j=0;j<matrix_size;j++){
            double_errors_3[i][j] = abs(ansd_3(i,j) - matrix_answer_3.at(i,j).get_double_value());
            float_errors_3[i][j] = abs(ansf_3(i,j) - matrix_answer_3.at(i,j).get_float_value());
        }
    }
    double max_double_error_1=0;
    double max_double_error_2=0;
    double max_double_error_3=0;

    double avg_double_error_1=0;
    double avg_double_error_2=0;
    double avg_double_error_3=0;

    float max_float_error_1=0;
    float max_float_error_2=0;
    float max_float_error_3=0;

    float avg_float_error_1=0;
    float avg_float_error_2=0;
    float avg_float_error_3=0;

    for(int i=0;i<matrix_size;i++){
        if(double_errors_1[i]>max_double_error_1) max_double_error_1 = double_errors_1[i];
        if(double_errors_2[i]>max_double_error_2) max_double_error_2 = double_errors_2[i];
        if(float_errors_1[i]>max_float_error_1) max_float_error_1 = float_errors_1[i];
        if(float_errors_2[i]>max_float_error_2) max_float_error_2 = float_errors_2[i];
        avg_double_error_1+=double_errors_1[i];
        avg_double_error_2+=double_errors_2[i];
        avg_float_error_1+=float_errors_1[i];
        avg_float_error_2+=float_errors_2[i];
        for(int j=0;j<matrix_size;j++){
            if(double_errors_3[i][j]>max_double_error_3) max_double_error_3 = double_errors_3[i][j];
            if(float_errors_3[i][j]>max_float_error_3) max_float_error_3 = float_errors_3[i][j];
            avg_double_error_3+=double_errors_3[i][j];
            avg_float_error_3+=float_errors_3[i][j];
        }
    }
    avg_double_error_1/=matrix_size;
    avg_double_error_2/=matrix_size;
    avg_double_error_3/=(matrix_size*matrix_size);
    avg_float_error_1/=matrix_size;
    avg_float_error_2/=matrix_size;
    avg_float_error_3/=(matrix_size*matrix_size);
    error_file << fixed << "Errors:" << endl;
    error_file << "Max double errors:" << endl;
    error_file << fixed << max_double_error_1 << endl << max_double_error_2 << endl << max_double_error_3 << endl;
    error_file << "Average double errors:" << endl;
    error_file << fixed << avg_double_error_1 << endl << avg_double_error_2 << endl << avg_double_error_3 << endl;
    error_file << "Max float errors:" << endl;
    error_file << fixed << max_float_error_1 << endl << max_float_error_2 << endl << max_float_error_3 << endl;
    error_file << "Average float errors" << endl;
    error_file << fixed << avg_float_error_1 << endl << avg_float_error_2 << endl <<avg_float_error_3 << endl;
    error_file.close();
    cout << "> TESTS" << endl;
    return 0;
}
