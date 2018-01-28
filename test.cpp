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

int main(){
    cout.precision(2);
    int matrix_size = 4;
    fstream myfile("data_my.txt", ios_base::in);
    vector<vector<double> > double_data;
    double_data.resize(matrix_size);
      for (unsigned i=0; i<double_data.size(); i++) {
        double_data[i].resize(matrix_size+1, 0);
      }
    double a;
    for(int i=0;i<matrix_size;i++){
    for(int j=0;j<matrix_size+1;j++){
            myfile >> a;
            double_data[i][j]=a;
            cout << double_data[i][j] << " ";

    }cout << endl;
    }
//--------------------------------------------------------------

    MyMatrix matrix_test(matrix_size,matrix_size+1);
    MyMatrix matrix_test1(matrix_size,matrix_size+1);
    MyMatrix matrix_test2(matrix_size,matrix_size+1);
    matrix_test = double_data;
    matrix_test1 = double_data;
    matrix_test2 = double_data;
    MatrixXd eigen_test(matrix_size,matrix_size);

    for(int i=0;i<matrix_size;i++){
        for(int j=0;j<matrix_size;j++){
            eigen_test(i,j) = double_data[i][j];
        }
    }

    VectorXd eigen_vars(matrix_size);
    for(int i=0;i<matrix_size;i++)
        eigen_vars[i]=double_data[i][matrix_size];


    VectorXd eigen_partial_result =  eigen_test.partialPivLu().solve(eigen_vars);
    cout << "Eigen_partial : ";
    for(int i=0;i<matrix_size;i++)
        cout << fixed<< eigen_partial_result(i) << " ";//prawidlowy wynik eigenowy
    cout <<endl;


    vector<double>  gauss_partial_result = matrix_test.Gauss_partial();
    cout << "GausPartial   : ";
    for(int i=0;i<matrix_size;i++)
        cout << fixed<< gauss_partial_result[i] << " ";//implementacja partial gaussa
    cout <<endl;

    vector<double>  jacob_result = matrix_test1.Jacob(200);
    cout << "Jacob         : ";
    for(int i=0;i<matrix_size;i++)
        cout << fixed<< jacob_result[i] << " ";//implementacja partial gaussa
    cout <<endl;

    vector<double>  gauss_seidel_result = matrix_test2.Gauss_seidel(200);
    cout << "GausSeidel    : ";
    for(int i=0;i<matrix_size;i++)
        cout << fixed<< gauss_seidel_result[i] << " ";//implementacja partial gaussa
    cout <<endl;

    myfile.close();
return 0;
}
