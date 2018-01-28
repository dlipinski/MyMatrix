#include <iostream>
#include <ctime>
#include <stdio.h>
#include <cstdlib>
#include <Eigen/Dense>
#include <iomanip>
#include <vector>

using namespace std;


class MyMatrix {
 private:
  unsigned height;
  unsigned width;

 double abs(double x) {
        return (x >= 0) ? x : -x;
    }



 public:
vector<vector <double> > matrix;
MyMatrix(unsigned _rows, unsigned _cols) {
  matrix.resize(_rows);
  for (unsigned i=0; i<matrix.size(); i++) {
    matrix[i].resize(_cols, 0);
  }
  height = _rows;
  width = _cols;
}

MyMatrix(const MyMatrix& arg) {
  matrix = arg.matrix;
  height = arg.get_rows();
  width = arg.get_cols();
}



MyMatrix& operator=(const MyMatrix& arg) {
 /* if (&arg == this)
    return *this;

  unsigned new_rows = arg.get_rows();
  unsigned new_cols = arg.get_cols();

  matrix.resize(new_rows);
  for (unsigned i=0; i<matrix.size(); i++) {
    matrix[i].resize(new_cols);
  }
*/
  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      matrix[i][j] = arg(i, j);
    }
  }
 // height = new_rows;
 // width = new_cols;

  return *this;
}

void print_matrix(){
    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++)
            cout << matrix[i][j] << " ";
        cout << endl;
    }
}

MyMatrix& operator=(const double tab[]) {
    for(unsigned i=0; i< height;i++)
        for(unsigned j=0; j<width;j++)
            matrix[i][j]=tab[i*width+j];
    return *this;
}
MyMatrix& operator=( vector<vector<double> > double_data) {
    for(unsigned i=0; i< height;i++)
        for(unsigned j=0; j<width;j++)
            matrix[i][j]=double_data[i][j];
    return *this;
}



MyMatrix operator+(const MyMatrix& arg) {
  MyMatrix result(height, width);

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      result(i,j) = matrix[i][j] + arg(i,j);
    }
  }

  return result;
}


MyMatrix& operator+=(const MyMatrix& arg) {
  unsigned height = arg.get_rows();
  unsigned width = arg.get_cols();

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      this->matrix[i][j] += arg(i,j);
    }
  }

  return *this;
}

MyMatrix operator-(const MyMatrix& arg) {
  unsigned height = arg.get_rows();
  unsigned width = arg.get_cols();
  MyMatrix result(height, width);

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      result(i,j) = this->matrix[i][j] - arg(i,j);
    }
  }

  return result;
}

MyMatrix& operator-=(const MyMatrix& arg) {
  unsigned height = arg.get_rows();
  unsigned width = arg.get_cols();

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      this->matrix[i][j] -= arg(i,j);
    }
  }

  return *this;
}


MyMatrix operator*(const MyMatrix& arg) {
  MyMatrix result(height, width);

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      for (unsigned k=0; k<height; k++) {
        result(i,j) += matrix[i][k] * arg(k,j);
      }
    }
  }

  return result;
}


MyMatrix operator*(const int& arg) {
  MyMatrix result(height, width);

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      result(i,j)=matrix[i][j] *  arg;
    }
  }

  return result;
}


MyMatrix& operator*=(const MyMatrix& arg) {
  MyMatrix result = (*this) * arg;
  (*this) = result;
  return *this;
}


MyMatrix transpose() {
  MyMatrix result(height, width);

  for (unsigned i=0; i<height; i++) {
    for (unsigned j=0; j<width; j++) {
      result(i,j) = this->matrix[j][i];
    }
  }

  return result;
}

double determinant(){
    double sum=0;
    double temp[2*height-1][width];
    double num=0;
    for(int i=0;i<height;i++)
        for(int j=0;j<width;i++)
            temp[i][j]=matrix[i][j];

    for(int i=0;i<height-1;i++)
        for(int j=0;j<width;i++)
            temp[i+height][j]=matrix[i][j];


    for(int i=0;i<height;i++){
        for(int j=0;j<width;j++){
            for(int y=0;y<width;y++){
                sum*=temp[i-y][j-y];
            }
            sum+=num;
        }
    }

    for(int i=0;i<height;i++){
        for(int j=width-1;j>=0;j++){
            for(int y=0;y<width;y++){
                sum*=temp[i-y][j-y];
            }
            sum-=num;
        }

    }
    return sum;
}

MyMatrix inverse(){
    double temp[height][width];
    MyMatrix woo(height,width);

    MyMatrix d(height-1,width-1);

    for(int i=0;i<height;i++)
        for(int j=0;j<width;j++){
            temp[i][j]=matrix[i][j];
            woo(i,j)=matrix[i][j];
        }
    int l=0,r=0;
    for(int i=0;i<height;i++)
        for(int j=0;j<width;j++){
            l=0;r=0;
            for(int in=0;in<height;in++)
                for(int jn=0;jn<width;jn++){
                    if(in!=i && jn!=j){
                        d(l,r)=temp[i][j];
                        l++;
                        if(l==width-1){ l=0;r++;}
                    }
                }
            matrix[i][j] = ((-1)^(i+j) ) *d.determinant();
        }

    double deto = 1/woo.determinant();
}

MyMatrix inverse_diag(){
    MyMatrix result(height,width);
    for(int i=0;i<height;i++)
        for(int j=0;j<width;j++)
            result(i,j)=matrix[i][j];
    for(int i=0;i<height;i++)
        for(int j=0;j<width;j++)
        if(i==j){
                double temp = result[i][j];
                result(i,j) = 1/temp;
        }
    return result;
}
// Macierz * wektor
vector<double> operator[](std::size_t idx)
{return matrix[idx];}

vector<double> operator*(const vector<double>& arg) {
            vector<double> result(arg.size(),0);
            for (unsigned i=0; i<height; i++) {
                for (unsigned j=0; j<width; j++) {
                    result[i] += this->matrix[i][j] * arg[j];
                }
            }
            return result;
        }



vector<double> diag_vec() {
  vector<double> result(height, 0);

  for (unsigned i=0; i<height; i++) {
    result[i] = this->matrix[i][i];
  }

  return result;
}


double& operator()(const unsigned& row, const unsigned& col) {
  return this->matrix[row][col];
}


const double& operator()(const unsigned& row, const unsigned& col) const {
  return this->matrix[row][col];
}

unsigned get_rows() const {
  return this->height;
}

unsigned get_cols() const {
  return this->width;
}

//-------------------

        void print_matrix(MyMatrix& mat, int prec){
        	cout << "\n";
        	for (int i=0; i<mat.get_rows(); i++) {
    			for (int j=0; j<mat.get_cols(); j++) {
    		cout << setprecision(prec) << mat(i,j) << "\double";
   			}
   		 	cout << endl;
  			}
		}

		 void print_matrix(MyMatrix& mat){
        	cout << "\n";
        	for (int i=0; i<mat.get_rows(); i++) {
    			for (int j=0; j<mat.get_cols(); j++) {
    		cout << mat(i,j) << "\double";
   			}
   		 	cout << endl;
  			}
		}

		void add_rows(int a, int b){
            for(int i=0;i<width;i++)
                matrix[a][i]+=matrix[b][i];
		}

        void delete_zeros(){
            for(unsigned i=0;i<height;i++){
                for(unsigned j=0;j<width;j++){
                    if(i==j && matrix[i][j]==0){
                        unsigned b=j;
                            for(unsigned y=0;y<height;y++){
                                if(matrix[y][b]!=0){
                                    add_rows(i,y);
                                }
                            }
                    }

                }
            }
        }

//-----------------------------------------------------------Gauss-Jordan po³owicznego wyboru
        vector<double> Gauss_partial() {

            int n = get_rows();
            for (int i = 0; i < n; i++) {
                // znajd wiersz z maksymalnym elementem
                double maxEl = abs(matrix[i][i]);
                int maxRow = i;
                for (int k = i+1; k < n; k++) {
                    if (abs(matrix[k][i]) > maxEl) {
                        maxEl = abs(matrix[k][i]);
                        maxRow = k;
                    }
                }
                // zamieñ maksymalny wiersz z obecnym
                for (int k = i; k < n+1; k++) {
                    double pom = matrix[maxRow][k];
                    matrix[maxRow][k] = matrix[i][k];
                    matrix[i][k] = pom;
                }
                // wyprowad zera przed obecnym wierszem
                for (int k = i+1; k < n; k++) {
                    double c = -matrix[k][i] / matrix[i][i];
                    for (int j = i; j < n+1; j++) {
                        if (i == j) {
                            matrix[k][j] = 0;
                        } else {
                            matrix[k][j] += c * matrix[i][j];
                        }
                    }
                }
            }
            // rozwi¹¿ Ax = B za pomoc¹ powsta³ej macierzy trójk¹tnej
            vector<double> x(n);
            for (int i=n-1; i>=0; i--) {
                x[i] = matrix[i][n] / matrix[i][i];
                for (int k=i-1;k>=0; k--) {
                    matrix[k][n] -= matrix[k][i] * x[i];
                }
            }
            return x;
        }




        vector<double> Jacob(int iter){
        //cout << "---------------------"<<endl<<"my_jacob:"<<endl<<endl;
        MyMatrix A(height,width-1),L(height,width-1),D(height,width-1),U(height,width-1),N(height,width-1),M(height,width-1);
        vector<double> B;
        B.resize(height);

       // cout <<"A:"<<endl;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                A(i,j)=matrix[i][j];
                //cout<<A[i][j]<<" ";
            }
           // cout << endl;
        }
        //cout <<"L:"<<endl;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                if(j<i)
                    L(i,j)=matrix[i][j];
               // cout<<L[i][j]<<" ";
            }
            //cout << endl;
        }
         //cout <<"D:"<<endl;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                if(j==i)
                    D(i,j)=matrix[i][j];
                //cout<<D[i][j]<<" ";
            }
            //cout << endl;
        }
        //cout <<"U:"<<endl;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                if(j>i)
                    U(i,j)=matrix[i][j];
              //  cout<<U[i][j]<<" ";
            }
           // cout << endl;
        }

        //cout <<"N:(D^-1):"<<endl;
        N = D.inverse_diag();
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
             //   cout<<N[i][j]<<" ";
            }
           // cout << endl;
        }
       // cout << "-N" << endl;
        MyMatrix Nminus(height,width);
        Nminus = N * (-1);
         for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
               // cout<<Nminus[i][j]<<" ";
            }
            //cout << endl;
        }
        //cout << "L+U" << endl;
        MyMatrix tempp(height,width);
        tempp = L+U;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
            //    cout<<tempp[i][j]<<" ";
            }
         //   cout << endl;
        }
       // cout <<"M(-N*(L+U)):"<<endl;

        M = ( Nminus * tempp );
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
               // cout<<M[i][j]<<" ";
            }
           // cout << endl;
        }
        //cout << endl<<"B:"<<endl;
        for(int i=0;i<height;i++){
            B[i]=matrix[i][width-1];
           // cout << B[i] << " ";
        }

        vector<double> ans;
        ans.resize(height,0);

        vector<double> _ans;
        _ans.resize(height,0);

       // cout << endl;
        for(int y=0;y<iter;y++){
            for(int i=0;i<ans.size();i++)
                ans[i]=0;
            for(int i=0;i<ans.size();i++){
                    for(int j=0;j<ans.size();j++){
                        ans[i]+=(M[i][j]*_ans[j] + N[i][j]*B[j]);
                    }

            }
            for(int i=0;i<height;i++)
             //   cout << _ans[i] << " ";
          //  cout << endl;
            for(int i=0;i<height;i++)
            //    cout << ans[i] << " ";
           // cout << endl;
           // cout << endl;
            _ans=ans;
        }
        return ans;
        }

        vector<double> Gauss_seidel(int iter){
          //cout << "---------------------"<<endl<<"my_jacob:"<<endl<<endl;
        MyMatrix A(height,width-1),L(height,width-1),D(height,width-1),U(height,width-1),N(height,width-1),M(height,width-1),NL(height,width-1),NU(height,width-1);
        vector<double> B,NB;
        B.resize(height,0);
        NB.resize(height,0);
       //create A (left matrix)
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                A(i,j)=matrix[i][j];
                //cout<<A[i][j]<<" ";
            }
           // cout << endl;
        }
        //create B ( right matrix)
        for(int i=0;i<height;i++){
            B[i]=matrix[i][width-1];
            //cout << B[i] << " ";
        }
        //create L (lower A)
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                if(j<i)
                    L(i,j)=matrix[i][j];
               // cout<<L[i][j]<<" ";
            }
            //cout << endl;
        }
        //create D (diag A)
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                if(j==i)
                    D(i,j)=matrix[i][j];
                //cout<<D[i][j]<<" ";
            }
            //cout << endl;
        }
        //create U (upper A)
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
                if(j>i)
                    U(i,j)=matrix[i][j];
                //cout<<U[i][j]<<" ";
            }
           // cout << endl;
        }

        //create N ( D^-1)
        N = D.inverse_diag();
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
               // cout<<N[i][j]<<" ";
            }
            //cout << endl;
        }

        NB = N*B;
        NL = N*L;
        NU = N*U;

        //cout <<endl<<"NB:"<<endl;
        for(int i=0;i<height;i++){
            //cout << NB[i]<<" ";
        }
        //cout << endl;
        //cout << "NL:"<<endl;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
               // cout<<NL[i][j]<<" ";
            }
           // cout << endl;
        }
        //cout <<endl<< "NU:"<<endl;
        for(int i=0;i<height;i++){
            for(int j=0;j<width-1;j++){
          //     cout<<NU[i][j]<<" ";
            }
            //cout << endl;
        }
        vector<double> ans;
        ans.resize(height,0);

        vector<double> _ans;
        _ans.resize(height,0);

       // cout << endl;
        for(int y=0;y<iter;y++){
            for(int i=0;i<ans.size();i++)
                ans[i]=0;
            for(int i=0;i<ans.size();i++){
                    for(int j=0;j<ans.size();j++){
                        ans[i]=(NB[j] - NU[i][j]*_ans[j])/(NL[i][j]+1);

                    }

            }

            _ans=ans;
        }
        return ans;
        }

};
