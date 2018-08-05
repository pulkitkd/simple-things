#include <iostream>
#include <fstream>
#include <math.h>
#include "square_matrix.h"
#define N 20
#define PI 3.14159265359

using namespace std;

void create_grid(vector<double, N + 1> &y)
{
    for(int i = 0; i <= N; i++)
        y(i) = cos(PI * i / N);
}

void create_D(vector<double, N + 1> &y, square_matrix<double, N + 1> &D)
{
    double c[N + 1];
    for(int i = 0; i <= N; i++)
        c[i] = 1.0;
    c[0] = 2.0;
    c[N] = 2.0;

    for(int i = 0; i <= N; i++)
        for(int j = 0; j <= N; j++)
        {
            if(i == j && i != 0 && i != N)
                D(i, j) = -y(i) / (2 * (1.0 - y(i) * y(i)));
            else
                D(i, j) = (c[i] / c [j]) * (pow(-1, (i + j)) / (y(i) - y(j)));
        }
    D(0, 0) = (2.0 * N * N + 1.0) / 6.0;
    D(N, N) = -(2.0 * N * N + 1.0) / 6.0;
}

void boundary_conditions(square_matrix<double, N + 1> &D2, square_matrix<double, N + 1> &I)
{
    for(int j = 0; j <= N; j++)
    {
        D2(0, j) = 0.0;
        D2(N, j) = 0.0;
        I(0, j) = 0.0;
        I(N, j) = 0.0;
    }
    D2(0, 0) = 1.0;
    D2(N, N) = 1.0;

}

void write_matrix(square_matrix<double, N + 1> &A, std::ofstream &file)
{
  for(int j = 0; j < N + 1; j++)
  {
      for(int i = 0; i < N + 1; i++)
    {
        file<<A(j ,i)<<" \t ";
    }
     file<<std::endl;
   }
}

int main()
{
    vector<double, N + 1> y;
    square_matrix<double, N + 1> I;
    square_matrix<double, N + 1> D;
    square_matrix<double, N + 1> D2;
    square_matrix<double, N + 1> D4;

    identity(I);
    create_grid(y);

    //cout<<y;
    create_D(y,D);
    //cout<<D;
    D2 = D * D;
    D4 = D2 * D2;
    I = -0.0625 * I;
    boundary_conditions(D2, I);
    cout<<I;

    ofstream fout;
    fout.open("A.txt");
    write_matrix(D2, fout);
    fout.close();

    ofstream fout2;
    fout.open("I.txt");
    write_matrix(I, fout);
    fout2.close();

    ofstream fout3;
    fout3.open("y.txt");
    for(int i = 0; i <= N; i++)
    fout3<<y(i)<<endl;
    fout3.close();

    return 0;
}
