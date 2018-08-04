//A code to invert a 4X4 matrix - can be generalized to larger matrices

#include<iostream>
#include"matrix.h"
#include <time.h>
#include<stdlib.h>

using namespace std;

matrix_class<double, 2, 2> invert2X2(matrix_class<double, 2, 2> temp1)
{
    matrix_class<double,2,2> temp;
    double t;
    double det;
    temp = temp1;

    det = temp(0,0)*temp(1,1) - temp(1,0)*temp(0,1);

    //swap leading diagonal elements
    t = temp(0,0);
    temp(0,0) = temp(1,1);
    temp(1,1) = t;
    //multiply by -1 other diagonal elements
    temp(1,0) = -temp(1,0);
    temp(0,1) = -temp(0,1);

    //calculate determinant and divide
    det = 1 /det;
    temp(0,0) = temp(0,0) * det;
    temp(0,1) = temp(0,1) * det;
    temp(1,0) = temp(1,0) * det;
    temp(1,1) = temp(1,1) * det;

    return temp;
}

void split(matrix_class<double,4,4> B, matrix_class<double,2,2> &b11,
           matrix_class<double,2,2> &b12, matrix_class<double,2,2> &b21,
           matrix_class<double,2,2> &b22)
{
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
        {
            b11(i,j) = B(i,j);
            b12(i,j) = B(i,j+2);
            b21(i,j) = B(i+2,j);
            b22(i,j) = B(i+2,j+2);
        }
}

void combine(matrix_class<double,4,4> &B, matrix_class<double,2,2> b11,
           matrix_class<double,2,2> b12, matrix_class<double,2,2> b21,
           matrix_class<double,2,2> b22)
{
    for(int i = 0; i < 2; i++)
        for(int j = 0; j < 2; j++)
        {
            B(i,j) = b11(i,j);
            B(i,j+2) = b12(i,j);
            B(i+2,j) = b21(i,j);
            B(i+2,j+2) = b22(i,j);
        }
}


int main()
{
    srand(time(NULL));
    matrix_class<double,4,4> A;
    matrix_class<double,4,4> B;

    matrix_class<double, 2, 2> a11;
    matrix_class<double, 2, 2> a12;
    matrix_class<double, 2, 2> a21;
    matrix_class<double, 2, 2> a22;
    matrix_class<double, 2, 2> a11inverse;

    matrix_class<double, 2, 2> b11;
    matrix_class<double, 2, 2> b12;
    matrix_class<double, 2, 2> b21;
    matrix_class<double, 2, 2> b22;
    matrix_class<double, 2, 2> buff;


    matrix_class<double,2,2> I;
    I(0,0) = 1.0;
    I(1,0) = 0.0;
    I(0,1) = 0.0;
    I(1,1) = 1.0;

    matrix_class<double,2,2> Null;
    Null(0,0) = 0.0;
    Null(1,0) = 0.0;
    Null(0,1) = 0.0;
    Null(1,1) = 0.0;

    cout<<"Randomly generating a matrix..."<<endl;

    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
                A(i,j) = rand() % 11;


    cout<<A;
    split(A,a11,a12,a21,a22);

    a11inverse = invert2X2(a11);
    b22 = a11inverse * a12;
    b22 = a21 * b22;
    b22 = a22 - b22;
    b22 = invert2X2(b22);

    b12 = a12 * b22;
    b12 = a11inverse * b12;
    b12 = Null - b12;

    b11 = b12 * a21;
    b11 = I - b11;
    b11 = b11 * a11inverse;

    b21 = b22 * a21;
    b21 = b21 * a11inverse;
    b21 = Null - b21;

    combine(B, b11, b12, b21, b22);
    cout<<endl;
    cout<<"The inverted matrix is"<<endl;
    cout<<B;
    cout<<endl;
    cout<<endl;
    cout<<"Multiplying by the original matrix. We expect Identity here"<<endl;
    cout<<A * B;

    return 0;
}
