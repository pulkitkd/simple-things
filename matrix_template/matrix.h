#ifndef MATRIX_H
#define MATRIX_H
#include <stdlib.h>
#include <fstream>
#include "vector.h"

template<typename T, int R, int C>

class matrix //: public vector
{
public:
    int rows;
    int columns;
    T *mat;
    matrix()
    {
        mat = new T[R * C];
        for(int i = 0; i < R * C; i++)
            mat[i] = 0.0;
     }

    //overload index access operator () - for cout<< purposes
    T operator()(int i, int j) const {return mat[C * i + j];}

    //overload index access operator () - for cin>> purposes
    T &operator()(int i, int j){return mat[C * i + j];}

    //overload output stream <<
    friend std::ostream& operator<<(std::ostream& os, const matrix<T,R,C> M)
    {
        std::cout<<std::endl;
        for(int i = 0; i < R; i++)
        {
            for(int j = 0; j < C; j++)
            {
                os<<M(i,j)<<"\t";
            }
            std::cout<<std::endl;
        }
        return os;
    }

    //overload input stream >>
    friend std::istream & operator>>(std::istream &is, matrix<T, R, C> &M)
    {
        std::cout<<"Enter a "<<R<<" X "<<C<<" matrix"<<std::endl;
        std::cout<<std::endl;
        for(int i = 0; i < R; i++)
            for(int j = 0; j < C; j++)
                is>>M(i,j);

        return is;
    }

    matrix<T, C, R> transpose()
    {
        matrix<T, C, R> P;
        int k = 0;
        for(int i = 0; i < C; i++)
            for(int j = 0; j < R; j++)
            {
                P.mat[k] = this->mat[C * j + i];
                k++;
            }
        return P;
    }
};

//overload + operator
template<typename T, int R, int C>
matrix<T,R,C> operator+(matrix<T, R, C> M1, matrix<T, R, C> M2)
{
    matrix<T, R, C> SUM;
    for(int i = 0; i < R; i++)
        for(int j = 0; j < C; j++)
            SUM(i,j) = M1(i,j) + M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload - operator
template<typename T, int R, int C>
matrix<T,R,C> operator-(matrix<T, R, C> M1, matrix<T, R, C> M2)
{
    matrix<T, R, C> DIFF;
    for(int i = 0; i < R; i++)
        for(int j = 0; j < C; j++)
            DIFF(i,j) = M1(i,j) - M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DIFF;
}

//overload * operator
template<typename T, int R1, int C1, int R2, int C2>
matrix<T,R1,C2> operator*(matrix <T,R1,C1> &M1, matrix<T,R2,C2> &M2)
{
    matrix<T,R1,C2> P;
    if(C1 == R2)
    {
        for(int i = 0; i < R1; i++)
            for(int k = 0; k < C2; k++)
            {
                double sum = 0.0;
                for(int j = 0; j < R2; j++)
                    sum = sum + M1(i,j) * M2(j,k);
                P(i,k) = sum;
            }
    }
    else
        std::cout<<"Matrices do not have correct dimension for multiplication" <<std::endl;
    return P;
}

//overload * operator for multiplication by scalar
template<typename T, int R, int C>
matrix<T, R, C> operator*(double a, matrix<T,R,C> M)
{
    matrix<T, R, C> P;
    for(int i = 0; i < R; i++)
        for(int j = 0; j < C; j++)
            P(i,j) = a * M(i,j);
    return P;
}

//overload % operator to calculate power of a matrix
template<typename T, int R, int C>
matrix<T, R, C> operator % (matrix<T,R,C> A, int N)
{
    if(R != C)
    {
        std::cout<<"\n Not a square matrix. Process will not affect the argument. \n";
        return A;
    }
    matrix<T, R, C> EA;
    identity(EA);
    while(N > 0)
    {
        if(N % 2 == 1)
            EA = EA * A;
        A = A * A;
        N = N / 2;
    }
    return EA;
}

//overload * operator for multiplication by vector
template<typename T, int R, int C>
vector<T, C> operator*(matrix<T, R, C> M, vector<T, C> V)
{
    vector<T, C> Res;
    vector<T, C> temp;
    for(int i = 0; i < R; i++)
        for(int j = 0; j < C; j++)
            Res(i) = Res(i) + M(i,j) * V(j);

    return Res;
}

//function to fill the matrix with random integers up to 9
template<typename T, int R, int C>
void random(matrix<T,R,C> &A)
{
    for(int i = 0; i < R; i++)
        for(int j = 0; j < C; j++)
            A(i,j) = rand() % 10;
}

template<typename T, int R, int C>
void identity(matrix<T,R,C> &A)
{
    for(int i = 0; i < R; i++)
        for(int j = 0; j < C; j++)
        {
            if(i == j)
            A(i,j) = 1.0;
            else
            A(i,j) = 0.0;
        }
}


#endif // MATRIX_H
