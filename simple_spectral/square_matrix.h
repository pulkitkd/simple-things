#ifndef SQUARE_MATRIX_H
#define SQUARE_MATRIX_H
#include "vector.h"
#include "complex_class.h"
#include <stdlib.h>

template<typename T, int D>

class square_matrix
{
public:
    int dim;
    T data[D][D];
    square_matrix()
    {
        dim = D;
        for(int i = 0; i < dim; i++)
            for(int j = 0; j < dim; j++)
                data[i][j] = 0.0;
    }
    //overload index access operator () - for cout<< purposes
    T operator()(int i, int j) const {return data[i][j];}

    //overload index access operator () - for cin>> purposes
    T &operator()(int i, int j){return data[i][j];}

    //overload output stream <<
    friend std::ostream& operator<<(std::ostream& os, const square_matrix<T,D> &M)
    {
        std::cout<<std::endl;
        for(int i = 0; i < D; i++)
        {
            for(int j = 0; j < D; j++)
            {
                os<<M(i,j)<<"\t";
            }
            std::cout<<std::endl;
        }
        return os;
    }

    //overload input stream >>
    friend std::istream & operator>>(std::istream& is, square_matrix<T, D> &M)
    {
        std::cout<<"Enter a "<<D<<" X "<<D<<" matrix"<<std::endl;
        std::cout<<std::endl;
        for(int i = 0; i < D; i++)
            for(int j = 0; j < D; j++)
                is>>M(i,j);

        return is;
    }

    square_matrix<T, D> transpose()
    {
        square_matrix<T,D> P;

        for(int i = 0; i < D; i++)
            for(int j = 0; j < D; j++)
                P(i,j) = this->data[j][i];

        return P;
    }

    //overload = operator
    square_matrix<T, D>& operator=(const square_matrix<T, D> &M)
    {
        for(int i = 0; i < D; i++)
            for(int j = 0; j < D; j++)
                this->data[i][j] = M(i,j);
            //SUM.vect[i] = v1.vect[i] + v2.vect[i];
        return *this;
    }

    square_matrix<T, D> operator+=(square_matrix<T, D> &M)
    {
        for(int i = 0; i < D; i++)
            for(int j = 0; j < D; j++)
                data[i][j] += M(i,j);
            //SUM.vect[i] = v1.vect[i] + v2.vect[i];
        return *this;
    }
};

//overload + operator
template<typename T, int D>
square_matrix<T, D> operator+(square_matrix<T, D> &M1, square_matrix<T, D> &M2)
{
    square_matrix<T, D> SUM;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            SUM(i,j) = M1(i,j) + M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload + operator
template<typename T, int D>
square_matrix<complex<T>, D> operator+(square_matrix<complex<T>, D> &M1, square_matrix<T, D> &M2)
{
    square_matrix<complex<T>, D> SUM;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            SUM(i,j) = M1(i,j) + M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload + operator
template<typename T, int D>
square_matrix<complex<T>, D> operator+(square_matrix<T, D> &M1, square_matrix<complex<T>, D> &M2)
{
    square_matrix<complex<T>, D> SUM;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            SUM(i,j) = M1(i,j) + M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload + operator
template<typename T, int D>
square_matrix<complex<T>, D> operator+(square_matrix<complex<T>, D> &M1, square_matrix<complex<T>, D> &M2)
{
    square_matrix<complex<T>, D> SUM;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            SUM(i,j) = M1(i,j) + M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload - operator
template<typename T, int D>
square_matrix<T, D> operator-(square_matrix<T, D> &M1, square_matrix<T, D> &M2)
{
    square_matrix<T, D> DIFF;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            DIFF(i,j) = M1(i,j) - M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DIFF;
}

//overload - operator
template<typename T, int D>
square_matrix<complex<T>, D> operator-(square_matrix<complex<T>, D> &M1, square_matrix<T, D> &M2)
{
    square_matrix<complex<T>, D> DIFF;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            DIFF(i,j) = M1(i,j) - M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DIFF;
}

//overload - operator
template<typename T, int D>
square_matrix<complex<T>, D> operator-(square_matrix<T, D> &M1, square_matrix<complex<T>, D> &M2)
{
    square_matrix<complex<T>, D> DIFF;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            DIFF(i,j) = M1(i,j) - M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DIFF;
}


//overload * operator
template<typename T, int D>
square_matrix<T, D> operator*(square_matrix <T, D> &M1, square_matrix<T, D> &M2)
{
    square_matrix<T, D> P;
    for(int i = 0; i < D; i++)
        for(int k = 0; k < D; k++)
        {
            double sum = 0.0;
            for(int j = 0; j < D; j++)
                sum = sum + M1(i,j) * M2(j,k);
            P(i,k) = sum;
        }

    return P;
}

//overload * operator for multiplication by vector
template<typename T, int D>
vector<T, D> operator*(square_matrix<T, D> &M, vector<T, D> &V)
{
    vector<T, D> Res;

    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            Res(i) = Res(i) + M(i,j) * V(j);

    return Res;
}

//overload * operator for multiplication by scalar
template<typename T, int D>
square_matrix<T, D> operator*(double a, square_matrix<T, D> &M)
{
    square_matrix<T, D> P;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            P(i,j) = a * M(i,j);
    return P;
}

//overload * operator for multiplication by complex
template<typename T, int D>
square_matrix<complex<T>, D> operator*(complex<T> &a, square_matrix<complex<T>, D> &M)
{
    square_matrix<complex<T>, D> P;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            P(i,j) = a * M(i,j);
    return P;
}

//overload * operator for multiplication of complex with real matrix
template<typename T, int D>
square_matrix<complex<T>, D> operator*(complex<T> &a, square_matrix<T, D> &M)
{
    square_matrix<complex<T>, D> P;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            P(i,j) = a * M(i,j);
    return P;
}

//overload % operator to calculate power of a square_matrix
template<typename T, int D>
square_matrix<T, D> operator % (square_matrix<T, D> A, int N)
{
    square_matrix<T, D> EA;
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

//function to fill the square_matrix with random integers up to 9
template<typename T, int D>
void random(square_matrix<T, D> &A, int n)
{
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            A(i,j) = rand() % n;
}

//function to fill the vector with zeroes
template<typename T, int D>
void null(square_matrix<T, D> &A)
{
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            A(i,j) = 0.0;
}

//function to fill the vector with zeroes
template<typename T, int D>
void complexnull(square_matrix<complex<T>, D> &A)
{
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            {
                A(i,j).re = 0.0;
                A(i,j).im = 0.0;
            }
}

//function to create identity matrix of required dimensions
template<typename T, int D>
void identity(square_matrix<T, D> &A)
{
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
        {
            if(i == j)
            A(i,j) = 1.0;
            else
            A(i,j) = 0.0;
        }
}

//function to create complex identity matrix of required dimensions
template<typename T, int D>
void complexidentity(square_matrix<complex<T>, D> &A)
{
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
        {
            if(i == j)
            {
                A(i,j).re = 1.0;
                A(i,j).im = 0.0;
            }
            else
            {
                A(i,j).re = 0.0;
                A(i,j).im = 0.0;
            }
        }
}


#endif // square_matrix_H
