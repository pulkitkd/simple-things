#ifndef LITTLE_MATRIX_H
#define LITTLE_MATRIX_H
#include "vector.h"

template<typename T, int D>

class little_matrix
{
public:
    int dim;
    T data[D][D];
    little_matrix()
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
    friend std::ostream& operator<<(std::ostream& os, const little_matrix<T,D> M)
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
    friend std::istream & operator>>(std::istream& is, little_matrix<T, D> &M)
    {
        std::cout<<"Enter a "<<D<<" X "<<D<<" matrix"<<std::endl;
        std::cout<<std::endl;
        for(int i = 0; i < D; i++)
            for(int j = 0; j < D; j++)
                is>>M(i,j);

        return is;
    }

    little_matrix<T, D> transpose()
    {
        little_matrix<T,D> P;

        for(int i = 0; i < D; i++)
            for(int j = 0; j < D; j++)
                P(i,j) = this->data[j][i];

        return P;
    }

};

//overload + operator
template<typename T, int D>
little_matrix<T, D> operator+(little_matrix<T, D> M1, little_matrix<T, D> M2)
{
    little_matrix<T, D> SUM;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            SUM(i,j) = M1(i,j) + M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload - operator
template<typename T, int D>
little_matrix<T, D> operator-(little_matrix<T, D> M1, little_matrix<T, D> M2)
{
    little_matrix<T, D> DIFF;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            DIFF(i,j) = M1(i,j) - M2(i,j);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DIFF;
}

//overload * operator
template<typename T, int D>
little_matrix<T, D> operator*(little_matrix <T, D> &M1, little_matrix<T, D> &M2)
{
    little_matrix<T, D> P;
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

//overload * operator for multiplication by scalar
template<typename T, int D>
little_matrix<T, D> operator*(double a, little_matrix<T, D> M)
{
    little_matrix<T, D> P;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            P(i,j) = a * M(i,j);
    return P;
}

//overload % operator to calculate power of a little_matrix
template<typename T, int D>
little_matrix<T, D> operator % (little_matrix<T, D> A, int N)
{
    little_matrix<T, D> EA;
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
template<typename T, int D>
vector<T, D> operator*(little_matrix<T, D> M, vector<T, D> V)
{
    vector<T, D> Res;
    vector<T, D> temp;
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            Res(i) = Res(i) + M(i,j) * V(j);

    return Res;
}

//function to fill the little_matrix with random integers up to 9
template<typename T, int D>
void random(little_matrix<T, D> &A)
{
    for(int i = 0; i < D; i++)
        for(int j = 0; j < D; j++)
            A(i,j) = rand() % 10;
}

//function to create identity matrix of required dimensions
template<typename T, int D>
void identity(little_matrix<T, D> &A)
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


#endif // LITTLE_MATRIX_H
