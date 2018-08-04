#ifndef VECTOR_H
#define VECTOR_H
#include<math.h>

template<typename T, int D>
class vector
{
public:
    int dim;
    T *vect;
    //constructor
    vector();
    //~vector();

    //overload output stream <<
    friend std::ostream& operator<<(std::ostream& os, const vector<T,D> v)
    {
        std::cout<<std::endl;
        for(int i = 0; i < D; i++)
            os<<v.vect[i]<<std::endl;
        std::cout<<std::endl;
        return os;
    }

    //overload input stream >>
    friend std::istream & operator>>(std::istream &is, vector<T, D> &v)
    {
        std::cout<<"Enter a "<<D<<" dimensional vector"<<std::endl;
        std::cout<<std::endl;
        for(int i = 0; i < D; i++)
            is>>v.vect[i];
        return is;
    }

    //overload index access operator () - for cout<< purposes
    T operator()(int i) const
    {
        return vect[i];
    }

    //overload index access operator () - for cin>> purposes
    T &operator()(int i)
    {
        return vect[i];
    }

    //overload = operator
    void operator=(const vector<T,D> &v1)
    {
        for(int i = 0; i < dim; i++)
            this->vect[i] = v1.vect[i];
    }

    //obtain the norm of the vector
    double norm()
    {
        double mag = 0.0;
        for(int i = 0; i < dim; i++)
             mag += this->vect[i] * this->vect[i];
        mag = sqrt(mag);
        return mag;
    }

    //obtain the unit vector corresponding to this vector
    vector<T, D> unit_vector()
    {
        vector<T, D> U;
        for(int i = 0; i < dim; i++)
            U.vect[i] = this->vect[i] / this->norm();

        return U;
    }

};

template<typename T, int D>
vector<T, D>::vector()
{
        dim = D;
        vect = new T[dim];
        for(int i = 0; i < dim; i++)
            vect[i] = 0.0;
}


//overload + operator
template<typename T, int D>
vector<T,D> operator+(vector<T,D> v1, vector<T,D> v2)
{
    vector<T,D> SUM;
    for(int i = 0; i < D; i++)
        SUM(i) = v1(i) + v2(i);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return SUM;
}

//overload - operator
template<typename T, int D>
vector<T,D> operator-(vector<T,D> v1, vector<T,D> v2)
{
    vector<T,D> DIFF;
    for(int i = 0; i < D; i++)
        DIFF(i) = v1(i) - v2(i);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DIFF;
}

//overload * operator for dot product
template<typename T, int D>
double operator*(vector<T,D> v1, vector<T,D> v2)
{
    double DOT;
    for(int i = 0; i < D; i++)
        DOT += v1(i) * v2(i);
        //SUM.vect[i] = v1.vect[i] + v2.vect[i];
    return DOT;
}

//overload * operator for multiplication by scalar
template<typename T, int D>
vector<T, D> operator*(vector<T, D> v1, double n)
{
    vector<T, D> V;
    for(int i = 0; i < D; i++)
        V(i) = n * v1(i);
    return V;
}

//overload * operator for multiplication by scalar
template<typename T, int D>
vector<T, D> operator*(double n, vector<T, D> v1)
{
    vector<T, D> V;
    for(int i = 0; i < D; i++)
        V(i) = n * v1(i);
    return V;
}

#endif // VECTOR_H
//.
//.
//.
//This work of art is due to Pulkit Kumar Dubey
