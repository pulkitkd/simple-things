#ifndef COMPLEX_H
#define COMPLEX_H
#include <math.h>

template<typename T>
class complex
{
public:
    T re;
    T im;
    //constructor
    complex()
        {
            re = 0.0;
            im = 0.0;
        }
    //overload index access operator () - for cout<< purposes
    T operator()(int i) const
        {
            if(i == 0)
            return re;

            else
            return im;
        }

    //overload index access operator () - for cin>> purposes
    T &operator()(int i)
        {
            if(i == 0)
            return re;

            else
            return im;
        }

    complex<T>& operator = (const complex<T> &rhs)
        {
            this->re = rhs.re;
            this->im = rhs.im;

            return *this;
        }

    complex<T>& operator = (const double &a)
        {
            this->re = a;
            this->im = 0.0;
            return *this;
        }


};

//template<typename T>
//complex<T> operator = (complex<T> &c , double &a)
//{
//    c.re = a;
//    c.im = 0.0;
//    return c;
//}


template<typename T>
std::ostream& operator<< (std::ostream& os, const complex<T> &num){
    os<<num.re<<"  "<<num.im<<"\t";
    return os;
}

template<typename T>
std::istream& operator>> (std::istream& os, complex<T> &num){
    os>>num.re>>num.im;
    return os;
}

template<typename T>
const complex<T> operator+(complex<T> first, complex<T> second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first.re + second.re;
    sum.im = first.im + second.im;
    return sum;
}

template<typename T>
const complex<T> operator+(complex<T> first, T second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first.re + second;
    sum.im = first.im ;
    return sum;
}

template<typename T>
const complex<T> operator+(T first, complex<T> second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first + second.re;
    sum.im = second.im ;
    return sum;
}

template<typename T>
const complex<T> operator-(complex<T> first, complex<T> second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first.re - second.re;
    sum.im = first.im - second.im;
    return sum;
}

template<typename T>
const complex<T> operator-(T first, complex<T> second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first - second.re;
    sum.im = -1 * second.im;
    return sum;
}

template<typename T>
const complex<T> operator-(complex<T> first, T second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first.re - second;
    sum.im = first.im;
    return sum;
}

template<typename T>
const complex<T>  operator*(const complex<T>  &first, const  complex<T>  &second) {
    complex<T> result;
    result.re = first.re * second.re- first.im *second.im;
    result.im = first.re * second.im + first.im *second.re;
    return result;
}

template<typename T>
const complex<T>  operator*(const T &first, const  complex<T>  &second) {
    complex<T> result;
    result.re = first * second.re;
    result.im = first * second.im;
    return result;
}

template<typename T>
const complex<T>  operator*(const complex<T>  &first, const T &second) {
    complex<T> result;
    result.re = first.re * second;
    result.im = first.im * second;
    return result;
}

template<typename T>
const complex<T>  operator/(const complex<T>  &first, const  complex<T>  &second) {
    complex<T> result;
    result.re = (first.re * second.re + first.im *second.im)
              / (second.re * second.re + second.im * second.im );
    result.im = first.im * second.re - first.re *second.im
              / (second.re * second.re + second.im * second.im );
    return result;
}

template<typename T>
const complex<T> csqrt(complex<T> &A)
{
    complex<T> S;
    double r = sqrt(A.re * A.re + A.im * A.im);
    double theta = atan(A.im / A.re);

    S.re = sqrt(r) * cos(theta / 2);
    S.im = sqrt(r) * sin(theta / 2);

    return S;
}


#endif // COMPLEX_H
