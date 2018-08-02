#ifndef COMPLEX_H
#define COMPLEX_H

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

    complex<T> operator = (const complex<T> &rhs)
        {
            this->re = rhs.re;
            this->im = rhs.im;

            return *this;
        }

};


template<typename T>
const complex<T> operator+(complex<T> first, complex<T> second)  {
    //complex<T>   tmp(first.re + second.re, first.im + second.im);
    complex<T> sum;
    sum.re = first.re + second.re;
    sum.im = first.im + second.im;
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
const complex<T>  operator*(const complex<T>  &first, const  complex<T>  &second) {
    complex<T> result;
    result.re = first.re * second.re- first.im *second.im;
    result.im = first.re * second.im + first.im *second.re;
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
std::ostream& operator<< (std::ostream& os, const complex<T> &num){
    if(num.im >= 0)
    os<<num.re<<" + i"<<num.im<<"";

    else
    os<<num.re<<" - i"<<-1 * num.im<<"";

    return os;

}

template<typename T>
std::istream& operator>> (std::istream& os, complex<T> &num){

    os>>num.re>>num.im;
    return os;

}

#endif // COMPLEX_H
