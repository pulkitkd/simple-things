#ifndef _MATRIX_H_
#define _MATRIX_H_

template<typename T, int R, int C>
struct matrix_class
{
    T mat_elem[R][C];
    matrix_class();   //class constructor


    //overload ostream <<
    friend std::ostream& operator<<(std::ostream& os, const matrix_class<T,R,C> &os_matrix)
    {
    std::cout<<" Matrix is " <<std::endl;
    for(int i=0; i<R; i++)
    {
        for(int j=0;j<C;j++)
        {
            os<<os_matrix.mat_elem[i][j] <<"\t";
        }
        std::cout<<std::endl;
    }
    return os;
    }


    //overload subscript type operator
    //T operator()(const int i, const int j) const {return &mat_elem[i][j];}
    //T & operator()(const int i, const int j) {return mat_elem[i][j];}

    T operator()( int i, int j)const {return mat_elem[i][j];} //operator is pointing to an address of matrix element
    T &operator()( int i, int j) {return mat_elem[i][j];}




    //overload >> operator
    friend std::istream& operator>>(std::istream& is, matrix_class<T,R,C> &is_matrix)
    {
    std::cout<<"Please enter " <<R<<"x"<<C <<" elements of matrix" <<std::endl;
        for(int i=0; i<R; i++) //(sizeof(os_matrix)/sizeof(os_matrix[0]))
        {
            for(int j=0;j<C;j++)
            {
                //is>>is_matrix(i,j) ;
                is>>is_matrix.mat_elem[i][j] ;
            }
        }
        //check for symmetry
        int temp=1;
        for(int i=0; i<R; i++) //(sizeof(os_matrix)/sizeof(os_matrix[0]))
        {
            for(int j=0;j<i;j++)
            {
                if(is_matrix(i,j)== is_matrix(j,i)) {temp=temp*1;}
                if(is_matrix(i,j)!= is_matrix(j,i)) {temp=temp*0;}
            }
        }
        if(temp==1)
        {std::cout<<"Symmetric matrix found!" <<std::endl;}
    return is;
    }
};

//constructor
template<typename T, int R, int C>
matrix_class<T, R,C>::matrix_class()
    {
        mat_elem[R][C];
    }



// overload + operator
template<typename T, int R1, int C1, int R2, int C2>
matrix_class<T,R1,C1> operator+(matrix_class<T,R1,C1> &first, matrix_class<T,R2,C2> &second)
{
    matrix_class<T,R1,C1> matrix_result;
    if(R1==R2 && C1==C2)
    {
        for(int k=0;k<R1;k++)
        {
            for(int l=0;l<C1;l++)
            {
                matrix_result.mat_elem[k][l]=first.mat_elem[k][l]+second.mat_elem[k][l];
            }
        }
    }
    else
        std::cout<<"Unequal size matrices can't be added." <<std::endl;
    return matrix_result;
}



// overload + operator
template<typename T, int R1, int C1, int R2, int C2>
matrix_class<T,R1,C1> operator-(matrix_class<T,R1,C1> &first, matrix_class<T,R2,C2> &second)
{
    matrix_class<T,R1,C1> matrix_result;
    if(R1==R2 && C1==C2)
    {
        for(int k=0;k<R1;k++)
        {
            for(int l=0;l<C1;l++)
            {
                matrix_result.mat_elem[k][l]=first.mat_elem[k][l]-second.mat_elem[k][l];
            }
        }
    }
    else
        std::cout<<"Unequal size matrices can't be added." <<std::endl;
    return matrix_result;
}



//overload * operator
template<typename T, int R1, int C1, int R2, int C2>
matrix_class<T,R1,C2> operator*(matrix_class<T,R1,C1> &first, matrix_class<T,R2,C2> &second)
{
    matrix_class<T,R1,C2> matrix_result;
    if(C1==R2)
    {
        for(int i=0;i<R1;i++)
        {
            for(int k=0;k<C2;k++)
            {
                double sum2=0.0;
                for(int j=0;j<R2;j++)
                {
                    sum2=sum2+first.mat_elem[i][j]*second.mat_elem[j][k];
                }
                matrix_result.mat_elem[i][k]=sum2;
            }
        }
    }
    else
        std::cout<<"Given matrices can't be multiplied." <<std::endl;
    return matrix_result;
}


#endif
