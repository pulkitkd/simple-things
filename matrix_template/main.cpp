#include <iostream>
#include <fstream>
#include "matrix.h"
#include "little_matrix.h"
//#include "vector.h" - vector.h is included in matrix.h

using namespace std;

int main()
{
    ofstream out;
    out.open ("test.txt");
    little_matrix<double, 3> M1;
    little_matrix<double, 3> M2;
    little_matrix<double, 3> M4;
    little_matrix<double, 3> M8;

    vector<double, 3> V1;

    cin>>V1;
    cout<<"\n Norm of V1 \n"<<V1.norm();
    cout<<"\n Unit vector in direction of V1 \n"<<V1.unit_vector();

    //fill M1 with random entries
    random(M1);

    out<<M1;
    cout<<"\n Randomly generated matrix M1 \n"<<M1;
    //matrix vector multiplication
    cout<<"\n M1 * V1 = \n"<<M1 * V1;
    //transpose
    cout<<"\n Transpose of M1 \n"<<M1.transpose();

    //matrix multiplication followed by same result from exponentiation
    M2 = M1 * M1;
    M4 = M2 * M2;
    M8 = M4 * M4;
    cout<<"\n M1 \n"<<M1;
    cout<<"\n M2 \n"<<M2;
    cout<<"\n M4 \n"<<M4;
    cout<<"\n M8 \n"<<M8;
    cout<<"\n Same result from exponentiation"<<endl;
    cout<<M1 % 8;

    out.close();

    return 0;
}

/*
int main()
{
    matrix<double, 3, 3> M1;
    matrix<double, 3, 3> M2;
    matrix<double, 3, 3> M4;
    matrix<double, 3, 3> M8;

    vector<double, 3> V1;

    cin>>V1;
    cout<<"\n Norm of V1 \n"<<V1.norm();
    cout<<"\n Unit vector in direction of V1 \n"<<V1.unit_vector();

    //fill M1 with random entries
    random(M1);

    cout<<"\n Randomly generated matrix M1 \n"<<M1;
    //matrix vector multiplication
    cout<<"\n M1 * V1 = \n"<<M1 * V1;
    //transpose
    cout<<"\n Transpose of M1 \n"<<M1.transpose();

    //matrix multiplication followed by same result from exponentiation
    M2 = M1 * M1;
    M4 = M2 * M2;
    M8 = M4 * M4;
    cout<<"\n M1 \n"<<M1;
    cout<<"\n M2 \n"<<M2;
    cout<<"\n M4 \n"<<M4;
    cout<<"\n M8 \n"<<M8;
    cout<<"\n Same result from exponentiation"<<endl;
    cout<<M1 % 8;

    return 0;
}
*/
