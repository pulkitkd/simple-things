#include <iostream>
#include <fstream>
#include <conio.h>
#include <cmath>

#define PI 3.14159265359

using namespace std;

int main()
{
    double xn, xn1, pn, pnhalf, pn1;

    double omega = 0.1;
    double dt = 0.5;
    double T = 2 * PI / omega;
    int itr = (int)(T / dt);

    ofstream fout;
    fout.open("leapingfrog.txt");


    //initial conditions
    xn = 1.0;
    pn = 0.0;

    for(int i = 0; i <= itr; i++)
    {
        pnhalf = pn + dt * (-omega * omega * xn) / 2;
        xn1 = xn + dt * pnhalf;
        pn1 = pnhalf + dt * (-omega * omega * xn1) / 2;

        xn = xn1;
        pn = pn1;

        fout<<i<<"\t"<<xn1<<"\t"<<pn1<<endl;
    }

    cout<<"new momentum \t"<<pn1<<endl;
    cout<<"new position \t"<<xn1<<endl;

    fout.close();
    return 0;
}
