#include <iostream>
#include <cmath>
#include <fstream>
#include <conio.h>
#define PI 3.14159
#define N 128
using namespace std;

void create_grid(double X[N + 1], double dx)
{
    X[0] = 0;
    for(int i = 1; i <= N; i++)
        X[i] = X[i-1] + dx;
}

void write(double X[N + 1], double u0[N + 1], double u1[N + 1], int a)
{
    char filename[150];
    sprintf(filename, "results/solution_%d.txt", a);
    ofstream file1;
    file1.open(filename);
    for(int i = 0; i <= N; i++)
        file1<<i<<"\t"<<X[i]<<"\t"<<u0[i]<<"\t"<<u1[i]<<endl;

    file1.close();
}

void print(double X[N + 1], double u0[N + 1], double u1[N + 1])
{
    for(int i = 0; i <= N; i++)
        cout<<i<<"\t"<<X[i]<<"\t"<<u0[i]<<"\t"<<u1[i]<<endl;
}

int main()
{
    double X[N + 1];
    double u0[N + 1];
    double u1[N + 1];
    double L = 1;
    double T = 0;
    double t = 0;
    double dx = L / (N);
    double dt = 0.5 * dx * dx;
    int a = 0;

    cout<<"Enter final time \n";
    cin>>T;

    double dtbydxsq = dt/(dx * dx);

    create_grid(X, dx);

    //set initial condition
    for(int i = 0; i <= N; i++)
    {
        //u0[i] = sin(2 * PI * X[i]);
        //u0[i] = 10 * X[i] * X[i] * X[i] * (X[i] - 1);
        //u0[i] = exp(-(5*X[i] - 2.5) * (5*X[i] - 2.5));
        u0[i] = (sin(PI * X[i]) + sin(6 * PI * X[i]) + sin(32 * PI * X[i])) / 3;
    }

    while(t < T)
    {
        for(int i = 1; i < N; i++)
        {
            u1[0] = 0; //Boundary conditions
            u1[N] = 0;
            u1[i] = u0[i] + dtbydxsq * (u0[i + 1] + u0[i - 1] - 2 * u0[i]);
        }

        write(X, u0, u1, a);
        for(int j = 0; j <= N; j++)
            u0[j] = u1[j];

        t += dt;
        a++;
    }

  return 0;
}
