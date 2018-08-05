/*Steady state heat transfer with one dirichilet (const. temperature)
and one Neumann (constant heat flux) boundary condition where flux
 is specified using convection coefficient*/

#include <iostream>
#include <conio.h>
#include <cmath>
#include <fstream>

#define N 200
#define PI 3.14159

using namespace std;

void solve(double* a, double* b, double* c, double* d, int n) {
     n--; // since we start from x0 (not x1)
     c[0] /= b[0];
     d[0] /= b[0];

     for (int i = 1; i < n; i++) {
         c[i] /= b[i] - a[i]*c[i-1];
         d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
     }

     d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

     for (int i = n; i-- > 0;) {
         d[i] -= c[i]*d[i+1];
     }
}

void create_grid(double X[N + 1], double dx)
{
    X[0] = 0.0;
    for(int i = 0; i <= N; i++)
        X[i + 1] = X[i] + dx;
}

void set_initial_condition(double U0[N + 1], double X[N + 1])
{
    for(int i = 0; i <= N; i++)
        U0[i] = exp(-50 * (X[i] - 0.5) * (X[i] - 0.5));
}

void set_boundary_condition(double U0[N + 1])
{
    U0[0] = 0;
    U0[N] = 0;
}

void create_RHS(double C[N + 1], double U0[N + 1], double lambda)
{
    C[0] = 0;
    C[N] = 0;
    for(int i = 1; i < N; i++)
        C[i] = lambda * (U0[i - 1] + U0[i + 1]) + (1 - 2 * lambda) * U0[i];

}

void create_Tri_Diagonal_Matrix(double A[N + 1][N + 1],double U1[N + 1], double lambda)
{
    for(int i = 1; i <= N - 1; i++)
        {
            A[i][i - 1] = -lambda;
            A[i][i] = 1 + 2 * lambda;
            A[i][i + 1] = -lambda;
        }

    A[0][0] = 1.0;
    A[N][N] = 1.0;

}

void put_diagonals_in_arrays(double A[N + 1][N + 1], double D1[N + 1], double D2[N + 1], double D3[N + 1])
{
    for(int i = 0; i <= N; i++)
    {
        D1[i] = 0.0;
        D2[i] = 0.0;
        D3[i] = 0.0;
    }

    for(int i = 0; i < N; i++) //assign elements in the diagonal arrays
	D1[i + 1] = A[i + 1][i];

	for(int i = 0; i <= N; i++) //assign elements in the diagonal arrays
	D2[i] = A[i][i];

	for(int i = 0; i < N; i++) //assign elements in the diagonal arrays
	D3[i] = A[i][i + 1];

}

void write(double X[N + 1], double U0[N + 1], int a)
{
    char filename[150];
    sprintf(filename, "results/solution_%d.txt", a);
    ofstream file1;
    file1.open(filename);
    for(int i = 0; i <= N; i++)
        file1<<i<<"\t"<<X[i]<<"\t"<<U0[i]<<endl;

    file1.close();
}

void display_vector(double U[N + 1])
{
    for(int i = 0; i <= N; i++)
        cout<<U[i]<<endl;
}

void display_matrix(double A[N + 1][N + 1])
{
    for(int j = 0; j <= N; j++){
        for(int i = 0; i <= N; i++)
            cout<<A[j][i]<<"\t";
            cout<<endl;
    }
}

int main()
{
    double U0[N + 1];
    double U1[N + 1];
    double X[N + 1];
    double A[N + 1][N + 1];
    double D1[N + 1], D2[N + 1], D3[N + 1];
    double C[N + 1];

    int L = 1; //length of domain
    double T = 0.1; //total time
    double f = 2.0; //diffusion factor, e.g. thermal conductivity
    double dx = 1.0 / (N); //grid size
    double dt = dx * dx;

    double lambda = f * dt / (dx * dx);

    for(int i = 0; i <= N; i++)
        for(int j = 0; j <= N; j++)
            A[i][j] = 0.0;


    create_grid(X, dx);

    set_initial_condition(U0, X);
    set_boundary_condition(U0);

    int iterations = (int)(T / dt);
    cout<<"Total iterations to be done = "<<iterations<<endl;
    cout<<"Press Enter to continue"<<endl;
    getch();

    for(int itr = 0; itr <= iterations; itr++)
    {

        create_RHS(C, U0, lambda);

        create_Tri_Diagonal_Matrix(A, U1, lambda);

        put_diagonals_in_arrays(A, D1, D2, D3);

        solve(D1, D2, D3, C, N + 1);

        for(int i = 0; i <= N; i++)
        U0[i] = C[i];

        write(X, U0, itr);
    }

    display_vector(C);




    return 0;
}
