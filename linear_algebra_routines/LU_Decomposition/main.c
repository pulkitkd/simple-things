#include <stdio.h>
#include <stdlib.h>
#define N 5

void create_matrix(double m[N][N])
{
    int i, j;
    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
        scanf("%lf", &m[i][j]);
}

void display_matrix(double m[N][N])
{
    int i, j;
    for(i=0; i<N; i++)
    {
        for(j=0; j<N; j++)
        printf("%lf \t", m[i][j]);
        printf("\n");
    }
}

void initialize_matrix(double m[N][N])
{
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            m[i][j]=0;

}


void factorise_matrix(double A[N][N], double L[N][N], double U[N][N])
{
    int i=0, j=0, k=0;
    initialize_matrix(L);
    initialize_matrix(U);

    for(j=0; j<N; j++)
    {
        for(i=0; i<N; i++)
        {
            if(i<=j)
            {
                U[i][j] = A[i][j];
                for(k = 0; k < i; k++)
                    U[i][j]-= L[i][k]*U[k][j];
                if(i==j)
                    L[i][j]=1;
                else
                    L[i][j]=0;
            }
            else
            {
                L[i][j] = A[i][j];
                for(k = 0; k <= j-1; k++)
                    L[i][j]-= L[i][k]*U[k][j];
                L[i][j]/= U[j][j];
                U[i][j] = 0;
            }
        }
    }

}

int main()
{
    double A[N][N], L[N][N], U[N][N];
    printf("Enter the matrix  \n");
    create_matrix(A);
    printf("The matrix entered is \n");
    display_matrix(A);
    factorise_matrix(A, L, U);
    printf("The matrix L is \n");
    display_matrix(L);
    printf("The matrix U is \n");
    display_matrix(U);
    return 0;
}
