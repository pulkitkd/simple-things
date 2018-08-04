#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 3

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

void display_vector(double m[N][1])
{
    int j;
    for(j=0; j<N; j++)
    printf("%lf \n", m[j][0]);
}

void matrix_vector_multiply(double m[N][N], double v_in[N][1], double v_out[N][1])
{
    int i, j;
    for(i = 0; i < N; i++)
    v_out[i][0] = 0;

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
                v_out[i][0] = v_out[i][0] + m[i][j]*v_in[j][0];
}

void initialize_vector(double v[N][1])
{
    int i;
    for(i = 0; i < N; i++)
    v[i][0] = 0.0;
}

double norm(double v[N][1])
{
    int i;
    double temp = 0;
    for(i = 0; i < N; i++)
    temp += v[i][0] * v[i][0];

    return sqrt(temp);
}

void normalize_vector(double v[N][1])
{
    double L;
    L = norm(v);
    int i;
    for(i = 0; i < N; i++)
    v[i][0] = (v[i][0] / L);
}

void equate_vectors(double v_in[N][1], double v_out[N][1])
{
    int i;
    for(i = 0; i < N; i++)
        v_in[i][0] = v_out[i][0];
}

double dot_product(double v1[N][1], double v2[N][1])
{
    int i;
    double dot = 0;
    for(i = 0; i < N; i++)
        dot = dot + v1[i][0] * v2[i][0];

    return dot;
}

//evaluates projection of v in the direction e and stores it in r
void projection(double v[N][1], double e[N][1], double r[N][1])
{
    int i;
    double a, b, c;
    b = dot_product(v, e);
    c = dot_product(e, e);
    a = b / c;

    for(i = 0; i < N; i++)
        r[i][0] = a * e[i][0];
}

void extract_column_vector(int c, double m[N][N], double v[N][1])
{
    int i;
    for(i = 0; i < N; i++)
        v[i][0] = m[i][c];
}

void equate_matrices(double in[N][N], double out[N][N])
{
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        out[i][j] = in[i][j];
}


void GS_orthogonalization(double A[N][N], double Q[N][N])
{
    int k, p;
    double temp[N][1];
    double project[N][1];
    double column_of_A[N][1];
    double column_of_Q[N][1];

        extract_column_vector(0, A, temp);
        normalize_vector(temp);
        for(p = 0; p < 3; p++)
            Q[p][0] = temp[p][0];

        extract_column_vector(1, A, column_of_A);
        extract_column_vector(0, Q, column_of_Q);
        projection(column_of_A, column_of_Q, project);
        Q[0][1] = A[0][1] - project[0][0];
        Q[1][1] = A[1][1] - project[1][0];
        Q[2][1] = A[2][1] - project[2][0];
        extract_column_vector(1, Q, column_of_Q);
        normalize_vector(column_of_Q);
        for(p = 0; p < 3; p++)
            Q[p][1] = column_of_Q[p][0];


        extract_column_vector(2, A, column_of_A);
        extract_column_vector(0, Q, column_of_Q);
        projection(column_of_A, column_of_Q, project);
        Q[0][2] = A[0][2] - project[0][0];
        Q[1][2] = A[1][2] - project[1][0];
        Q[2][2] = A[2][2] - project[2][0];

        extract_column_vector(1, Q, column_of_Q);
        projection(column_of_A, column_of_Q, project);
        Q[0][2] = Q[0][2] - project[0][0];
        Q[1][2] = Q[1][2] - project[1][0];
        Q[2][2] = Q[2][2] - project[2][0];

        extract_column_vector(2, Q, column_of_Q);
        normalize_vector(column_of_Q);
        for(p = 0; p < 3; p++)
            Q[p][2] = column_of_Q[p][0];

}

void transpose_matrix(double A[N][N], double Atranspose[N][N])
{
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
        Atranspose[i][j] = A[j][i];
}

void multiply_matrices(double A[N][N], double B[N][N], double result[N][N])
{
    int i, j, k;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            for(k = 0; k < N; k++)
            {
                result[i][j] = result[i][j] + A[i][k]*B[k][j];
                //printf("R[%d][%d] = %lf \n", i, j, result[i][j]);
            }

}


int main()
{
    double A[N][N];
    double Q[N][N];
    double QT[N][N];
    double R[N][N];

    printf("Enter matrix \n");
    create_matrix(A);
    GS_orthogonalization(A, Q);
    printf("The matrix Q is \n");
    display_matrix(Q);
    transpose_matrix(Q, QT);
    printf("The matrix Q Transpose is \n");
    display_matrix(QT);
    printf("The matrix A is \n");
    display_matrix(A);
    multiply_matrices(QT, A, R);
    printf("The matrix R is \n");
    display_matrix(R);

    return 0;

}
