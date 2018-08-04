
#include <stdio.h>
#include <stdlib.h>
# define N 2

int add(int a, int b)
{
    int c = a + b;
    return c;
}

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
    printf("resulting vector is \n");

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
    v[i][0] = 1;
}

double largest_element(double v[N][1])
{
    int i;
    double temp = v[0][0];
    for(i = 0; i < N; i++)
    {
        if(v[i][0] > temp)
            temp = v[i][0];
    }
    return temp;
}

void normalize(double v[N][1], double largest)
{
    int i;
    for(i = 0; i < N; i++)
    v[i][0] = v[i][0] / largest;
}

void equate_vectors(double v_in[N][1], double v_out[N][1])
{
    int i;
    for(i = 0; i < N; i++)
        v_in[i][1] = v_out[i][1];
}

double dot_product(double v1[N][1], double v2[N][1])
{
    int i;
    double dot = 0;
    for(i = 0; i < N; i++)
        dot = dot + v1[i][1] * v2[i][1];

    return dot;
}

double eigenvalue(double m[N][N], double v[N][1])
{
    double v_out[N][1];
    double eig;
    matrix_vector_multiply(m, v, v_out);
    eig = dot_product(v_out, v) / dot_product(v, v);

    return eig;

}

int main()
{
    int i;
    int iterations = 5;
    double L, eigenval;
    double m[N][N];
    double v_in[N][1];
    double v_out[N][1];
    printf("Enter matrix row wise\n");
    create_matrix(m);
    display_matrix(m);
    initialize_vector(v_in);

    for(i = 0; i < iterations; i++)
    {
        matrix_vector_multiply(m, v_in, v_out);
        display_vector(v_out);
        display_vector(v_in);
        L = largest_element(v_out);
        normalize(v_out, L);
        equate_vectors(v_out, v_in);

    }
    display_vector(v_out);
    eigenval = eigenvalue(m, v_out);
    printf("eigenvalue %lf \n", eigenval);

    //scanf("%d", &a);
    //scanf("%d", &b);
    return 0;
}
