#include<stdio.h>
#include<math.h>
#define N 3
void main()
{
    int i, j, itr, m;
    double x[N],a[N][N],b[N],c;
    printf("Enter number of iterations \n");
    scanf("%d",&itr);
    printf("Enter the right hand side constants \n");
    for(i = 0; i < N; i++)
        scanf("%lf",&b[i]);

    printf("Enter the coefficients row wise \n");
    for(i=0;i<N;i++)
    {
        x[i]=0;
        for(j=0;j<N;j++)
            scanf("%lf",&a[i][j]);

    }
    for(m = 0; m < itr; m++)
    for(i = 0; i < N; i++)
    {
        c = b[i];
        for(j = 0; j < N; j++)
            if(i != j)
                c = c - a[i][j] * x[j];
                x[i] = c / a[i][i];
    }
        printf("The Solution is \n");
        for(i = 0; i < N; i++)
            printf("x(%d) = %lf \n", i, x[i]);


}
