#include<stdio.h>
#include<math.h>
#include<limits.h>


double Max(double u[], int size) //function to return maximum value from an array
{
  int i=0;
  double maximum=0;
  for(i=1;i<size-1;i++)
    if(u[i]>maximum)
      maximum=fabs(u[i]);

  return maximum;
}

double Max2(double a ,double b) //function to return the mod of greater of two supplied values
{
  if (fabs(a) > fabs(b))
    return fabs(a);
  else
    return fabs(b);
}

double Min2(double a ,double b) //function to return the greater of two supplied values
{
  if (a < b)
    return a;
  else
    return b;
}

double avgflux(double a, double b) //arithmetic mean flux - blows up
{
  double c=(a*a+b*b)/4;
  return c;
}

double hmflux(double a,double b) //harmonic mean flux - bad results
{
  return a*a*b*b/2*(a*a+b*b) ;
}

double laxflux(double a, double b, double h, double dt) //lax-friedrichs flux
{
  double c=(a*a+b*b)/4 - h*(b-a)/(2*dt);
  return c;
}


double gudflux(double a, double b) //gudonov's flux
{
  double ml, mr;
  ml=Max2(0,a);
  mr=Min2(b,0);
  return Max2(ml*ml/2, mr*mr/2);
}

double roeflux(double a, double b) //roe flux
{
  double lambda=fabs((b+a)/2);
  return (a*a + b*b)/4-lambda*(b-a)/2;
}

double rusanflux(double a, double b) //rusanov flux
{
  double lambda=Max2(fabs(b),fabs(a));
  return (a*a + b*b)/4-lambda*(b-a)/2;
}

//the main function can be altered to put entire data into a single file or generate seperate data files for each time step
void main()
{
  char filename[PATH_MAX];
  FILE *fp,*fp2;
  int N;
  double CFL;
  printf("enter no. of divisions \n");
  scanf("%d",&N);
  printf("enter CFL no. \n");
  scanf("%lf",&CFL);

  int i,j,k=0,init;
  double u1[N],u0[N],ue[N],Tf=0.25,dt,h,fo,fi,t;
  h=1.0/N;

  fp=fopen("burgersnum.dat","w");
  fp2=fopen("burgersexact.dat","w");

  printf("choose initial condition \n 1-> down step \n 2-> up step \n 3-> peak \n") ;
  scanf("%d",&init);

  //initial conditions
  //////////////////////
  if(init == 1)
    {
      for(i=0;i<N/2;i++)
	u0[i]=1;
      for(i=N/2;i<N;i++)
	u0[i]=0;
    }

  else if (init == 2)
    {
      for(i=0;i<N/2;i++)
	u0[i]=0;
      for(i=N/2;i<N;i++)
	u0[i]=1;
    }
  else
    for(i=0;i<N;i++)
      u0[i]=pow(10,(-1*(10*i*h-5)*(10*i*h-5)));
  //////////////////////

  dt=CFL*h/Max(u0,N); //time step

  //FVM Implementation
  //////////////////////
  for(t=0;t<Tf;t=t+dt)
    {

      for(j=1;j<N-1;j++)
	{
	  fi=gudflux(u0[j-1],u0[j]);
	  fo=gudflux(u0[j],u0[j+1]);
	  u1[j]=u0[j]-(dt)*(fo-fi)/h;
	  fi=fo;
	}

      for(i=1;i<N-1;i++)
	u0[i]=u1[i];
	
      //calculate time step for next iteration
      dt=CFL*h/Max(u0,N);  
    }
  //sprintf(filename,"%d_burgers.dat",k);
  //fp=fopen(filename,"w");
  for(j=0;j<N;j++)
    fprintf(fp,"%lf \t  %lf \n",j*h,u0[j]);

  //fclose(fp);
  //k++;

  //plotting with analytical solution for comparison
  //analytical solution is only available for step-jump initial conditions
  if (init==2)
    {
  for(j=0;j<N/2;j++)
    ue[j]=0;
  for(j=N/2;j<0.75*N;j++)
    ue[j]=(j*h-0.5)/t;
  for(j=0.75*N;j<N;j++)
    ue[j]=1;
  for(j=0;j<N;j++)
    fprintf(fp2,"%lf \t  %lf \n",j*h,ue[j]);
    }

  else if(init==1)
    {
      for(j=0;j<0.625*N;j++)
	ue[j]=1;
      for(j=0.625*N;j<N;j++)
	ue[j]=0;
      for(j=0;j<N;j++)
	fprintf(fp2,"%lf \t  %lf \n",j*h,ue[j]);
    }
  fclose(fp);
  fclose(fp2);
}
