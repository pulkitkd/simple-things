#include <iostream>
#include <fstream>
#include <cmath>

enum direction{MINUS,ZERO,PLUS};

#define SIZE 800
#define N_DV 3


//Define struct for the model
struct ModelD1Q3
{
 double ci[N_DV];
 double wt[N_DV];
 double theta0;
 double oneBytheta0;
 double rho;
 double u1;
 double fEq[N_DV];
};

//function declaration
void setModelParameters(ModelD1Q3 &myModel);

// definition, (Update the structure values)
void setModelParameters(ModelD1Q3 &myModel)
{
 myModel.ci[ZERO ]    =  0.0;
 myModel.ci[PLUS]     =  1.0;
 myModel.ci[MINUS]    = -1.0;
 myModel.wt[ZERO ]    = 4.0/6.0;
 myModel.wt[PLUS]     = 1.0/6.0;
 myModel.wt[MINUS]    = 1.0/6.0;
 myModel.theta0       = 1.0/3.0;
 myModel.oneBytheta0  = 3.0;
}


void checkModelParameters(ModelD1Q3 myModel)
{
double sum(0.0);
for(int i=0;i<N_DV;i++)
  sum = sum + myModel.wt[i];
std::cout<<"sum of weights = "<<sum<<std::endl;

sum = 0.0;
for(int i=0;i<N_DV;i++)
  sum = sum + myModel.ci[i]*myModel.wt[i];
std::cout<<"sum of wi*ci = "<<sum<<std::endl;

sum = 0.0;
for(int i=0;i<N_DV;i++)
  sum = sum + myModel.ci[i]*myModel.ci[i]*myModel.wt[i];
std::cout<<"sum of wi*ci^2 = "<<sum/myModel.theta0<<std::endl;
}

// Grid Class

class myGrid
{
public:
    myGrid()
    {
     G1   =  0;
     nB1  =  1;
     nE1  = SIZE;
     G2   = SIZE+1;
    }

  void initialiseIntegers()
  {
   for(int i=0;i<SIZE+2;i++)
    for(int j=0;j<N_DV;j++)
       data[i][j] = i*N_DV + j;
  }

 void printGrid(int population)
 {
   for(int i=nB1;i<=nE1;i++)
     std::cout<<data[i][population]<<"  ";

   std::cout<<std::endl;
 }

public:
  double data[SIZE+2][N_DV];
  int G1,G2;
  int nB1,nE1;
};

//periodic boundary conditions
void prepareBoundaryConditions(myGrid &gridLB)
{
 gridLB.data[gridLB.G1][PLUS]  =  gridLB.data[gridLB.nE1][PLUS] ;
 gridLB.data[gridLB.G2][MINUS] =  gridLB.data[gridLB.nB1][MINUS] ;
}


void advection(myGrid &gridLB)
{
// MINUS
 for(int i=gridLB.nB1;i<=gridLB.nE1;i++)
  gridLB.data[i][MINUS] = gridLB.data[i+1][MINUS];
// PLUS
 for(int i=gridLB.nE1;i>=gridLB.nB1;i--)
  gridLB.data[i][PLUS] = gridLB.data[i-1][PLUS];
}

// Get moments at a given Node point
void getMoments(ModelD1Q3 &myModel,myGrid &gridLB, int point)
{
 myModel.rho = 0.0;
 for(int i=0;i<N_DV;i++)
     myModel.rho += gridLB.data[point][i];
 myModel.u1 = 0.0;
 for(int i=0;i<N_DV;i++)
     myModel.u1 = myModel.u1 + (gridLB.data[point][i]*myModel.ci[i]);

 myModel.u1 = myModel.u1/ myModel.rho;
}


// Get Feq
void getFeq(ModelD1Q3 &myModel)
{
 for(int i=0;i<N_DV;i++)
  myModel.fEq[i] = myModel.wt[i]*myModel.rho * (1.0 + myModel.u1*myModel.ci[i]*myModel.oneBytheta0);

 double tmp, secTerm, MaSq;

 MaSq = myModel.u1*myModel.u1*myModel.oneBytheta0;
 secTerm = sqrt(1.0 + MaSq);

 myModel.fEq[ZERO]= myModel.wt[ZERO]*myModel.rho * (2.0-secTerm);

 tmp = myModel.rho * myModel.wt[PLUS];

 myModel.fEq[MINUS]=tmp*(myModel.u1*myModel.ci[MINUS]*myModel.oneBytheta0-1.0+2.0*secTerm);
 myModel.fEq[PLUS] =tmp*(myModel.u1*myModel.ci[PLUS ]*myModel.oneBytheta0-1.0+2.0*secTerm);

}


void copyToGrid(ModelD1Q3 &myModel,myGrid &gridLB, int point)
{
 for(int i=0;i<N_DV;i++)
  gridLB.data[point][i] =  myModel.fEq[i];
}


// Collision
// f = f + 2*beta*(feq-f)
void collide(ModelD1Q3 &myModel,myGrid &gridLB,double beta)
{
 for(int i=gridLB.nB1;i<=gridLB.nE1;i++)
 {
  getMoments(myModel,gridLB,i);
  getFeq(myModel);
  for(int j =0;j<N_DV;j++)
   gridLB.data[i][j] = gridLB.data[i][j] + 2.0*beta*(myModel.fEq[j]-gridLB.data[i][j]);
 }
}

// Wall BC
void prepareWallBC(myGrid &gridLB)
{
 gridLB.data[gridLB.G1][PLUS]  =  gridLB.data[gridLB.nB1][MINUS] ;
 gridLB.data[gridLB.G2][MINUS] =  gridLB.data[gridLB.nE1][PLUS] ;
}

void initialConditions(ModelD1Q3 &myModel,myGrid &gridLB)
{
 myModel.rho = 1.5;
 myModel.u1  = 0.0;
 getFeq(myModel);

 int halfDomain = (int)(0.5*SIZE);
 for(int i=gridLB.nB1;i<=halfDomain;i++)
  copyToGrid(myModel,gridLB,i);

 myModel.rho = 0.75;
 myModel.u1  = 0.0;
 getFeq(myModel);

 for(int i=halfDomain+1;i<=gridLB.nE1;i++)
  copyToGrid(myModel,gridLB,i);
}


void printdensity(ModelD1Q3 &myModel,myGrid &gridLB,int step)
{
 char filename[150];
 sprintf(filename, "results/density_%d.txt",step);
 std::ofstream file1;
 file1.open(filename);
 for(int i=gridLB.nB1;i<=gridLB.nE1;i++)
 {
  getMoments(myModel,gridLB,i);
  file1<<i<<"   "<<myModel.rho<<std::endl;
  }
 file1.close();
}



int main()
{

ModelD1Q3 lbModel;
setModelParameters(lbModel);
checkModelParameters(lbModel);
myGrid gridLBM;

// Simulation Parameters
 double deltaT = 1.0/SIZE;
 double nu     = 0.00001;
 double tau    = nu * lbModel.oneBytheta0;
 double beta   = deltaT/(2.0*tau + deltaT);

 int simulationTime = 500;

 initialConditions(lbModel,gridLBM);
 printdensity(lbModel,gridLBM,0);

 for(int time=0;time<simulationTime;time++)
 {
   collide(lbModel,gridLBM,beta);
   prepareWallBC(gridLBM);
   advection(gridLBM);
   printdensity(lbModel,gridLBM,time);
 }

// std::cout<<"executed"<<std::endl;
return 0;
}
