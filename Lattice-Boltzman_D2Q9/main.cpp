#include <iostream>
#include <conio.h>
#include <fstream>
#include <iomanip>
#include <cmath>

#define N 50 //number of nodes in x and y direction
#define Nc 9 //number of velocity directions
/*
//initialize the f, rho and u variables
//->assign initial conditions at all points
//->use initial conditions to obtain f's for all points using eq for feq's
//->ADVECT <------------------------------------------------
//->get new f's at all nodes                                |
//->calculate new rho and u                                 |
//->use new rho and u to get new feq's                      |
//->using f's and feq's, perform COLLIDE; obtain new f's    |
//->use BOUNDARY CONDITIONS                                 |
//->ADVECT -------------------------------------------------
*/
using namespace std;

enum direction {
     center = 0, east = 1, north, west, south,
     neast, nwest, swest, seast};

//Define structure for the D2Q9 model
//This much information is stored for each node
//These variables must be calculated for each node per iteration
struct modelD2Q9
{
 double ci[Nc][2];  //given
 double wt[Nc];     //given
 double fEq[Nc];
 double rho;
 double u1x;
 double u1y;
 double theta0;     //given
 double onebytheta0;//given
};


void checkModelParameters(modelD2Q9 myModel)
{
double sum = 0.0;
for(int i = 0; i < Nc; i++)
  sum = sum + myModel.wt[i];
std::cout<<"sum of weights = "<<sum<<std::endl;

sum = 0.0;
for(int i = 0; i < Nc; i++)
  sum = sum + myModel.ci[i][0]*myModel.wt[i];
std::cout<<"sum of wi*ci in x direction= "<<sum<<std::endl;

sum = 0.0;
for(int i = 0; i < Nc; i++)
  sum = sum + myModel.ci[i][1]*myModel.wt[i];
std::cout<<"sum of wi*ci in y direction= "<<sum<<std::endl;

sum = 0.0;
for(int i=0;i<Nc;i++)
  sum = sum + (myModel.ci[i][0]*myModel.ci[i][0] + myModel.ci[i][1]*myModel.ci[i][1]) * myModel.wt[i];
std::cout<<"sum of wi*ci^2 = "<<sum / myModel.theta0<<std::endl;
}

void setModelParameters(modelD2Q9 &myModel)
{
 myModel.ci[center][0] =  0.0;
 myModel.ci[center][1] =  0.0;

 myModel.ci[east][0]   =  1.0;
 myModel.ci[east][1]   =  0.0;

 myModel.ci[north][0]  =  0.0;
 myModel.ci[north][1]  =  1.0;

 myModel.ci[west][0]   =  -1.0;
 myModel.ci[west][1]   =  0.0;

 myModel.ci[south][0]  =  0.0;
 myModel.ci[south][1]  =  -1.0;

 myModel.ci[neast][0]  =  1.0;
 myModel.ci[neast][1]  =  1.0;

 myModel.ci[nwest][0]  =  -1.0;
 myModel.ci[nwest][1]  =  1.0;

 myModel.ci[seast][0]  = 1.0;
 myModel.ci[seast][1]  =  -1.0;

 myModel.ci[swest][0]  =  -1.0;
 myModel.ci[swest][1]  =  -1.0;


 myModel.wt[center] = 4.0/9.0;
 myModel.wt[north]  = 1.0/9.0;
 myModel.wt[east]   = 1.0/9.0;
 myModel.wt[west]   = 1.0/9.0;
 myModel.wt[south]  = 1.0/9.0;
 myModel.wt[neast]  = 1.0/36.0;
 myModel.wt[nwest]  = 1.0/36.0;
 myModel.wt[seast]  = 1.0/36.0;
 myModel.wt[swest]  = 1.0/36.0;


 myModel.theta0       = 1.0/3.0;
 myModel.onebytheta0  = 3.0;
}

class grid{
public:
grid();

 int NB, NE;
 int G1, G2;
 int a;
 int data[(N + 2) * (N + 2) * Nc];
 int operator()(int i,int j, int k)const {return data[(i * (N+2) + j) * 9 + k];}
 int & operator()(int i,int j, int k) {return data[(i * (N+2) + j) * 9 + k];}



void initializeIntegers()
{
    for(int i = 0; i <= N + 1 ; i++)
        for(int j = 0; j <= N + 1 ; j++)
            for(int k = 0; k < Nc; k++)
                data[((i * N + j) * 9 + k)] = 0.0;
}

void printGrid(int population)
{
    for(int i = 0; i < N ; i++){
        for(int j = 0; j < N ; j++)
            cout<<data[(i * N + j) * 9 + population]<<" ";
        cout<<endl;
    }
}
//Row number 0 and Nx + 1 are ghost locations
//Column number 0 and Ny + 1 are ghost locations

};

grid::grid()
{
    NB = 1;
    NE = N;
    G1 = 0;
    G2 = NE + 1;
}

// Advection
void advect(grid &gridLB)
{
//East
   for(int i = gridLB.NB; i <= gridLB.NE; i++)
       for(int j = gridLB.NE; j >= gridLB.NB; j--)
           gridLB(i,j,east) = gridLB(i,j-1,east);

//North
   for(int i = gridLB.NE; i >= gridLB.NB; i--)
       for(int j = gridLB.NB; j <= gridLB.NE; j++)
           gridLB(i,j,north) = gridLB(i-1,j,north);

//West
   for(int i = gridLB.NB; i <= gridLB.NE; i++)
       for(int j = gridLB.NB; j <= gridLB.NE; j++)
           gridLB(i,j,west) = gridLB(i,j+1,west);

//South
   for(int i = gridLB.NB; i <= gridLB.NE ; i++)
       for(int j = gridLB.NB; j <= gridLB.NE; j++)
           gridLB(i,j,south) = gridLB(i+1,j,south);

//North-East
   for(int i = gridLB.NE ; i >= gridLB.NB; i--)
       for(int j = gridLB.NE; j >= gridLB.NB; j--)
           gridLB(i,j,neast) = gridLB(i-1,j-1,neast);

//North-West
   for(int i = gridLB.NE ; i >= gridLB.NB; i--)
       for(int j = gridLB.NB; j <= gridLB.NE; j++)
           gridLB(i,j,nwest) = gridLB(i-1,j+1,nwest);

//South-West
   for(int i = gridLB.NB ; i <= gridLB.NE ; i++)
       for(int j = gridLB.NB; j <= gridLB.NE; j++)
           gridLB(i,j,swest) = gridLB(i+1,j+1,swest);

//South-East
   for(int i = gridLB.NB ; i <= gridLB.NE; i++)
       for(int j = gridLB.NE; j >= gridLB.NB; j--)
           gridLB(i,j,seast) = gridLB(i+1,j-1,seast);
}

//converts f's to rho and u1
//final values get stored in modelD2Q9
void getMoments(modelD2Q9 &model/*out*/, grid &gridLB/*in*/, int row, int col, double &rho, double &u1x, double &u1y)
{
    //density = f0 + f1 + ...+ f8
    rho = 0.0;
    for(int k = 0; k < Nc; k++)
        rho += gridLB(row,col,k) ;

    //x-momentum = f0cx0 + f1cx1 ... + f8cx8
    u1x = 0.0;
    for(int k = 0; k < Nc; k++)
        u1x += (gridLB(row,col,k) * model.ci[k][0]);

    u1x = u1x / rho ;
    //y-momentum = f0cy0 + f1cy1 ... + f8cy8
    u1y = 0.0;
    for(int k = 0; k < Nc; k++)
        u1y += (gridLB(row,col,k) * model.ci[k][1]);

    u1y = u1y / rho ;
}

//converts rho and u1 to fEq's at a node
//Feq equation is using only linear terms
void getFeq(modelD2Q9 &mymodel)
{
 for(int i = 0; i < Nc; i++) {
     double sum = 0.0;
     sum = (mymodel.ci[i][0]*mymodel.u1x +  mymodel.ci[i][1]*mymodel.u1y);
     mymodel.fEq[i] = mymodel.wt[i] * mymodel.rho * (1.0 + 3 * sum + 4.5 * sum * sum
                      -1.5 * (mymodel.u1x * mymodel.u1x + mymodel.u1y*mymodel.u1y));
 }
}

//copy feq value to data
void copyToGrid(modelD2Q9 &mymodel/*in*/, grid &mygrid/*out*/, int row, int col)
{
 for(int i = 0; i < Nc; i++)
    mygrid(row,col,i) = mymodel.fEq[i];
}

//Collision
//f = f + 2*beta*(fEq - f)

/*With gravity
gridLB(i,j,k]+=2.0*beta*(myModel.fEq[k]-gridLB(i,j,k])+
(2.0*beta*tow*myModel.rho[i,j]*myModel.wt[k]*3.0*(myModel.g_y*myModel.ci_y[k]));
*/
void collide(modelD2Q9 &model, grid &gridLB, double beta)
{
    for (int i = 1; i <= N; i++)
    {
        for(int j = 1; j <= N; j++)
        {
             getMoments(model, gridLB, i, j, model.rho, model.u1x, model.u1y);
             getFeq(model);
            for(int k = 0; k < Nc; k++)
            {
                gridLB(i,j,k) = gridLB(i,j,k) + 2.0*beta*(model.fEq[k] - gridLB(i,j,k));
            }
        }
    }
}



//The function simply copies the velocities at the edge nodes to corresponding ghost nodes
void prepareBC(modelD2Q9 &myModel, grid &gridLB)
{
//    //EAST WEST PERIODIC
//    for(int i = 0; i <= N+1; i++)
//      for(int dv = 0; dv < Nc ; dv++)
//    {
//	gridLB(i,0,dv] = gridLB(i,1,dv];
//    }
//    for(int i = 0; i <= N+1; i++)
//      for(int dv = 0; dv < Nc ; dv++)
//    {
//	gridLB(i,N+1,dv] = gridLB(i,N,dv];
//    }
////    //TOP BOTTOM PERIODIC
//    for(int j = 0; j <= N+1; j++)
//      for(int dv = 0; dv < Nc ; dv++)
//    {
//	gridLB(0,j,dv] = gridLB(1,j,dv];
//    }
//    for(int j = 0; j <= N+1; j++)
//      for(int dv = 0; dv < Nc ; dv++)
//    {
//	gridLB(N+1,j,dv] = gridLB(N,j,dv];
//    }

//EAST-WEST AND TOP-BOTTOM PERIODIC
//    for(int i = 0; i <= N+1; i++)
//      for(int dv = 0; dv < Nc ; dv++)
//    {
//	gridLB(i,0,dv] = gridLB(i,1,dv];
//	gridLB(i,N+1,dv] = gridLB(i,N,dv];
//	gridLB(0,i,dv] = gridLB(1,i,dv];
//	gridLB(N+1,i,dv] = gridLB(N,i,dv];
//    }

////Bottom wall - copy to ghost
// for (int x = 1; x <= N; x++) {
//     gridLB(0,x,south] = gridLB(1,x,south];
//     gridLB(0,x,seast] = gridLB(1,x,seast];
//     gridLB(0,x,swest] = gridLB(1,x,swest];
// }
//
// //Left wall - copy to ghost
// for (int x = 1; x <= N; x++) {
//     gridLB(x,0,west] = gridLB(x,1,west];
//     gridLB(x,0,nwest] = gridLB(x,1,nwest];
//     gridLB(x,0,swest] = gridLB(x,1,swest];
// }
//
// //Right wall - copy to ghost
// for (int x = 1; x <= N; x++) {
//     gridLB(x,N + 1,east] = gridLB(x,N,east];
//     gridLB(x,N + 1,seast] = gridLB(x,N,seast];
//     gridLB(x,N + 1,neast] = gridLB(x,N,neast];
// }
//
////Top wall - copy to ghost
// for (int x = 1; x <= N; x++) {
// gridLB(N + 1,x,north] = gridLB(N,x,north];
// gridLB(N + 1,x,nwest] = gridLB(N,x,nwest];
// gridLB(N + 1,x,neast] = gridLB(N,x,neast];
// }

//Bottom wall - copy to ghost
 for (int x = 1; x <= N; x++) {
     gridLB(0,x,south) = gridLB(1,x,south);
     gridLB(0,x,seast) = gridLB(1,x,seast);
     gridLB(0,x,swest) = gridLB(1,x,swest);
//Left wall - copy to ghost
     gridLB(x,0,west) = gridLB(x,1,west);
     gridLB(x,0,nwest) = gridLB(x,1,nwest);
     gridLB(x,0,swest) = gridLB(x,1,swest);
//Right wall - copy to ghost
     gridLB(x,N + 1,east) = gridLB(x,N,east);
     gridLB(x,N + 1,seast) = gridLB(x,N,seast);
     gridLB(x,N + 1,neast) = gridLB(x,N,neast);
//Top wall - copy to ghost
     gridLB(N + 1,x,north) = gridLB(N,x,north);
     gridLB(N + 1,x,nwest) = gridLB(N,x,nwest);
     gridLB(N + 1,x,neast) = gridLB(N,x,neast);
     }


}


void applyBC(modelD2Q9 &myModel, grid &gridLB)
{

//Bounce back at lower wall
    for(int j = 1; j <= N; j++){
        gridLB(gridLB.NB,j,north) = gridLB(gridLB.G1,j,south);
        gridLB(gridLB.NB,j,nwest) = gridLB(gridLB.G1,j,seast);
        gridLB(gridLB.NB,j,neast) = gridLB(gridLB.G1,j,swest);
    }

//Bounce back at Left wall,
    for(int i = 1; i <= N; i++){
        gridLB(i,gridLB.NB,east) = gridLB(i,gridLB.G1,west);
        gridLB(i,gridLB.NB,neast) = gridLB(i,gridLB.G1,swest);
        gridLB(i,gridLB.NB,seast) = gridLB(i,gridLB.G1,nwest);
    }

//Bounce back at Right wall
    for(int i = 1; i <= N; i++){
        gridLB(i,gridLB.NE,west) = gridLB(i,gridLB.G2,east);
        gridLB(i,gridLB.NE,nwest) = gridLB(i,gridLB.G2,seast);
        gridLB(i,gridLB.NE,swest) = gridLB(i,gridLB.G2,neast);
    }


//Bounce back at upper wall
//    for(int j = 1; j <= N; j++){
//        gridLB(gridLB.NE,j,south) = gridLB(gridLB.G2,j,north];
//        gridLB(gridLB.NE,j,swest] = gridLB(gridLB.G2,j,neast];
//        gridLB(gridLB.NE,j,seast] = gridLB(gridLB.G2,j,nwest];
//    }

    double factor = 0.0;

//Moving wall at top wall
     for(int j = gridLB.NB; j <= gridLB.NE; j++){
         myModel.rho = 1.0;
         myModel.u1x = 0.01;
         myModel.u1y = 0.0;
         getFeq(myModel);
         factor = (gridLB(gridLB.G2,j,north) +
                   gridLB(gridLB.G2,j,neast) + gridLB(gridLB.G2,j,nwest)) /
                  (myModel.fEq[north] + myModel.fEq[neast] + myModel.fEq[nwest]);
         gridLB(gridLB.NE,j,south) = myModel.fEq[south] * factor ;
         gridLB(gridLB.NE,j,seast) = myModel.fEq[seast] * factor ;
         gridLB(gridLB.NE,j,swest) = myModel.fEq[swest] * factor ;

     }
}

void initialConditions(modelD2Q9 &mymodel, grid &gridLB)
{
mymodel.rho = 1.0;
mymodel.u1x = 0.01;
mymodel.u1y = 0.0;
getFeq(mymodel);
 for(int i = 1; i <= N; i++)
    for(int j = 1; j <= N; j++)
        copyToGrid(mymodel, gridLB, i, j);

//mymodel.rho = 0.0001;
//mymodel.u1x = 0.0;
//mymodel.u1y = 0.0;
//getFeq(mymodel);
//for(int i = N/2; i <= N + 1; i++)
// for(int j = 0; j <= N + 1; j++)
//     copyToGrid(mymodel, gridLB, i, j);
}


//These functions are just for checking that moments are
//being calculated correctly. They don't print values being
//calculated in the program.
void print_density(modelD2Q9 &myModel, grid &gridLB)
{
      double rho, ux, uy;
    for(int i = 0; i <= N + 1; i++)
    {
        for(int j = 0; j <= N + 1; j++)
        {
            getMoments(myModel, gridLB, i, j, rho, ux,uy);
            std::cout<<rho<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print_u1x(modelD2Q9 &myModel, grid &gridLB)
{
    double rho, ux, uy;
  for(int i = 1; i <= N; i++)
    {
        //for(int j = 0; j <= N + 1; j++)
        //{
            getMoments(myModel, gridLB, i, (int)(N/2) ,rho, ux,uy);
            cout<<setprecision(15)<<ux<<"\n";
        //}
        //std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print_u1y(modelD2Q9 &myModel, grid &gridLB)
{
  double rho, ux, uy;
  for(int i = 0; i <= N + 1; i++)
    {
        for(int j = 0; j <= N + 1; j++)
        {
            getMoments(myModel, gridLB, i,j, rho, ux,uy);
            std::cout<<uy<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print_feq(modelD2Q9 &myModel, grid &gridLB, int a)
{
    for(int i = 0; i <= N + 1; i++)
    {
        for(int j = 0; j <= N + 1; j++)
        {
//             getFeq(myModel);
            std::cout<<myModel.fEq[a]<<"\t";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void check_mass(grid gridlbm)
{
    double sum = 0.0;
    for(int i = 1; i <= N ; i++)
    {
        for(int j = 1; j <= N ; j++)
        {
            for(int k = 0; k < Nc; k++)
            {
                sum += gridlbm(i,j,k);
            }
        }
    }
    cout<<setprecision(15) <<sum<<"\n";
}


void write_density(modelD2Q9 &myModel, grid &gridLB,int step)
{
 //cout<<"Write_density function called"<<endl;
 char filename[150];
 sprintf(filename, "results/density/density_%d.txt",step);
 std::ofstream file1;
 file1.open(filename);
 for(int i = gridLB.NB;i <= gridLB.NE; i++)
    for(int j = gridLB.NB; j <= gridLB.NE; j++)
 {
  getMoments(myModel,gridLB,i,j,myModel.rho, myModel.u1x, myModel.u1y);
  file1<<i<<" "<<j<<"   "<<myModel.rho<<std::endl;
  }
 file1.close();
}

void write_velocity(modelD2Q9 &myModel, grid &gridLB,int step)
{
 //cout<<"Write_density function called"<<endl;
 char filename[150];
 sprintf(filename, "results/velocity/velocity_%d.txt",step);
 std::ofstream file1;
 file1.open(filename);
 for(int i = gridLB.NB;i <= gridLB.NE; i++)
    for(int j = gridLB.NB; j <= gridLB.NE; j++)
 {
  getMoments(myModel,gridLB,i,j,myModel.rho, myModel.u1x, myModel.u1y);
  file1<<j<<" "<<i<<"   "<<myModel.u1x * 500<<" "<<myModel.u1y * 500<<
  " "<<sqrt(myModel.u1x * myModel.u1x + myModel.u1y * myModel.u1y) <<std::endl;
  }
 file1.close();
}

int main()
{
modelD2Q9 lbmodel;
setModelParameters(lbmodel);
grid gridLBM ;
gridLBM.initializeIntegers();
initialConditions(lbmodel, gridLBM);

double cs2 = 1.0/3.0;
double cs =  sqrt(cs2);
double refLen = 1;
double Ma = 0.01;
double uBoundary = Ma * cs;

double Re = 100;
double Kn = Ma/Re;
double kinVisc = uBoundary*refLen/Re;
double tau = kinVisc/cs2;
double dt = refLen/N;
double tauNdim = tau/dt;
double beta = 1.0 / (2.0*tauNdim + 1.0);
int convectionTime = N/uBoundary;
int diffusionTime  = 1.0/kinVisc;
int simulationTime = 3 * diffusionTime;


cout<<"Simulation time = "<<simulationTime<<endl;
cout<<"Press enter to continue"<<endl;
getch();
cout<<"Process started"<<endl;
for(int time = 0; time < simulationTime; time++)
    {
        collide(lbmodel, gridLBM, beta);
        prepareBC(lbmodel, gridLBM);
        advect(gridLBM);
        applyBC(lbmodel, gridLBM);
        if(time % 5000 == 0){
            cout<<time<<endl;
            write_velocity(lbmodel, gridLBM, time);
        }

    }
    print_u1x(lbmodel, gridLBM);
    write_velocity(lbmodel, gridLBM, 999);
//    write_density(lbmodel, gridLBM, 999);
return 0;
}
