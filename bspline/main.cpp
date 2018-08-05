#include <iostream>
#include <fstream>


using namespace std;

 void bspline( double p0[2], double p1[2], double p2[2], double p3[2], double a[25][2])
 {
     int n = 25;

     double dt = 1.0 / n ;
     double t = dt;
     int k = 0;

     cout<<dt<<" " <<t<< " "<<k;

     while(k < 25)
     {
         a[k][0] = (1-t)*(1-t)*(1-t)*p0[0] + 3*t*(1-t)*(1-t)*p1[0] + 3*t*t*(1-t)*p2[0] + t*t*t*p3[0];
         a[k][1] = (1-t)*(1-t)*(1-t)*p0[1] + 3*t*(1-t)*(1-t)*p1[1] + 3*t*t*(1-t)*p2[1] + t*t*t*p3[1];
         t = t + dt;
         k++;
     }
 }

int main()
{
    cout << "Hello world!" << endl;

    double a[25][2];
    double p0[2];
    double p1[2];
    double p2[2];
    double p3[2];

    p0[0] = 0.1;
    p0[1] = 3.2;

     p1[0] = 1.1;
    p1[1] = 7.8;

     p2[0] = 5.6;
    p2[1] = 9.0;

     p3[0] = 9.0;
    p3[1] = 1.1;

    bspline(p0, p1, p2, p3, a);

    for(int i = 0; i < 25; i++)
        cout<<a[i][0]<<"  "<<a[i][1]<<endl;

    ofstream fout;
     fout.open("curve.txt");

     for(int i = 0; i < 25; i++)
        fout<<a[i][0]<<"  "<<a[i][1]<<endl;

     fout.close();

     ofstream fout2;
     fout2.open("points.txt");
        fout2<<p0[0]<<" "<<p0[1]<<endl;
        fout2<<p1[0]<<" "<<p1[1]<<endl;
        fout2<<p2[0]<<" "<<p2[1]<<endl;
        fout2<<p3[0]<<" "<<p3[1]<<endl;

     fout2.close();


    return 0;
}
