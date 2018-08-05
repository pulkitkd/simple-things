#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
    double dt = 0.001;
    double t = 0;
    double T = 0;
    double y0 = 10;
    double y1 = 0;
    double lambda = 1;
    double exact = 0.0;
    double absolute_error = 0.0;
    double relative_error = 0.0;
    ofstream fout;
    fout.open("solution.txt");
    cout<<"Enter final time"<<endl;
    cin>>T;

    while(t < T){
    y1 = y0 - lambda * y0 * dt;
    exact = y0 * exp(-lambda * t);
    absolute_error = y1 - exact;
    relative_error = absolute_error / exact;
    //cout<<t<<" "<<y1<<" "<<exact<<endl;
    fout<<t<<" "<<y1<<" "<<exact<<"  "<< absolute_error<<" "<<relative_error<<endl;
    y0 = y1;
    t = t + dt;
    }
    fout.close();
    return 0;


}
