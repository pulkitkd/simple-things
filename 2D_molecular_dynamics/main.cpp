#include <iostream>
#include "matrix.h"
#include <fstream>
#include <time.h>
#include <math.h>
#include <conio.h>

#define N 10 //no. of particles
#define D 2 //no. of spatial dimensions
#define PI 3.14159265359

using namespace std;

void get_distances(matrix<double, N, D> &pos, matrix<double, N, N> &dist, double &L)
{
    Null(dist);

    for(int i = 0; i < N - 1; i++)
        for(int j = i + 1; j < N; j++)
            {
                for(int k = 0; k < D; k++)
                    {
                        dist(i, j) += (pos(i,k)-pos(j,k)) * (pos(i,k)-pos(j,k));
                    }
                dist(j, i) = dist(i, j);
            }

}

void get_distances_2(matrix<double, N, D> &pos, matrix<double, N, N> &dist, double &L)
{
    Null(dist);

    for(int i = 0; i < N - 1; i++)
        for(int j = i + 1; j < N; j++)
            {
                for(int k = 0; k < D; k++)
                    {
                        if(abs(pos(i,k)-pos(j,k)) > L/2)
                            dist(i, j) += (L - (pos(i,k)-pos(j,k))) * (L - (pos(i,k)-pos(j,k)));
                        else
                            dist(i, j) += (pos(i,k)-pos(j,k)) * (pos(i,k)-pos(j,k));
                    }
                dist(j, i) = dist(i, j);
            }

}


void get_acceleration(matrix<double, N, D> &pos, matrix<double, N, D> &acc,
                      matrix<double, N, N> &dist, double &L)
{
    //matrix<double, N, D> force;
    Null(acc);
    double d8 = 0;
    double d14 = 0;
    double factor = 0;
    double Rc = 9; /* Square of Rc, since we are comparing it with dist^2 */

    get_distances_2(pos, dist, L);

    for(int i = 0; i < N - 1; i++)
    {
        for(int k = i + 1; k < N; k++)
        {
            if(dist(i, k) < Rc)
            {
                d8 = pow(dist(i, k), -4);
                d14 = pow(dist(i, k), -7);
                factor = 48*(d14 - 0.5 * d8);
                for(int j = 0; j < D; j++)
                {
                    acc(i, j) += factor * (pos(i,j)-pos(k,j));
                    acc(k, j) -= factor * (pos(i,j)-pos(k,j));
                }
             }
        }
    }
}

void get_acceleration_2(matrix<double, N, D> &pos, matrix<double, N, D> &acc,
                      matrix<double, N, N> &dist, double &L)
{
    //matrix<double, N, D> force;
    Null(acc); //set all entries to zero
    double d8 = 0;
    double d14 = 0;
    double factor = 0;
    double Rc = 9; /* Square of Rc, since we are comparing it with dist^2 */

    get_distances_2(pos, dist, L);

    for(int i = 0; i < N; i++)
    {
        for(int k = 0; k < N; k++)
        {
            if(dist(i,k) < Rc)
            {
                if(i!=k)
                {
                d8 = pow(dist(i, k), -4);
                d14 = pow(dist(i, k), -7);
                factor = 48*(d14 - 0.5 * d8);
                for(int j = 0; j < D; j++)
                {
                    if(abs(pos(i,j)-pos(k,j)) > L/2)
                        acc(i, j) += factor * (pos(k,j) - pos(i,j));
                    else
                        acc(i, j) -= factor * (pos(k,j) - pos(i,j));

                    //acc(k, j) -= factor * (pos(i,j)-pos(k,j));
                }
                }
                //cout<<"***"<<i<<" "<<k<<" "<<acc(i,k)<<endl;
             }

        }

    }
        //cout<<acc;
}


void write_position(matrix<double, N, D> position, int a)
{
    char filename[150];
    sprintf(filename, "results/positions/position_%d.txt",a);
    std::ofstream file1;
    file1.open(filename);
    for(int i = 0; i < N; i++)
    {
        for(int j = 0; j < D; j++)
        {
        file1<<position(i,j)<<"\t";
        }
        file1<<i<<endl;
    }
    file1.close();
}

void periodic_boundary_condition(matrix<double, N, D> &position, double L)
{
    for(int i = 0; i < N; i++)
        for(int j = 0; j < D; j++)
            {
                if(position(i,j) > 0.99*L)
                    while(position(i,j) > 0.99*L)
                        position(i,j) = position(i,j) - L;
                if(position(i,j) < 0.01)
                    while(position(i,j) < 0.01)
                        position(i,j) = position(i,j) + L;
            }
}

double min_dist(matrix<double, N, N> &dist)
{
    double a = dist(0,1);
    for(int i = 0; i < N; i++)
        for(int j = i + 1; j < N; j++)
            if(dist(i,j) < a)
                a = dist(i,j);
    return a;
}

double get_temperature(matrix<double, N, D> velocity, double temp)
{
    double T;
    for(int i = 0; i < N; i++)
        for(int j = 0; j < D; j++)
    {
        T += velocity(i,j) * velocity(i,j);
    }

    if(T < 1e4)
    return T/(3 * (N - 1));

    else
        return temp;
}


int main()
{
    srand(time(NULL));

    matrix<double, N, D> position_new;
    matrix<double, N, D> position_old;
    matrix<double, N, D> position;
    matrix<double, N, D> velocity;
    matrix<double, N, D> acceleration;
    matrix<double, N, N> dist;

    //Potential :
    // _A_  - _B_
    // r^12   r^6

    double L = 20;
    int Np = (int)(1 * L);
    int Nv = (int)(2 * L);
    double dt = 0.0001; //time step
    double dt1 = 0.01; //image printing step
    double dtsq = dt * dt;
    double onebydt = 1 / dt;
    double T = 10;
    double t = 0;
    double t1 = 0;
    int a = 0;
    double dmin = 0;
    double tempcutoff = 500;
    double temp = 0;
    double s = 0;

    //ensure that initial configuration is not too dense
    while(dmin < 3)
        {
            Null(dist);
            cout<<"minimum distance is "<<dmin<<" Reassigning positions."<<endl;
            random(position_old, Np);
            //Put particles roughly in the center
            for(int i = 0; i < N; i++)
                for(int j = 0; j < D; j++)
                    position_old(i, j) += 0.5 * (L - Np);

            get_distances_2(position_old, dist, L);
            dmin = min_dist(dist);
        }

    cout<<dist;
    cout<<"The minimum distance is "<<min_dist(dist)<<endl;
    cout<<"The positions are given by "<<endl<<position_old;

    //initialize with random velocities
    random(velocity, Nv);
    for(int i = 0; i < N; i++)
        for(int j = 0; j < D; j++)
            velocity(i, j) = velocity(i, j) - Nv / 2;

    cout<<velocity;

    ofstream fout;
    fout.open("results/temperature/temperature.txt");

    ofstream fout2;
    fout2.open("results/velocity/velocity.txt");

    get_acceleration_2(position_old, acceleration, dist, L);
    velocity = velocity + dt * acceleration;
    position = position_old + dt * velocity + 0.5 * dtsq * acceleration;
    periodic_boundary_condition(position, L);

    cout<<"position_old "<<endl<<position_old;
    cout<<"position "<<endl<<position;
    cout<<"velocity "<<endl<<velocity;
    cout<<"acceleration "<<endl<<acceleration;
    cout<<"distances "<<endl<<dist;
    cout<<"The minimum distance is "<<min_dist(dist);
    getch();

    //velocity verlet scheme
    while(t < T)
    {
        if((t1 - t) < 1e-8)
        {
            cout<<"Write file at time = "<<t<<endl;
            write_position(position, a);
            a++;
            t1 = t1 + dt1;
        }
//         velocity is not needed for iterations of verlet scheme
//         we compute velocity to calculate temperature
//         this calculation can be omitted to save memory
        position_new = 2 * position - position_old + dtsq * acceleration;
        velocity = 0.5 * onebydt * (position_new - position_old);
        get_acceleration_2(position_new, acceleration, dist, L);
        periodic_boundary_condition(position_new, L);

        position_old = position;
        position = position_new;

        t = t + dt;

        temp = get_temperature(velocity, tempcutoff);

        fout<<t<<" "<<temp<<endl;
    }
    for(int i = 0; i < N; i++)
        {
            s = 0;
            for(int j = 0; j < D; j++)
                s += velocity(i,j)*velocity(i,j);
            fout2<<i<<" "<<s<<endl;
        }
    fout.close();
    fout2.close();

    return 0;
}
