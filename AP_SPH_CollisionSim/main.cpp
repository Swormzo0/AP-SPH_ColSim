#include <iostream>
#include <stdio.h>
//#include <windows.h>
//#include <iomanip>
#include <cmath>

using namespace std;

// Constants
// Timestep
const double dt = 0.0001;

// Convenience variables for coordianates
const int X = 0;
const int Y = 1;

// Convenience variables for identifying balls
const int A = 0;
const int B = 1;

// Global variables
double M, R; // Mass and radius
bool isStressBall = false; // Whether to use mass and radius of stress ball

class Ball
{
public:
    // Set up kinematics variables
    double s[2], v[2], a[2];
    // Declare net force on ball
    double F[2];
    // Individual mass and radius
    double m, r;

// Update displacement, velocity and acceleration
    void updateKM()
    {
        // Update acceleration
        for(int i=0; i < 2; i++)
            a[i] = F[i] / m;
        // Update position
        for(int i=0; i < 2; i++)
            // Factored version of the position equation for speed
            s[i] = s[i] + (v[i]+.5*a[i]*dt)*dt;
        // Update velocity
        for(int i=0; i < 2; i++)
            v[i] = v[i] + a[i]*dt;
    }

    // Constructor
    Ball(double m, double r, double[] v, double[] s)
    {
        // Set variables
        this->m = abs(m);
        this->r = abs(r);
        this->v = v;
        this->s = s;
    }
};

class BSystem
{
    Ball b[2];

    bool displacement(double &xA[])
    {
        // Prematurely calculate displacement from equilibrium wrt ball 1
        double ds[2];
        for(int i=0; i < 2; i++)
            ds[i] = b[0].s[i] - b[0].s[1];

        if (sqrt(ds[0]^2+ds[1]^2) < 2*R)
            xA[X] = xA[Y] = 0;
        else
            for(int i=0; i < 2; i++)
                xA[i] = .5 * ds[i] - abs(ds[i]) * R/ds[i];
    }

    void ballForce(double x[], double &F[])
    {
        for(int i=0; i < 2; i++)
            F[i] = x * (x* (2534468*x-20795.8)+2127.939)
    }

public:
    void updateForce(Ball b)
    {

    }
};

int main()
{
    double m1, m2;
    // Angle of first ball relative to second ball
    // Must be in (-90, 90)in degrees
    double theta;
    double v1[2];
    double v2[2];

    // Prompt for variables
    cout<<"Mass 1:";
    cin>>m1;
    cout<<"Mass 2:";
    cin>>m2;

    cout<<"Co-ordinates of ball 2 relative to ball 1:"<<endl<<endl;
    cin>>theta;
    cout<<endl;

    cout<<"Components of initial velocity of mass 1:"<<endl<<endl;
    for (int m=0; m<2; m++)
    {
        cin>>v1[m];
        cout<<endl;
    }

    cout<<"Components of initial velocity of mass 2:"<<endl<<endl;
    for (int m=0; m<2; m++)
    {
        cin>>v2[m];
        cout<<endl;
    }
}
