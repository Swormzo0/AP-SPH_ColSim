#include <iostream>
#include <stdio.h>
#include <cmath>

// Program controls
#define USER_MASS // Allow the user to define mass for each ball?
#define DEBUG // Provide a debug table with t, s, v, a and F
//#define DEBUG_POSN // Debug the position from equilibrium
//#define CUBIC_F // Using the cubic regression? (and not quadratic?)

using namespace std;

// Configuration
const bool isStressBall = false; // Whether to use mass and radius of stress ball
double dt = 0.0000001; // Timestep
const double stopDist = 0.0001; // Stop when the balls are this far from each other
const double timeOut = 21; // How long before giving up on simulation
const double minTSPerc = .001; // Largest percentage of radius length permitted for first timestep
const double dispPrec = 15; // How many decimal places to show in debug table
// displacement from equilibrium

// Debug config
bool debugBall = 0; // Debug ball A
const unsigned int updateTableT = 500; // Period for table updates

// Constants
const double M = .25, R = 3.5, TWO_R = R * 2; // Default mass and radius of ball
const double PI = 3.14159265358979323846264338;
const double F_CUBIC[] = {2534468.26698494, -20795.82059, 2127.93853}; // Coefficients of cubic equation
const double F_QUAD[] = {48234.5591, 1699.472766}; // Coefficients of quadratic equation

// Convenience variables for coordinates
const int X = 0;
const int Y = 1;

// Convenience variables for identifying balls
const int A = 0;
const int B = 1;


// Helper functions
// Return the magnitude of a 2D vector
double norm(double (&a)[2])
{
    return sqrt(a[0]*a[0]+a[1]*a[1]);
}


// Calculate angle in degrees with the domain [0, 360)
double angleDeg(double y, double x)
{
    double theta = atan2(y, x) * 180 / PI;
    return theta < 0? theta + 360: theta;
}


/*// Copy the sign of one variable to another
double signCopy(double orig, double dest)
{
    if(orig < 0)
        return -dest;
    if(orig == 0)
        return 0;
    return dest;
}*/


// A stress ball class
class Ball
{
public:
    // Set up kinematics variables
    double s[2], v[2], a[2];
    // Declare net force on ball
    double F[2];
    // Individual mass
    double m;

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

    // Constructors
    // Default Constructor
    Ball()
    {
        // Set the mass to the default one and leave everything else at 0
        m = M;
        for(int i=0; i < 2; i++)
        {
            s[i] = 0;
            v[i] = 0;
            a[i] = 0;
            F[i] = 0;
        }
    }

    // Construct a ball in a given circumstance
    Ball(double m, double (&v)[2], double (&s)[2])
    {
        // Set variables
        this->m = abs(m);
        for(int i=0; i < 2; i++)
        {
            this->v[i] = v[i];
            this->s[i] = s[i];
        }

        // Set other variables to zero
        for(int i=0; i < 2; i++)
        {
            a[i] = 0;
            F[i] = 0;
        }
    }
};


// A class for the system of two balls
class BSystem
{
    Ball b[2]; // Two balls
    double dist, t=0, vi[2][2]; // Distance, time and initial velocity, respectively
    bool hasEverCollided = false; // Record: Was there any collision?
#ifdef DEBUG
    unsigned int iterNum = 0; // Simulation iteration index
#endif // DEBUG


    // Calculate displacement from equilibrium and update distance from each other with respect to the first ball
    void displacement(double (&xA)[2])
    {
        double ds[2];
        bool isInCollision;

        // Calculate difference in position
        for(int i=0; i < 2; i++)
            ds[i] = b[0].s[i] - b[1].s[i];

        // Calculate distance and see if the balls are in collision
        dist = norm(ds);
        isInCollision = dist < TWO_R;

        // Update whether the ball has ever collided
        hasEverCollided |= isInCollision;

        // If there is no collision, displacement from equilibrium is 0
        if (!isInCollision)
            xA[X] = xA[Y] = 0;
        // Otherwise
        // Calculate distance from equilibrium
        else
            for(int i=0; i < 2; i++)
                xA[i] = ds[i] * (-.5 + R / dist);

#if defined(DEBUG) && defined(DEBUG_POSN)
        // Calculate x using Igi's formula
        double x = R - .5 * dist;

        // If time to update
        if(iterNum % updateTableT == 0)
        {
            // Display the displacement from equilibrium wrt ball 1
            printf("\t%.5f\t%.5f", xA[X], xA[Y]);
            // Display difference in magnitudes of x
            printf("\t%.5f", x-norm(xA));
        }
#endif // DEBUG
    }


    // Calculate the force for a given distance from equilibrium
    void ballForce(double (&x)[2], double (&F)[2])
    {
        double newX;

        for(int i=0; i < 2; i++)
        {
            // The regressed force requires a negative x
            newX = -abs(x[i]);

            // Negate force as to oppose direction of displacement
#ifdef CUBIC_F
            F[i] = -x[i] * (newX*(F_CUBIC[0]*newX+F_CUBIC[1])+F_CUBIC[2]);
#else
            F[i] = -x[i] * (F_QUAD[0]*newX+F_QUAD[1]);
#endif
        }
    }


    // Update forces on balls for the time step
    void updateForces()
    {
        // Calculate displacement from equilibrium for ball A
        double xA[] = {0, 0};
        displacement(xA);

        // Update force on ball A
        ballForce(xA, b[A].F);

        // Effect equal and opposite force on ball B
        for(int i=0; i < 2; i++)
            b[B].F[i] = -b[A].F[i];
    }


    // Update the state and increment time step
    void update()
    {
        // Debug stuff
#ifdef DEBUG
        // If time to update, display next row of debug table
        if(iterNum % updateTableT == 0)
        {
            printf("%d\t%.*f", iterNum, dispPrec, t);
            printf("\t%.*f\t%.*f", dispPrec, b[debugBall].s[X], dispPrec, b[debugBall].s[Y]);
            printf("\t%.*f\t%.*f", dispPrec, b[debugBall].v[X], dispPrec, b[debugBall].v[Y]);
            printf("\t%.*f\t%.*f", dispPrec, b[debugBall].a[X], dispPrec, b[debugBall].a[Y]);
            printf("\t%.*f\t%.*f", dispPrec, b[debugBall].F[X], dispPrec, b[debugBall].F[Y]);
        }
#endif // DEBUG

        // Evaluate forces
        updateForces();

        // Calculate new a, v and s
        for(int i=0; i < 2; i++)
            b[i].updateKM();

        // Debug table
#ifdef DEBUG
        // End row for debug table when needed
        if(iterNum % updateTableT == 0)
            printf("\n");
        // Increment iteration number
        iterNum++;
#endif // DEBUG

        // Increment time step
        t += dt;
    }


    // Return whether to continue simulation
    bool getContinueState()
    {
        // Is there no (small) gap between the balls?
        // Or was there a timeout?
        return t < timeOut && (dist <= TWO_R + stopDist);
    }


    // Helper function for displaying velocities at the end
    void getDisplayVelocities(double (&a)[2][2][2])
    {
        // Get initial velocities
        for(int i=0; i < 2; i++)
            for(int j=0; j < 2; j++)
                a[0][i][j] = vi[i][j];

        // Get final velocities
        for(int i=0; i < 2; i++)
        {
            a[1][A][i] = b[A].v[i];
            a[1][B][i] = b[B].v[i];
        }
    }


public:
    // Run simulation
    void run()
    {
#ifdef DEBUG
        // Display table debug header for the ball to debug
        printf("b%c\tt\tsx\tsy\tvx\tvy\tax\tay\tFx\tFy", debugBall? 'B': 'A');
#ifdef DEBUG_POSN
        printf("\txAx\txAy\txCheck");
#endif // DEBUG_POSN
        printf("\n--------------------------------------------------------------------------------\n");
#endif // DEBUG

        // Decrease dt if it seems too high
        double minTSX = R * minTSPerc;
        dt = fmin(fmin(dt, minTSX/norm(b[A].v)), minTSX/norm(b[B].v));

        do
        {
            // Update the state for this time step
            update();
        }
        // Run until the balls are moving away from each other or there is a timeout
        while(getContinueState());

        #ifdef DEBUG
        printf("\n");
        #endif // DEBUG
    }


    // Display the variables
    void print()
    {
        char state[] = {'i', 'f'};
        double dispV[2][2][2];
        double vx, vy, normV, theta;

        // Display whether the balls collided or if the simulation timed out
        if(!hasEverCollided)
            printf("The balls did not collide.\n\n");
        if(timeOut <= t)
            printf("Simulation timed out.\n\n");

        // Display table of velocity data
        printf("\tx\ty\tnorm\ttheta\n");
        printf("--------------------------------------\n");
        getDisplayVelocities(dispV);
        for(int stateIndex = 0; stateIndex < 2; stateIndex++)
            for(int ball=0; ball < 2; ball++)
            {
                vx = dispV[stateIndex][ball][X];
                vy = dispV[stateIndex][ball][Y];
                normV = norm(dispV[stateIndex][ball]);
                theta = angleDeg(dispV[stateIndex][ball][Y], dispV[stateIndex][ball][X]);
                printf("v%d%c\t%.4f\t%.4f\t%.4f\t%.2f\n", ball+1, state[stateIndex],
                       vx, vy, normV, theta);
            }

#ifdef DEBUG
        // Calculate momentum and energy errors
        double dpX = b[A].m * vi[A][X] + b[B].m * vi[B][X] - b[A].m * b[A].v[X] - b[B].m * b[B].v[X];
        double dpY = b[A].m * vi[A][Y] + b[B].m * vi[B][Y] - b[A].m * b[A].v[Y] - b[B].m * b[B].v[Y];
        double dT = .5 * (b[A].m*norm(vi[A])*norm(vi[A]) + b[B].m*norm(vi[B])*norm(vi[B])
                          - b[A].m*norm(b[A].v)*norm(b[A].v) - b[B].m*norm(b[B].v)*norm(b[B].v));

        // Display momentum and energy errors
        printf("Difference in px: %.7f\n", dpX);
        printf("Difference in py: %.7f\n", dpY);
        printf("Difference in T: %.7f\n", dT);
#endif // DEBUG
    }


    // Default constructor
    BSystem()
    {
        // Angle of first ball relative to second ball
        // Must be in (-90, 90)in degrees
        double theta;
        double m[2], v[2][2], s[2][2];

#ifdef DEBUG
        printf("Debug which ball? ");
        scanf("%d", &debugBall);
#endif // DEBUG

        // Prompt for variables
        cout<<"Mass 1 in kg:";
        cin>>m[A];
        cout<<"Mass 2 in kg:";
        cin>>m[B];

        cout<<"Angle of ball 2 relative to ball 1 in degrees:"<<endl<<endl;
        cin>>theta;
        cout<<endl;

        cout<<"Components of initial velocity of ball 1 in ms^-1:"<<endl<<endl;
        for (int m=0; m<2; m++)
        {
            cin>>v[A][m];
        }
        cout<<endl;

        cout<<"Components of initial velocity of ball 2 in ms^-1:"<<endl<<endl;
        for (int m=0; m<2; m++)
        {
            cin>>v[B][m];
        }
        cout<<endl;

        // Calculate positions
        // Set ball 1 at the origin
        // Set ball 2 where the balls are just touching at the desired angle from ball 1
        theta *= PI / 180;
        s[A][X] = 0;
        s[B][X] = TWO_R*cos(theta);
        s[A][Y] = 0;
        s[B][Y] = TWO_R*sin(theta);

        // Remember initial velocity
        for(int i=0; i < 2; i++)
            for(int j=0; j < 2; j++)
                vi[i][j] = v[i][j];

        // Create balls
        for(int i=0; i < 2; i++)
            b[i] = Ball(m[i], v[i], s[i]);
    }
};


// Run the simulation
int main()
{
    // Create a system of balls
    BSystem colSys = BSystem();

    // Run the system
    colSys.run();

    // Dump variables
    colSys.print();
}
