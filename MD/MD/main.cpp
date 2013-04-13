//  This program was written to simulate the intermolecular interactions of a lennard-jones fluid (namely argon initialised in a lattice in a vacuum)
//  Changing dt changes the time interval of the loop, changing N changes number of particles in lattice, 
//  'animate' outputs coordinates of all atoms. Designed to be imported into jmol (an animation software which plots 3d coordinates of particles at time intervals)
//  Created by Scott Adams and Anna Jordan - 2012


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <math.h>
using namespace std;


double dt = 0.01;
const int N = 125;
double r[N][3], v[N][3], a[N][3];   // distances, velocities and accelerations of particles
double f;


int main()
{
	ofstream animate("XYZ_.txt");
	
	int n = int(ceil(pow(N, 1.0/3)));    // calculates number of atoms in each direction
    double L = 10;
	double d = L / n;                  // lattice spacing
	int p = 0;                         // particles placed so far
    for (int x = 0; x < n; x++)         // for x,y,z dimensions 
		for (int y = 0; y < n; y++) 
			for (int z = 0; z < n; z++) 
			{            
				if (p <= N)       
				{
					r[p][0] = (x) * d;  // updates position for x direction
					r[p][1] = (y) * d;  // updates position for y direction 
					r[p][2] = (z) * d;  // updates position for z direction
				}
				++p; 
			}
    
	// initialise velocity and acceleration of particles to zero
	for (int p = 0; p < N; p++)
		for (int k = 0; k < 3; k++)
		{
			v[p][k] = 0;
			a[p][k] = 0;
		}
	
    
	for (double time = 0; time <= 10; time += dt)
	{
		for (int i = 0; i < N; i++) 
		{
			for (int k = 0; k < 3; k++)
			{
				v[i][k] += 0.5 * a[i][k] * dt;  // updates velocity
				r[i][k] += v[i][k] * dt;    // updates position
			}
		}
		
		for (int p = 0; p < N; p++)
			for (int k = 0; k < 3; k++)
			{
				a[p][k] = 0;        // initialise accelerations to zero
			}
		
		
		for (int i = 0; i < N-1; i++) 
		{
			for (int j = i+1; j < N; j++) 
			{
				double rij[3];      // initialise array for relative distances in each direction
				double rSqd = 0; 
				for (int k = 0; k < 3; k++) 
				{
					rij[k] = r[i][k] - r[j][k];    // distance between two particles
					rSqd += rij[k] * rij[k];       // square distance
				} 
				
				double f = 24 * (2 * pow(rSqd, -7) - pow(rSqd, -4));  // update forces
				
				for (int k = 0; k < 3; k++)
				{
					a[i][k] += rij[k] * f;  
					a[j][k] -= rij[k] * f;  // updates accelerations based on forces
				}
			}
		}
		
		for (int i = 0; i < N; i++)
			for (int k = 0; k < 3; k++)
			{
				v[i][k] += 0.5 * a[i][k] * dt;     // updates velocities
			}
		
		// outputs all particle coordinates in jmol format
		animate << N << endl << endl ;
		for (int i = 0; i < N; i++)
        {
            animate << "particle \t" << r[i][0] << "\t" << r[i][1] << "\t" << r[i][2] << endl;
        }
		
	}
	return 0; 
}
