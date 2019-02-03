//FrisbeeModel.cpp
#include <cstdlib>
#include <cmath>
#include <iomanip>
#import <iostream>
using namespace std;

int main() {
	//INITIAL CONDITION - Change variables to model specific flight
	double x = 0; 					//x position of the frisbee [m]
	double y = 1;					//y position of the frisbee [m]
	double v = 16;					//velocity magnitude [m/s]
	double dvx;						//change of velocity in x direction [m/s]
	double dvy;						//change of velocity in y direction [m/s]
	double dt = 0.001;				//time increment [s]
	double t = 0;					//total time elapsed [s]

	double g = 9.81;				//gravity [m/s^2]
	double m = 0.175;				//mass of frisbee [kg]
	double r = 0.135;
	double area = M_PI*pow(r, 2);//area of frisbee [m^2]
	double rho = 1.23; 				//density of air in [kg/m^3]

	double fd;						//drag force
	double cd;						//drag coefficient
	double cd0 = 0.18; 			//drag coefficent at alpha0
	double cda = 0.69;			//drag coefficient for alpha

	double fl;						//lift force
	double cl;						//lift coefficient
	double cl0 = 0.33;			//lift coefficient at alpha0
	double cla = 1.91;			//lift coefficient for alpha

	double alpha0 = 15;				//initial angle between frisbee and velocity
	double alpha = alpha0;			//angle between frisbee and velocity
	double theta = 5;			   //angle between velocity and x-axis
	double beta = alpha + theta;	//angle between x-axis and frsibee

	double vx = v*cos(theta*M_PI / 180);//velocity of frisbee in x-direction [m/s]
	double vy = v*sin(theta*M_PI / 180);//velocity of frisbee in y-direction [m/s]

										//Display Header Titles (for future data display every 10 loops)
	cout << "t\tx\ty\tv\ttheta\n";

	//loop index
	int i = 0;

	//loop until frisbee hits the ground
	while (y>0) {

		//Display data every 100 iterations
		if (i % 100 == 0) {
			cout << fixed << setprecision(3) << t << "\t" << x << "\t" << y << "\t" << v << "\t" << theta << "\t" << alpha << "\n";
		}


		t += dt; //increment time

				 //Calculate drag and lift coefficients
		cd = cd0 + cda*pow((alpha - alpha0)*M_PI / 180, 2);
		cl = cl0 + cla*alpha*M_PI / 180;

		//Calaculate drag and lift forces
		fd = rho*area*pow(v, 2)*cd / 2;
		fl = rho*area*pow(v, 2)*cl / 2;

		//Calculate the change of velocity in x and y direction
		dvx = (abs(fl)*sin(-theta*M_PI / 180) - abs(fd)*cos(theta*M_PI / 180))*(dt / m);
		dvy = (abs(fl)*cos(theta*M_PI / 180) + abs(fd)*sin(-theta*M_PI / 180) - m*g)*(dt / m);

		//Calculate new velocity value in each direction
		vx += dvx;
		vy += dvy;

		//Calculate total velocity
		v = sqrt(pow(vx, 2) + pow(vy, 2));

		//Calculate the angles theta and alpha
		theta = atan(vy / vx) * 180 / M_PI;
		alpha = beta - theta;

		//Calculate the current x and y position
		x += vx*dt;
		y += vy*dt;

		i++;
	}
}