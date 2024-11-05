//Riley Wind
//Computational physics midterm
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//choose units such that g=1
//pendulum such that m=1, l=1
//therefore no need for parameter for any constants
//legrangian is time-independent so no need for parameter for time
//evolution in this case refers to time derivative not path integral
//o is angle, w is angular velocity

double o1evolution(double o1, double o2, double w1, double w2){
double evo=6*((2*w1-3*w2*cos(o1-o2))/(16-9*pow(cos(o1-o2), 2));
return evo;
}

double o2evolution(double o1, double o2, double w1, double w2){
double evo=6*((8*w2-3*w1*cos(o1-o2))/(16-9*pow(cos(o1-o2), 2));
return evo;
}


double w1evolution(double o1, double o2, double w1, double w2){
double evo=(-.5)*(w1*w2*sin(o1-o2)+3*sin(o1));
return evo;
}
double w2evolution(double o1, double o2, double w1, double w2){
double evo=(-.5)*((-1)*w1*w2*sin(o1-o2)+*sin(o2));
return evo;
}

double step(double dt, double f, double *x){
double newval=*x+f*dt;
return newval;
}

void rkO2(double dt, int steps, double *t, double *o1, double *o2, double *w1, double *w2, 
			double (*do1dt)(double , double , double , double), double (*do2dt)(double , double , double , double),
			double (*dw1dt)(double , double , double , double), double (*dw2dt)(double , double , double , double)){
double o1temp, o2temp, w1temp, w3temp, fo1, fo2, fw1, fw2, fo1temp, fo2temp, fw1temp, fw2temp;
int i;
//calculate what the values would be after half a time step, then calculate the derivatives based on these values, then step
//from the original values using the new derivative
for(i=0; i<steps-1; i++){
fo1temp=o1evolution(*(o1+i), *(o2+i), *(w1+i), *(w2+i));
fo2temp=o2evolution(*(o1+i), *(o2+i), *(w1+i), *(w2+i));
fw1temp=w1evolution(*(o1+i), *(o2+i), *(w1+i), *(w2+i));
fw2temp=w2evolution(*(o1+i), *(o2+i), *(w1+i), *(w2+i));
o1temp=step(0.5*dt, fo1temp, *(o1+i));
o2temp=step(0.5*dt, fo2temp, *(o2+i));
w1temp=step(0.5*dt, fw1temp, *(w1+i));
w2temp=step(0.5*dt, fw2temp, *(w2+i));
fo1=o1evolution(o1temp, o2temp, w1temp, w2temp);
fo2=o2evolution(o1temp, o2temp, w1temp, w2temp);
fw1=w1evolution(o1temp, o2temp, w1temp, w2temp);
fw2=w2evolution(o1temp, o2temp, w1temp, w2temp);
*(o1+i+1)=step(dt, fo1temp, *(o1+i));
*(o2+i+1)=step(dt, fo2temp, *(o2+i));
*(w1+i+1)=step(dt, fw1temp, *(w1+i));
*(w2+i+1)=step(dt, fw2temp, *(w2+i));
*(t+i+1)=*(t+i)+dt;
}
}

main(){
//time interval of .05 seconds, 50 seconds
double *o1=(double *)calloc(1000, sizeof(double));
double *o2=(double *)calloc(1000, sizeof(double));
double *w1=(double *)calloc(1000, sizeof(double));
double *w2=(double *)calloc(1000, sizeof(double));
double *t=(double *)calloc(1000, sizeof(double));
//I think the equations in the literature only work in radians so I'm using radians
//otherwise, define a new cos and sin that converts argument to degrees before calculating
*o1=1.7;
*o2=1.7;
*w1=0;
*w2=0;
*t=0;
double dt=.05;
int steps=1000;
int j;
rkO2(dt, steps, t, o1, o2, w1, w2, o1evolution, o2evolution, w1evolution, w2evolution);
FILE *nf=fopen("DoublePendulumTrajectory.dat", "w");
for(j=0; j<steps; j++){
fprintf(nf, "t: %f\to1: %f\to2: %f\tw1: %f\tw2: %f\n", *(t+j), *(o1+j), *(o2+j), *(w1+j), *(w2+j));
}
}