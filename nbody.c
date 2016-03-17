/*
 ============================================================================
 Name        : nbody.c
 Author      : Eric Stone
 Version     : 3
 Copyright   :
 Description : n-body gravity simulation

A n-body simulation modeling the (simplified) affect  the interacting 
gravitational pull of the object on each other in 2-D space.

The program times the difference, for varying numbers of objects, between the 
simulation using the results of Newtons third law to halve the number of 
calculations and the time needed without that optimization.


See: https://en.wikipedia.org/wiki/N-body_problem
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <time.h>

//program config settings
const int NUM_OBJECTS = 5000;
int thirdLaw = 1;
int printProgress = 0;
int printEveryTimeSteps = 10;

//simulation config
const unsigned long TIME_STEPS = 4;
const double TIME_STEP_SIZE = 0.001;
const double GRAV_CONST = 6.673E-11;

//stores the mass and position in space of each object
enum object {
	MASS, POS_X, POS_Y
};
//stores each objects current velocity
enum objectV {
	VEL_X, VEL_Y
};
enum twoSpace {
	X, Y
};

//default config is to position all the objects along a line
//the default makes the simulation reproducable
void makeRandomObjects(int totalObjects, double objects[][3], double objectsV[][2]) {

	const int MAX_X_DIM = 1.0;

	int i = 0;
	for (i = 0; i < totalObjects; i++) {

		objects[i][MASS] = 2.0;
		objects[i][POS_X] = i * (double) MAX_X_DIM / (double) totalObjects;
		objects[i][POS_Y] = 0;
		objectsV[i][VEL_X] = 0.0;
		objectsV[i][VEL_Y] = 0.0;
	}
}

//prints the status of all the objects to console
void printObjects(int timeStep, double objects[][3], double objectsV[][2], int totalObjects, int printEveryTimeSteps) {
	if(timeStep % printEveryTimeSteps == 0) {
		int i = 0;
		printf("\nOBJ_ID,    MASS,     POS_X,    POS_Y,    VEL_X,    VEL_Y,    TIME  \n");

		for (i = 0; i < totalObjects; i++) {

			printf("%d,       %f, %f, %f, %f, %f, ", i, objects[i][MASS],
					objects[i][POS_X], objects[i][POS_Y], objectsV[i][VEL_X],
					objectsV[i][VEL_Y]);
			printf("%d     \n", timeStep);
		}
	}
}

void runSimulation(int totalObjects, int tL) {

	thirdLaw = tL;

	long time = 0;
	int obj = 0;
	int otherObj = 0;

	double postition_X;
	double postition_Y;

	double r_length;

	double objects[NUM_OBJECTS][3];
	double objectsV[NUM_OBJECTS][2];
	double accelerationVector[NUM_OBJECTS][2];
	// the static double with hard-coded max values is needed
	// otherwise we overflow the stack on creation and seg fault
	static double forceVector[5000][5000][2];
	double velocityVector[NUM_OBJECTS][2];
	double forceTotal[NUM_OBJECTS][2];

	makeRandomObjects(NUM_OBJECTS, objects,objectsV);

	// setup to time the run
	clock_t begin, end;
	double time_spent;
	begin = clock();

	for (time = 0; time < TIME_STEPS; time++) {

		if(printProgress)
			printObjects(time, objects, objectsV, totalObjects, printEveryTimeSteps);

		for (obj = 0; obj < totalObjects; obj++) {
			forceTotal[obj][X] = 0.0;
			forceTotal[obj][Y] = 0.0;
		}



		for (obj = 0; obj < totalObjects; obj++) {

			if (thirdLaw) {
				for (otherObj = obj + 1; otherObj < totalObjects; otherObj++) {
					if (obj != otherObj) {
						postition_X = objects[otherObj][POS_X] - objects[obj][POS_X];
						postition_Y = objects[otherObj][POS_Y] - objects[obj][POS_Y];

						r_length = pow(postition_X, 2) + pow(postition_Y, 2);
						r_length = sqrt(r_length);

						forceVector[obj][otherObj][X] = GRAV_CONST * 1 * objects[otherObj][MASS] * postition_X / pow(r_length, 3);
						forceVector[obj][otherObj][Y] = GRAV_CONST * 1 * objects[otherObj][MASS] * postition_Y / pow(r_length, 3);

						forceTotal[obj][X] = forceTotal[obj][X] + forceVector[obj][otherObj][X];
						forceTotal[obj][Y] = forceTotal[obj][Y] + forceVector[obj][otherObj][Y];

						//implements Newton's third law symmetry
						forceTotal[otherObj][X] = forceTotal[otherObj][X] - forceVector[obj][otherObj][X];
						forceTotal[otherObj][Y] = forceTotal[otherObj][Y] - forceVector[obj][otherObj][Y];
					}
				}
			}

			if (!thirdLaw) {
				for (otherObj = 0; otherObj < totalObjects; otherObj++) {
					if (obj != otherObj) {
						postition_X = objects[otherObj][POS_X] - objects[obj][POS_X];
						postition_Y = objects[otherObj][POS_Y] - objects[obj][POS_Y];

						r_length = pow(postition_X, 2) + pow(postition_Y, 2);
						r_length = sqrt(r_length);

						forceVector[obj][otherObj][X] = GRAV_CONST * 1 * objects[otherObj][MASS] * postition_X / pow(r_length, 3);
						forceVector[obj][otherObj][Y] = GRAV_CONST * 1 * objects[otherObj][MASS] * postition_Y / pow(r_length, 3);

						forceTotal[obj][X] = forceTotal[obj][X] + forceVector[obj][otherObj][X];
						forceTotal[obj][Y] = forceTotal[obj][Y] + forceVector[obj][otherObj][Y];
					}
				}
			}


		}

		for (obj = 0; obj < totalObjects; obj++) {

			accelerationVector[obj][X] = forceTotal[obj][X];
			accelerationVector[obj][Y] = forceTotal[obj][Y];

			velocityVector[obj][X] = objectsV[obj][VEL_X];
			velocityVector[obj][Y] = objectsV[obj][VEL_Y];

			objectsV[obj][VEL_X] = objectsV[obj][VEL_X] + TIME_STEP_SIZE * accelerationVector[obj][X];
			objectsV[obj][VEL_Y] = objectsV[obj][VEL_Y] + TIME_STEP_SIZE * accelerationVector[obj][Y];

			objectsV[obj][POS_X] = objectsV[obj][POS_X] + TIME_STEP_SIZE * velocityVector[obj][X];
			objectsV[obj][POS_Y] = objectsV[obj][POS_Y] + TIME_STEP_SIZE * velocityVector[obj][Y];
		}
	}

	// timing code
	end = clock();
	time_spent = (double) (end - begin) / CLOCKS_PER_SEC;

	if(printProgress)
		printObjects(time + 1, objects, objectsV, totalObjects, printEveryTimeSteps);

	printf("  %d,   %d,  ", totalObjects, thirdLaw);
	printf(" %f \n", time_spent);
}

int main(void) {
	printf("NUM_OBJECTS, thirdLaw, time \n");


	runSimulation(10, 0);
	runSimulation(10, 0);
	runSimulation(10, 0);
	runSimulation(10, 0);
	printf("\n");
	runSimulation(10, 1);
	runSimulation(10, 1);
	runSimulation(10, 1);
	runSimulation(10, 1);
	printf("\n");
	runSimulation(100, 0);
	runSimulation(100, 0);
	runSimulation(100, 0);
	runSimulation(100, 0);
	printf("\n");
	runSimulation(100, 1);
	runSimulation(100, 1);
	runSimulation(100, 1);
	runSimulation(100, 1);
	printf("\n");
	runSimulation(1000, 0);
	runSimulation(1000, 0);
	runSimulation(1000, 0);
	runSimulation(1000, 0);
	printf("\n");
	runSimulation(1000, 1);
	runSimulation(1000, 1);
	runSimulation(1000, 1);
	runSimulation(1000, 1);
	printf("\n");
	runSimulation(2000, 0);
	runSimulation(2000, 0);
	runSimulation(2000, 0);
	runSimulation(2000, 0);
	printf("\n");
	runSimulation(2000, 1);
	runSimulation(2000, 1);
	runSimulation(2000, 1);
	runSimulation(2000, 1);
	printf("\n");
	runSimulation(3000, 0);
	runSimulation(3000, 0);
	runSimulation(3000, 0);
	runSimulation(3000, 0);
	printf("\n");
	runSimulation(3000, 1);
	runSimulation(3000, 1);
	runSimulation(3000, 1);
	runSimulation(3000, 1);
	printf("\n");
	runSimulation(4000, 0);
	runSimulation(4000, 0);
	runSimulation(4000, 0);
	runSimulation(4000, 0);
	printf("\n");
	runSimulation(4000, 1);
	runSimulation(4000, 1);
	runSimulation(4000, 1);
	runSimulation(4000, 1);
	printf("\n");
	runSimulation(5000, 0);
	runSimulation(5000, 0);
	runSimulation(5000, 0);
	runSimulation(5000, 0);
	printf("\n");
	runSimulation(5000, 1);
	runSimulation(5000, 1);
	runSimulation(5000, 1);
	runSimulation(5000, 1);
	printf("\n");
	return EXIT_SUCCESS;
}
