#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <iostream>
#include "omp.h"
#include <cmath>
#include <iomanip>
#include <cstdlib>
//#include <string>


using namespace std;

struct SphericalBody {
	double coordX;
	double coordY;
	double coordZ;

	double velocityX;
	double velocityY;
	double velocityZ;

	double forceX;
	double forceY;
	double forceZ;

	double radius;
	double mass;
};

struct OnlyVec
{
	double vecX;
	double vecY;
	double vecZ;
};

double EPSILON = 0.1;
int NUMBER_OF_ITERATIONS = 10000;
int NUMBER_OF_BODIES;
int MAX_BUFFER;

const double MAX_FORCE = 1000000000;
const double DELTA_TIME = 0.1;
const double GRAVITATION = 0.0067;
int numberOfThreads;

OnlyVec** forceSupplements;
OnlyVec** velocitySupplements;
OnlyVec* velocitySupplementsForSingle;

SphericalBody* alterBodiesArray;
SphericalBody* bodiesArray;
bool* isDivided;

bool* probabilityOfDivide;
int offset;



bool getRandomNumber(double probability)
{
	static const double fraction = 1.0 / ((double)(RAND_MAX)+1.0);
	return (((double)rand() * fraction) < probability) ? true : false;
}

double getRandomNumber(double min, double max)
{
	static const double fraction = 1.0 / ((double)(RAND_MAX)+1.0);
	return rand() * fraction * (max - min + 1.0) + min;
}

void generatePlanets(SphericalBody* generatingArray, char* fileName) {
	FILE* planetsArray;
	int bodiesNumber;
	double coordMin, coordMax;
	double velosityMin, velosityMax;
	double massMin, massMax;
	double radiusMin, radiusMax;
	double distanceX, distanceY, distanceZ;
	double rootOfInverseSquareDistance;

	cout << "how much bodies do yo want?" << endl;
	cin >> bodiesNumber;
	if (bodiesNumber < 1) {
		bodiesNumber = 1;
	}
	generatingArray = new SphericalBody[bodiesNumber * 4];
	NUMBER_OF_BODIES = bodiesNumber;
	MAX_BUFFER = NUMBER_OF_BODIES * 4;
	cout << "choose min and max values for coordinates" << endl;
	cin >> coordMin >> coordMax;
	cout << "choose min and max values for velosity" << endl;
	cin >> velosityMin >> velosityMax;
	cout << "choose min and max values for mass" << endl;
	cin >> massMin >> massMax;
	cout << "choose min and max values for radius" << endl;
	cin >> radiusMin >> radiusMax;

	generatingArray[0].coordX = getRandomNumber(coordMin, coordMax);
	generatingArray[0].coordY = getRandomNumber(coordMin, coordMax);
	generatingArray[0].coordZ = getRandomNumber(coordMin, coordMax);

	generatingArray[0].velocityX = getRandomNumber(velosityMin, velosityMax);
	generatingArray[0].velocityY = getRandomNumber(velosityMin, velosityMax);
	generatingArray[0].velocityZ = getRandomNumber(velosityMin, velosityMax);

	generatingArray[0].forceX = 0;
	generatingArray[0].forceY = 0;
	generatingArray[0].forceZ = 0;

	generatingArray[0].mass = getRandomNumber(massMin, massMax);
	generatingArray[0].radius = getRandomNumber(radiusMin, radiusMax);

	for (int i = 1; i < NUMBER_OF_BODIES; i++) {
		generatingArray[i].coordX = getRandomNumber(coordMin, coordMax);
		generatingArray[i].coordY = getRandomNumber(coordMin, coordMax);
		generatingArray[i].coordZ = getRandomNumber(coordMin, coordMax);

		generatingArray[i].velocityX = getRandomNumber(velosityMin, velosityMax);
		generatingArray[i].velocityY = getRandomNumber(velosityMin, velosityMax);
		generatingArray[i].velocityZ = getRandomNumber(velosityMin, velosityMax);

		generatingArray[i].forceX = 0;
		generatingArray[i].forceY = 0;
		generatingArray[i].forceZ = 0;

		generatingArray[i].mass = getRandomNumber(massMin, massMax);
		generatingArray[i].radius = getRandomNumber(radiusMin, radiusMax);

		for (int j = i - 1; j >= 0; j--) {
			distanceX = generatingArray[j].coordX - generatingArray[i].coordX;
			distanceY = generatingArray[j].coordY - generatingArray[i].coordY;
			distanceZ = generatingArray[j].coordZ - generatingArray[i].coordZ;

			rootOfInverseSquareDistance = sqrt(1 / (distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ));
			if ((1 / rootOfInverseSquareDistance) < (generatingArray[i].radius + generatingArray[j].radius)) {
				i--;
				NUMBER_OF_BODIES--;
				break;
			}
		}

	}

	getchar();
	cout << "Number of bodies: " << NUMBER_OF_BODIES << endl;
	cout << "Name the file:" << endl;
	scanf("%[^\n]", fileName);
	planetsArray = fopen(fileName, "wb");
	fwrite(&NUMBER_OF_BODIES, sizeof(int), 1, planetsArray);
	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		fwrite(&generatingArray[i], sizeof(SphericalBody), 1, planetsArray);
	}
	fclose(planetsArray);

	delete[] generatingArray;


}

/*void preprocess() {
	SphericalBody temp;
	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		temp.coordX = (i % 8) * 150000;
		temp.coordY = ((i / 8) % 8) * 150000;
		temp.coordZ = ((i / 64) % 8) * 150000;
		temp.forceX = 0;
		temp.forceY = 0;
		temp.forceZ = 0;
		temp.velocityX = 0;
		temp.velocityY = 0;
		temp.velocityZ = 0;
		temp.mass = 100000000000 + (double)(i % 100) * 1000000000;
		temp.radius = 10000 + (double)(i % 100) * 100;
		bodiesArray[i] = temp;
	}
}*/

SphericalBody revectorize(int indexOfDivide) {
	double oldX;
	double oldY;
	double oldZ;

	oldX = bodiesArray[indexOfDivide].velocityX;
	oldY = bodiesArray[indexOfDivide].velocityY;
	oldZ = bodiesArray[indexOfDivide].velocityZ;

	bodiesArray[indexOfDivide].mass /= 2;

	bodiesArray[indexOfDivide].velocityX = (1.0 - sqrt(2.0)) / 2.0 * oldX + (1.0 + sqrt(2.0)) / 2.0 * oldY + oldZ / sqrt(2.0);
	bodiesArray[indexOfDivide].velocityY = (1.0 + sqrt(2.0)) / 2.0 * oldX + (1.0 - sqrt(2.0)) / 2.0 * oldY + oldZ / sqrt(2.0);
	bodiesArray[indexOfDivide].velocityZ = oldX / sqrt(2.0) + oldY / sqrt(2.0) - oldZ;

	double magnificationFirst;

	magnificationFirst = bodiesArray[indexOfDivide].radius / 2 /
		sqrt(bodiesArray[indexOfDivide].velocityX * bodiesArray[indexOfDivide].velocityX +
			bodiesArray[indexOfDivide].velocityY * bodiesArray[indexOfDivide].velocityY +
			bodiesArray[indexOfDivide].velocityZ * bodiesArray[indexOfDivide].velocityZ);

	SphericalBody newBody;

	newBody.mass = bodiesArray[indexOfDivide].mass;

	newBody.velocityX = oldX * 2 - bodiesArray[indexOfDivide].velocityX;
	newBody.velocityY = oldY * 2 - bodiesArray[indexOfDivide].velocityY;
	newBody.velocityZ = oldZ * 2 - bodiesArray[indexOfDivide].velocityZ;

	newBody.coordX = bodiesArray[indexOfDivide].coordX;
	newBody.coordY = bodiesArray[indexOfDivide].coordY;
	newBody.coordZ = bodiesArray[indexOfDivide].coordZ;

	bodiesArray[indexOfDivide].coordX += magnificationFirst * bodiesArray[indexOfDivide].velocityX;
	bodiesArray[indexOfDivide].coordY += magnificationFirst * bodiesArray[indexOfDivide].velocityY;
	bodiesArray[indexOfDivide].coordZ += magnificationFirst * bodiesArray[indexOfDivide].velocityZ;

	newBody.coordX -= magnificationFirst * bodiesArray[indexOfDivide].velocityX;
	newBody.coordY -= magnificationFirst * bodiesArray[indexOfDivide].velocityY;
	newBody.coordZ -= magnificationFirst * bodiesArray[indexOfDivide].velocityZ;

	bodiesArray[indexOfDivide].radius /= 2.1;
	newBody.radius = bodiesArray[indexOfDivide].radius;

	newBody.forceX = 0;
	newBody.forceY = 0;
	newBody.forceZ = 0;

	return newBody;
}

void collisionSolve(int first, int second, double distanceX, double distanceY, double distanceZ, double inverseSquareRoot) {

	double scalar;
	double massForFirst;
	double massForSecond;

	scalar = (bodiesArray[second].velocityX - bodiesArray[first].velocityX) * distanceX +
		(bodiesArray[second].velocityY - bodiesArray[first].velocityY) * distanceY +
		(bodiesArray[second].velocityZ - bodiesArray[first].velocityZ) * distanceZ;

	massForFirst = 2 * bodiesArray[second].mass / (bodiesArray[first].mass + bodiesArray[second].mass);
	massForSecond = 2 * bodiesArray[first].mass / (bodiesArray[first].mass + bodiesArray[second].mass);

	velocitySupplementsForSingle[first].vecX += massForFirst * scalar * inverseSquareRoot * (-distanceX);
	velocitySupplementsForSingle[first].vecY += massForFirst * scalar * inverseSquareRoot * (-distanceY);
	velocitySupplementsForSingle[first].vecZ += massForFirst * scalar * inverseSquareRoot * (-distanceZ);

	velocitySupplementsForSingle[second].vecX += massForSecond * scalar * inverseSquareRoot * distanceX;
	velocitySupplementsForSingle[second].vecY += massForSecond * scalar * inverseSquareRoot * distanceY;
	velocitySupplementsForSingle[second].vecZ += massForSecond * scalar * inverseSquareRoot * distanceZ;
}

void collisionSolveMultiple(int first, int second, double distanceX, double distanceY, double distanceZ, double inverseSquareRoot, int threadNumber) {

	double scalar;
	double massForFirst;
	double massForSecond;

	scalar = (bodiesArray[second].velocityX - bodiesArray[first].velocityX) * distanceX +
		(bodiesArray[second].velocityY - bodiesArray[first].velocityY) * distanceY +
		(bodiesArray[second].velocityZ - bodiesArray[first].velocityZ) * distanceZ;

	massForFirst = 2 * bodiesArray[second].mass / (bodiesArray[first].mass + bodiesArray[second].mass);
	massForSecond = 2 * bodiesArray[first].mass / (bodiesArray[first].mass + bodiesArray[second].mass);

	velocitySupplements[first][threadNumber].vecX -= massForFirst * scalar * inverseSquareRoot * (-distanceX);
	velocitySupplements[first][threadNumber].vecY -= massForFirst * scalar * inverseSquareRoot * (-distanceY);
	velocitySupplements[first][threadNumber].vecZ -= massForFirst * scalar * inverseSquareRoot * (-distanceZ);

	velocitySupplements[second][threadNumber].vecX -= massForSecond * scalar * inverseSquareRoot * distanceX;
	velocitySupplements[second][threadNumber].vecY -= massForSecond * scalar * inverseSquareRoot * distanceY;
	velocitySupplements[second][threadNumber].vecZ -= massForSecond * scalar * inverseSquareRoot * distanceZ;
}

void forceCalculation() {
	double distanceX;
	double distanceY;
	double distanceZ;

	double inverseSquareDistance;
	double rootOfInverseSquareDistance;

	double gravityForce;

	double forceX;
	double forceY;
	double forceZ;

	for (int i = 0; i < NUMBER_OF_BODIES - 1; i++) {
		for (int j = i + 1; j < NUMBER_OF_BODIES; j++) {
			distanceX = bodiesArray[j].coordX - bodiesArray[i].coordX;
			distanceY = bodiesArray[j].coordY - bodiesArray[i].coordY;
			distanceZ = bodiesArray[j].coordZ - bodiesArray[i].coordZ;

			inverseSquareDistance = 1 / (distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ);
			rootOfInverseSquareDistance = sqrt(inverseSquareDistance);

			if ((1 / rootOfInverseSquareDistance) < (bodiesArray[i].radius + bodiesArray[j].radius)) {
				collisionSolve(i, j, distanceX, distanceY, distanceZ, inverseSquareDistance);
				isDivided[i] = true;
			}

		}
	}


	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		bodiesArray[i].velocityX -= velocitySupplementsForSingle[i].vecX;
		bodiesArray[i].velocityY -= velocitySupplementsForSingle[i].vecY;
		bodiesArray[i].velocityZ -= velocitySupplementsForSingle[i].vecZ;

		velocitySupplementsForSingle[i].vecX = 0;
		velocitySupplementsForSingle[i].vecY = 0;
		velocitySupplementsForSingle[i].vecZ = 0;

	}

	for (int i = 0; i < NUMBER_OF_BODIES - 1; i++) {
		int threadNumber = omp_get_thread_num();
		for (int j = i + 1; j < NUMBER_OF_BODIES; j++) {

			distanceX = bodiesArray[j].coordX - bodiesArray[i].coordX;
			distanceY = bodiesArray[j].coordY - bodiesArray[i].coordY;
			distanceZ = bodiesArray[j].coordZ - bodiesArray[i].coordZ;

			inverseSquareDistance = 1 / (distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ);
			rootOfInverseSquareDistance = sqrt(inverseSquareDistance);

			gravityForce = GRAVITATION * bodiesArray[i].mass * bodiesArray[j].mass * inverseSquareDistance;
			if (gravityForce > MAX_FORCE) {
				gravityForce = MAX_FORCE;
			}

			forceX = gravityForce * distanceX * rootOfInverseSquareDistance;
			forceY = gravityForce * distanceY * rootOfInverseSquareDistance;
			forceZ = gravityForce * distanceZ * rootOfInverseSquareDistance;

			bodiesArray[i].forceX += forceX;
			bodiesArray[i].forceY += forceY;
			bodiesArray[i].forceZ += forceZ;

			bodiesArray[j].forceX -= forceX;
			bodiesArray[j].forceY -= forceY;
			bodiesArray[j].forceZ -= forceZ;

		}
	}

}

void forceCalculationMultiple() {
	double distanceX;
	double distanceY;
	double distanceZ;

	double inverseSquareDistance;
	double rootOfInverseSquareDistance;

	double gravityForce;

	double forceX;
	double forceY;
	double forceZ;




#pragma omp parallel for num_threads(numberOfThreads) shared(bodiesArray, velocitySupplements, isDivided) private(distanceX, distanceY, distanceZ, inverseSquareDistance, rootOfInverseSquareDistance)
	for (int i = 0; i < NUMBER_OF_BODIES - 1; i++) {
		int threadNumber = omp_get_thread_num();
		for (int j = i + 1; j < NUMBER_OF_BODIES; j++) {
			distanceX = bodiesArray[j].coordX - bodiesArray[i].coordX;
			distanceY = bodiesArray[j].coordY - bodiesArray[i].coordY;
			distanceZ = bodiesArray[j].coordZ - bodiesArray[i].coordZ;

			inverseSquareDistance = 1 / (distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ);
			rootOfInverseSquareDistance = sqrt(inverseSquareDistance);

			if ((1 / rootOfInverseSquareDistance) < (bodiesArray[i].radius + bodiesArray[j].radius)) {
				collisionSolveMultiple(i, j, distanceX, distanceY, distanceZ, inverseSquareDistance, threadNumber);
				isDivided[i] = true;
			}

		}
	}


#pragma omp parallel for num_threads(numberOfThreads) shared(bodiesArray, velocitySupplements)
	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		for (int j = 0; j < numberOfThreads; j++) {
			bodiesArray[i].velocityX += velocitySupplements[i][j].vecX;
			bodiesArray[i].velocityY += velocitySupplements[i][j].vecY;
			bodiesArray[i].velocityZ += velocitySupplements[i][j].vecZ;

			velocitySupplements[i][j].vecX = 0;
			velocitySupplements[i][j].vecY = 0;
			velocitySupplements[i][j].vecZ = 0;
		}
	}

#pragma omp parallel for num_threads(numberOfThreads) shared(bodiesArray, forceSupplements) private(distanceX, distanceY, distanceZ, inverseSquareDistance, rootOfInverseSquareDistance, gravityForce, forceX, forceY, forceZ)
	for (int i = 0; i < NUMBER_OF_BODIES - 1; i++) {
		int threadNumber = omp_get_thread_num();
		//threadNumber = 0;
		for (int j = i + 1; j < NUMBER_OF_BODIES; j++) {

			distanceX = bodiesArray[j].coordX - bodiesArray[i].coordX;
			distanceY = bodiesArray[j].coordY - bodiesArray[i].coordY;
			distanceZ = bodiesArray[j].coordZ - bodiesArray[i].coordZ;

			inverseSquareDistance = 1 / (distanceX * distanceX + distanceY * distanceY + distanceZ * distanceZ);
			rootOfInverseSquareDistance = sqrt(inverseSquareDistance);

			gravityForce = GRAVITATION * bodiesArray[i].mass * bodiesArray[j].mass * inverseSquareDistance;
			if (gravityForce > MAX_FORCE) {
				gravityForce = MAX_FORCE;
			}

			forceX = gravityForce * distanceX * rootOfInverseSquareDistance;
			forceY = gravityForce * distanceY * rootOfInverseSquareDistance;
			forceZ = gravityForce * distanceZ * rootOfInverseSquareDistance;

			forceSupplements[i][threadNumber].vecX += forceX;
			forceSupplements[i][threadNumber].vecY += forceY;
			forceSupplements[i][threadNumber].vecZ += forceZ;

			forceSupplements[j][threadNumber].vecX -= forceX;
			forceSupplements[j][threadNumber].vecY -= forceY;
			forceSupplements[j][threadNumber].vecZ -= forceZ;
		}
	}

#pragma omp parallel for num_threads(numberOfThreads) shared(bodiesArray, forceSupplements)
	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		for (int j = 0; j < numberOfThreads; j++) {
			bodiesArray[i].forceX += forceSupplements[i][j].vecX;
			bodiesArray[i].forceY += forceSupplements[i][j].vecY;
			bodiesArray[i].forceZ += forceSupplements[i][j].vecZ;

			forceSupplements[i][j].vecX = 0;
			forceSupplements[i][j].vecY = 0;
			forceSupplements[i][j].vecZ = 0;
		}
	}


}

void moveBodiesAndNullifyForces() {
	int realCountOfSupplements = 0;
	int countOfSupplements = 0;
	int index = 0;

	double accelerationPerTimeX;
	double accelerationPerTimeY;
	double accelerationPerTimeZ;

	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		accelerationPerTimeX = bodiesArray[i].forceX * DELTA_TIME / bodiesArray[i].mass;
		accelerationPerTimeY = bodiesArray[i].forceY * DELTA_TIME / bodiesArray[i].mass;
		accelerationPerTimeZ = bodiesArray[i].forceZ * DELTA_TIME / bodiesArray[i].mass;

		bodiesArray[i].coordX += bodiesArray[i].velocityX * DELTA_TIME + accelerationPerTimeX * DELTA_TIME / 2;
		bodiesArray[i].coordY += bodiesArray[i].velocityY * DELTA_TIME + accelerationPerTimeY * DELTA_TIME / 2;
		bodiesArray[i].coordZ += bodiesArray[i].velocityZ * DELTA_TIME + accelerationPerTimeZ * DELTA_TIME / 2;

		bodiesArray[i].velocityX += accelerationPerTimeX;
		bodiesArray[i].velocityY += accelerationPerTimeY;
		bodiesArray[i].velocityZ += accelerationPerTimeZ;

		bodiesArray[i].forceX = 0;
		bodiesArray[i].forceY = 0;
		bodiesArray[i].forceZ = 0;
	}

	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		if (isDivided[i]) {
			countOfSupplements++;
		}
	}

	if (((NUMBER_OF_BODIES + countOfSupplements) < MAX_BUFFER) && (countOfSupplements != 0)) {
		int* newElements = new int[countOfSupplements];

		for (int i = 0; i < NUMBER_OF_BODIES; i++) {
			if (isDivided[i]) {
				newElements[index] = i;
				isDivided[i] = false;
				index++;
			}
		}

		for (int i = 0; i < countOfSupplements; i++) {
			if (probabilityOfDivide[(offset + i) % MAX_BUFFER]) {
				bodiesArray[NUMBER_OF_BODIES + i] = revectorize(newElements[i]);
				realCountOfSupplements++;
			}
		}

		offset += realCountOfSupplements;

		NUMBER_OF_BODIES = NUMBER_OF_BODIES + realCountOfSupplements;
	}



}

void moveBodiesAndNullifyForcesMultiple() {
	int realCountOfSupplements = 0;
	int countOfSupplements = 0;
	int index = 0;

	double accelerationPerTimeX;
	double accelerationPerTimeY;
	double accelerationPerTimeZ;

#pragma omp parallel for shared(bodiesArray) private(accelerationPerTimeX, accelerationPerTimeY, accelerationPerTimeZ) 
	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		accelerationPerTimeX = bodiesArray[i].forceX * DELTA_TIME / bodiesArray[i].mass;
		accelerationPerTimeY = bodiesArray[i].forceY * DELTA_TIME / bodiesArray[i].mass;
		accelerationPerTimeZ = bodiesArray[i].forceZ * DELTA_TIME / bodiesArray[i].mass;

		bodiesArray[i].coordX += bodiesArray[i].velocityX * DELTA_TIME + accelerationPerTimeX * DELTA_TIME / 2;
		bodiesArray[i].coordY += bodiesArray[i].velocityY * DELTA_TIME + accelerationPerTimeY * DELTA_TIME / 2;
		bodiesArray[i].coordZ += bodiesArray[i].velocityZ * DELTA_TIME + accelerationPerTimeZ * DELTA_TIME / 2;

		bodiesArray[i].velocityX += accelerationPerTimeX;
		bodiesArray[i].velocityY += accelerationPerTimeY;
		bodiesArray[i].velocityZ += accelerationPerTimeZ;

		bodiesArray[i].forceX = 0;
		bodiesArray[i].forceY = 0;
		bodiesArray[i].forceZ = 0;
	}

#pragma omp parallel for reduction(+: countOfSupplements)
	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		if (isDivided[i]) {
			countOfSupplements++;
		}
	}

	if (((NUMBER_OF_BODIES + countOfSupplements) < MAX_BUFFER) && (countOfSupplements != 0)) {
		int* newElements = new int[countOfSupplements];

		for (int i = 0; i < NUMBER_OF_BODIES; i++) {
			if (isDivided[i]) {
				newElements[index] = i;
				isDivided[i] = false;
				index++;
			}
		}

		for (int i = 0; i < countOfSupplements; i++) {
			if (probabilityOfDivide[(offset + i) % MAX_BUFFER]) {
				bodiesArray[NUMBER_OF_BODIES + i] = revectorize(newElements[i]);
				realCountOfSupplements++;
			}
		}

		offset += realCountOfSupplements;

		NUMBER_OF_BODIES = NUMBER_OF_BODIES + realCountOfSupplements;
	}

}



void physicsCalculation() {
	forceCalculation();
	moveBodiesAndNullifyForces();
	//cout << bodiesArray[0].coordX;
}

void physicsCalculationMultiple() {
	forceCalculationMultiple();
	moveBodiesAndNullifyForcesMultiple();
	//cout << bodiesArray[0].coordX;
}

void subMenu() {

	int tempNUMBOFBOD = NUMBER_OF_BODIES;
	SphericalBody* tempBodiesArray = new SphericalBody[MAX_BUFFER];

	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		tempBodiesArray[i] = bodiesArray[i];
	}

	numberOfThreads = omp_get_max_threads();

	forceSupplements = new OnlyVec * [MAX_BUFFER];
	for (int i = 0; i < MAX_BUFFER; i++) {
		forceSupplements[i] = new OnlyVec[numberOfThreads];
	}

	for (int i = 0; i < MAX_BUFFER; i++) {
		for (int j = 0; j < numberOfThreads; j++) {
			forceSupplements[i][j].vecX = 0;
			forceSupplements[i][j].vecY = 0;
			forceSupplements[i][j].vecZ = 0;
		}
	}

	velocitySupplements = new OnlyVec * [MAX_BUFFER];
	for (int i = 0; i < MAX_BUFFER; i++) {
		velocitySupplements[i] = new OnlyVec[numberOfThreads];
	}

	for (int i = 0; i < MAX_BUFFER; i++) {
		for (int j = 0; j < numberOfThreads; j++) {
			velocitySupplements[i][j].vecX = 0;
			velocitySupplements[i][j].vecY = 0;
			velocitySupplements[i][j].vecZ = 0;
		}
	}

	velocitySupplementsForSingle = new OnlyVec[MAX_BUFFER];
	for (int i = 0; i < MAX_BUFFER; i++) {
		velocitySupplementsForSingle[i].vecX = 0;
		velocitySupplementsForSingle[i].vecY = 0;
		velocitySupplementsForSingle[i].vecZ = 0;
	}

	isDivided = new bool[MAX_BUFFER];

	for (int i = 0; i < MAX_BUFFER; i++) {
		isDivided[i] = false;
	}

	alterBodiesArray = new SphericalBody[MAX_BUFFER];

	double time = omp_get_wtime();

	for (int i = 0; i < NUMBER_OF_ITERATIONS; i++) {
		physicsCalculationMultiple();
	}

	time = omp_get_wtime() - time;

	cout << "time for multiple: " << time << endl;
	cout << NUMBER_OF_BODIES << endl;

	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		alterBodiesArray[i] = bodiesArray[i];
	}

	delete[] bodiesArray;

	offset = 0;

	NUMBER_OF_BODIES = tempNUMBOFBOD;

	bodiesArray = new SphericalBody[MAX_BUFFER];

	for (int i = 0; i < NUMBER_OF_BODIES; i++) {
		bodiesArray[i] = tempBodiesArray[i];
	}


	for (int i = 0; i < MAX_BUFFER; i++) {
		isDivided[i] = false;
	}

	time = omp_get_wtime();

	for (int i = 0; i < NUMBER_OF_ITERATIONS; i++) {
		physicsCalculation();
	}

	time = omp_get_wtime() - time;

	cout << "time for single: " << time << endl;
	cout << NUMBER_OF_BODIES << endl;
}


int main(int argc, char** argv) {

	setlocale(LC_ALL, "Russian");

	char* fileName = new char[256];

	FILE* planetsArray;

	int menuKey = -1;

	srand(time(0));

	cout << "Choose number of operation:" << endl;
	cout << "[1] - Choose file from disk" << endl;
	cout << "[2] - Generate data" << endl;
	cout << "[3] - Calculation from a single- and multi-threaded program" << endl;
	cout << "[4] - Compare results" << endl;
	cout << "[0] - Exit" << endl;

	cin >> menuKey;
	getchar();

	while (menuKey != 0) {
		switch (menuKey)
		{
		case 1:
		{
			cout << "Choose file:" << endl;
			scanf("%[^\n]", fileName);
			if ((planetsArray = fopen(fileName, "rb")) == NULL) {
				printf("He удается открыть файл.\n");
				exit(1);
			}
			fread(&NUMBER_OF_BODIES, sizeof(int), 1, planetsArray);
			MAX_BUFFER = NUMBER_OF_BODIES * 4;
			bodiesArray = new SphericalBody[MAX_BUFFER];
			for (int i = 0; i < NUMBER_OF_BODIES; i++) {
				fread(&bodiesArray[i], sizeof(SphericalBody), 1, planetsArray);
			}

			fclose(planetsArray);
			break;
		}
		case 2:
		{
			generatePlanets(bodiesArray, fileName);
			break;
		}
		case 3:
		{
			double probability;
			cout << "Choose N, where T = N * deltaT , where deltaT = " << DELTA_TIME << endl;
			cin >> NUMBER_OF_ITERATIONS;
			cout << "Choose probability of divide" << endl;
			cin >> probability;
			getchar();
			probabilityOfDivide = new bool[MAX_BUFFER];
			for (int i = 0; i < MAX_BUFFER; i++) {
				probabilityOfDivide[i] = getRandomNumber(probability);
			}
			offset = 0;
			subMenu();
			break;
		}
		case 4:
		{
			cout << "Select the EPSILON you want" << endl;
			cin >> EPSILON;
			getchar();
			cout << "If the program finds a mismatch for the given EPSILON = " << EPSILON << " , it will print them " << endl;
			cout.fixed;
			cout.precision(10);

			for (int i = 0; i < NUMBER_OF_BODIES; i++) {
				if ((abs(alterBodiesArray[i].coordX / bodiesArray[i].coordX - 1) < (EPSILON)) &&
					(abs(alterBodiesArray[i].coordY / bodiesArray[i].coordY - 1) < (EPSILON)) &&
					(abs(alterBodiesArray[i].coordZ / bodiesArray[i].coordZ - 1) < (EPSILON)) &&
					(abs(alterBodiesArray[i].velocityX / bodiesArray[i].velocityX - 1) < (EPSILON)) &&
					(abs(alterBodiesArray[i].velocityY / bodiesArray[i].velocityY - 1) < (EPSILON)) &&
					(abs(alterBodiesArray[i].velocityZ / bodiesArray[i].velocityZ - 1) < (EPSILON))) {

				}
				else {
					cout << "---------------------------------" << endl;
					cout << "mult body number: " << i << endl;
					cout << "coords " << endl;
					cout << alterBodiesArray[i].coordX << "  " << alterBodiesArray[i].coordY << "  " << alterBodiesArray[i].coordZ << endl;
					cout << "velocities " << endl;
					cout << alterBodiesArray[i].velocityX << "  " << alterBodiesArray[i].velocityY << "  " << alterBodiesArray[i].velocityZ << endl;
					cout << "---------------------------------" << endl;
					cout << "single body number: " << i << endl;
					cout << "coords " << endl;
					cout << alterBodiesArray[i].coordX << "  " << alterBodiesArray[i].coordY << "  " << alterBodiesArray[i].coordZ << endl;
					cout << "velocities " << endl;
					cout << alterBodiesArray[i].velocityX << "  " << alterBodiesArray[i].velocityY << "  " << alterBodiesArray[i].velocityZ << endl;
					cout << "---------------------------------" << endl;
					cout << "differense between bodies number: " << i << endl;
					cout << "coords " << endl;
					cout << abs(alterBodiesArray[i].coordX / bodiesArray[i].coordX - 1) << "  " << abs(alterBodiesArray[i].coordY / bodiesArray[i].coordY - 1) << "  " << abs(alterBodiesArray[i].coordZ / bodiesArray[i].coordZ - 1) << endl;
					cout << "velocities " << endl;
					cout << abs(alterBodiesArray[i].velocityX / bodiesArray[i].velocityX - 1) << "  " << abs(alterBodiesArray[i].velocityY / bodiesArray[i].velocityY - 1) << "  " << abs(alterBodiesArray[i].velocityZ / bodiesArray[i].velocityZ - 1) << endl;
					cout << "---------------------------------" << endl;

				}
			}

			cout << "Choose output file for single calculation: " << endl;
			scanf("%[^\n]", fileName);
			planetsArray = fopen(fileName, "wb");
			fwrite(&NUMBER_OF_BODIES, sizeof(int), 1, planetsArray);
			for (int i = 0; i < NUMBER_OF_BODIES; i++) {
				fwrite(&bodiesArray[i], sizeof(SphericalBody), 1, planetsArray);
			}
			fclose(planetsArray);

			getchar();
			cout << "Choose output file for multiple calculation: " << endl;
			scanf("%[^\n]", fileName);
			planetsArray = fopen(fileName, "wb");
			fwrite(&NUMBER_OF_BODIES, sizeof(int), 1, planetsArray);
			for (int i = 0; i < NUMBER_OF_BODIES; i++) {
				fwrite(&bodiesArray[i], sizeof(SphericalBody), 1, planetsArray);
			}
			fclose(planetsArray);

			break;
		}
		default:
			break;
		}

		cout << "Choose number of operation:" << endl;
		cout << "[1] - Choose file from disk" << endl;
		cout << "[2] - Generate data" << endl;
		cout << "[3] - Calculation from a single- and multi-threaded program" << endl;
		cout << "[4] - Compare results" << endl;
		cout << "[0] - Exit" << endl;

		cin >> menuKey;
		getchar();
	}

}

