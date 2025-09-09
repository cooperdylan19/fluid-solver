/*
This C script solves Euler's equations in 1D, using the two methods covered in lectures 4, 5, 6
    (x) achieved using the numerical Lax Friedrichs method.
    () achieved using the HLL method.
Soltuions used to produce shock tubes in terms of density, velocity, pressure, internal energy.

[] denotes a reference/ source of information 
// is used to talk through code
*/

// CODE BEGINS HERE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // import relevent libraries

typedef struct // Euler's equations (1D) can be expressed (dQ / dt) + (dF / dx) = 0 [Lecture 4]
{
    double *Q1;
    double *Q2;
    double *Q3;
} Q;

Q currentVectorQ; // this will hold the current values at each spatial grid point Q = [Q1, Q2, Q3]

typedef struct // F(Q)
{
    double *F1;
    double *F2;
    double *F3;
} F;

F currentVectorF; // "..." F = [F1, F2, F3]

/*
function below determines initial values for quantities relevent to producing shock tubes.
    ie velocity simply found from v = rho * v / rho
called in main before calling numerical solver function (LF or HLL)
*/
void initialize(double adiabaticIndex, // gamma() is defined in math.h
            double rhoLEFT, double rhoRIGHT, double vLEFT, double vRIGHT, double pLEFT, double pRIGHT, // BCs
            int N, double L, // no. points N in interval [0,L]
            double *e, double *p, // specific internal energy, pressure
            double shockPosition)  // different for each figure
{
    double currentPosition; // an integer no. spatial steps
    double spatialStep = L / N; // calculate stepsize which covers full domain

    int i;

    for (int i = 0; i < N; i++)
    {
        currentPosition = i * spatialStep;

        if (currentPosition < shockPosition) 
        {
            currentVectorQ.Q1[i] = rhoLEFT; // rho = mass density 
            currentVectorQ.Q2[i] = rhoLEFT * vLEFT; // rho * v = momentum density

            p[i] = pLEFT; // pressure
            e[i] = pLEFT / ((adiabaticIndex - 1) * rhoLEFT); // internal energy. known equation [Lecture 4]
        } else 
        {
            currentVectorQ.Q1[i] = rhoRIGHT;
            currentVectorQ.Q2[i] = rhoRIGHT * vRIGHT;

            p[i] = pRIGHT;
            e[i] = pRIGHT / ((adiabaticIndex - 1) * rhoRIGHT);
        }

        double q = currentVectorQ.Q2[i] / currentVectorQ.Q1[i]; // defined for convenience. Line 68 was long

        currentVectorQ.Q3[i] = (p[i] / (adiabaticIndex - 1)) + 0.5 * currentVectorQ.Q1[i] * pow(q, 2);
        // Q3 calculated using Q2 and Q1, so not in 'if, else' loop


        currentVectorF.F1[i] = currentVectorQ.Q2[i];
        currentVectorF.F2[i] = (currentVectorQ.Q2[i] * currentVectorQ.Q2[i]) / currentVectorQ.Q1[i] + p[i];
        currentVectorF.F3[i] = currentVectorQ.Q2[i] * ((currentVectorQ.Q3[i] + p[i]) / currentVectorQ.Q1[i]);
        // calculate corresponding F(Q) values. known equations [Lecture 4]
    }
}

/*
Lax Friedrichs scheme implementation.
    1. allocate memory for new arrays to store updated state variables and fluxes
    2. copy initial flux and state variable values to temporary structures
    3. calculate half-step flux values at each cell boundary using averaged values and differences from neighboring cells
    4. update state variables Q based on the computed flux differences across cell boundaries
    5. calculate internal energy and pressure (e and p) from the updated state variables
    6. update the main flux and state variables with the newly computed values
    7. free allocated memory used for new arrays

    8. loop through time steps (performed in runSimulationLF defined below)
*/
void laxFriedrichsScheme(double adiabaticIndex, double timeStep, int N, double L, double *e, double *p) 
{
    double spatialStep = L / N;

    double *pNew = malloc(N * sizeof(double));
    double *eNew = malloc(N * sizeof(double)); // allocate memory for e and p

    Q copyQ; F copyF;

    copyQ.Q1 = malloc(N * sizeof(double));
    copyQ.Q2 = malloc(N * sizeof(double));
    copyQ.Q3 = malloc(N * sizeof(double)); // copy stores updated 

    for (int i = 0; i < N; ++i) 
    {
        copyQ.Q1[i] = currentVectorQ.Q1[i];
        copyQ.Q2[i] = currentVectorQ.Q2[i];
        copyQ.Q3[i] = currentVectorQ.Q3[i];
    }

    Q newCurrentVectorQ = copyQ;


    copyF.F1 = malloc(N * sizeof(double));
    copyF.F2 = malloc(N * sizeof(double));
    copyF.F3 = malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i) 
    {
        copyF.F1[i] = currentVectorF.F1[i];
        copyF.F2[i] = currentVectorF.F2[i];
        copyF.F3[i] = currentVectorF.F3[i];
    }

    F newCurrentVectorF = copyF;


    double stepRatio = spatialStep / timeStep; 
    double inverseStepRatio = 1 / stepRatio; // convenience definitions

    newCurrentVectorQ.Q1[0] = newCurrentVectorQ.Q1[1]; //outflow BCs
    newCurrentVectorQ.Q2[0] = newCurrentVectorQ.Q2[1];
    newCurrentVectorQ.Q3[0] = newCurrentVectorQ.Q3[1]; // left

    // For the right boundary:
    newCurrentVectorQ.Q1[N-1] = newCurrentVectorQ.Q1[N-2];
    newCurrentVectorQ.Q2[N-1] = newCurrentVectorQ.Q2[N-2];
    newCurrentVectorQ.Q3[N-1] = newCurrentVectorQ.Q3[N-2]; // right

    int i;
    
    for (i = 1; i < N - 1; i++) // consider grid boundaries
    {
        // a conditionally stable Lax-Friedrichs method [Lecture 5]
        double F1LEFT = 0.5 * (currentVectorF.F1[i] + currentVectorF.F1[i-1]) - 0.5 * (stepRatio * (currentVectorQ.Q1[i] - currentVectorQ.Q1[i-1]));
        double F1RIGHT = 0.5 * (currentVectorF.F1[i] + currentVectorF.F1[i+1]) + 0.5 * (stepRatio * (currentVectorQ.Q1[i] - currentVectorQ.Q1[i+1]));

        double F2LEFT = 0.5 * (currentVectorF.F2[i] + currentVectorF.F2[i-1]) - 0.5 * (stepRatio * (currentVectorQ.Q2[i] - currentVectorQ.Q2[i-1]));
        double F2RIGHT = 0.5 * (currentVectorF.F2[i] + currentVectorF.F2[i+1]) + 0.5 * (stepRatio * (currentVectorQ.Q2[i] - currentVectorQ.Q2[i+1]));

        double F3LEFT = 0.5 * (currentVectorF.F3[i] + currentVectorF.F3[i-1]) - 0.5 * (stepRatio * (currentVectorQ.Q3[i] - currentVectorQ.Q3[i-1]));
        double F3RIGHT = 0.5 * (currentVectorF.F3[i] + currentVectorF.F3[i+1]) + 0.5 * (stepRatio * (currentVectorQ.Q3[i] - currentVectorQ.Q3[i+1]));

        // combine fluxes at cell boundaries with state values that are averages over cell volumes [Lecture 5]
        newCurrentVectorQ.Q1[i] = newCurrentVectorQ.Q1[i] + (inverseStepRatio) * (F1LEFT - F1RIGHT);
        newCurrentVectorQ.Q2[i] = newCurrentVectorQ.Q2[i] + (inverseStepRatio) * (F2LEFT - F2RIGHT);
        newCurrentVectorQ.Q3[i] = newCurrentVectorQ.Q3[i] + (inverseStepRatio) * (F3LEFT - F3RIGHT);

        double velocity = newCurrentVectorQ.Q2[i] / newCurrentVectorQ.Q1[i];  
        double kineticEnergy = 0.5 * newCurrentVectorQ.Q1[i] * velocity * velocity;
        double internalEnergyUnitVolume = newCurrentVectorQ.Q3[i] - kineticEnergy; // step by step, required for internal energy

        eNew[i] = internalEnergyUnitVolume / newCurrentVectorQ.Q1[i];
        pNew[i] = (adiabaticIndex - 1) * newCurrentVectorQ.Q1[i] * eNew[i];

        p[i] = pNew[i];
        e[i] = eNew[i];

        // update flux storage with the new calculated fluxes
        newCurrentVectorF.F1[i] = newCurrentVectorQ.Q2[i];
        newCurrentVectorF.F2[i] = (newCurrentVectorQ.Q2[i] * newCurrentVectorQ.Q2[i]) / newCurrentVectorQ.Q1[i] + p[i];
        newCurrentVectorF.F3[i] = newCurrentVectorQ.Q2[i] * ((newCurrentVectorQ.Q3[i] + p[i]) / newCurrentVectorQ.Q1[i]);
    }

    currentVectorQ = newCurrentVectorQ;
    currentVectorF = newCurrentVectorF; // update F and Q for next iteration

    free(pNew);
    free(eNew); // free memory
}

/*
THIS IS THE EXACT SAME FUNCTION AS ABOVE IN SPC
    considers an extra source term
*/
void laxFriedrichsSchemeSPC(double adiabaticIndex, double timeStep, int N, double L, double *e, double *p) 
{
    double spatialStep = L / N;

    double *pNew = malloc(N * sizeof(double));
    double *eNew = malloc(N * sizeof(double)); // allocate memory for e and p

    Q copyQ; F copyF;

    copyQ.Q1 = malloc(N * sizeof(double));
    copyQ.Q2 = malloc(N * sizeof(double));
    copyQ.Q3 = malloc(N * sizeof(double)); // copy stores updated 

    for (int i = 0; i < N; ++i) 
    {
        copyQ.Q1[i] = currentVectorQ.Q1[i];
        copyQ.Q2[i] = currentVectorQ.Q2[i];
        copyQ.Q3[i] = currentVectorQ.Q3[i];
    }

    Q newCurrentVectorQ = copyQ;


    copyF.F1 = malloc(N * sizeof(double));
    copyF.F2 = malloc(N * sizeof(double));
    copyF.F3 = malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i) 
    {
        copyF.F1[i] = currentVectorF.F1[i];
        copyF.F2[i] = currentVectorF.F2[i];
        copyF.F3[i] = currentVectorF.F3[i];
    }

    F newCurrentVectorF = copyF;


    double stepRatio = spatialStep / timeStep; 
    double inverseStepRatio = 1 / stepRatio; // convenience definitions

    newCurrentVectorQ.Q1[0] = newCurrentVectorQ.Q1[1]; //outflow BCs
    newCurrentVectorQ.Q2[0] = newCurrentVectorQ.Q2[1];
    newCurrentVectorQ.Q3[0] = newCurrentVectorQ.Q3[1]; // left

    // For the right boundary:
    newCurrentVectorQ.Q1[N-1] = newCurrentVectorQ.Q1[N-2];
    newCurrentVectorQ.Q2[N-1] = newCurrentVectorQ.Q2[N-2];
    newCurrentVectorQ.Q3[N-1] = newCurrentVectorQ.Q3[N-2]; // right

    int i;
    
    for (i = 1; i < N - 1; i++) // consider grid boundaries
    {
        double S1 = (2 * timeStep * currentVectorF.F1[i]) / (i * spatialStep);
        double S2 = (2 * timeStep * (currentVectorF.F2[i]-p[i])) / (i * spatialStep);
        double S3 = (2 * timeStep * currentVectorF.F3[i]) / (i * spatialStep); // THIS IS NEW - SOURCE TERMS [reference in report]

    

        // a conditionally stable Lax-Friedrichs method [Lecture 5]
        double F1LEFT = 0.5 * (currentVectorF.F1[i] + currentVectorF.F1[i-1]) - 0.5 * (stepRatio * (currentVectorQ.Q1[i] - currentVectorQ.Q1[i-1]));
        double F1RIGHT = 0.5 * (currentVectorF.F1[i] + currentVectorF.F1[i+1]) + 0.5 * (stepRatio * (currentVectorQ.Q1[i] - currentVectorQ.Q1[i+1]));

        double F2LEFT = 0.5 * (currentVectorF.F2[i] + currentVectorF.F2[i-1]) - 0.5 * (stepRatio * (currentVectorQ.Q2[i] - currentVectorQ.Q2[i-1]));
        double F2RIGHT = 0.5 * (currentVectorF.F2[i] + currentVectorF.F2[i+1]) + 0.5 * (stepRatio * (currentVectorQ.Q2[i] - currentVectorQ.Q2[i+1]));

        double F3LEFT = 0.5 * (currentVectorF.F3[i] + currentVectorF.F3[i-1]) - 0.5 * (stepRatio * (currentVectorQ.Q3[i] - currentVectorQ.Q3[i-1]));
        double F3RIGHT = 0.5 * (currentVectorF.F3[i] + currentVectorF.F3[i+1]) + 0.5 * (stepRatio * (currentVectorQ.Q3[i] - currentVectorQ.Q3[i+1]));

        // combine fluxes at cell boundaries with state values that are averages over cell volumes [Lecture 5]
        newCurrentVectorQ.Q1[i] = newCurrentVectorQ.Q1[i] + (inverseStepRatio) * (F1LEFT - F1RIGHT) - S1;
        newCurrentVectorQ.Q2[i] = newCurrentVectorQ.Q2[i] + (inverseStepRatio) * (F2LEFT - F2RIGHT) - S2;
        newCurrentVectorQ.Q3[i] = newCurrentVectorQ.Q3[i] + (inverseStepRatio) * (F3LEFT - F3RIGHT) - S3; // SOURCE TERMS SUBTRACTED FROM Q HERE

        double velocity = newCurrentVectorQ.Q2[i] / newCurrentVectorQ.Q1[i];  
        double kineticEnergy = 0.5 * newCurrentVectorQ.Q1[i] * velocity * velocity;
        double internalEnergyUnitVolume = newCurrentVectorQ.Q3[i] - kineticEnergy; // step by step, required for internal energy

        eNew[i] = internalEnergyUnitVolume / newCurrentVectorQ.Q1[i];
        pNew[i] = (adiabaticIndex - 1) * newCurrentVectorQ.Q1[i] * eNew[i];

        p[i] = pNew[i];
        e[i] = eNew[i];

        // update flux storage with the new calculated fluxes
        newCurrentVectorF.F1[i] = newCurrentVectorQ.Q2[i];
        newCurrentVectorF.F2[i] = (newCurrentVectorQ.Q2[i] * newCurrentVectorQ.Q2[i]) / newCurrentVectorQ.Q1[i] + p[i];
        newCurrentVectorF.F3[i] = newCurrentVectorQ.Q2[i] * ((newCurrentVectorQ.Q3[i] + p[i]) / newCurrentVectorQ.Q1[i]);
    }

    currentVectorQ = newCurrentVectorQ;
    currentVectorF = newCurrentVectorF; // update F and Q for next iteration

    free(pNew);
    free(eNew); // free memory
}

/*
HLL scheme implementation
follows same procedures as LF Scheme with new flux formula + defined variables
*/
void HLLScheme(double adiabaticIndex, double timeStep, int N, double L, double *e, double *p) 
{
    double spatialStep = L / N;

    double *pNew = malloc(N * sizeof(double));
    double *eNew = malloc(N * sizeof(double));

    Q copyQ;
    F copyF;

    copyQ.Q1 = malloc(N * sizeof(double));
    copyQ.Q2 = malloc(N * sizeof(double));
    copyQ.Q3 = malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i) 
    {
        copyQ.Q1[i] = currentVectorQ.Q1[i];
        copyQ.Q2[i] = currentVectorQ.Q2[i];
        copyQ.Q3[i] = currentVectorQ.Q3[i];
    }

    Q newCurrentVectorQ = copyQ;

    copyF.F1 = malloc(N * sizeof(double));
    copyF.F2 = malloc(N * sizeof(double));
    copyF.F3 = malloc(N * sizeof(double));

    for (int i = 0; i < N; ++i) 
    {
        copyF.F1[i] = currentVectorF.F1[i];
        copyF.F2[i] = currentVectorF.F2[i];
        copyF.F3[i] = currentVectorF.F3[i];
    }

    F newCurrentVectorF = copyF;

    double stepRatio = spatialStep / timeStep; 
    double inverseStepRatio = 1 / stepRatio;

    for (int i = 1; i < N - 1; i++) // again, consider grid
    {
        double vL = currentVectorQ.Q2[i-1] / currentVectorQ.Q1[i-1];
        double vR = currentVectorQ.Q2[i] / currentVectorQ.Q1[i]; // define relevent variables for guess

        double pL = p[i-1];
        double pR = p[i];

        double rhoL = currentVectorQ.Q1[i-1];
        double rhoR = currentVectorQ.Q1[i];

        double csL = sqrt(adiabaticIndex * pL / rhoL);
        double csR = sqrt(adiabaticIndex * pR / rhoR);

        double pGuess = 0.5 * (pL + pR) - 0.5 * (vR - vL) * 0.5 * (rhoL + rhoR) * 0.5 * (csL + csR);

        double qL = (pGuess > pL) ? sqrt(1 + ((adiabaticIndex + 1) / (2 * adiabaticIndex)) * (pGuess / pL - 1)) : 1;
        double qR = (pGuess > pR) ? sqrt(1 + ((adiabaticIndex + 1) / (2 * adiabaticIndex)) * (pGuess / pR - 1)) : 1;

        double SL = vL - (csL * qL);
        double SR = vR + (csR * qR);

        double F1, F2, F3;
        if (SL >= 0) 
        {
            F1 = currentVectorF.F1[i-1];
            F2 = currentVectorF.F2[i-1];
            F3 = currentVectorF.F3[i-1];
        } else if (SR <= 0) 
        {
            F1 = currentVectorF.F1[i];
            F2 = currentVectorF.F2[i];
            F3 = currentVectorF.F3[i];
        } else 
        {
            F1 = (SR * currentVectorF.F1[i-1] - SL * currentVectorF.F1[i] + SL * SR * (currentVectorQ.Q1[i] - currentVectorQ.Q1[i-1])) / (SR - SL);
            F2 = (SR * currentVectorF.F2[i-1] - SL * currentVectorF.F2[i] + SL * SR * (currentVectorQ.Q2[i] - currentVectorQ.Q2[i-1])) / (SR - SL);
            F3 = (SR * currentVectorF.F3[i-1] - SL * currentVectorF.F3[i] + SL * SR * (currentVectorQ.Q3[i] - currentVectorQ.Q3[i-1])) / (SR - SL);
        }

        newCurrentVectorQ.Q1[i] = newCurrentVectorQ.Q1[i] - inverseStepRatio * (newCurrentVectorF.F1[i-1] - newCurrentVectorF.F1[i]);
        newCurrentVectorQ.Q2[i] = newCurrentVectorQ.Q2[i] - inverseStepRatio * (newCurrentVectorF.F2[i-1] - newCurrentVectorF.F2[i]);
        newCurrentVectorQ.Q3[i] = newCurrentVectorQ.Q3[i] - inverseStepRatio * (newCurrentVectorF.F3[i-1] - newCurrentVectorF.F3[i]);

        double velocity = newCurrentVectorQ.Q2[i] / newCurrentVectorQ.Q1[i];  
        double kineticEnergy = 0.5 * newCurrentVectorQ.Q1[i] * velocity * velocity;
        double internalEnergyUnitVolume = newCurrentVectorQ.Q3[i] - kineticEnergy;

        eNew[i] = internalEnergyUnitVolume / newCurrentVectorQ.Q1[i];
        pNew[i] = (adiabaticIndex - 1) * newCurrentVectorQ.Q1[i] * eNew[i];

        p[i] = pNew[i];
        e[i] = eNew[i];
    }

    currentVectorQ = newCurrentVectorQ;
    currentVectorF = newCurrentVectorF;

    free(pNew);
    free(eNew);
}

/*
functions below merely for convenience
runSim: uses numerical scheme functions over full time period
ResultsToFile: trivial
*/
typedef void (*SchemeFunction)(double, double, int, double, double *, double *); // enables function to be argument for runSimulation

void runSimulation(double adiabaticIndex, double t1, int N, double L, double *e, double *p, double initialTimeStep, SchemeFunction scheme) 
{
    double t = 0.0; // start time
    double timeStep = initialTimeStep;

    while (t < t1) 
    {
        if (t + timeStep > t1) 
        {
            timeStep = t1 - t;
        }

        // call the scheme function passed as an argument
        scheme(adiabaticIndex, timeStep, N, L, e, p);

        t += timeStep;
    }
}

int writeResultsToFile(const char *filename, double *Q1, double *Q2, double *p, double *e, int N) 
{
    FILE *outputFile = fopen(filename, "w");
    if (outputFile == NULL) 
    {
        fprintf(stderr, "Error opening output file.\n");
        return 1; // return an error code indicating failure
    }

    for (int i = 0; i < N; i++) 
    {
        double velocity = Q2[i] / Q1[i]; // calculate velocity at point i
        fprintf(outputFile, "%f, %f, %f, %f, %f\n", (double) i / N, Q1[i], velocity, p[i], e[i]);
    }

    fclose(outputFile);
    return 0; // success
}




int main() // main is used for definining arguments for function calling, allocating memory, exporting data to file
{
    int N = 100; // set N in main

    currentVectorQ.Q1 = malloc(N * sizeof(double));
    currentVectorQ.Q2 = malloc(N * sizeof(double));
    currentVectorQ.Q3 = malloc(N * sizeof(double));

    currentVectorF.F1 = malloc(N * sizeof(double));
    currentVectorF.F2 = malloc(N * sizeof(double));
    currentVectorF.F3 = malloc(N * sizeof(double));

    double *e = malloc(N * sizeof(double));
    double *p = malloc(N * sizeof(double)); // allocate memory

    if (!currentVectorQ.Q1 || !currentVectorQ.Q2 || !currentVectorQ.Q3 || !currentVectorF.F1 || !currentVectorF.F2 || !currentVectorF.F3) 
    {
        fprintf(stderr, "Failed to allocate memory for storage structures.\n"); 
        return 1;
    }

    if (!e || !p) 
    {
        fprintf(stderr, "Error allocating memory for e or p\n"); // memory allocation error messages
        return 1;
    }


    // now we will begin to call numerical solver functions...

    double L = 1.0; // domain length
    double spatialStep = L / N;
    double preFactor =  0.3; // conventional value for CFL condition [lecture 5]
    double adiabaticIndex = 1.4;


    /*
    calling functions and definings values to produce figure 1 (LF)
    data -> CSV file LaxFreidrichsData1
    */
    double t1_1 = 0.2; // time of snapshot
    double rhoLEFT_1 = 1.0, rhoRIGHT_1 = 0.125, vLEFT_1 = 0.75, vRIGHT_1 = 0.0, pLEFT_1 = 1.0, pRIGHT_1 = 0.1; // figure 1

    double soundSpeed_1 = sqrt(adiabaticIndex * pLEFT_1 / rhoLEFT_1);
    double maxSound_1 = fabs(vLEFT_1) + soundSpeed_1; // consider flow velocity
    double timeStep_1 = preFactor * spatialStep / maxSound_1; // CFL-based time step calculation

    double shockPosition_1 = 0.3;


    initialize(adiabaticIndex, rhoLEFT_1, rhoRIGHT_1, vLEFT_1, vRIGHT_1, pLEFT_1, pRIGHT_1, N, L, e, p, shockPosition_1);

    runSimulation(adiabaticIndex, t1_1, N, L, e, p, timeStep_1, laxFriedrichsScheme);

    if (writeResultsToFile("LaxFriedrichsData1.csv", currentVectorQ.Q1, currentVectorQ.Q2, p, e, N)) 
    {
        return 1; // Exit if writing to file fails
    }

    /*
    calling functions and definings values to produce figure 2 (LF)
    data -> CSV file LaxFreidrichsData2
    */
    double t1_2 = 0.15; // new time of snapshot
    double rhoLEFT_2 = 1.0, rhoRIGHT_2 = 1.0, vLEFT_2 = -2.0, vRIGHT_2 = 2.0, pLEFT_2 = 0.4, pRIGHT_2 = 0.4; // new BCs for figure 2


    double soundSpeed_2 = sqrt(adiabaticIndex * pLEFT_2 / rhoLEFT_2);
    double maxSound_2 = fabs(vLEFT_2) + soundSpeed_2; // accounting for the flow velocity
    double timeStep_2 = preFactor * spatialStep / maxSound_2; // CFL time step calculation

    double shockPosition_2 = 0.5;


    initialize(adiabaticIndex, rhoLEFT_2, rhoRIGHT_2, vLEFT_2, vRIGHT_2, pLEFT_2, pRIGHT_2, N, L, e, p, shockPosition_2);

    runSimulation(adiabaticIndex, t1_2, N, L, e, p, timeStep_2, laxFriedrichsScheme);

    if (writeResultsToFile("LaxFriedrichsData2.csv", currentVectorQ.Q1, currentVectorQ.Q2, p, e, N)) 
    {
        return 1; // exit if writing to file fails
    }

    /*
    calling functions and definings values to produce figure 1 (HLL)
    data -> CSV file HLLData1
    */
    double t1_1H = 0.2; // time of snapshot
    double rhoLEFT_1H = 1.0, rhoRIGHT_1H = 0.125, vLEFT_1H = 0.75, vRIGHT_1H = 0.0, pLEFT_1H = 1.0, pRIGHT_1H = 0.1; // figure 1

    double soundSpeed_1H = sqrt(adiabaticIndex * pLEFT_1H / rhoLEFT_1H);
    double maxSound_1H = fabs(vLEFT_1H) + soundSpeed_1H; // accounting for flow velocity 
    double timeStep_1H = preFactor * spatialStep / maxSound_1H; // CFL time step calculation

    double shockPosition_1H = 0.3;


    initialize(adiabaticIndex, rhoLEFT_1H, rhoRIGHT_1H, vLEFT_1H, vRIGHT_1H, pLEFT_1H, pRIGHT_1H, N, L, e, p, shockPosition_1H);

    runSimulation(adiabaticIndex, t1_1, N, L, e, p, timeStep_1, HLLScheme);

    if (writeResultsToFile("HLLData1.csv", currentVectorQ.Q1, currentVectorQ.Q2, p, e, N)) 
    {
        return 1; // exit if writing to file fails
    }

    /*
    calling functions and definings values to produce figure 2 (HLL)
    data -> CSV file HLLData2
    */
    double t1_2H = 0.15; // new time of snapshot
    double rhoLEFT_2H = 1.0, rhoRIGHT_2H = 1.0, vLEFT_2H = -2.0, vRIGHT_2H = 2.0, pLEFT_2H = 0.4, pRIGHT_2H = 0.4; // new BCs for figure 2


    double soundSpeed_2H = sqrt(adiabaticIndex * pLEFT_2H / rhoLEFT_2H);
    double maxSound_2H = fabs(vLEFT_2H) + soundSpeed_2; // accounting for the flow velocity
    double timeStep_2H = preFactor * spatialStep / maxSound_2H; // CFL time step calculation

    double shockPosition_2H = 0.5;


    initialize(adiabaticIndex, rhoLEFT_2H, rhoRIGHT_2H, vLEFT_2H, vRIGHT_2H, pLEFT_2H, pRIGHT_2H, N, L, e, p, shockPosition_2H);

    runSimulation(adiabaticIndex, t1_2, N, L, e, p, timeStep_2, HLLScheme);

    if (writeResultsToFile("HLLData2.csv", currentVectorQ.Q1, currentVectorQ.Q2, p, e, N)) 
    {
        return 1; // exit if writing to file fails
    }

    /*
    calling functions and definings values to produce figure 3 (LF, SPC)
    data -> CSV file LaxFreidrichsData3
    */
    double t1_3 = 0.25; // new time of snapshot
    double rhoLEFT_3 = 1.0, rhoRIGHT_3 = 0.125, vLEFT_3 = 0.0, vRIGHT_3 = 0.0, pLEFT_3 = 1.0, pRIGHT_3 = 0.1; // new BCs for figure 2


    double soundSpeed_3 = sqrt(adiabaticIndex * pLEFT_3 / rhoLEFT_3);
    double maxSound_3 = fabs(vLEFT_3) + soundSpeed_3; // accounting for the flow velocity
    double timeStep_3 = preFactor * spatialStep / maxSound_3; // CFL time step calculation

    double shockPosition_3 = 0.4;


    initialize(adiabaticIndex, rhoLEFT_3, rhoRIGHT_3, vLEFT_3, vRIGHT_3, pLEFT_3, pRIGHT_3, N, L, e, p, shockPosition_3);

    runSimulation(adiabaticIndex, t1_3, N, L, e, p, timeStep_3, laxFriedrichsSchemeSPC);

    if (writeResultsToFile("LaxFriedrichsData3.csv", currentVectorQ.Q1, currentVectorQ.Q2, p, e, N)) 
    {
        return 1; // exit if writing to file fails
    }


    free(currentVectorQ.Q1);
    free(currentVectorQ.Q2);
    free(currentVectorQ.Q3);

    free(currentVectorF.F1);
    free(currentVectorF.F2);
    free(currentVectorF.F3);

    free(e);
    free(p); // free memory 

    return 0;
}