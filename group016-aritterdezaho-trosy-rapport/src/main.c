
 
#include "fem.h"
#include <time.h>

double fun(double x, double y) 
{
    return 1;
}

double computeTime(clock_t start, clock_t end) {
    return (double)(end - start) / CLOCKS_PER_SEC;
}


int main(int argc, char *argv[])
{  
    double start, end;

    char *meshFile = "../data/mesh_precis.txt"; // Maillage par défaut
    char *problemFile = "../data/problem_inconel600.txt"; // Problème par défaut : inconel 600
    femSolverType solver = SOLVEUR_BANDE;
    femRenumType renumType = FEM_XNUM;

    // Lecture des arguments
    if (argc > 6) {
        Error("Unexpected argument");
    }
    switch (argc) {
        case 5:
            if (*argv[4] == '0')
                renumType = FEM_NO;
            else if (*argv[4] == 'X')
                renumType = FEM_XNUM;
            else if (*argv[4] == 'Y')
                renumType = FEM_YNUM;
            else
                Error("Unexpected argument");
        case 4:
            if (*argv[3] == 'B')
                solver = SOLVEUR_BANDE;
            else if (*argv[3] == 'F')
                solver = SOLVEUR_PLEIN;
            else if (*argv[3] == 'G')
                solver = GRADIENTS_CONJUGUES;
            else
                Error("Unexpected argument");
        case 3:
            problemFile = argv[2];
        case 2:
            meshFile = argv[1];
            if (argc == 2) problemFile = "../data/problem.txt";
            break;
        default:
            printf("\nValeurs par défaut (inconel600) : \n");
    }

    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    geoMeshRead(meshFile); 
    
//
//  -2- Creation probleme 
//  

    femProblem* theProblem = femElasticityRead(theGeometry, problemFile); 
    theProblem->solver = solver;
    theProblem->renumType = renumType;
    
//
//  -3- Resolution du probleme et calcul des forces
//

    start = clock();
    double *theSoluce = femElasticitySolve(theProblem);
    double *theForces = femElasticityForces(theProblem);
    double area = femElasticityIntegrate(theProblem, fun); 
    end = clock();
   
//
//  -4- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    
    femNodes *theNodes = theGeometry->theNodes;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    double deformationFactor = 1;
    
    for (int i=0; i<theNodes->nNodes; i++){
        
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor; 
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1]; }
  
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

//
//  -5- Calcul de la force globaleresultante
//

    double theGlobalForce[2] = {0, 0};
    for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1]; }
    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * theProblem->rho * theProblem->g);

//
//  -6- Visualisation du maillage
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
    printf("\nComputing solution takes %.6f seconds\n", computeTime(start, end));

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    
    exit(EXIT_SUCCESS);
    
    return 0;  
}
