/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"
#include "fem.h"
#include <time.h>

double fun(double x, double y) 
{
    return 1;
}

double computeTime(clock_t start, clock_t end) {
    return (double)(end - start) / CLOCKS_PER_SEC;
}

void createNewMesh(const char *filename, double (*sizeCallback)(double x, double y)) {
    femGeo* theGeometry = geoGetGeometry();

    // Set the size callback
    geoSetSizeCallback(sizeCallback);

    // Set geometry parameters
    double rInner = 1;
    double rOuter = 1.2;

    theGeometry->xCenter = 0.0;
    theGeometry->yCenter = 0.0;
    theGeometry->rInner  = rInner;
    theGeometry->rOuter  = rOuter;
    theGeometry->dInner  = rOuter - rInner;
    theGeometry->dOuter  = 0.66 * rInner;
    theGeometry->hInner  = 0.05;
    theGeometry->hOuter  = 0.05;
    theGeometry->h       = 0.15;
    theGeometry->elementType = FEM_TRIANGLE;

    // Generate and import the mesh
    geoMeshGenerate();
    geoMeshImport();
    geoSetDomainName(1, "ExternalBoundary");
    geoSetDomainName(0, "InternalBoundary");
    geoMeshWrite(filename);
}

femProblem* createNewProblem(femGeo* theGeometry, const char *filename) {
    double vMass = 1200; // mass of the vehicle in kg
    double wMass = vMass / 4; // mass distributed on one wheel
    double E = 2e6; // Young's modulus in Pa
    double nu = 0.49; // Poisson's ratio
    double rho = 9.02e3; // Density in kg/m³
    double g = 9.81;

    femProblem* theProblem = femElasticityCreate(theGeometry, E, nu, rho, g, PLANAR_STRESS);

    // Weight of the car (force applied downward)
    double F_car = g * wMass;
    // Reaction force from the ground (force applied upward)
    double F_reaction = F_car;
    // Internal air pressure
    double P_internal = 200000.0; // 200 kPa (2 bars)

    femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", NEUMANN_Y, -F_car); // Force vers le bas aucune diff entre les 2 axes
    femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", NEUMANN_Y, F_reaction); // Force vers le haut
    femElasticityAddBoundaryCondition(theProblem, "InternalBoundary", NEUMANN_X, P_internal); // Pression interne (x) aucune diff entre les 2 axes
    femElasticityAddBoundaryCondition(theProblem, "InternalBoundary", NEUMANN_Y, P_internal); // Pression interne (y)

    // for (int i = 0; i < theGeometry->theNodes->nNodes; i++) {
    //     double x = theGeometry->theNodes->X[i];
    //     double y = theGeometry->theNodes->Y[i];
    //     double xCenter = theGeometry->xCenter;
    //     double EPSILON = 1e-2; // Tolérance pour comparer
    
    //     // Si le nœud est à gauche du centre, impose un déplacement vers la gauche
    //     if (x < xCenter - EPSILON) {
    //         femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", DIRICHLET_X, -1e-2); // Déplacement vers la gauche
    //     }
    
    //     // Si le nœud est à droite du centre, impose un déplacement vers la droite
    //     if (x > xCenter + EPSILON) {
    //         femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", DIRICHLET_X, 1e-2); // Déplacement vers la droite
    //     }
    
    //     // Bloquer les déplacements verticaux sur l'axe vertical (x = xCenter)
    //     if (fabs(x - xCenter) < EPSILON) {
    //         femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", DIRICHLET_Y, 0.0);
    //     }
    // }
    //femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", DIRICHLET_X, 1e-2); // Bloque déplacement horizontal
    //femElasticityAddBoundaryCondition(theProblem, "ExternalBoundary", DIRICHLET_Y, 0.0);
    // Déplacement symétrique selon X à gauche et à droite du centre
    //femElasticityAddBoundaryCondition(theProblem, "LeftBoundary", DIRICHLET_X, -1e-2);  // Vers la gauche
    //femElasticityAddBoundaryCondition(theProblem, "RightBoundary", DIRICHLET_X, 1e-2);  // Vers la droite

    // Bloquer les déplacements verticaux pour conserver la symétrie
    //femElasticityAddBoundaryCondition(theProblem, "LeftBoundary", DIRICHLET_Y, 0.0); 
    //femElasticityAddBoundaryCondition(theProblem, "RightBoundary", DIRICHLET_Y, 0.0);

    // condition domaine interieur
    femElasticityAddBoundaryCondition(theProblem, "InternalBoundary", DIRICHLET_X, 0.0); // Bloque déplacement horizontal
    femElasticityAddBoundaryCondition(theProblem, "InternalBoundary", DIRICHLET_Y, 0.0); // Bloque déplacement vertical
    femElasticityWrite(theProblem, filename);
    return theProblem;
}

int main(int argc, char *argv[])
{  
    double start, end;

    char *meshFile = "../data/mesh_precis.txt"; // Maillage par défaut
    char *problemFile = "../data/problem.txt"; // Problème par défaut
    femSolverType solver = SOLVEUR_BANDE;
    femRenumType renumType = FEM_XNUM;
    double (*sizeCallback)(double x, double y) = geoSize; // Choisir entre geoSize et geoSizeDefault

    // Lecture des arguments
    if (argc > 5) {
        Error("Unexpected argument");
    }
    switch(argc) {
        case 5:
            if (*argv[4] == '0')
                renumType  = FEM_NO;
            else if (*argv[4] == 'X')
                renumType  = FEM_XNUM;
            else if (*argv[4] == 'Y')
                renumType = FEM_YNUM;
            else
                Error("Unexpected argument");
        case 4:
            if (*argv[3] == 'B')
                solver  = SOLVEUR_BANDE;
            else if (*argv[3] == 'F')
                solver  = SOLVEUR_PLEIN;
            else if (*argv[3] == 'G')
                solver = GRADIENTS_CONJUGUES;
            else
                Error("Unexpected argument");
        case 3:
            problemFile = argv[2];
        case 2:
            meshFile = argv[1];
            if (argc==2) problemFile = "../data/problem.txt";
            break;
        default: printf("\nValeurs par défaut : ");
    }

    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n\n\n");

    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    //createNewMesh("../data/mesh.txt", sizeCallback); //  Décommenter cette ligne pour créer un nouveau maillage
    geoMeshRead(meshFile); // Commenter cette ligne si tu décommente celle au dessus
    
        
//
//  -2- Creation probleme 
//  

    femProblem* theProblem = createNewProblem(theGeometry, "../data/problem.txt");   // Décommenter cette ligne pour créer un nouveau problème
    //femProblem* theProblem = femElasticityRead(theGeometry, problemFile); // Commenter cette ligne si tu décommente celle au dessus
    femElasticityPrint(theProblem);
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
    double deformationFactor = 10;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    //double *normForces = malloc(theNodes->nNodes * sizeof(double));

    
    for (int i=0; i<theNodes->nNodes; i++){
        //double fx = theForces[2*i+0];
        //double fy = theForces[2*i+1];

        //normForces[i] = sqrt(fx*fx + fy*fy);
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor; //partie que j ai modif
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
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        if (t-told > 0.5) {freezingButton = FALSE; }
        
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        // if (mode == 4) {  
        //     // Calculer et afficher le champ normForces (si non déjà calculé hors boucle) (ca j ai rajouté)
        //     glfemPlotField(theGeometry->theElements, normForces);
        //     sprintf(theMessage, "Visualisation des efforts internes (norme) ");
        //     glColor3f(1.0,0.0,0.0);
        //     glfemMessage(theMessage); }
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    printf("\nComputing solution takes %.6f seconds\n", computeTime(start, end));

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    //free(normForces);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    
    return 0;  
}
