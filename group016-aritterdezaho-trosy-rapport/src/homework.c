#include "fem.h"

double* val;


// Strip : BEGIN
double **A_copy = NULL;
double *B_copy  = NULL;
// Strip : END

int compare(const void *a, const void *b) {
    if (val[*(int*)a] < val[*(int*)b]) return 1;
    if (val[*(int*)a] > val[*(int*)b]) return -1;
    return 0;
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int nNodes = theMesh->nodes->nNodes;
    int *tab = malloc(sizeof(int) * nNodes);

    if (tab == NULL) {
        Error("Memory allocation failed for renumbering array");
    }

    for (int i = 0; i < nNodes; i++) {
        tab[i] = i;
    }

    if (renumType == FEM_XNUM || renumType == FEM_YNUM) {
        val = (renumType == FEM_XNUM) ? theMesh->nodes->X : theMesh->nodes->Y;
        if (val == NULL) {
            free(tab);
            return;
        }
        qsort(tab, nNodes, sizeof(int), compare);
    }

    for (int i = 0; i < nNodes; i++) {
        theMesh->number[tab[i]] = i;
    }

    free(tab);
}

int femMeshComputeBand(femMesh *theMesh) {
    int nLocal = theMesh->nLocalNode;
    int myBand = 0;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        int map[4];
        int myMin, myMax;

        for (int j = 0; j < nLocal; j++) {
            map[j] = theMesh->number[theMesh->elem[iElem * nLocal + j]];
        }

        myMin = map[0];
        myMax = map[0];
        for (int j = 1; j < nLocal; j++) {
            myMax = fmax(map[j], myMax);
            myMin = fmin(map[j], myMin);
        }

        myBand = fmax(myBand, myMax - myMin);
    }

    return (myBand + 1) * 2;
}

void conjugateGradient(double **A, double *b, double *x, int size) {
    double *r = malloc(size * sizeof(double));
    double *d = malloc(size * sizeof(double));
    double *s = malloc(size * sizeof(double));
    if (r == NULL || d == NULL || s == NULL) {
        Error("Memory allocation failed for conjugate gradient vectors");
    }

    double alpha, beta, delta, deltaNew;
    int nMax = 2000;

    // Initialization
    for (int i = 0; i < size; i++) {
        r[i] = b[i];
        d[i] = r[i];
        x[i] = 0.0;
    }

    delta = 0.0;
    for (int i = 0; i < size; i++) {
        delta += r[i] * r[i];
    }

    int k = 0;
    while (delta > 1e-8 && k < nMax) {  
        for (int i = 0; i < size; i++) {
            s[i] = 0.0;
            for (int j = 0; j < size; j++) {
                s[i] += A[i][j] * d[j];
            }
        }

        double sd = 0.0;
        for (int i = 0; i < size; i++) {
            sd += d[i] * s[i];
        }
        alpha = delta / sd;

        for (int i = 0; i < size; i++) {
            x[i] += alpha * d[i];
            r[i] -= alpha * s[i];
        }

        deltaNew = 0.0;
        for (int i = 0; i < size; i++) {
            deltaNew += r[i] * r[i];
        }

        beta = deltaNew / delta;

        for (int i = 0; i < size; i++) {
            d[i] = r[i] + beta * d[i];
        }

        delta = deltaNew;
        k++;
    }

    printf("Convergence after %d iterations.\n", k);

    free(r);
    free(d);
    free(s);
}

void femElasticityAssembleElements(femProblem *theProblem)
{
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femGeo *theGeometry = theProblem->geometry;
    femNodes *theNodes = theGeometry->theNodes;
    femMesh *theMesh = theGeometry->theElements;

    double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    int nLocal = theMesh->nLocalNode;
    double a = theProblem->A, b = theProblem->B, c = theProblem->C;
    double rho = theProblem->rho, g = theProblem->g;
    double **A = theSystem->A;
    double *B = theSystem->B;

    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        int map[4], mapX[4], mapY[4];

        for (int j = 0; j < nLocal; j++) {
            map[j] = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
        }

        for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (int i = 0; i < theSpace->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (int i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            double weightedJac = jac * weight;

            for (int i = 0; i < theSpace->n; i++) {
                for (int j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * weightedJac;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * weightedJac;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * weightedJac;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * weightedJac;
                }
                B[mapY[i]] -= phi[i] * g * rho * weightedJac;
            }
        }
    }
}

void femElasticityAssembleNeumann(femProblem *theProblem)
{
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete *theSpace = theProblem->spaceEdge;
    femGeo *theGeometry = theProblem->geometry;
    femNodes *theNodes = theGeometry->theNodes;
    femMesh *theEdges = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int nLocal = 2;
    double *B = theSystem->B;

    for (int iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;

        if (type == DIRICHLET_X || type == DIRICHLET_Y) {
            continue;
        }

        femDomain *domain = theCondition->domain;
        double value = theCondition->value;
        int shift = (type == NEUMANN_X) ? 0 : 1;

        for (int iEdge = 0; iEdge < domain->nElem; iEdge++) {
            int iElem = domain->elem[iEdge];
            int map[2], mapU[2];

            for (int j = 0; j < nLocal; j++) {
                map[j] = theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
            }

            double dx = x[1] - x[0];
            double dy = y[1] - y[0];
            double length = sqrt(dx * dx + dy * dy);
            double jac = length / 2;

            for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
                double xsi = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace, xsi, phi);

                for (int i = 0; i < theSpace->n; i++) {
                    B[mapU[i]] += phi[i] * value * jac * weight;
                }
            }
        }
    }
}

void femApplyBoundaryConditions(femProblem *theProblem) {

    femFullSystem *theSystem = theProblem->system;
    femMesh *theMesh         = theProblem->geometry->theElements;
    double **A  = theSystem->A;
    double *B   = theSystem->B;
    int iCase   = theProblem->planarStrainStress;

    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {

        femBoundaryCondition* cnd = theProblem->conditions[i];
        femMesh * bndMesh = cnd->domain->mesh;
        femBoundaryType bndType = cnd->type;
        double * X = bndMesh->nodes->X;
        double * Y = bndMesh->nodes->Y;
        int * bndElem = cnd->domain->elem;
        int nElem = cnd->domain->nElem;
        int node0, node1;
                   

        if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y) {
            for (int j = 0; j < nElem; j++) {
                node0 = bndMesh->elem[2 * bndElem[j]];
                node1 = bndMesh->elem[2 * bndElem[j] + 1];                
                double jac = 0.5 * sqrt( (X[node0] - X[node1]) *  (X[node0] - X[node1]) +
                                    (Y[node0] - Y[node1]) * (Y[node0] - Y[node1]));

                if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y) {
                    if (iCase == AXISYM) { // Prise en compte du cas axisymétrique
                        
                        // ce sont les deux valeurs prises par les fcts de formes, mais pas les fcts de formes
                        double phi[2] = { (1+sqrt(3))/2, (1-sqrt(3))/2}; 
                        double x[2] = {X[node0], X[node1]};
                        double xLoc[2] = {x[0] * phi[1] + x[1] * phi[0], x[0] * phi[0] + x[1] * phi[1]};

                        int shift = (cnd->type == NEUMANN_X ) ? 0 : 1;    
                        node0 = theMesh->number[node0]; // renumérotation
                        node1 = theMesh->number[node1];       
                        B[2 * node0 + shift] += jac * cnd->value * phi[1] * xLoc[0] ; // premier point d'intégration
                        B[2 * node0 + shift] += jac * cnd->value * phi[0] * xLoc[1] ; // deuxième point d'intégration

                        B[2 * node1 + shift] += jac * cnd->value * phi[0] * xLoc[0] ; // premier point d'intégration
                        B[2 * node1 + shift] += jac * cnd->value * phi[1] * xLoc[1] ; // deuxième point d'intégration                                                                             
                    }

                    else {
                        node0 = theMesh->number[node0]; // renumérotation
                        node1 = theMesh->number[node1];
                        int shift = (cnd->type == NEUMANN_X ) ? 0 : 1;
                        B[2 * node0 + shift] += jac * cnd->value;
                        B[2 * node1 + shift] += jac * cnd->value;                      
                    }
                }

            }
        }
    }
    
    // Conditions DIRICHLET X-Y
    int *theConstrainedNodes = theProblem->constrainedNodes; // NB : ici nodes ne correspond pas à un noeud, mais à la composante x ou y
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femBoundaryType cndType = theProblem->conditions[theConstrainedNodes[i]]->type;
            if (cndType == DIRICHLET_X || cndType == DIRICHLET_Y) {
                int shift = (cndType ==DIRICHLET_X) ? 0 : 1;
                int node = (i - shift) / 2;
                int index = theMesh->number[node] * 2 + shift;  // renumérotation
                femFullSystemConstrain(theSystem, index,value); }
            }
        }
}

double* femBandSystemEliminate(femBandSystem* mySystem)
{
    double **A = mySystem->A;
    double *B = mySystem->B;
    int size = mySystem->size;
    int band = mySystem->band;

    for (int k = 0; k < size; k++) {
        if (fabs(A[k][k]) <= 1e-16) {
            printf("Pivot index %d  ", k);
            printf("Pivot value %e  ", A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }
        int jmax = fmin(size, k + band / 2 + 1);
        for (int i = k + 1; i < jmax; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k + 1; j < jmax; j++) {
                A[i][j] -= A[k][j] * factor;
            }
            B[i] -= B[k] * factor;
        }
    }

    for (int i = size - 1; i >= 0; i--) {
        double factor = 0.0;
        int jmax = fmin(size, i + band / 2 + 1);
        for (int j = i + 1; j < jmax; j++) {
            factor += A[i][j] * B[j];
        }
        B[i] = (B[i] - factor) / A[i][i];
    }

    return mySystem->B;
}

void femBandSystemInit(femBandSystem *myBandSystem)
{
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    memset(myBandSystem->B, 0, sizeof(double) * size * (band + 1));
}

femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    if (myBandSystem == NULL) {
        Error("Memory allocation failed for femBandSystem");
    }

    myBandSystem->B = malloc(sizeof(double) * size * (band + 1));
    myBandSystem->A = malloc(sizeof(double *) * size);
    if (myBandSystem->B == NULL || myBandSystem->A == NULL) {
        Error("Memory allocation failed for femBandSystem arrays");
    }

    myBandSystem->size = size;
    myBandSystem->band = band;

    myBandSystem->A[0] = myBandSystem->B + size;
    for (int i = 1; i < size; i++) {
        myBandSystem->A[i] = myBandSystem->A[i - 1] + band - 1;
    }

    femBandSystemInit(myBandSystem);
    return myBandSystem;
}

double *femElasticitySolve(femProblem *theProblem)
{
    femFullSystem *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femGeo *theGeometry = theProblem->geometry;
    femNodes *theNodes = theGeometry->theNodes;
    femMesh *theMesh = theGeometry->theElements;

    femFullSystemInit(system);

    // femElasticityAssembleElements(theProblem);
    // femElasticityAssembleNeumann(theProblem);

    int nLocal = theMesh->nLocalNode;
    double x[nLocal], y[nLocal], phi[nLocal], dphidxsi[nLocal], dphideta[nLocal], dphidx[nLocal], dphidy[nLocal];
    int map[nLocal], mapX[nLocal], mapY[nLocal];

    double a = theProblem->A, b = theProblem->B, c = theProblem->C;
    double rho = theProblem->rho, g = theProblem->g;
    double **A = theSystem->A;
    double *B = theSystem->B;

    femMeshRenumber(theMesh, theProblem->renumType);
    int myBand = femMeshComputeBand(theMesh);

    for (int iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (int j = 0; j < nLocal; j++) {
            map[j] = theMesh->elem[iElem * nLocal + j];
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
            map[j] = theMesh->number[map[j]];  
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
        }

        for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0, xLoc = 0.0;
            for (int i = 0; i < theSpace->n; i++) {
                xLoc += x[i] * phi[i];
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            for (int i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            if (jac <= 0) {
                printf("Error: Jacobian is non-positive for element %d.\n", iElem);
                return NULL;
            }

            for (int i = 0; i < theSpace->n; i++) {
                for (int j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight;
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight;
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight;
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight;
                }
                B[mapY[i]] -= phi[i] * g * rho * jac * weight;
            }
        }
    }

    femApplyBoundaryConditions(theProblem);

    double *sol = malloc(sizeof(double) * theSystem->size);
    if (theProblem->solver == SOLVEUR_PLEIN) {
        B = femFullSystemEliminate(theSystem);
        memcpy(sol, B, sizeof(double) * theSystem->size);
    } else if (theProblem->solver == SOLVEUR_BANDE) {
        myBand = myBand * 2 + 1;
        femBandSystem *theBandSystem = femBandSystemCreate(theSystem->size, myBand);
        for (int i = 0; i < theSystem->size; i++) {
            int jmin = fmax(0, i - myBand / 2);
            int jmax = fmin(theSystem->size, i + myBand / 2 + 1);
            for (int j = jmin; j < jmax; j++) {
                theBandSystem->A[i][j] = A[i][j];
            }
            theBandSystem->B[i] = B[i];
        }
        theBandSystem->B = femBandSystemEliminate(theBandSystem);
        memcpy(sol, theBandSystem->B, sizeof(double) * theSystem->size);
        free(theBandSystem->A);
        free(theBandSystem->B);
        free(theBandSystem);
    } else if (theProblem->solver == GRADIENTS_CONJUGUES) {
        conjugateGradient(A, B, sol, theSystem->size);
    }

    for (int i = 0; i < theMesh->nodes->nNodes; i++) {
        B[2 * i] = sol[2 * theMesh->number[i]];
        B[2 * i + 1] = sol[2 * theMesh->number[i] + 1];
    }
    free(sol);

    return B;
}

double *femElasticityForces(femProblem *theProblem)
{
    int systemSize = theProblem->system->size;
    double *solution = theProblem->soluce;
    double *residuals = theProblem->residuals;

    A_copy = malloc(systemSize * sizeof(double *));
    B_copy = malloc(systemSize * sizeof(double));
    if (A_copy == NULL || B_copy == NULL) {
        Error("Memory allocation failed for A_copy or B_copy");
    }

    memcpy(B_copy, theProblem->system->B, sizeof(double) * systemSize);
    for (int i = 0; i < systemSize; i++) {
        B_copy[i] = -B_copy[i];
    }

    memcpy(A_copy, theProblem->system->A, sizeof(double *) * systemSize);
    for (int i = 0; i < systemSize; i++) {
        A_copy[i] = malloc(systemSize * sizeof(double));
        if (A_copy[i] == NULL) {
            Error("Memory allocation failed for A_copy row");
        }
        for (int j = 0; j < systemSize; j++) {
            A_copy[i][j] = -theProblem->system->A[i][j];
        }
    }

    if (residuals == NULL) {
        residuals = malloc(systemSize * sizeof(double));
        if (residuals == NULL) {
            Error("Memory allocation failed for residuals");
        }
        theProblem->residuals = residuals;
    }

    for (int i = 0; i < systemSize; i++) {
        residuals[i] = 0.0;
        for (int j = 0; j < systemSize; j++) {
            residuals[i] += A_copy[i][j] * solution[j];
        }
        residuals[i] -= B_copy[i];
    }

    for (int i = 0; i < systemSize; i++) {
        free(A_copy[i]);
    }
    free(A_copy);
    free(B_copy);

    return residuals;
}
