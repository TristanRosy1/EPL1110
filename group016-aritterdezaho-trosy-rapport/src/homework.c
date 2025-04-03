#include "fem.h"

double* VAL;


// Strip : BEGIN
double **A_copy = NULL;
double *B_copy  = NULL;
// Strip : END

int compare(const void *a, const void *b) {
    if (VAL[*(int*)a] < VAL[*(int*)b]) return 1;
    if (VAL[*(int*)a] > VAL[*(int*)b]) return -1;
    return 0;
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i;
    int *tab = malloc(sizeof(int) * theMesh->nodes->nNodes);

    for (i = 0; i < theMesh->nodes->nNodes; i++) 
        tab[i] = i;

    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM : 
            VAL = theMesh->nodes->X;
            qsort(tab, theMesh->nodes->nNodes, sizeof(int), compare);
            break;
        case FEM_YNUM : 
            VAL = theMesh->nodes->Y; 
            qsort(tab, theMesh->nodes->nNodes, sizeof(int), compare);
            break;            
        default : Error("Unexpected renumbering option"); }
    
    for (i = 0; i < theMesh->nodes->nNodes; i++) {
        theMesh->number[tab[i]] = i;
    }
    free(tab);
}

int femMeshComputeBand(femMesh *theMesh)
{   
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theMesh->number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; 
    }        
    myBand *=2 ;
    return(++myBand); // Pour chaque noeud, on a une composante U et une composante V
}

void conjugateGradient(double** A, double *b, double *x, int size) {
    double *r = malloc(size*sizeof(double));
    double *d = malloc(size*sizeof(double));
    double *s = malloc(size*sizeof(double));
    double alpha, beta, delta, deltaNew;
    int nMax = 2000;

    // Initialisation
    for (int i = 0; i < size; i++) {
        r[i] = b[i];
        d[i] = r[i];
        x[i] = 0.0;
    }

    delta = 0.0;
    for (int i = 0; i < size; i++)
        delta += r[i] * r[i];

    int k = 0;
    while (delta > 1e-8 & k < nMax) {  // Critère de convergence
        // Calcul de s
        for (int i = 0; i < size; i++) {
            s[i] = 0.0;
            for (int j = 0; j < size; j++)
                s[i] += A[i][j] * d[j];
        }

        // Calcul de alpha
        double sd = 0.0;
        for (int i = 0; i < size; i++)
            sd += d[i] * s[i];

        alpha = delta / sd;

        // Mise à jour de x et r
        for (int i = 0; i < size; i++) {
            x[i] += alpha * d[i];
            r[i] -= alpha * s[i];
        }

        deltaNew = 0.0;
        for (int i = 0; i < size; i++)
            deltaNew += r[i] * r[i];

        // Calcul de beta
        beta = deltaNew / delta;

        // Mise à jour de d
        for (int i = 0; i < size; i++)
            d[i] = r[i] + beta * d[i];

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
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->rule;
    femDiscrete    *theSpace    = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theMesh     = theGeometry->theElements;
    femMesh        *theEdges    = theGeometry->theEdges;

    double x[4], y[4], phi[4], dphidxsi[4] ,dphideta[4], dphidx[4], dphidy[4];
    int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++)
    {
        for (j = 0; j < nLocal; j++)
        {
            map[j]  = theMesh->elem[iElem * nLocal + j];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];
        } 
        
        for (iInteg = 0; iInteg < theRule->n; iInteg++)
        {   
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];

            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            
            double dxdxsi = 0.0; double dxdeta = 0.0;
            double dydxsi = 0.0; double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++)
            {  
                dxdxsi += x[i] * dphidxsi[i];       
                dxdeta += x[i] * dphideta[i];   
                dydxsi += y[i] * dphidxsi[i];   
                dydeta += y[i] * dphideta[i];
            }

            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++)
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }

            double weightedJac = jac * weight;

            for (i = 0; i < theSpace->n; i++)
            { 
                for (j = 0; j < theSpace->n; j++)
                {
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
    femFullSystem  *theSystem   = theProblem->system;
    femIntegration *theRule     = theProblem->ruleEdge;
    femDiscrete    *theSpace    = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes    = theGeometry->theNodes;
    femMesh        *theEdges    = theGeometry->theEdges;

    double x[2], y[2], phi[2];
    int iBnd, iElem, iInteg, iEdge, i, j, d, map[2], mapU[2];
    
    int nLocal = 2;
    double *B  = theSystem->B;

    for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++)
    {
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        femDomain *domain = theCondition->domain;
        double value = theCondition->value;

        if (type == DIRICHLET_X || type == DIRICHLET_Y) { continue; }

        int shift = (type == NEUMANN_X) ? 0 : 1;
        
        for (iEdge = 0; iEdge < domain->nElem; iEdge++)
        {
            iElem = domain->elem[iEdge];

            for (j = 0; j < nLocal; j++)
            {
                map[j] = theEdges->elem[iElem * nLocal + j];
                mapU[j] = 2 * map[j] + shift;
                x[j] = theNodes->X[map[j]];
                y[j] = theNodes->Y[map[j]];
            }
            
            double dx = x[1] - x[0];
            double dy = y[1] - y[0];
            double length = sqrt(dx * dx + dy * dy);
            double jac = length / 2;

            for (iInteg = 0; iInteg < theRule->n; iInteg++)
            {
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace, xsi, phi);

                for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value * jac * weight; }
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem)
{
    femFullSystem *system = theProblem->system;
    femFullSystemInit(system);

    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    int systemSize = system->size;

    if (A_copy == NULL) {
        A_copy = (double **) malloc(sizeof(double *) * systemSize);
        for (int row = 0; row < systemSize; row++) {
            A_copy[row] = (double *) malloc(sizeof(double) * systemSize);
        }
    }
    if (B_copy == NULL) {
        B_copy = (double *) malloc(sizeof(double) * systemSize);
    }

    for (int row = 0; row < systemSize; row++) {
        for (int col = 0; col < systemSize; col++) {
            A_copy[row][col] = system->A[row][col];
        }
        B_copy[row] = system->B[row];
    }

    int *constrainedNodes = theProblem->constrainedNodes;
    for (int node = 0; node < systemSize; node++) {
        if (constrainedNodes[node] != -1) {
            double constraintValue = theProblem->conditions[constrainedNodes[node]]->value;
            femFullSystemConstrain(system, node, constraintValue);
        }
    }

    femFullSystemEliminate(system);
    memcpy(theProblem->soluce, system->B, systemSize * sizeof(double));
    return theProblem->soluce;
}

double *femElasticityForces(femProblem *theProblem)
{
    double *residuals = theProblem->residuals;
    double *solution = theProblem->soluce;
    int systemSize = theProblem->system->size;

    if (residuals == NULL) {
        residuals = (double *) malloc(sizeof(double) * systemSize);
    }

    for (int i = 0; i < systemSize; i++) {
        residuals[i] = 0.0;
    }

    for (int row = 0; row < systemSize; row++) {
        for (int col = 0; col < systemSize; col++) {
            residuals[row] += A_copy[row][col] * solution[col];
        }
        residuals[row] -= B_copy[row];
    }

    for (int row = 0; row < systemSize; row++) {
        free(A_copy[row]);
        A_copy[row] = NULL;
    }
    free(A_copy);
    free(B_copy);
    A_copy = NULL;
    B_copy = NULL;

    return residuals;
}
