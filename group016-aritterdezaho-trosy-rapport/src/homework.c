#include "fem.h"

// Strip : BEGIN
double **A_copy = NULL;
double *B_copy  = NULL;
// Strip : END

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
    // Résumé de la fonction
    // But : Intégrer les conditions de Neumann dans le vecteur des forces ( B ).
    // Processus :
    // 1: Parcourt toutes les conditions de frontières.
    // 2: Ignore les conditions de Dirichlet (pas pertinentes ici).
    // 3: Pour chaque condition de Neumann : 1) Parcourt les arêtes du domaine.
    //                                       2)Calcule et intègre les contributions aux forces globales.
    // Sortie : Mise à jour du vecteur ( B ), prêt pour la résolution du système global.



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

                for (i = 0; i < theSpace->n; i++) { B[mapU[i]] += phi[i] * value * jac * weight; } //À la fin, toutes les contributions des conditions de Neumann sur chaque arête sont ajoutées au vecteur global ( B ) du système.
            }
        }
    }
}

double *femElasticitySolve(femProblem *theProblem) //Assemble et résout le système d'équations global pour obtenir la solution (déplacements ou forces).
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
            double constraintValue = theProblem->conditions[constrainedNodes[node]]->value; //contrainte de dirichlet sont applique ici 
            femFullSystemConstrain(system, node, constraintValue);
        }
    }

    femFullSystemEliminate(system);
    memcpy(theProblem->soluce, system->B, systemSize * sizeof(double));
    return theProblem->soluce;
}

double *femElasticityForces(femProblem *theProblem) //calcule des residus, Elle calcule la différence entre les forces internes et externes après résolution. Cela permet de vérifier l'équilibre des forces dans le système et de détecter d'éventuelles erreurs ou instabilités.
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
