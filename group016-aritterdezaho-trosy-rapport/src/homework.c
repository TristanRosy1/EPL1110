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

void femApplyBoundaryConditions(femProblem *theProblem) {
    femFullSystem *theSystem = theProblem->system;
    femMesh *theMesh = theProblem->geometry->theElements;
    double **A = theSystem->A;
    double *B = theSystem->B;
    int iCase = theProblem->planarStrainStress;

    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition *cnd = theProblem->conditions[i];
        femMesh *bndMesh = cnd->domain->mesh;
        double *X = bndMesh->nodes->X, *Y = bndMesh->nodes->Y;
        int *bndElem = cnd->domain->elem, nElem = cnd->domain->nElem;

        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T || 
            cnd->type == NEUMANN_N || cnd->type == NEUMANN_T) {

            if (!cnd->domain->n_t_malloced) {
                cnd->domain->n_t_malloced = 1;
                cnd->domain->normales = calloc(2 * (bndMesh->nElem + 1), sizeof(double));
                cnd->domain->tangentes = calloc(2 * (bndMesh->nElem + 1), sizeof(double));
                double *normales = cnd->domain->normales, *tangentes = cnd->domain->tangentes;

                for (int j = 0; j < nElem; j++) {
                    int node0 = bndMesh->elem[2 * bndElem[j]];
                    int node1 = bndMesh->elem[2 * bndElem[j] + 1];
                    double dx = X[node1] - X[node0], dy = Y[node1] - Y[node0];

                    tangentes[2 * j] += dx; tangentes[2 * j + 1] += dy;
                    normales[2 * j] -= dy; normales[2 * j + 1] += dx;
                    tangentes[2 * (j + 1)] += dx; tangentes[2 * (j + 1) + 1] += dy;
                    normales[2 * (j + 1)] -= dy; normales[2 * (j + 1) + 1] += dx;
                }

                for (int j = 0; j < nElem + 1; j++) {
                    double tNorm = sqrt(tangentes[2 * j] * tangentes[2 * j] + tangentes[2 * j + 1] * tangentes[2 * j + 1]);
                    double nNorm = sqrt(normales[2 * j] * normales[2 * j] + normales[2 * j + 1] * normales[2 * j + 1]);
                    tangentes[2 * j] /= tNorm; tangentes[2 * j + 1] /= tNorm;
                    normales[2 * j] /= nNorm; normales[2 * j + 1] /= nNorm;
                }
            }

            if ((cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) && !cnd->domain->n_t_matrix) {
                cnd->domain->n_t_matrix = 1;
                double *normales = cnd->domain->normales, *tangentes = cnd->domain->tangentes;

                for (int j = 0; j < nElem + 1; j++) {
                    int node0 = (j == nElem) ? bndMesh->elem[2 * bndElem[j - 1] + 1] : bndMesh->elem[2 * bndElem[j]];
                    node0 = theMesh->number[node0];
                    double nx = normales[2 * j], ny = normales[2 * j + 1];
                    double tx = tangentes[2 * j], ty = tangentes[2 * j + 1];

                    double B_U = B[2 * node0], B_V = B[2 * node0 + 1];
                    B[2 * node0] = nx * B_U + ny * B_V;
                    B[2 * node0 + 1] = tx * B_U + ty * B_V;

                    for (int k = 0; k < theSystem->size; k++) {
                        double A_U = A[2 * node0][k], A_V = A[2 * node0 + 1][k];
                        A[2 * node0][k] = nx * A_U + ny * A_V;
                        A[2 * node0 + 1][k] = tx * A_U + ty * A_V;

                        A_U = A[k][2 * node0]; A_V = A[k][2 * node0 + 1];
                        A[k][2 * node0] = nx * A_U + ny * A_V;
                        A[k][2 * node0 + 1] = tx * A_U + ty * A_V;
                    }
                }
            }
        }

        if (cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) {
            for (int j = 0; j < nElem + 1; j++) {
                int node0 = (j == nElem) ? bndMesh->elem[2 * bndElem[j - 1] + 1] : bndMesh->elem[2 * bndElem[j]];
                int shift = (cnd->type == DIRICHLET_N) ? 0 : 1;
                femFullSystemConstrain(theSystem, 2 * node0 + shift, cnd->value);
            }
        }

        if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y || cnd->type == NEUMANN_N || cnd->type == NEUMANN_T) {
            for (int j = 0; j < nElem; j++) {
                int node0 = bndMesh->elem[2 * bndElem[j]], node1 = bndMesh->elem[2 * bndElem[j] + 1];
                double jac = 0.5 * sqrt((X[node0] - X[node1]) * (X[node0] - X[node1]) + (Y[node0] - Y[node1]) * (Y[node0] - Y[node1]));

                if (cnd->type == NEUMANN_X || cnd->type == NEUMANN_Y) {
                    int shift = (cnd->type == NEUMANN_X) ? 0 : 1;
                    node0 = theMesh->number[node0];
                    node1 = theMesh->number[node1];
                    B[2 * node0 + shift] += jac * cnd->value;
                    B[2 * node1 + shift] += jac * cnd->value;
                } else {
                    double *n_or_t = (cnd->type == NEUMANN_N) ? cnd->domain->normales : cnd->domain->tangentes;
                    node0 = theMesh->number[node0];
                    node1 = theMesh->number[node1];
                    B[2 * node0] += jac * cnd->value * n_or_t[2 * j];
                    B[2 * node0 + 1] += jac * cnd->value * n_or_t[2 * j + 1];
                    B[2 * node1] += jac * cnd->value * n_or_t[2 * (j + 1)];
                    B[2 * node1 + 1] += jac * cnd->value * n_or_t[2 * (j + 1) + 1];
                }
            }
        }
    }

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femBoundaryType cndType = theProblem->conditions[theConstrainedNodes[i]]->type;
            int shift = (cndType == DIRICHLET_X) ? 0 : 1;
            int node = (i - shift) / 2;
            int index = theMesh->number[node] * 2 + shift;
            femFullSystemConstrain(theSystem, index, value);
        }
    }
}

void femDomChange(femProblem *theProblem, double *sol) {
    femFullSystem *theSystem = theProblem->system;
    femMesh *theMesh = theProblem->geometry->theElements;
    double *B = sol;

    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition *cnd = theProblem->conditions[i];
        if ((cnd->type == DIRICHLET_N || cnd->type == DIRICHLET_T) && cnd->domain->n_t_matrix == 1) {
            cnd->domain->n_t_malloced = 0;
            femMesh *bndMesh = cnd->domain->mesh;
            int *bndElem = cnd->domain->elem;
            double *normales = cnd->domain->normales, *tangentes = cnd->domain->tangentes;

            for (int j = 0; j <= cnd->domain->nElem; j++) {
                int node0 = (j == cnd->domain->nElem) 
                            ? bndMesh->elem[2 * bndElem[j - 1] + 1] 
                            : bndMesh->elem[2 * bndElem[j]];
                node0 = theMesh->number[node0];

                double nx = normales[2 * j], ny = normales[2 * j + 1];
                double tx = tangentes[2 * j], ty = tangentes[2 * j + 1];
                double B_U = B[2 * node0], B_V = B[2 * node0 + 1];

                B[2 * node0] = nx * B_U + ny * B_V;
                B[2 * node0 + 1] = tx * B_U + ty * B_V;
            }
        }
    }
}

void femFreeNT(femProblem *theProblem) {
    for (int i = 0; i < theProblem->nBoundaryConditions; i++) {
        femBoundaryCondition *cnd = theProblem->conditions[i];
        if (cnd->domain->n_t_malloced) {
            free(cnd->domain->normales);
            free(cnd->domain->tangentes);
            cnd->domain->n_t_malloced = 0; // Reset flag after freeing
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
            map[j] = theMesh->number[theMesh->elem[iElem * nLocal + j]];
            x[j] = theNodes->X[map[j]];
            y[j] = theNodes->Y[map[j]];
            mapX[j] = 2 * map[j];
            mapY[j] = 2 * map[j] + 1;
        }

        for (int iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg], eta = theRule->eta[iInteg], weight = theRule->weight[iInteg];
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
            if (jac <= 0) {
                printf("Error: Jacobian is non-positive for element %d.\n", iElem);
                return NULL;
            }

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

    femApplyBoundaryConditions(theProblem);

    double *sol = malloc(sizeof(double) * theSystem->size);
    if (theProblem->solver == SOLVEUR_PLEIN) {
        memcpy(sol, femFullSystemEliminate(theSystem), sizeof(double) * theSystem->size);
    } else if (theProblem->solver == SOLVEUR_BANDE) {
        myBand = myBand * 2 + 1;
        femBandSystem *theBandSystem = femBandSystemCreate(theSystem->size, myBand);
        for (int i = 0; i < theSystem->size; i++) {
            int jmin = fmax(0, i - myBand / 2), jmax = fmin(theSystem->size, i + myBand / 2 + 1);
            for (int j = jmin; j < jmax; j++) {
                theBandSystem->A[i][j] = A[i][j];
            }
            theBandSystem->B[i] = B[i];
        }
        memcpy(sol, femBandSystemEliminate(theBandSystem), sizeof(double) * theSystem->size);
        free(theBandSystem->A);
        free(theBandSystem->B);
        free(theBandSystem);
    } else if (theProblem->solver == GRADIENTS_CONJUGUES) {
        conjugateGradient(A, B, sol, theSystem->size);
    }

    femDomChange(theProblem, sol);
    femFreeNT(theProblem);

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

    // Allocate and initialize A_copy and B_copy
    double **A_copy = malloc(systemSize * sizeof(double *));
    double *B_copy = malloc(systemSize * sizeof(double));
    if (A_copy == NULL || B_copy == NULL) {
        Error("Memory allocation failed for A_copy or B_copy");
    }

    for (int i = 0; i < systemSize; i++) {
        B_copy[i] = -theProblem->system->B[i];
        A_copy[i] = malloc(systemSize * sizeof(double));
        if (A_copy[i] == NULL) {
            Error("Memory allocation failed for A_copy row");
        }
        for (int j = 0; j < systemSize; j++) {
            A_copy[i][j] = -theProblem->system->A[i][j];
        }
    }

    // Allocate residuals if not already allocated
    if (residuals == NULL) {
        residuals = malloc(systemSize * sizeof(double));
        if (residuals == NULL) {
            Error("Memory allocation failed for residuals");
        }
        theProblem->residuals = residuals;
    }

    // Compute residuals
    for (int i = 0; i < systemSize; i++) {
        residuals[i] = -B_copy[i];
        for (int j = 0; j < systemSize; j++) {
            residuals[i] += A_copy[i][j] * solution[j];
        }
    }

    // Free allocated memory
    for (int i = 0; i < systemSize; i++) {
        free(A_copy[i]);
    }
    free(A_copy);
    free(B_copy);

    return residuals;
}
