
#include "fem.h"
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Assemble les matrices A et B propres aà notre probleme et à nos éléments.
 * 
 * Cette fonction prend notre problème et génèrent les matrices associées.
 * 
 * @param theProblem le problème à résoudre.
 * @return void.
 */

void femElasticityAssembleElements(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theMesh = theGeometry->theElements;
  double x[4], y[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
  int iElem, iInteg, iEdge, i, j, d, map[4], mapX[4], mapY[4];
  int nLocal = theMesh->nLocalNode;
  //données de notre problème
  double a = theProblem->A;
  double b = theProblem->B;
  double c = theProblem->C;
  double rho = theProblem->rho;
  double gx = theProblem->gx;
  double gy = theProblem->gy;
  double **A = theSystem->A;
  double *B = theSystem->B;

  for (iElem = 0; iElem < theMesh->nElem; iElem++) {
    for (j = 0; j < nLocal; j++) {
      map[j] = theMesh->elem[iElem * nLocal + j];
      mapX[j] = 2 * map[j];
      mapY[j] = 2 * map[j] + 1;
      x[j] = theNodes->X[map[j]];
      y[j] = theNodes->Y[map[j]];
    }

    for (iInteg = 0; iInteg < theRule->n; iInteg++) {
      double xsi = theRule->xsi[iInteg];
      double eta = theRule->eta[iInteg];
      double weight = theRule->weight[iInteg];
      femDiscretePhi2(theSpace, xsi, eta, phi);
      femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

      double dxdxsi = 0.0;
      double dxdeta = 0.0;
      double dydxsi = 0.0;
      double dydeta = 0.0;
      //correspond au r du cas axisymétrique (coordonnées polaires)
      double r = 0.0;
      for (i = 0; i < theSpace->n; i++) {
        dxdxsi += x[i] * dphidxsi[i];
        dxdeta += x[i] * dphideta[i];
        dydxsi += y[i] * dphidxsi[i];
        dydeta += y[i] * dphideta[i];
        r += x[i]*phi[i];
      }
      double jac = dxdxsi * dydeta - dxdeta * dydxsi;
      if (jac < 0.0)
        printf("Negative jacobian! Your mesh is oriented in reverse. The normals will be wrong\n");
      jac = fabs(jac);

      for (i = 0; i < theSpace->n; i++) {
        dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
        dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
      }
      if (theProblem->planarStrainStress == AXISYM){
        for (i = 0; i < theSpace->n; i++) {
          for (j = 0; j < theSpace->n; j++) {
            // on a besoin de r pour le cas axisymétrique 
            A[mapX[i]][mapX[j]] += (dphidx[i] * a * r * dphidx[j] + dphidy[i] * c * r * dphidy[j] + phi[i] * ((b * dphidx[j]) + (a * phi[j] / r)) + dphidx[i] * b * phi[j]) * jac * weight ;
            A[mapX[i]][mapY[j]] += (dphidx[i] * b * r * dphidy[j] + dphidy[i] * c * r * dphidx[j] + phi[i] * b * dphidy[j]) * jac * weight ;
            A[mapY[i]][mapX[j]] += (dphidy[i] * b * r * dphidx[j] + dphidx[i] * c * r * dphidy[j] + dphidy[i] * b * phi[j]) * jac * weight ;
            A[mapY[i]][mapY[j]] += (dphidy[i] * a * r * dphidy[j] + dphidx[i] * c * r * dphidx[j]) * jac * weight ;
          }
        }
      }
      else{
        for (i = 0; i < theSpace->n; i++) {
          for (j = 0; j < theSpace->n; j++) { 
            A[mapX[i]][mapX[j]] +=  (dphidx[i] * a * dphidx[j] + dphidy[i] * c * dphidy[j]) * jac * weight ;
            A[mapX[i]][mapY[j]] +=  (dphidx[i] * b * dphidy[j] + dphidy[i] * c * dphidx[j]) * jac * weight ;
            A[mapY[i]][mapX[j]] +=  (dphidy[i] * b * dphidx[j] + dphidx[i] * c * dphidy[j]) * jac * weight ;
            A[mapY[i]][mapY[j]] +=  (dphidy[i] * a * dphidy[j] + dphidx[i] * c * dphidx[j]) * jac * weight ;
          }
        }
      }
      
      for (i = 0; i < theSpace->n; i++) {
        if(theProblem->planarStrainStress == AXISYM){
          // on a besoin de r pour le cas axisymétrique
          B[mapX[i]] += r* phi[i] * gx * rho * jac * weight;
          B[mapY[i]] += r *  phi[i] * gy * rho * jac * weight;
        }
        else{
          // on assemble la partie droite de l'équation
          B[mapX[i]] += phi[i] * gx * rho * jac * weight;
          B[mapY[i]] += phi[i] * gy * rho * jac * weight;
        }
      }
    }
  }
}
/**
 * @brief Impose nos conditions de Neumann et on modifie notre matrice B.
 * 
 * Cette fonction prend notre problème en entrée et pour chaque frontière impose les conditions de neumann demandées en modifiant la matrice B .
 * 
 * @param theProblem le problème à résoudre.
 * @return void.
 */
void femElasticityAssembleNeumann(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T && type != NEUMANN_WIND ){
      continue;
    }

    int POS ; // position du noeud par rapport à 0

    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];
      }
      // pour trouver une position moyenne, nous faisons donc une approximation
      POS = y[0]+y[1]/2;

      double tx = x[1] - x[0];
      double ty = y[1] - y[0];
      double length = hypot(tx, ty);
      double jac = length / 2.0;
      
      double f_x = 0.0;
      double f_y = 0.0;
      if (type == NEUMANN_X) {
        f_x = value;
      }
      if (type == NEUMANN_Y) {
        f_y = value;
      }
      if (type == NEUMANN_N) {
        double nx =  ty / length;
        double ny = -tx / length;
        f_x = value * nx;
        f_y = value * ny;
      }
      if (type == NEUMANN_T) {
        f_x = value * tx / length;
        f_y = value * ty / length;
      }

      if (type ==NEUMANN_WIND){
        // on veut imposer une force qui varie selon la hauteur de l'élément étudié.
        // cette force varie de manière sinusoidale avec la position de l'élément
        f_x = value * sin(POS)* POS;
        
      }



      
      

      // Attention, pour le normal tangent on calcule la normale (sortante) au SEGMENT, surtout PAS celle de constrainedNodes
      // Une petite aide pour le calcul de la normale :-)
      // double nx =  ty / length;
      // double ny = -tx / length;

      // on ne modifie que la partie B du système car elles représentent les forces nodales, les autres sont à 0 
      for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double weight = theRule->weight[iInteg];
        femDiscretePhi(theSpace, xsi, phi);
        if (theProblem->planarStrainStress==AXISYM){
          double dphidxsi[4], dphideta[4];
          double eta = theRule->eta[iInteg];
          femDiscretePhi2(theSpace, xsi, eta, phi);
          
          femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);

          double dxdxsi = 0.0;
          double dxdeta = 0.0;
          double dydxsi = 0.0;
          double dydeta = 0.0;
          //correspond au r du cas axisymétrique
          double r = 0.0;
          for (i = 0; i < theSpace->n; i++) {
            dxdxsi += x[i] * dphidxsi[i];
            dxdeta += x[i] * dphideta[i];
            dydxsi += y[i] * dphidxsi[i];
            dydeta += y[i] * dphideta[i];
            r += x[i]*phi[i];
          }
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x *r;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y *r;
          }

        }
        else {
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y;
          }
        }
      }
    }
  }
}

/**
 * @brief Impose nos conditions de Dirichket et on modifie nos matrices à l'aide de femSystemConstrain.
 * 
 * Cette fonction prend notre problème en entrée et pour chaque frontière impose les conditions de Dirichlet demandées en modifiant les matrices .
 * 
 * @param theProblem le problème à résoudre.
 * @return void.
 */
void femElasticityApplyDirichlet(femProblem *theProblem) {
  femBandSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      // on récupère la valeur de la contrainte
      double value = theConstrainedNode->value1;
      // on l'impose dans le système au noeud correspondant
      // on le fait avce une procédure d'assemblage de matrice
      femSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femSystemConstrain(theSystem, 2 * node + 0, value_x);
      femSystemConstrain(theSystem, 2 * node + 1, value_y);
    }

    if (type == DIRICHLET_N) {
      double value = theConstrainedNode->value1;
      // donc valeurs normalisées
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      femSystemConstrain(theSystem, 2 * node + 0, value*nx);
      femSystemConstrain(theSystem, 2 * node + 1, value*ny);
    }
    if (type == DIRICHLET_T) {
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double tx = ny ;
      double ty = -nx ;
      femSystemConstrain(theSystem, 2 * node + 0, value*tx);
      femSystemConstrain(theSystem, 2 * node + 1, value*ty);

    }
    if (type == DIRICHLET_NT) {
      double value_n = theConstrainedNode->value1;
      double value_t = theConstrainedNode->value2;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double tx = ny ;
      double ty = -nx ;
      
      femSystemConstrain(theSystem, 2 * node + 0, value_n*nx + value_t*tx);
      femSystemConstrain(theSystem, 2 * node + 1, value_n*ny + value_t*ty);
    }
  }
}
  
/**
 * @brief Résout un système linéaire.
 * 
 * Cette fonction prend une matrice  A, un vecteur B et une taille en entrée et résout le système linéaire .
 * 
 * @param A la matrice .
 * @param B le vecteur.
 * @param size la taille de la matrice.
 * @return Le vecteur solution.
 */
double  *femBandSystemEliminate(double **A, double *B, int size)
{   
    double factor;
    int     i, j, k, jend;
    
    // On récupère la largeur de bande avant la renumérotation
    int band = matrixComputeBand(A, size);
    printf("Bande initiale: %d\n", band);
    // On trouve le vecteur de renumeration avec reverse cuthillmackee
    int *permutation = reverseCuthillMackee(A, size);
    // On renumérote avec reverse cuthillmackee avec la permutation obtenue
    problem *newProblem = renumberMesh(A, B, permutation, size);
    // on trouve la largeur de bande après la renumérotation
    band = matrixComputeBand(newProblem->A, size);
    printf("Bande renumerotee: %d\n", band);
    
    A = newProblem->A;
    B = newProblem->B;
    // On va résoudre notre problème linéaire avec la méthode de Gauss
    // Gauss elimination
    for (k = 0; k < size; ++k) {
        if (fabs(A[k][k]) < 1e-8) Error("Pivot nul !");
        jend = (k+band < size) ? k+band : size;
        for (j = k+1; j < jend; ++j) {
            factor = A[k][j]/A[k][k];
            for (i = j; i < jend; ++i) {
                A[j][i] -= factor * A[k][i];
            }
            B[j] -= factor * B[k];
        }
    }
    // Back-substitution
    for (k = size-1; k >= 0; --k) {
        factor = 0;
        jend = (k+band < size) ? k+band : size;
        for (j = k+1; j < jend; ++j) {
            factor += A[k][j]*B[j];
        }
        B[k] = (B[k] - factor)/A[k][k];
    }
    // on fait les permutations inverses pour retrouver les valeurs dans notre matrice initiale
    double *sol = malloc(size*sizeof(double));
    sol = inversePermutation(B, permutation, size);
    free(permutation);
    freeProblem(newProblem, size);
    return(sol);
}

/**
 * @brief Résout un problème d'élasticité.
 * 
 * Cette fonction prend un problème en entrée et résout le problème d'élasticité.
 * 
 * @param theProblem le problème à résoudre.
 * @return Le vecteur solution.
 */
double *femElasticitySolve(femProblem *theProblem) {
  // on assemble les éléments de notre problème et les conitions frontières imposées
  femElasticityAssembleElements(theProblem);
  femElasticityAssembleNeumann(theProblem);
  femElasticityApplyDirichlet(theProblem);
  // on rédoud notre système linéaire après une renumérotation des noeuds avec reverse cuthillmackee
  double *soluce = femBandSystemEliminate(theProblem->system->A, theProblem->system->B, theProblem->system->size) ;
  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}
