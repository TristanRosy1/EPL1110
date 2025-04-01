#include "fem.h"
#include <stdlib.h>
#include <stddef.h>
#include <stdbool.h>

int *DEGREE;
// Il faut un fifrelin generaliser ce code.....
//  (1) Ajouter l'axisymétrique !    (mandatory)
//  (2) Ajouter les conditions de Neumann !   (mandatory)
//  (3) Ajouter les conditions en normal et tangentiel !   (strongly advised)
//  (4) Et remplacer le solveur plein par un truc un fifrelin plus subtil  (mandatory)

void femElasticityAssembleElements(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
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


void femElasticityAssembleNeumann(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femIntegration *theRule = theProblem->ruleEdge;
  femDiscrete *theSpace = theProblem->spaceEdge;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;
  femMesh *theEdges = theGeometry->theEdges;
  double x[2], y[2], phi[2];
  int iBnd, iElem, iInteg, iEdge, i, j, d, map[2];
  int nLocal = 2;
  double *B = theSystem->B;
  double pos;
  double y_moy;
  double hauteur;

  for (iBnd = 0; iBnd < theProblem->nBoundaryConditions; iBnd++) {
    femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
    femBoundaryType type = theCondition->type;
    double value = theCondition->value1;

    if(type != NEUMANN_X && type != NEUMANN_Y && type != NEUMANN_N && type != NEUMANN_T && type != NEUMANN_HYDROSTAT){
      continue;
    }
    

////// Maintenant on applique la condition de Neumann sur chaque morceaux de chaque domaine 
    int iEl;//utiliser d'autres variables pour ne pas ecraser les variables de boucle
    double hauteur;// calcule de la hauteur du domaine en question
    double y_max = -1;//
    double y_min = -1;//
    double y_moy;//distance au centre d'une arete depuis l'origine
    double POS; // position correcte avec 0 au dessus du barrage, donc négative. Correspond à y dans rho*g*y
    //double length; //longueur de l'arete
    int number_elem = theCondition->domain->nElem;
    for (int iEd = 0; iEd < theCondition->domain->nElem; iEd++) {
      iEl = theCondition->domain->elem[iEd];
      for (int k = 0; k < nLocal; k++) {
        map[k] = theEdges->elem[iEl * nLocal + k];
        y[k] = theNodes->Y[map[k]];
        if(y_max == -1 && y_min == -1){
          y_max = y[k];
          y_min = y[k];
        }
        if(y[k] > y_max){
          y_max = y[k];
        }
        if(y[k] < y_min){
          y_min = y[k];
        }
      }
    }

    hauteur = y_max - y_min;//
////////////////////////////////////////////
////// Maintenant on applique la condition de Neumann sur chaque morceaux de chaque domaine 

    for (iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++) {
      iElem = theCondition->domain->elem[iEdge];
      for (j = 0; j < nLocal; j++) {
        map[j] = theEdges->elem[iElem * nLocal + j];
        x[j] = theNodes->X[map[j]];
        y[j] = theNodes->Y[map[j]];

      }

      double y_moy = (y[0]+y[1])/2; //
      POS = y_moy - hauteur;//
      //length = fabs(y[1] - y[0]);

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
        f_x = value * tx/length;
        f_y = value * ty/length;
      }
      
      if (type == NEUMANN_HYDROSTAT){//
        f_x = value * POS;//
        
      }//

      double r = 0.0;
      for (iInteg = 0; iInteg < theRule->n; iInteg++) {
        double xsi = theRule->xsi[iInteg];
        double weight = theRule->weight[iInteg];
        femDiscretePhi(theSpace, xsi, phi);
        for(i = 0; i < theSpace->n; i++){
          r = x[i]*phi[i];
          }
        if (theProblem->planarStrainStress == AXISYM){
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x*r;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y*r;
          }
        }else{
          for (i = 0; i < theSpace->n; i++) {
            B[2*map[i] + 0] += jac * weight * phi[i] * f_x;
            B[2*map[i] + 1] += jac * weight * phi[i] * f_y;
        }

        }
        }
      }
    }
  }


void femElasticityApplyDirichlet(femProblem *theProblem) {
  femFullSystem *theSystem = theProblem->system;
  femGeo *theGeometry = theProblem->geometry;
  femNodes *theNodes = theGeometry->theNodes;

  for (int node = 0; node < theNodes->nNodes; node++) {
    femConstrainedNode *theConstrainedNode = &theProblem->constrainedNodes[node];
    if (theConstrainedNode->type == UNDEFINED)
      continue;
    femBoundaryType type = theConstrainedNode->type;

    if (type == DIRICHLET_X) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 0, value);
    }
    if (type == DIRICHLET_Y) {
      double value = theConstrainedNode->value1;
      femFullSystemConstrain(theSystem, 2 * node + 1, value);
    }
    if (type == DIRICHLET_XY) {
      double value_x = theConstrainedNode->value1;
      double value_y = theConstrainedNode->value2;
      femFullSystemConstrain(theSystem, 2 * node + 0, value_x);
      femFullSystemConstrain(theSystem, 2 * node + 1, value_y);
    }


    if (type == DIRICHLET_N) {
      
      double value = theConstrainedNode->value1;
      double nx = theConstrainedNode->nx;
      double ny = theConstrainedNode->ny;
      double tx = -ny;
      double ty = nx;

      double **a, *b;
      int i, size, band;
      double cx, cy, lx, ly; //attention si t'as mis int à la place de double ici tu perds instantanément 3h de ta vie à débeug

      a    = theSystem->A;
      b    = theSystem->B;
      size = theSystem->size;
      
      double att =  (tx * (tx * a[2*node][2*node] + ty * a[2*node+1][2*node]) + ty * (tx * a[2*node][2*node+1] + ty * a[2*node+1][2*node+1]));
      double atn =  (nx * (tx * a[2*node][2*node] + ty * a[2*node+1][2*node]) + ny * (tx * a[2*node][2*node+1] + ty * a[2*node+1][2*node+1]));
      double bt = tx * b[2*node] + ty * b[2*node+1];

      for(int i = 0; i < size; i++) { //on mets les autres lignes et collones à jour
        
        cx=a[i][2*node];
        cy=a[i][2*node+1];
        lx=a[2*node][i];
        ly=a[2*node+1][i];

        a[i][2*node] = tx*(cx*tx + cy*ty);
        a[i][2*node+1] = ty*(cx*tx + cy*ty);
        a[2*node][i] = tx*(lx*tx + ly*ty);
        a[2*node+1][i] = ty*(lx*tx + ly*ty);
        b[i] -= value*(nx*cx + ny*cy);
      }

      // On effectue les opérations sur la petite matrice
      a[2*node][2*node] = nx*nx + tx*tx*att;
      a[2*node][2*node+1] = nx*ny + tx*ty*att;
      a[2*node+1][2*node] = nx*ny + ty*tx*att;
      a[2*node+1][2*node+1] = ny*ny + ty*ty*att;

      b[2*node] = value*nx + tx*(bt - value*atn);
      b[2*node+1] = value*ny + ty*(bt - value*atn); 
    }



      
      if (type == DIRICHLET_T) {

        double value = theConstrainedNode->value1;
        double nx = theConstrainedNode->nx;
        double ny = theConstrainedNode->ny;
        double tx = ny;
        double ty = -nx;

        double **a, *b;
        int i, size, band;
        double cx, cy, lx, ly;

        a    = theSystem->A;
        b    = theSystem->B;
        size = theSystem->size;
        
        double att =  (nx * (nx * a[2*node][2*node] + ny * a[2*node+1][2*node]) + ny * (nx * a[2*node][2*node+1] + ny * a[2*node+1][2*node+1]));
        double atn =  (tx * (nx * a[2*node][2*node] + ny * a[2*node+1][2*node]) + ty * (nx * a[2*node][2*node+1] + ny * a[2*node+1][2*node+1]));
        double bt = nx * b[2*node] + ny * b[2*node+1];

        for(int i = 0; i < size; i++) { //on mets les autres lignes et collones à jour
          
          cx=a[i][2*node];
          cy=a[i][2*node+1];
          lx=a[2*node][i];
          ly=a[2*node+1][i];

          a[i][2*node] = nx*(cx*nx + cy*ny);
          a[i][2*node+1] = ny*(cx*nx + cy*ny);
          a[2*node][i] = nx*(lx*nx + ly*ny);
          a[2*node+1][i] = ny*(lx*nx + ly*ny);
          b[i] -= value*(tx*cx + ty*cy);
        }

        // On effectue les opérations sur la petite matrice 
        a[2*node][2*node] = tx*tx + nx*nx*att;
        a[2*node][2*node+1] = tx*ty + nx*ny*att;
        a[2*node+1][2*node] = tx*ty + ny*nx*att;
        a[2*node+1][2*node+1] = ty*ty + ny*ny*att;

        b[2*node] = value*tx + nx*(bt - value*atn);
        b[2*node+1] = value*ty + ny*(bt - value*atn);
      }


      if (type == DIRICHLET_NT) {
        double value_n = theConstrainedNode->value1;
        double value_t = theConstrainedNode->value2;
        double nx = theConstrainedNode->nx;
        double ny = theConstrainedNode->ny;
  
        double tx = -ny;
        double ty = nx;
        femFullSystemConstrain(theSystem, 2 * node + 0, value_n * nx + value_t * tx);
        femFullSystemConstrain(theSystem, 2 * node + 1, value_n * ny + value_t * ty);
      }
  }
}


int matrixComputeBand(double **A, int size) {
    int myBand = 0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            if (A[i][j] != 0 && abs(i-j) > myBand){
                myBand = abs(i-j);
            }
        }
    }
    myBand++;
    return myBand;
}



  void *matrixLocalCreate(femFullSystem * theSystem, int size, double **A, double *B, int *r){
    double **Arenum = malloc(sizeof(double*) * size);
    for(int i = 0; i<size ;++i){
      Arenum[i] = malloc(sizeof(double) * size);
    }
    for(int i = 0; i<size ;++i){
      for(int j = 0; j<size ;++j){ 
        Arenum[i][j] = A[r[i]] [r[j]];
      }
    }
    double *Brenum = malloc(sizeof(double) * size);
    for(int i = 0; i<size ;++i){
      Brenum[i] = B[r[i]];
    }
   
    theSystem->A = Arenum;
    theSystem->B = Brenum;
  }


  double* InverseMatrix(double* B, int* r, int size){ 
    int *invX = malloc(sizeof(int) * size);
    for(int i = 0; i<size; ++i){
      invX[r[i]] = i;
    }
    double* Xsol = malloc(sizeof(double) * size);
    for(int i = 0; i<size; ++i){
      Xsol[i] = B[invX[i]];
    }
    free(invX);
    return Xsol;
  }

 


double *femBandSystemEliminate(femProblem *theProblem) {
    femGeo *theGeometry = theProblem->geometry;
    femFullSystem *theSystem = theProblem->system;
    double** A = theSystem->A;
    double * B = theSystem->B;
    int size = theProblem->system->size;
    int i,j,k,jend;
    double factor;
    int bandA, band;
    
    bandA = matrixComputeBand(A, size);
    int* renumber = RenumberCuthill(A, size);

    matrixLocalCreate(theSystem, size, A, B, renumber);

    A=theSystem->A;
    B=theSystem->B;
 
  
    band = matrixComputeBand(A, size);

    printf("Bande renumerotee: %d\n", band);
    
    for ( k = 0; k < size; ++k){
        if (fabs(A[k][k]) <= 1e-8){ Error("Cannot eliminate with such a pivot.");}
        jend = (k+band < size) ? k+band : size;
        for ( i = k+1; i < jend; ++i) {
            factor = A[k][i] / A[k][k];
            for ( j = i; j < jend; ++j) {
                A[i][j] -= A[k][j] * factor;
            }
            B[i] -= B[k] * factor;
        }
    }

    for (i = size-1; i >= 0; --i) {
        factor = 0;
        jend = (i+band < size) ? i+band : size;
        for ( j = i+1; j < jend; ++j) {
            factor += A[i][j] * B[j];
        }
        B[i] = (B[i] - factor)/A[i][i];
    }
    
    double *X = malloc(sizeof(double) * size);
    
    X = InverseMatrix(B, renumber, size);
    free(renumber);
    return (X);
}

  int compare(const void *a, const void *b){

    int *x = *(int**)a;
    int *y = *(int**)b;

    if(x[0] == y[0])
      return (x[1] - y[1]);
    return (x[0] - y[0]);
  }

  int compare2(const void *a, const void *b){
      int x = *(int*)a;
      int y = *(int*)b;
  return (DEGREE[x] - DEGREE[y]);
  }

int* RenumberCuthill(double** matrice, int matrice_size){

  //liste permettant de connaittre les voisins de chaque noeud de la matrice
  //list[i][j] nous donne le j-ème voisin du noeud i
  //on parle bien ici des points de la matrice et non des noeuds du maillage
  int **list = malloc(sizeof(int*) * matrice_size); 
  
  //liste des degrés de chaque noeud
  int *degrees = malloc(sizeof(int) * matrice_size);

  for(int i = 0; i < matrice_size; i++){
    list[i] = malloc(sizeof(int)*50);
    degrees[i] = 0;
    for(int j = 0; j < matrice_size; j++){
      if(matrice[i][j] != 0 && i != j){
        list[i][degrees[i]++] = j;
        if(degrees[i]>=50){printf("Attention, le nombre de noeuds est trop grand pour la liste du noeud : %d\n",i);}
      }
    }
  }




  
  int *q = malloc(sizeof(int)*matrice_size*matrice_size);
  int *r = malloc(sizeof(int)*matrice_size);
  int next_empty_q = 0;
  int first_filled_q = 0;
  int next_empty_r = 0;
  bool *not_visited = malloc(sizeof(int)*matrice_size);

//S1:  Nous trouvons d'abord l'objet de degré minimum dont l'index n'a pas encore été ajouté à R.  
//     Disons que l'objet correspondant à la pième ligne a été identifié comme l'objet ayant un degré minimum. On ajoute p à R.

  for (int i = 0; i < matrice_size; i++) {
    not_visited[i] = true;
  }

  int actual_degree = degrees[0];
  int minimum_node = 0;
  bool change = true;

    for (int i = 0; i< matrice_size; i++) {
      if(degrees[i] < actual_degree && not_visited[i]){
        actual_degree = degrees[i];
        minimum_node = i;
      }
    }

  while(change){

      r[next_empty_r++] = minimum_node;
      not_visited[minimum_node] = false;


      // S2 : Quand un indice est ajouté à R, ajoutez tous les voisins de l'objet correspondant à l'indice,
      //      dans l'ordre croissant du degré, à Q. 
      
      for (int k = 0; k < degrees[minimum_node]; k++) {
        q[next_empty_q++] = list[minimum_node][k];
      }
      
      DEGREE = degrees;
      qsort(&q[first_filled_q], next_empty_q-first_filled_q, sizeof(int), compare2);


      // S3 : Extraire le premier nœud de Q, appelé C. Si C n'a pas été inséré dans R, l'ajouter à R, ajouter à Q les voisins
      // de C dans l'ordre croissant du degré.
      // S4 : Si Q n'est pas vide, répéter S3.
      int C;
        while(first_filled_q < next_empty_q){
          C = q[first_filled_q++];
          if(not_visited[C]){
            r[next_empty_r++] = C;
            not_visited[C] = false;
            int scission = next_empty_q;
            for (int j = 0; j < degrees[C]; j++) {
              q[next_empty_q++] = list[C][j];
            }
            qsort(&q[scission], next_empty_q-scission, sizeof(int), compare2);//il faut sort à partir de scission à mon avis ce que j'ai mis est faux 
        }
      
    }
      // S5 : Si Q est vide, mais qu'il reste des objets dans la matrice qui n'ont pas été inclus dans R, recommencez à partir de S1 une fois de plus. (Cela pourrait se produire s'il y a des graphes disjoints)
      // S6 : Terminez cet algorithme une fois que tous les objets sont inclus dans R.
    change = false;
    actual_degree = matrice_size+1;
    for (int i = 0; i< matrice_size; i++) {
      if(degrees[i] < actual_degree && not_visited[i]){
        actual_degree = degrees[i];
        minimum_node = i;
        change = true;
        }
    }
  }


    // S7 : Enfin, inverser les indices dans R, c'est-à-dire (échanger R[i] avec R[P-i+1]).
  
    int temp;
    for (int i = 0; i < matrice_size / 2; i++) {
        temp = r[i];
        r[i] = r[matrice_size - i - 1];
        r[matrice_size - i - 1] = temp;
    }

 
    // section free (j'ai été un peu vite la dessus)   
    for (int i = 0; i < matrice_size; i++) {
        free(list[i]);
    }
    free(list);
    free(degrees);
    free(q);
    free(not_visited);
    return r; //il faut free r après l'avoir utilisé 
}

  double *femElasticitySolve(femProblem *theProblem) {
  
  femElasticityAssembleElements(theProblem);

  femElasticityAssembleNeumann(theProblem);

  femElasticityApplyDirichlet(theProblem);



  double *soluce = femBandSystemEliminate(theProblem);

  memcpy(theProblem->soluce, soluce, theProblem->system->size * sizeof(double));
  return theProblem->soluce;
}