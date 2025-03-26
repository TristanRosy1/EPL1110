#include "fem.h"




double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xCenter;
    double y0 = theGeometry->yCenter;

    double rInner = theGeometry->rInner;
    double dInner = theGeometry->dInner;
    double hInner = theGeometry->hInner;

    double rOuter = theGeometry->rOuter;
    double dOuter = theGeometry->dOuter;
    double hOuter = theGeometry->hOuter;

    double hfinal = h;

    // Distance de la frontière interne (influence de la pression interne)
    double d = sqrt((x-x0)*(x-x0) + (y-y0)*(y-y0)) - rInner;
    if (d < dInner) {
        double a = (-2*h + 2*hInner)/(dInner*dInner*dInner);
        double b = (3*h  - 3*hInner)/(dInner*dInner);
        double c = 0;
        hfinal = a*d*d*d + b*d*d + c*d + hInner; 
    }

    double tempX = x*x*x*x;

    if (y < 0){
        // Distance du point de contact de la réaction du sol (frontière externe)
        d = sqrt((tempX-x0)*(tempX-x0) + (y+rOuter)*(y+rOuter));
        if (d < dOuter) {
            double a = (-2*h + 2*hOuter)/(dOuter*dOuter*dOuter);
            double b = (3*h  - 3*hOuter)/(dOuter*dOuter);
            double c = 0;
            hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + hOuter); }
    } else {
        // Distance du point de contact du poids de la voiture (frontière externe)
        d = sqrt((tempX-x0)*(tempX-x0) + (y-rOuter)*(y-rOuter));
        if (d < dOuter) {
            double a = (-2*h + 2*hOuter)/(dOuter*dOuter*dOuter);
            double b = (3*h  - 3*hOuter)/(dOuter*dOuter);
            double c = 0;
            hfinal = fmin(hfinal,a*d*d*d + b*d*d + c*d + hOuter); }
    }

    return hfinal;

}


void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();
     
    double x = theGeometry->xCenter;
    double y = theGeometry->yCenter;
    
    double rInner = theGeometry->rInner;
    double rOuter = theGeometry->rOuter;
 
//
//  -1- Construction de la géométrie avec OpenCascade
//      On crée le rectangle
//      On crée les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idInner = gmshModelOccAddDisk(x,y,0.0,rInner,rInner,-1,NULL,0,NULL,0,&ierr); 
    int idOuter = gmshModelOccAddDisk(x,y,0.0,rOuter,rOuter,-1,NULL,0,NULL,0,&ierr); 
    
    int inner[] = {2,idInner};
    int outer[] = {2,idOuter};
    gmshModelOccCut(outer,2,inner ,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
 
//
//  -2- Définition de la fonction callback pour la taille de référence
//      Synchronisation de OpenCascade avec gmsh
//      Génération du maillage (avec l'option Mesh.SaveAll :-)
                  
    geoSetSizeCallback(geoSize);   
    gmshModelOccSynchronize(&ierr);  
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads (avec quelques triangles...) :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 11, &ierr);     // Auparavant = 8 (=11 dixit Alexandre)
//    gmshOptionSetNumber("Mesh.SmoothRatio", 21.5, &ierr);    // Nouvelle option magique
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
//    gmshModelMeshGenerate(2, &ierr);   
   
 
//
//  Plot avec Fltk et vous avez accès au graphique de gmsh :-)
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  
//   chk(ierr);

    
}