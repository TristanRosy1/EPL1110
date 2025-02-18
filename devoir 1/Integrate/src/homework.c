#include <stdio.h>
#include <math.h>
#include "glfem.h"


double interpolate ( double u [3] , double xi , double eta ){

    return u [0]* xi + u [1]* eta + u [2]*(1.0 - xi - eta );

}



double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;

    const double xi [3] = {0.166666666666667 ,0.666666666666667 ,0.166666666666667};
    const double eta [3] = {0.166666666666667 ,0.166666666666667 ,0.666666666666667};
    const double weight [3] = {0.166666666666667 ,0.166666666666667 ,0.166666666666667};

    double jacob = fabs (( x [0] - x [1]) * ( y [0] - y [2]) - ( x [0] - x [2]) * ( y [0] - y [1]) ) ;


    double xLoc ;
    double yLoc ;

    for ( int i =0 ; i <3 ; i ++) {
        double xiLoc = xi [ i ];
        double etaLoc = eta [ i ];
        xLoc = interpolate (x , xiLoc , etaLoc ) ;
        yLoc = interpolate (y , xiLoc , etaLoc ) ;
        I += f ( xLoc , yLoc ) * weight [ i ];
    }







//
// ... A modifier :-)
//
//
// Pour dessiner l'element, les sommets du triangle :-)
// Decommenter la ligne pour dessiner aussi les points d'integration
//

    glfemSetColor(GLFEM_BLACK); glfemDrawElement( x , y , 3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes( x , y ,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(&xLoc, &yLoc,3);



     return I * jacob;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

//
// ... A modifier :-)
// y-compris la ligne juste en dessous :-)
//
    const int nodes [4][3] = {{0 ,3 ,5} ,{3 ,1 ,4} ,{5 ,4 ,2} ,{3 ,4 ,5}};
    const double xi [6] = {0.0 ,1.0 ,0.0 ,0.5 ,0.5 ,0.0};
    const double eta [6] = {0.0 ,0.0 ,1.0 ,0.0 ,0.5 ,0.5};
    double xLoc [3];
    double yLoc [3];

    if ( n <= 0) return integrate (x , y , f  ) ;

    double I = 0.0;
    for (int i =0; i <4; i ++) {
        for ( int j =0; j <3; j ++) {
            double xiLoc = xi [ nodes [ i ][ j ]];
            double etaLoc = eta [ nodes [ i ][ j ]];
            xLoc [ j ] = interpolate (x , xiLoc , etaLoc ) ;
            yLoc [ j ] = interpolate (y , xiLoc , etaLoc ) ;
        }
        I += integrateRecursive ( xLoc , yLoc , f , n -1 ) ;
    }
    return I;
}
