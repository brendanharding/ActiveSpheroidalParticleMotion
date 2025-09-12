#include <math.h>
#include <stdio.h>

double Poiseuille_Rectangle_w(double x,double y,double W,int N,double AR) {
    // Centroid at origin, sides at x=+AR,-AR and y=+1,-1
    double w = 1.0-y*y;
    double k;
    double sign = -1;
    for (int n=0;n<N;n++) {
        k = 0.5*M_PI*(2*n+1);
        w += 4.0*sign*cosh(k*x)*cos(k*y)/(k*k*k*cosh(k*AR));
        sign *= -1;
    }
    return W*w;
}

double Poiseuille_Rectangle_dwdx(double x,double y,double W,int N,double AR) {
    double w = 0.0;
    double k;
    double sign = -1;
    for (int n=0;n<N;n++) {
        k = 0.5*M_PI*(2*n+1);
        w += 4.0*sign*sinh(k*x)*cos(k*y)/(k*k*cosh(k*AR));
        sign *= -1;
    }
    return W*w;
}

double Poiseuille_Rectangle_dwdy(double x,double y,double W,int N,double AR) {
    double w = -2.0*y;
    double k;
    double sign = +1; // note starting with opposite sign
    for (int n=0;n<N;n++) {
        k = 0.5*M_PI*(2*n+1);
        w += 4.0*sign*cosh(k*x)*sin(k*y)/(k*k*cosh(k*AR));
        sign *= -1;
    }
    return W*w;
}

double Poiseuille_Rectangle_d2wdx2(double x,double y,double W,int N,double AR) {
    double w = 0.0;
    double k;
    double sign = -1;
    for (int n=0;n<N;n++) {
        k = 0.5*M_PI*(2*n+1);
        w += 4.0*sign*cosh(k*x)*cos(k*y)/(k*cosh(k*AR));
        sign *= -1;
    }
    return W*w;
}

double Poiseuille_Rectangle_d2wdxdy(double x,double y,double W,int N,double AR) {
    double w = 0.0;
    double k;
    double sign = +1; // note starting with opposite sign
    for (int n=0;n<N;n++) {
        k = 0.5*M_PI*(2*n+1);
        w += 4.0*sign*sinh(k*x)*sin(k*y)/(k*cosh(k*AR));
        sign *= -1;
    }
    return W*w;
}

double Poiseuille_Rectangle_d2wdy2(double x,double y,double W,int N,double AR) {
    double w = -2.0;
    double k;
    double sign = +1; // note starting with opposite sign
    for (int n=0;n<N;n++) {
        k = 0.5*M_PI*(2*n+1);
        w += 4.0*sign*cosh(k*x)*cos(k*y)/(k*cosh(k*AR));
        sign *= -1;
    }
    return W*w;
}

double Poiseuille_Ellipse_w(double x,double y,double W,double AR) {
    // centroid at origin, major axis from x=-AR to x=+AR, minor axis from y=-1 to y=+1
    return W*(1.0-x*x/(AR*AR)-y*y);
}

double Poiseuille_Ellipse_dwdx(double x,double y,double W,double AR) {
    return -W*2.0*x/(AR*AR);
}

double Poiseuille_Ellipse_dwdy(double x,double y,double W,double AR) {
    return -W*2.0*y;
}

double Poiseuille_Ellipse_d2wdx2(double x,double y,double W,double AR) {
    return -W*2.0/(AR*AR);
}

double Poiseuille_Ellipse_d2wdxdy(double x,double y,double W,double AR) {
    return 0.0;
}

double Poiseuille_Ellipse_d2wdy2(double x,double y,double W,double AR) {
    return -W*2.0;
}

double Poiseuille_Equilateral_Triangle_w(double x,double y,double W,double L) {
    // Centroid at origin, side length L, base parallel to x-axis
    double c = L/sqrt(12.0);
    double w = (y+c)*((y-2.0*c)*(y-2.0*c)-3.0*x*x);
    return w*W/(4.0*c*c*c);
}

double Poiseuille_Equilateral_Triangle_dwdx(double x,double y,double W,double L) {
    // Centroid at origin, side length L, base parallel to x-axis
    double c = L/sqrt(12.0);
    double w = (y+c)*(-6.0*x);
    return w*W/(4.0*c*c*c);
}

double Poiseuille_Equilateral_Triangle_dwdy(double x,double y,double W,double L) {
    // Centroid at origin, side length L, base parallel to x-axis
    double c = L/sqrt(12.0);
    double w = ((y-2.0*c)*(y-2.0*c)-3.0*x*x)+(y+c)*(2.0*(y-2.0*c));
    return w*W/(4.0*c*c*c);
}

double Poiseuille_Equilateral_Triangle_d2wdx2(double x,double y,double W,double L) {
    // Centroid at origin, side length L, base parallel to x-axis
    double c = L/sqrt(12.0);
    double w = -6.0*(y+c);
    return w*W/(4.0*c*c*c);
}

double Poiseuille_Equilateral_Triangle_d2wdxdy(double x,double y,double W,double L) {
    // Centroid at origin, side length L, base parallel to x-axis
    double c = L/sqrt(12.0);
    double w = -6.0*x;
    return w*W/(4.0*c*c*c);
}

double Poiseuille_Equilateral_Triangle_d2wdy2(double x,double y,double W,double L) {
    // Centroid at origin, side length L, base parallel to x-axis
    double c = L/sqrt(12.0);
    double w = +6.0*(y-c);
    return w*W/(4.0*c*c*c);
}
