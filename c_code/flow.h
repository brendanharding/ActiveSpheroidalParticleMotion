#include <math.h>

double Poiseuille_Rectangle_w(double x,double y,double W,int N,double AR);
double Poiseuille_Rectangle_dwdx(double x,double y,double W,int N,double AR);
double Poiseuille_Rectangle_dwdy(double x,double y,double W,int N,double AR);
double Poiseuille_Rectangle_d2wdx2(double x,double y,double W,int N,double AR);
double Poiseuille_Rectangle_d2wdxdy(double x,double y,double W,int N,double AR);
double Poiseuille_Rectangle_d2wdy2(double x,double y,double W,int N,double AR);

double Poiseuille_Ellipse_w(double x,double y,double W,double AR);
double Poiseuille_Ellipse_dwdx(double x,double y,double W,double AR);
double Poiseuille_Ellipse_dwdy(double x,double y,double W,double AR);
double Poiseuille_Ellipse_d2wdx2(double x,double y,double W,double AR);
double Poiseuille_Ellipse_d2wdxdy(double x,double y,double W,double AR);
double Poiseuille_Ellipse_d2wdy2(double x,double y,double W,double AR);

double Poiseuille_Equilateral_Triangle_w(double x,double y,double W,double L);
double Poiseuille_Equilateral_Triangle_dwdx(double x,double y,double W,double L);
double Poiseuille_Equilateral_Triangle_dwdy(double x,double y,double W,double L);
double Poiseuille_Equilateral_Triangle_d2wdx2(double x,double y,double W,double L);
double Poiseuille_Equilateral_Triangle_d2wdxdy(double x,double y,double W,double L);
double Poiseuille_Equilateral_Triangle_d2wdy2(double x,double y,double W,double L);
