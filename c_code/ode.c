#include "flow.c"
//#include "flow_approx.c" // <-- uses the polynomial approximation for square duct flow
#include <math.h>
#include <stdio.h>

int ODE_N = 10;
double ODE_AR = 1.0;
double ODE_W = 1.0;
const int ODE_SYMMETRISE = 0;

double w(double x,double y) {
    if (ODE_SYMMETRISE==1) {
        return 0.5*( Poiseuille_Rectangle_w(x,y,ODE_W,ODE_N,ODE_AR)
                    +Poiseuille_Rectangle_w(y,x,ODE_W,ODE_N,ODE_AR));
    }
    return Poiseuille_Rectangle_w(x,y,ODE_W,ODE_N,ODE_AR);
}

double dwdx(double x,double y) {
    if (ODE_SYMMETRISE==1) {
        return 0.5*( Poiseuille_Rectangle_dwdx(x,y,ODE_W,ODE_N,ODE_AR)
                    +Poiseuille_Rectangle_dwdy(y,x,ODE_W,ODE_N,ODE_AR));
    }
    return Poiseuille_Rectangle_dwdx(x,y,ODE_W,ODE_N,ODE_AR);
}

double dwdy(double x,double y) {
    if (ODE_SYMMETRISE==1) {
        return 0.5*( Poiseuille_Rectangle_dwdy(x,y,ODE_W,ODE_N,ODE_AR)
                    +Poiseuille_Rectangle_dwdx(y,x,ODE_W,ODE_N,ODE_AR));
    }
    return Poiseuille_Rectangle_dwdy(x,y,ODE_W,ODE_N,ODE_AR);
}

double d2wdx2(double x,double y) {
    if (ODE_SYMMETRISE==1) {
        return 0.5*( Poiseuille_Rectangle_d2wdx2(x,y,ODE_W,ODE_N,ODE_AR)
                    +Poiseuille_Rectangle_d2wdy2(y,x,ODE_W,ODE_N,ODE_AR));
    }
    return Poiseuille_Rectangle_d2wdx2(x,y,ODE_W,ODE_N,ODE_AR);
}

double d2wdxdy(double x,double y) {
    if (ODE_SYMMETRISE==1) {
        return 0.5*( Poiseuille_Rectangle_d2wdxdy(x,y,ODE_W,ODE_N,ODE_AR)
                    +Poiseuille_Rectangle_d2wdxdy(y,x,ODE_W,ODE_N,ODE_AR));
    }
    return Poiseuille_Rectangle_d2wdxdy(x,y,ODE_W,ODE_N,ODE_AR);
}

double d2wdy2(double x,double y) {
    if (ODE_SYMMETRISE==1) {
        return 0.5*( Poiseuille_Rectangle_d2wdy2(x,y,ODE_W,ODE_N,ODE_AR)
                    +Poiseuille_Rectangle_d2wdx2(y,x,ODE_W,ODE_N,ODE_AR));
    }
    return Poiseuille_Rectangle_d2wdy2(x,y,ODE_W,ODE_N,ODE_AR);
}

int ode_s(double *dv,double t,double *v,double Cs) {
    // The caller is responsible for providing dv (pre-allocated)
    // t is ignored...
    dv[0] = v[2];
    dv[1] = v[3];
    double temp = 0.5*(Cs+0.5*w(v[0],v[1])); // half e_z
    dv[2] = -temp*dwdx(v[0],v[1]);
    dv[3] = -temp*dwdy(v[0],v[1]);
    //return dv;
    return 0;
}

int ode_variational_s(double *dv,double t,double *v,double Cs) {
    // The caller is responsible for providing dv (pre-allocated)
    // t is ignored...
    // Pre-calculate some things
    double w_ = w(v[0],v[1]);
    double temp = 0.5*(Cs+0.5*w_); // half e_z
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    double d2wdx2_ = d2wdx2(v[0],v[1]);
    double d2wdxdy_ = d2wdxdy(v[0],v[1]);
    double d2wdy2_ = d2wdy2(v[0],v[1]);
    // Calculate the usual ode component
    dv[0] = v[2];
    dv[1] = v[3];
    dv[2] = -temp*dwdx_;
    dv[3] = -temp*dwdy_;
    // Calculate the (matrix) product of the Jacobian with the variational component
    double a = -temp*d2wdx2_-0.25*dwdx_*dwdx_;
    double b = -temp*d2wdxdy_-0.25*dwdx_*dwdy_;
    double c = -temp*d2wdy2_-0.25*dwdy_*dwdy_;
    //double J[16] = {0.0,0.0,1.0,0.0,
    //                0.0,0.0,0.0,1.0,
    //                  a,  b,0.0,0.0,
    //                  b,  c,0.0,0.0};
    for (int i=0;i<4;i++) {
        dv[4+i] = v[12+i];
        dv[8+i] = v[16+i];
        dv[12+i] = a*v[4+i]+b*v[8+i];
        dv[16+i] = b*v[4+i]+c*v[8+i];
    }
    return 0;
}

int ode_p(double *dv,double t,double *v,double Cp,double Gp) {
    // The caller is responsible for providing dv (pre-allocated)
    // t is ignored...
    dv[0] = v[2];
    dv[1] = v[3];
    double w_ = w(v[0],v[1]);
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    double fp = sqrt(2.0*Gp/(1.0+Gp));
    double temp = 0.5*tanh(Cp+0.5*(1.0+Gp)*fp*w_)/fp; // half e_z
    double temp2 = 2.0*Gp*(v[2]*dwdx_+v[3]*dwdy_);
    dv[2] = -temp*((1.0-Gp)*dwdx_+v[2]*temp2);
    dv[3] = -temp*((1.0-Gp)*dwdy_+v[3]*temp2);
    //return dv;
    return 0;
}

int ode_variational_p(double *dv,double t,double *v,double Cp,double Gp) {
    // The caller is responsible for providing dv (pre-allocated)
    // t is ignored...
    // Pre-calculate some things
    double w_ = w(v[0],v[1]);
    double fp = sqrt(2.0*Gp/(1.0+Gp));
    double arg = Cp+0.5*(1+Gp)*fp*w_;
    double temp1 = 0.5*tanh(arg)/fp; // half e_z
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    double temp2 = 2.0*Gp*(v[2]*dwdx_+v[3]*dwdy_);
    double temp3 = ((1.0-Gp)*dwdx_+v[2]*temp2);
    double temp4 = ((1.0-Gp)*dwdy_+v[3]*temp2);
    // Calculate the usual ode component
    dv[0] = v[2];
    dv[1] = v[3];
    dv[2] = -temp1*temp3;
    dv[3] = -temp1*temp4;
    // Calculate the (matrix) product of the Jacobian with the variational component
    double d2wdx2_ = d2wdx2(v[0],v[1]);
    double d2wdxdy_ = d2wdxdy(v[0],v[1]);
    double d2wdy2_ = d2wdy2(v[0],v[1]);
    double temp = cosh(arg);
    //double temp5 = 0.5/(fp*temp*temp); // is this missing terms???
    double temp5 = (1.0+Gp)/(4.0*temp*temp);
    double temp6 = 2.0*Gp*(v[2]*d2wdx2_ +v[3]*d2wdxdy_);
    double temp7 = 2.0*Gp*(v[2]*d2wdxdy_+v[3]*d2wdy2_ );
    double temp8 = 2.0*Gp*temp1;
    double a1 = -temp5*dwdx_*temp3-temp1*(d2wdx2_ *(1.0-Gp)+v[2]*temp6);
    double b1 = -temp5*dwdy_*temp3-temp1*(d2wdxdy_*(1.0-Gp)+v[2]*temp7);
    double c1 = -temp8*(2.0*v[2]*dwdx_+v[3]*dwdy_);
    double d1 = -temp8*v[2]*dwdy_;
    double a2 = -temp5*dwdx_*temp4-temp1*(d2wdxdy_*(1.0-Gp)+v[3]*temp6);
    double b2 = -temp5*dwdy_*temp4-temp1*(d2wdy2_ *(1.0-Gp)+v[3]*temp7);
    double c2 = -temp8*v[3]*dwdx_;
    double d2 = -temp8*(v[2]*dwdx_+2.0*v[3]*dwdy_);
    //double J[16] = {0.0,0.0,1.0,0.0,
    //                0.0,0.0,0.0,1.0,
    //                 a1, b1, c1, d1,
    //                 a2, b2, c2, d2};
    for (int i=0;i<4;i++) {
        dv[4+i] = v[12+i];
        dv[8+i] = v[16+i];
        dv[12+i] = a1*v[4+i]+b1*v[8+i]+c1*v[12+i]+d1*v[16+i];
        dv[16+i] = a2*v[4+i]+b2*v[8+i]+c2*v[12+i]+d2*v[16+i];
    }
    return 0;
}

int ode_o(double *dv,double t,double *v,double Co,double Go) {
    // The caller is responsible for providing dv (pre-allocated)
    // t is ignored...
    dv[0] = v[2];
    dv[1] = v[3];
    double w_ = w(v[0],v[1]);
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    double fo = sqrt(-2.0*Go/(1.0+Go));
    double temp = 0.5*tan(Co+0.5*(1.0+Go)*fo*w_)/fo; // half e_z
    double temp2 = 2.0*Go*(v[2]*dwdx_+v[3]*dwdy_);
    dv[2] = -temp*((1.0-Go)*dwdx_+v[2]*temp2);
    dv[3] = -temp*((1.0-Go)*dwdy_+v[3]*temp2);
    //return dv;
    return 0;
}

int ode_variational_o(double *dv,double t,double *v,double Co,double Go) {
    // The caller is responsible for providing dv (pre-allocated)
    // t is ignored...
    // Pre-calculate some things
    double w_ = w(v[0],v[1]);
    double fo = sqrt(-2.0*Go/(1.0+Go));
    double arg = Co+0.5*(1+Go)*fo*w_;
    double temp1 = 0.5*tan(arg)/fo; // half e_z
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    double temp2 = 2.0*Go*(v[2]*dwdx_+v[3]*dwdy_);
    double temp3 = ((1.0-Go)*dwdx_+v[2]*temp2);
    double temp4 = ((1.0-Go)*dwdy_+v[3]*temp2);
    // Calculate the usual ode component
    dv[0] = v[2];
    dv[1] = v[3];
    dv[2] = -temp1*temp3;
    dv[3] = -temp1*temp4;
    // Calculate the (matrix) product of the Jacobian with the variational component
    double d2wdx2_ = d2wdx2(v[0],v[1]);
    double d2wdxdy_ = d2wdxdy(v[0],v[1]);
    double d2wdy2_ = d2wdy2(v[0],v[1]);
    double temp = cos(arg);
    //double temp5 = 0.5/(fo*temp*temp); // is this missing terms???
    double temp5 = (1.0+Go)/(4.0*temp*temp);
    double temp6 = 2.0*Go*(v[2]*d2wdx2_ +v[3]*d2wdxdy_);
    double temp7 = 2.0*Go*(v[2]*d2wdxdy_+v[3]*d2wdy2_ );
    double temp8 = 2.0*Go*temp1;
    double a1 = -temp5*dwdx_*temp3-temp1*(d2wdx2_ *(1.0-Go)+v[2]*temp6);
    double b1 = -temp5*dwdy_*temp3-temp1*(d2wdxdy_*(1.0-Go)+v[2]*temp7);
    double c1 = -temp8*(2.0*v[2]*dwdx_+v[3]*dwdy_);
    double d1 = -temp8*v[2]*dwdy_;
    double a2 = -temp5*dwdx_*temp4-temp1*(d2wdxdy_*(1.0-Go)+v[3]*temp6);
    double b2 = -temp5*dwdy_*temp4-temp1*(d2wdy2_ *(1.0-Go)+v[3]*temp7);
    double c2 = -temp8*v[3]*dwdx_;
    double d2 = -temp8*(v[2]*dwdx_+2.0*v[3]*dwdy_);
    //double J[16] = {0.0,0.0,1.0,0.0,
    //                0.0,0.0,0.0,1.0,
    //                 a1, b1, c1, d1,
    //                 a2, b2, c2, d2};
    for (int i=0;i<4;i++) {
        dv[4+i] = v[12+i];
        dv[8+i] = v[16+i];
        dv[12+i] = a1*v[4+i]+b1*v[8+i]+c1*v[12+i]+d1*v[16+i];
        dv[16+i] = a2*v[4+i]+b2*v[8+i]+c2*v[12+i]+d2*v[16+i];
    }
    return 0;
}

// The following is an implementation using the spherical angles as per our PRE paper
int ode_ss(double *dv,double t,double *v) {
    double w_ = w(v[0],v[1]);
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    dv[0] = -cos(v[2])*sin(v[3]);
    dv[1] =  sin(v[2]);
    dv[2] = 0.5*cos(v[3])*dwdy_;
    dv[3] = 0.5*tan(v[2])*sin(v[3])*dwdy_-0.5*dwdx_;
    return 0;
}

int ode_variational_ss(double *dv,double t,double *v) {
    // Calculate the usual ode component
    double w_ = w(v[0],v[1]);
    double dwdx_ = dwdx(v[0],v[1]);
    double dwdy_ = dwdy(v[0],v[1]);
    dv[0] = -cos(v[2])*sin(v[3]);
    dv[1] =  sin(v[2]);
    dv[2] =  0.5*dwdy_*cos(v[3]);
    dv[3] =  0.5*dwdy_*tan(v[2])*sin(v[3])-0.5*dwdx_;
    // Calculate the (matrix) product of the Jacobian with the variational component
    double d2wdx2_  = d2wdx2 (v[0],v[1]);
    double d2wdxdy_ = d2wdxdy(v[0],v[1]);
    double d2wdy2_  = d2wdy2 (v[0],v[1]);
    double c1 =  sin(v[2])*sin(v[3]);
    double d1 = -cos(v[2])*cos(v[3]);
    double c2 =  cos(v[2]);
    double a3 =  0.5*d2wdxdy_*cos(v[3]);
    double b3 =  0.5*d2wdy2_ *cos(v[3]);
    double d3 = -0.5*dwdy_   *sin(v[3]);
    double a4 =  0.5*d2wdxdy_*tan(v[2])*sin(v[3])           -0.5*d2wdx2_;
    double b4 =  0.5*d2wdy2_ *tan(v[2])*sin(v[3])           -0.5*d2wdxdy_;
    double c4 =  0.5*(dwdy_/(cos(v[2])*cos(v[2])))*sin(v[3]);
    double d4 =  0.5*dwdy_*tan(v[2])*cos(v[3]);
    //double J[16] = {0.0,0.0, c1, d1,
    //                0.0,0.0, c2,0.0,
    //                 a3, b3,0.0, d3,
    //                 a4, b4, c4, d4};
    for (int i=0;i<4;i++) {
        dv[4+i]  =                     c1*v[12+i]+d1*v[16+i];
        dv[8+i]  =                     c2*v[12+i]           ;
        dv[12+i] = a3*v[4+i]+b3*v[8+i]           +d3*v[16+i];
        dv[16+i] = a4*v[4+i]+b4*v[8+i]+c4*v[12+i]+d4*v[16+i];
    }
    return 0;
}
