#include "ode.c"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// Global
double IVP_Cs = 0.0;
double IVP_Cp = 0.0;
double IVP_Co = 0.0;
double IVP_Gp = +0.6;
double IVP_Go = -0.6;
const int IVP_CASE = 3;

int ode(double *dv,double t,double *v) {
    switch (IVP_CASE) {
        case 1:
            ode_variational_s(dv,t,v,IVP_Cs);
            break;
        case 2:
            ode_variational_p(dv,t,v,IVP_Cp,IVP_Gp);
            break;
        case 3:
            ode_variational_o(dv,t,v,IVP_Co,IVP_Go);
            break;
        case 4:
            ode_variational_ss(dv,t,v);
            break;
        default:
            ode_s(dv,t,v,IVP_Cs);
            break;
    }
    return 0;
}

int RungeKutta_Fixed_4(double* result,double t0,double* y0,int ny,double dt,double nt,int (*func)(double *,double,double*)) {
    // result should have length (nt+1)*ny
    double y[ny];
    for (int i=0;i<ny;i++) {
        y[i]=y0[i];
        result[i]=y[i];
    }
    double k1[ny],k2[ny],k3[ny],k4[ny],temp[ny];
    double t = t0;
    for (int k=0;k<nt;k++) {
        func(k1,t,y);
        for (int i=0;i<ny;i++) temp[i] = y[i]+0.5*dt*k1[i];
        func(k2,t+0.5*dt,temp);
        for (int i=0;i<ny;i++) temp[i] = y[i]+0.5*dt*k2[i];
        func(k3,t+0.5*dt,temp);
        for (int i=0;i<ny;i++) temp[i] = y[i]+dt*k3[i];
        func(k4,t+dt,temp);
        for (int i=0;i<ny;i++) {
            y[i] += (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*dt/6.0;
            result[(k+1)*ny+i] = y[i];
        }
    }
    return 0;
}

int QR_inplace(double* A,int n,double* R) {
    // Uses (classical) Gram-Schmidt,
    // so is known to be numerically unstable,
    // but will be okay for small matrices.
    double u_i[n];
    double inner;
    for (int i=0;i<n;i++) {
        for (int k=0;k<n;k++) u_i[k] = A[k*n+i]; // u_i = a_i
        if (i>0) {
            for (int j=0;j<i;j++) {
                inner = 0.0;
                for (int k=0;k<n;k++) inner += u_i[k]*A[k*n+j]; // <u_i,e_j> (= <a_i,e_j>)
                R[j*n+i] = inner;
                R[i*n+j] = 0.0; // could be avoided, but can't hurt...
                for (int k=0;k<n;k++) u_i[k] -= inner*A[k*n+j]; // u_i -= <u_i,e_j> e_j
            }
        }
        inner = 0.0;
        for (int k=0;k<n;k++) inner += u_i[k]*A[k*n+i]; // <u_i,a_i> (= <u_i,u_i>)
        //R[i*n+i] = inner;
        double norm = sqrt(inner);
        R[i*n+i] = norm;
        for (int k=0;k<n;k++) A[k*n+i] = u_i[k]/norm; // replace a_i with e_i
    }
    return 0;
}

int test_QR() {
    // Example from wikipedia
    double A[9] = {12.0,-51.0,  4.0,
        6.0,167.0,-68.0,
        -4.0, 24.0,-41.0};
    double R[9];
    QR_inplace(A,3,R);
    printf("Q:\n");
    for (int i=0;i<3;i++) printf("%.6f,%.6f,%.6f\n",A[3*i+0],A[3*i+1],A[3*i+2]);
    printf("R:\n");
    for (int i=0;i<3;i++) printf("%.6f,%.6f,%.6f\n",R[3*i+0],R[3*i+1],R[3*i+2]);
    return 0;
}

int RungeKutta_Fixed_Variational_4(double* result,double t0,double* y0,int ny,double dt,double nt,double* lyap,int no,int (*func)(double *,double,double*)) {
    // result should have length (nt+1)*ny
    // lyap should have length roughly (nt+1)*(ny/no)
    // simultaneously solves a variational problem
    int N = ny*(ny+1);
    double y[N];
    for (int i=0;i<N;i++) y[i] = y0[i];
    for (int i=0;i<ny;i++) {
        result[i] = y[i];
        lyap[i] = 1.0; // initialise this with 1.0?
    }
    double k1[N],k2[N],k3[N],k4[N],temp[N];
    double t = t0;
    int j = 1;
    for (int k=0;k<nt;k++) {
        func(k1,t,y);
        for (int i=0;i<N;i++) temp[i] = y[i]+0.5*dt*k1[i];
        func(k2,t+0.5*dt,temp);
        for (int i=0;i<N;i++) temp[i] = y[i]+0.5*dt*k2[i];
        func(k3,t+0.5*dt,temp);
        for (int i=0;i<N;i++) temp[i] = y[i]+dt*k3[i];
        func(k4,t+dt,temp);
        for (int i=0;i<N;i++) {
            y[i] += (k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])*dt/6.0;
        }
        // Orthonormalise the variational part as needed
        if ((k+1)%no==0) {
            QR_inplace(y+ny,ny,temp); // re-use temp here
            for (int i=0;i<ny;i++) lyap[j*ny+i] = temp[i*ny+i];
            j++;
        }
        for (int i=0;i<ny;i++) {
            result[(k+1)*ny+i] = y[i];
        }
    }
    return j;
}

int test_1() {
    // This will run the solver and output the trajectory.
    
    // Adjust the global parameter W (from ode.c), and optionally N
    //ODE_N = 100; // default global value is 10
    ODE_W = 10.0/Poiseuille_Rectangle_w(0.0,0.0,1.0,ODE_N,ODE_AR);
    printf("Flow velocity at origin: %.6f\n",w(0.0,0.0));
    // Other setup
    double y0[4] = {0.4,0.2,0.0,0.0};
    double ez0 = -1.0;
    double f_;
    switch (IVP_CASE) {
        case 1:
            IVP_Cs = ez0-0.5*w(y0[0],y0[1]); // <-- global
            printf("Computation will be inefficient!\n");
            break;
        case 2:
            f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
            IVP_Cp = atanh(f_*ez0)-0.5*(1.0+IVP_Gp)*f_*w(y0[0],y0[1]); // <-- global
            printf("Computation will be inefficient!\n");
            break;
        case 3:
            f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
            IVP_Co = atan(f_*ez0)-0.5*(1.0+IVP_Go)*f_*w(y0[0],y0[1]); // <-- global
            printf("Computation will be inefficient!\n");
            break;
        case 4:
            // nothing to do
            break;
        default:
            IVP_Cs = ez0-0.5*w(y0[0],y0[1]);
            break;
    }
    // Specify time stepping and run
    double t0 = 0.0;
    double dt = 0.1;
    int nt = 1000;
    double* result = malloc((nt+1)*4*sizeof(double));
    RungeKutta_Fixed_4(result,t0,y0,4,dt,nt,*ode);
    // Output
    for (int i=0;i<=nt;i++) {
        printf("%.6f: %.6f,%.6f,%.6f,%.6f\n",i*dt,result[4*i+0],result[4*i+1],result[4*i+2],result[4*i+3]);
    }
    // cleanup
    free(result);
    return 0;
}

int test_2() {
    // This will run the solver and output the (approximate) Lyapunov exponents.
    
    // Adjust the global parameter W (from ode.c), and optionally N
    //ODE_N = 100; // default global value is 10
    ODE_W = 10.0/Poiseuille_Rectangle_w(0.0,0.0,1.0,ODE_N,ODE_AR);
    printf("Flow velocity at origin: %.6f\n",w(0.0,0.0));
    // Other setup
    double y0[20] = {0.4,0.2,0.0,0.0,
                     1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0};
    double ez0 = -1.0;
    double f_;
    switch (IVP_CASE) {
        case 1:
            IVP_Cs = ez0-0.5*w(y0[0],y0[1]); // <-- global
            break;
        case 2:
            f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
            IVP_Cp = atanh(f_*ez0)-0.5*(1.0+IVP_Gp)*f_*w(y0[0],y0[1]); // <-- global
            break;
        case 3:
            f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
            IVP_Co = atan(f_*ez0)-0.5*(1.0+IVP_Go)*f_*w(y0[0],y0[1]); // <-- global
            break;
        case 4:
            printf("This test must be run with a variational solver! Exiting.\n");
            return 0;
        default:
            printf("This test must be run with a variational solver! Exiting.\n");
            return 0;
    }
    // Specify time stepping and run
    double t0 = 0.0;
    double dt = 0.1;
    int nt = 1000;
    int no = 10;
    double* result = malloc((nt+1)*4*sizeof(double));
    double* lyap = malloc((nt/no+1)*4*sizeof(double));
    int j_max = RungeKutta_Fixed_Variational_4(result,t0,y0,4,dt,nt,lyap,no,*ode);
    // post processing
    int j = 0;
    double lyap_exps[4] = {0.0,0.0,0.0,0.0};
    for (int i=0;i<=nt;i++) {
        if ((i%no)==0) {
            printf("%.6f: %.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",i*dt,
                   result[4*i+0],result[4*i+1],result[4*i+2],result[4*i+3],
                   lyap[4*j+0],lyap[4*j+1],lyap[4*j+2],lyap[4*j+3]);
            for (int k=0;k<4;k++) lyap_exps[k] += log(lyap[4*j+k]);
            j++;
        } else {
            printf("%.6f: %.6f,%.6f,%.6f,%.6f\n",i*dt,
                   result[4*i+0],result[4*i+1],result[4*i+2],result[4*i+3]);
        }
    }
    // Output
    for (int k=0;k<4;k++) lyap_exps[k] /= (dt*nt);
    printf("exponents: %.6f,%.6f,%.6f,%.6f\n",
           lyap_exps[0],lyap_exps[1],lyap_exps[2],lyap_exps[3]);
    // Cleanup
    free(result);
    free(lyap);
    return 0;
}

int test_3() {
    // This does a sweep over points in the cross-section outputing Lyapunov exponents
    
    // Adjust the global parameter W (from ode.c), and optionally N
    //ODE_N = 100; // default global value is 10
    ODE_W = 10.0/Poiseuille_Rectangle_w(0.0,0.0,1.0,ODE_N,ODE_AR);
    printf("Flow velocity at origin: %.6f\n",w(0.0,0.0));
    // Set various time stepping parameters...
    double t0 = 0.0;
    double dt = 1./64;
    int nt = 65536; // Total run time is T=256
    int no = 64;
    // Prepare intial condition and other working arrays
    double y0[20] = {0.4,0.2,0.0,0.0,
                     1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0};
    double ez0 = -1.0;
    double f_;
    double* result = malloc((nt+1)*4*sizeof(double));
    double* lyap = malloc((nt/no+1)*4*sizeof(double));
    // Loop over points in one quadrant
    int nx = 64;
    for (int ix=0;ix<=64;ix++) {
        y0[0] = -1.0+((double) ix)/nx;
        for (int iy=0;iy<=ix;iy++) {
            y0[1] = -1.0+((double) iy)/nx;
            
            // Carry out the integration
            switch (IVP_CASE) {
                case 1:
                    IVP_Cs = ez0-0.5*w(y0[0],y0[1]); // <-- global
                    break;
                case 2:
                    f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
                    IVP_Cp = atanh(f_*ez0)-0.5*(1.0+IVP_Gp)*f_*w(y0[0],y0[1]); // <-- global
                    break;
                case 3:
                    f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
                    IVP_Co = atan(f_*ez0)-0.5*(1.0+IVP_Go)*f_*w(y0[0],y0[1]); // <-- global global
                    break;
                case 4:
                    printf("This test must be run with a variational solver! Exiting.\n");
                    return 0; // neglects cleanup...
                default:
                    printf("This test must be run with a variational solver! Exiting.\n");
                    return 0; // neglects cleanup...
            }
            int j_max = RungeKutta_Fixed_Variational_4(result,t0,y0,4,dt,nt,lyap,no,*ode);
            
            // Calculate the lyapunov exponents
            double lyap_exps[4] = {0.0,0.0,0.0,0.0};
            for (int j=0;j<j_max;j++) {
                for (int k=0;k<4;k++) lyap_exps[k] += log(lyap[4*j+k]);
            }
            for (int k=0;k<4;k++) lyap_exps[k] /= (dt*nt);
            printf("%.15f,%.15f: %.6f,%.6f,%.6f,%.6f\n",y0[0],y0[1],
                   lyap_exps[0],lyap_exps[1],lyap_exps[2],lyap_exps[3]);
        }
    }
    // Cleanup
    free(result);
    free(lyap);
    return 0;
}

int classify_orbit(double *result,int nt) {
    // This is based on my *new* classification scheme
    // Check for escaping orbit
    int is_escaping = 0; // assume false
    for (int i=0;i<nt+1;i++) {
        if (result[i*4+0]<-ODE_AR || result[i*4+0]>+ODE_AR || result[i*4+1]<-1.0 || result[i*4+1]>+1.0) {is_escaping=1;break;};
    }
    if (is_escaping==1) return 0;
    // Determine whether tumbling
    int is_tumbling = 0; // assume false
    double W_max = w(0.0,0.0); // <-- Note: this is different from ODE_W in this implementation
    double f_;
    switch (IVP_CASE) {
        case 1:
            if (IVP_Cs>=1.0-0.5*W_max) is_tumbling = 1;
            break;
        case 2:
            f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
            if (IVP_Cp>=atanh(f_)-0.5*(1.0+IVP_Gp)*f_*W_max) is_tumbling = 1;
            break;
        case 3:
            f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
            if (IVP_Co>=atan(f_)-0.5*(1.0+IVP_Go)*f_*W_max) is_tumbling = 1;
            break;
        default:
            if (IVP_Cs>=1.0-0.5*W_max) is_tumbling = 1; // assume a spherical particle run
            break;
    }
    // Continue onto subsequent checks
    if (is_tumbling==0) {
        int is_simple = 1; // assume true
        double ez;
        for (int i=0;i<nt+1;i++) {
            switch (IVP_CASE) {
                case 1:
                    ez = IVP_Cs+0.5*w(result[i*4+0],result[i*4+1]);
                    break;
                case 2:
                    ez = tanh(IVP_Cp+0.5*(1.0+IVP_Gp)*f_*w(result[i*4+0],result[i*4+1]))/f_;
                    break;
                case 3:
                    ez = tan(IVP_Co+0.5*(1.0+IVP_Go)*f_*w(result[i*4+0],result[i*4+1]))/f_;
                    break;
                default:
                    ez = -1; // really should raise an error here
                    break;
            }
            if (ez>0.0) {is_simple=0;break;}
        }
        if (is_simple==1) return 1;
        // otherwise, check for other forms of confinement
        int x_confined = 1; // assume true
        int y_confined = 1; // assume true
        for (int i=0;i<nt+1;i++) {
            if (result[i*4+0]<-0.25 || result[i*4+0]>+0.25) {x_confined=0;break;};
        }
        for (int i=0;i<nt+1;i++) {
            if (result[i*4+1]<-0.25 || result[i*4+1]>+0.25) {y_confined=0;break;};
        }
        if (x_confined*y_confined==1) {
            return 1; // assume central swinging = simple swinging
        } else if (x_confined==1) {
            return 2;
        } else if (y_confined==1) {
            return 3;
        }
        return 6;
    } else {
        // check for off-centre trapping
        int x_pos = 0; // assume false
        int x_neg = 0; // assume false
        int y_pos = 0; // assume false
        int y_neg = 0; // assume false
        for (int i=0;i<nt+1;i++) {
            if (result[i*4+0]<0.0) x_neg = 1;
            if (result[i*4+0]>0.0) x_pos = 1;
            if (result[i*4+1]<0.0) y_neg = 1;
            if (result[i*4+1]>0.0) y_pos = 1;
        }
        if (x_neg*x_pos==0 || y_neg*y_pos==0) return 4;
        return 5;
    }
    return -1; // should never reach here
}

int classify_orbit_old(double *result,int nt) {
    // This is based on Rahil's PRE classification scheme
    // Check for escaping orbit (my addition)
    int is_escaping = 0; // assume false
    for (int i=0;i<nt+1;i++) {
        if (result[i*4+0]<-ODE_AR || result[i*4+0]>+ODE_AR || result[i*4+1]<-1.0 || result[i*4+1]>+1.0) {is_escaping=1;break;};
    }
    if (is_escaping==1) return 0;
    // Check for central confinement
    int x_confined = 1; // assume true
    int y_confined = 1; // assume true
    for (int i=0;i<nt+1;i++) {
        if (result[i*4+0]<-0.25 || result[i*4+0]>+0.25) {x_confined=0;break;};
    }
    for (int i=0;i<nt+1;i++) {
        if (result[i*4+1]<-0.25 || result[i*4+1]>+0.25) {y_confined=0;break;};
    }
    if (x_confined*y_confined==1) {
        return 1;
    } else if (x_confined==1) {
        return 2;
    } else if (y_confined==1) {
        return 3;
    }
    // Check for off-centre trapping
    int x_pos = 0; // assume false
    int x_neg = 0; // assume false
    int y_pos = 0; // assume false
    int y_neg = 0; // assume false
    for (int i=0;i<nt+1;i++) {
        if (result[i*4+0]<0.0) x_neg = 1;
        if (result[i*4+0]>0.0) x_pos = 1;
        if (result[i*4+1]<0.0) y_neg = 1;
        if (result[i*4+1]>0.0) y_pos = 1;
    }
    if (x_neg*x_pos==0 || y_neg*y_pos==0) return 4;
    // Check if it misses the centre (i.e. tumbling)
    int tumbling = 1; // assume true
    for (int i=0;i<nt+1;i++) {
        if (result[i*4+0]>=-0.25 && result[i*4+0]<=+0.25 && result[i*4+1]>=-0.25 && result[i*4+1]<=+0.25) {tumbling=0;break;};
    }
    if (tumbling==1) return 5;
    // If nothing else, classify as wandering
    return 6;
}

int test_4() {
    // This does a sweep over points in the cross-section
    // outputing Lyapunov exponents and a classification
    
    // Adjust the global parameter W (from ode.c), and optionally N
    //ODE_N = 100; // default global value is 10
    ODE_W = 5.0/Poiseuille_Rectangle_w(0.0,0.0,1.0,ODE_N,ODE_AR);
    printf("Flow velocity at origin: %.6f\n",w(0.0,0.0));
    // Set various time stepping parameters...
    // First run, mirroring what is done in my python code
    double t0 = 0.0;
    double dt = 1./16;
    int nt = 4096; // Total run time is T=256
    if (IVP_CASE==2) nt *= 2; // double the time in the prolate case
    int no = 64;
    // Prepare intial condition and other working arrays
    double y0[20] = {0.4,0.2,0.0,0.0,
                     1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0};
    double ez0 = -1.0;
    double f_;
    double* result = malloc((nt+1)*4*sizeof(double));
    double* lyap = malloc((nt/no+1)*4*sizeof(double));
    // Loop over quadrant
    int nx = 64;
    for (int ix=0;ix<=64;ix++) {
        y0[0] = -1.0+((double) ix)/nx;
        for (int iy=0;iy<=ix;iy++) {
            y0[1] = -1.0+((double) iy)/nx;
            // Carry out the integration
            switch (IVP_CASE) {
                case 1:
                    IVP_Cs = ez0-0.5*w(y0[0],y0[1]); // <-- global
                    break;
                case 2:
                    f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
                    IVP_Cp = atanh(f_*ez0)-0.5*(1.0+IVP_Gp)*f_*w(y0[0],y0[1]); // <-- global
                    break;
                case 3:
                    f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
                    IVP_Co = atan(f_*ez0)-0.5*(1.0+IVP_Go)*f_*w(y0[0],y0[1]); // <-- global
                    break;
                case 4:
                    printf("This test must be run with a variational solver! Exiting.\n");
                    return 0; // neglects cleanup...
                default:
                    printf("This test must be run with a variational solver! Exiting.\n");
                    return 0; // neglects cleanup...
            }
            int j_max = RungeKutta_Fixed_Variational_4(result,t0,y0,4,dt,nt,lyap,no,*ode);
            // Calculate the lyapunov exponents
            double lyap_exps[4] = {0.0,0.0,0.0,0.0};
            for (int j=0;j<j_max;j++) {
                for (int k=0;k<4;k++) lyap_exps[k] += log(lyap[4*j+k]);
            }
            for (int k=0;k<4;k++) lyap_exps[k] /= (dt*nt);
            // Calculate the classification
            int classification = classify_orbit(result,nt);
            //int classification = classify_orbit_old(result,nt);
            // Print the results
            printf("%.15f,%.15f: %.6f,%.6f,%.6f,%.6f,%d\n",y0[0],y0[1],
                   lyap_exps[0],lyap_exps[1],lyap_exps[2],lyap_exps[3],
                   classification);
        }
    }
    // Cleanup
    free(result);
    free(lyap);
    return 0;
}

double check_hamiltonian(double* result,int nt) {
    double max_error = 0.0;
    double f_;
    double x,y,e_x,e_y,e_z;
    for (int i=0;i<nt+1;i++) {
        x = result[i*4+0];
        y = result[i*4+1];
        e_x = result[i*4+2]; // is thie correct for this implementation?
        e_y = result[i*4+3];
        e_z = 0.0;
        switch (IVP_CASE) {
            case 1:
                e_z = IVP_Cs+0.5*w(x,y);
                break;
            case 2:
                f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
                e_z = tanh(IVP_Cp+0.5*(1.0+IVP_Gp)*f_*w(x,y))/f_;
                break;
            case 3:
                f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
                e_z = tan(IVP_Co+0.5*(1.0+IVP_Go)*f_*w(x,y))/f_;
                break;
            case 4:
                // nothing to do
                break;
            default:
                // nothing to do
                break;
        }
        double error = fabs(e_x*e_x+e_y*e_y+e_z*e_z-1);// /2;
        if (error>max_error) max_error = error;
    }
    return max_error;
}

int test_5() {
    // Some "production" runs, initial position sweep checking
    // Lyapunov exponents, classification and the Hamiltonian.
    
    // Adjust the global parameter W (from ode.c), and optionally N
    //ODE_N = 100; // default global value is 10
    ODE_W = 10.0/Poiseuille_Rectangle_w(0.0,0.0,1.0,ODE_N,ODE_AR);
    printf("Flow velocity at origin: %.6f\n",w(0.0,0.0));
    // Set various time stepping parameters...
    // First run, mirroring what is done in my python code
    double t0 = 0.0;
    double dt = 1./(1<<10);
    int nt = 1<<(10+10); // Total run time is T=1024
    if (IVP_CASE==2) nt *= 2; // double the time in the prolate case
    int no = (1<<8); // Frequency with which to orthonormalise
    // Prepare intial condition and other working arrays
    double y0[20] = {0.4,0.2,0.0,0.0,
                     1.0,0.0,0.0,0.0,
                     0.0,1.0,0.0,0.0,
                     0.0,0.0,1.0,0.0,
                     0.0,0.0,0.0,1.0};
    double ez0 = -1.0;
    double f_;
    double* result = malloc((nt+1)*4*sizeof(double));
    double* lyap = malloc((nt/no+1)*4*sizeof(double));
    // Loop over quadrant
    int nx = 64;
    for (int ix=0;ix<=64;ix++) {
        y0[0] = -1.0+((double) ix)/nx;
        for (int iy=0;iy<=ix;iy++) {
            y0[1] = -1.0+((double) iy)/nx;
            // Carry out the integration
            switch (IVP_CASE) {
                case 1:
                    IVP_Cs = ez0-0.5*w(y0[0],y0[1]); // <-- global
                    break;
                case 2:
                    f_ = sqrt( 2.0*IVP_Gp/(1.0+IVP_Gp));
                    IVP_Cp = atanh(f_*ez0)-0.5*(1.0+IVP_Gp)*f_*w(y0[0],y0[1]); // <-- global
                    break;
                case 3:
                    f_ = sqrt(-2.0*IVP_Go/(1.0+IVP_Go));
                    IVP_Co = atan(f_*ez0)-0.5*(1.0+IVP_Go)*f_*w(y0[0],y0[1]); // <-- global
                    break;
                case 4:
                    printf("This test must be run with a variational solver! Exiting.\n");
                    return 0; // neglects cleanup...
                default:
                    printf("This test must be run with a variational solver! Exiting.\n");
                    return 0; // neglects cleanup...
            }
            int j_max = RungeKutta_Fixed_Variational_4(result,t0,y0,4,dt,nt,lyap,no,*ode);
            // Calculate the lyapunov exponents
            double lyap_exps[4] = {0.0,0.0,0.0,0.0};
            for (int j=0;j<j_max;j++) {
                for (int k=0;k<4;k++) lyap_exps[k] += log(lyap[4*j+k]);
            }
            for (int k=0;k<4;k++) lyap_exps[k] /= (dt*nt);
            // Calculate the classification
            int classification = classify_orbit(result,nt);
            //int classification = classify_orbit_old(result,nt);
            // Check deviation in the hamiltonian (actually ||e||)
            double ham_error = check_hamiltonian(result,nt);
            // print the results
            printf("%.6f,%.6f: %.6f,%.6f,%.6f,%.6f,%d,%.1E\n",y0[0],y0[1],
                   lyap_exps[0],lyap_exps[1],lyap_exps[2],lyap_exps[3],
                   classification,ham_error);
        }
    }
    // Cleanup
    free(result);
    free(lyap);
    return 0;
}

int main() {
    
    //test_QR(); // tests the QR algorithm
    //test_1(); // computes and outputs a trajectory
    //test_2(); // computes a trajectory and outputs the Lyapunov exponents
    //test_3(); // computes Lyapunov exponents over a quadrant
    //test_4(); // computes Lyapunov exponents and classification over a quadrant
    test_5(); // as in test 4, production quality, also checks Hamiltonian
    
    return 0;
}
