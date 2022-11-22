#include <R.h>
#include <math.h>


void classo(double *y, double *X, double *b0, double *b, double *r, int *n, int *p, 
            double *lambda, int *pen_fac, double *tol, int *cd_maxit) 
{
    
    double b_curr[*p];
    for(int j = 0; j < *p; ++j) 
        b_curr[j] = b[j];

    double xtx[*p];
    for(int j = 0; j < *p; ++j)
    {
        xtx[j] = 0.0;
        for(int i = 0; i < *n; ++i) 
            xtx[j] += pow(X[j*(*n)+i], 2);
    }

    int it;
    double thresh = (*lambda)*(*n);
    for(it = 0; it < *cd_maxit; ++it) 
    {
        // b0 = mean(y - X * b_next)
        *b0 = 0.0;
        for(int i = 0; i < *n; ++i)
            *b0 += r[i];
        *b0 = *b0 / *n;

        for(int j = 0; j < *p; ++j)  // try random ordering
        {
            //b.next[j] <- sum(X[,j] * (y - b0 - X[,-j]%*%b.next[-j]))
            b[j] = b[j] * xtx[j];
            for(int i = 0; i < *n; ++i)
                b[j] += X[j*(*n)+i] * (r[i] - *b0);

            //b.next[j] <- soft(b[j],lambda *n / sum(X[,j]*X[,j]))
            if(pen_fac[j] == 1) 
            {
                if(b[j] > thresh) 
                    b[j] -= thresh;
                else if (b[j] < -thresh) 
                    b[j] += thresh;
                else 
                    b[j] = 0.0;
            }
            b[j] = b[j] / xtx[j];

            for(int i = 0; i < *n; ++i)
                r[i] += X[j*(*n)+i] * (b_curr[j] - b[j]);
        }

        // termination 
        double diff = 0.0;
        for(int j = 0; j < *p; ++j) 
            diff += fabs(b[j] - b_curr[j]);
        if (diff <= (*p) * (*tol))
            break;
        
        // update current b
        for(int j = 0; j < *p; ++j)
            b_curr[j] = b[j];
    }

    if(it == *cd_maxit) 
        Rprintf("Warning: the coordinate descent algorithm does not converge.");

}


void ctlpreg0(double *y, double *X, double *b0, double *b, double *r, int *n, int *p,
             double *tau, double *gamma, int *pen_fac, double *tol, int *dc_maxit, int *cd_maxit)
{
    double lambda = (*gamma) / (*tau);

    double b_curr[*p];
    for(int j = 0; j < *p; ++j) 
        b_curr[j] = b[j];

    int it;
    for(it = 0; it < *dc_maxit; ++it) 
    {
        int pen[*p];
        for(int j = 0; j < *p; ++j)
            pen[j] = (fabs(b[j]) < *tau ? pen_fac[j] : 0);
        
        classo(y, X, b0, b, r, n, p, &lambda, pen, tol, cd_maxit);

        // termination
        double diff = 0.0;
        for(int j = 0; j < *p; ++j) 
            diff += fabs(b[j] - b_curr[j]);
        if (diff <= (*p) * (*tol))
            break;
        
        // update current b
        for(int j = 0; j < *p; ++j)
            b_curr[j] = b[j];
    }

    if(it == *dc_maxit) 
        Rprintf("Warning: the difference of convex functions algorithm does not converge.");

}





