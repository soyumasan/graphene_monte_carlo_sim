/*
for dummy2 = 1:N                                                
    dummy1 = gridindex(dummy2);                                 
    gridSx(dummy1) = gridSx(dummy1) + Sx(dummy2);
    gridSy(dummy1) = gridSy(dummy1) + Sy(dummy2);
    gridSz(dummy1) = gridSz(dummy1) + Sz(dummy2);
end

calcgridS(gridpoints,gridindex,Sx,Sy,Sz)
*/

#include "mex.h"

void calcgridS(double gridpoints, double *gridindex, double *Sx, 
        double *Sy, double *Sz, double *gridSx, double *gridSy,
         double *gridSz, mwSize n)
{
    mwSize i;
    mwSize j;
    
    for (i=0; i<n; i++) {
        j = gridindex[i] - 1;
        gridSx[j] = gridSx[j] + Sx[i];
        gridSy[j] = gridSy[j] + Sy[i];
        gridSz[j] = gridSz[j] + Sz[i];
        
    }
    
    
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *gridindex;
    double *Sx;
    double *Sy;
    double *Sz;
    mwSize ncolsin;
    mwSize ncolsout;
    
    double *gridSx;
    double *gridSy;
    double *gridSz;
    
    gridindex = mxGetPr(prhs[1]);
    ncolsin = mxGetN(prhs[1]);
    ncolsout = mxGetScalar(prhs[0]);
    Sx = mxGetPr(prhs[2]);
    Sy = mxGetPr(prhs[3]);
    Sz = mxGetPr(prhs[4]);
    
    
    plhs[0] = mxCreateDoubleMatrix(1,ncolsout,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,ncolsout,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,ncolsout,mxREAL);
    
    gridSx = mxGetPr(plhs[0]);
    gridSy = mxGetPr(plhs[1]);
    gridSz = mxGetPr(plhs[2]);
    
    calcgridS(ncolsout,gridindex,Sx,Sy,Sz,gridSx,gridSy,gridSz,ncolsin);
    
}