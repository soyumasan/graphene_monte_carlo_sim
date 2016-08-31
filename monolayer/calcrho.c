/* for dummy2 = gridindex                                        
    gridrho(dummy2) = gridrho(dummy2) - 1;
end */

/*calcrho(gridpoints,gridindex)*/

#include "mex.h"

void calcrho(double gridpoints, double *gridindex, double *rhoout, mwSize n)
{
    mwSize i;
    mwSize j;
    
    for (i=0; i<n; i++) {
        j = gridindex[i] - 1;
        rhoout[j] = rhoout[j] - 1;
    }
    
    
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *inMat;
    mwSize ncolsin;
    mwSize ncolsout;
    double *outMat;
    
    inMat = mxGetPr(prhs[1]);
    ncolsin = mxGetN(prhs[1]);
    ncolsout = mxGetScalar(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(ncolsout,1,mxREAL);
    outMat = mxGetPr(plhs[0]);
    calcrho(ncolsout,inMat,outMat,ncolsin);
    
}