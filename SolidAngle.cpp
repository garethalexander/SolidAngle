#include <math.h>
#include <cmath>
#include <string.h>
#include "SolidAngle.h"
#include "Geometry.h"
#include "InputOutput.h"
#include "Constants.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main (int argc, char*argv[])
{
    if(argc<2)
    {
        std::cerr << "Please specify a folder name - see README" << endl;
        std::cerr << "e.g. SolidAngle Unknot if your executable is called SolidAngle and there's a folder within the knots folder called Unknot" << endl;
        return 1;
    }

    knot_filename = argv[1];

    if(0==InitialiseSystemParameters())
    {
        // omega is a 3D grid of solid angle values in space around our curve. It's stored in row-major order
        // with the origin of physical space in its centre.
        vector<double>omega(Nx*Ny*Nz);
        vector<double>Bx(0);
        vector<double>By(0);
        vector<double>Bz(0);
        // If the user wants the gradient of the solid angle computed directly too, arrays for its components
        if(Gradient)
        {
            Bx.resize(Nx*Ny*Nz);
            By.resize(Nx*Ny*Nz);
            Bz.resize(Nx*Ny*Nz);
        }
        // all data read in or computed about the curve itself is stored in this data structure.
        Link Curve;
        if(0==InitialiseFromFile(Curve))
        {
            cout << "Filling in the Geometry of the Input Curve \n";
            ComputeGeometry(Curve);

            cout << "Beginning the Solid Angle Calculations \n";
            ComputeSolidAngleAllPoints(Curve,omega);

            cout << "Printing the solid angle and the scaled knot \n";
            OutputSolidAngle(omega,knot_filename+"_0_4pi");
            for (int n=0; n<Nx*Ny*Nz; n++)
            {
                while(omega[n]>2*M_PI) omega[n] -= 4*M_PI;
                while(omega[n]<-2*M_PI) omega[n]+= 4*M_PI;
            }
            OutputSolidAngle(omega,knot_filename+"_-2pi_2pi");
            OutputScaledKnot(Curve);

            if(Gradient)
            {
                ComputeGradientAllPoints(Curve,Bx,By,Bz);
                OutputGradient(Bx,By,Bz,knot_filename+"_Gradient");
            }

            return 0;
        }
    }
    return 1;
}

double ComputeSolidAngleOnePoint(const Link& Curve, const viewpoint& View)
{
    double totalomega = 0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();


        // the users settings for n_infty
        double ninftyx = Initialninftyx;
        double ninftyy = Initialninftyy;
        double ninftyz = Initialninftyz;
        double mag = sqrt(ninftyx*ninftyx+ninftyy*ninftyy+ninftyz*ninftyz);
        if(std::abs(mag-1)>0.0001)
        {
            cout << "the magnitude of your n_infty vector isnt close to 1 - did you normalise? I've done it for you but best double check its really the direction you want \n";
        }
        ninftyx /= mag;
        ninftyy /= mag;
        ninftyz /= mag;

        // can we proceed to the calculation? Initially the answer is no
        bool thresholdokay = false;
        // how many times have we needed to pick a new vector?
        int numtimesthresholdexceeded = 0;
        while ( thresholdokay == false)
        {
            double ndotnmin = 1.0;
            double ndotnmax = -1.0;
            int smin;
            for (int s=0; s<NP; s++)
            {
                double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
                double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
                double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
                double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
                double ndotninfty = (viewx*ninftyx+viewy*ninftyy+viewz*ninftyz)/dist;
                if (ndotninfty<ndotnmin) {ndotnmin = ndotninfty; smin = s;}
                if (ndotninfty>ndotnmax) {ndotnmax = ndotninfty;}
            }

            thresholdokay = true;
            // is ninfty an okay choice ?
            if (ndotnmin < Threshold_ninf_dot_n)
            {
                // if i sent ninfty -> -ninfty would it be okay?
                if (ndotnmax > -Threshold_ninf_dot_n)
                {
                    // problem case - we need to choose a genuinely new direction. We do so from a list of hardcoded vectors.
                    // in practice I've never hit the end of this list, but if we do, just give up.
                    thresholdokay = false;
                    ninftyx = hardcodedvectors[3*numtimesthresholdexceeded];
                    ninftyy = hardcodedvectors[3*numtimesthresholdexceeded+1];
                    ninftyz = hardcodedvectors[3*numtimesthresholdexceeded+2];

                    numtimesthresholdexceeded ++;
                    if(numtimesthresholdexceeded == numhardcodedvectors) {thresholdokay = true;} // give up jack
                }
                else
                {
                    // flippings fine, do that
                    ninftyx = -ninftyx;
                    ninftyy = -ninftyy;
                    ninftyz = -ninftyz;
                }
            }
        }

        // okay we have an acceptable n_infty - now do the integration. Simple trapezium rule quadrature
        double Integral = 0;
        for (int s=0; s<NP; s++)
        {
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double dist = sqrt(viewx*viewx + viewy*viewy + viewz*viewz);
            double ndotninfty = viewx*ninftyx + viewy*ninftyy + viewz*ninftyz;
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            Integral += (ds/dist)*(ninftyz*(ty*viewx-tx*viewy)+ninftyx*(tz*viewy-ty*viewz)+ninftyy*(tx*viewz-tz*viewx))/(dist + ndotninfty);
        }

        totalomega += Integral;
    }
    while(totalomega>4*M_PI) totalomega -= 4*M_PI;
    while(totalomega<0) totalomega += 4*M_PI;
    return totalomega;
}

void ComputeSolidAngleAllPoints(const Link& Curve,vector<double>& omega)
{

#pragma omp parallel default(none) shared(omega,Curve,Nx,Ny,Nz,cout)
    {
        viewpoint Point;
        int progressbarcounter =0;

#pragma omp for
        for (int i=0; i<Nx; i++)
        {
            for (int j=0; j<Ny; j++)
            {
                for (int k=0; k<Nz; k++)
                {
                    Point.xcoord = x(i);
                    Point.ycoord = y(j);
                    Point.zcoord = z(k);

                    double SolidAngle = ComputeSolidAngleOnePoint(Curve,Point);

                    int n = pt(i,j,k);
                    omega[n]=SolidAngle;
                }
            }
            // progress bar!
#ifdef _OPENMP
            if(omp_get_thread_num() == 0)
            {
                if ((  ( ( ((double)i)*omp_get_num_threads() )/Nx ) * 100 ) > progressbarcounter)
                {
                    cout<< progressbarcounter << "%" << endl;
                    progressbarcounter += 10;
                }
            }
#endif
#ifndef _OPENMP
            if (( ((double)(i)/Nx)*100 )> progressbarcounter )
            {
                cout<< progressbarcounter << "%" << endl;
                progressbarcounter += 10;
            }
#endif
        }
    }

}
void ComputeGradientOnePoint(const Link& Curve, const viewpoint View, double& Bx, double& By, double& Bz)
{
    double internalBx=0;
    double internalBy=0;
    double internalBz=0;
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();

        for (int s=0; s<NP; s++)
        {
            double viewx = Curve.Components[i].knotcurve[s].xcoord - View.xcoord;
            double viewy = Curve.Components[i].knotcurve[s].ycoord - View.ycoord;
            double viewz = Curve.Components[i].knotcurve[s].zcoord - View.zcoord;
            double distcube = sqrt(  (viewx*viewx + viewy*viewy + viewz*viewz)*(viewx*viewx + viewy*viewy + viewz*viewz)*(viewx*viewx + viewy*viewy + viewz*viewz) );
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            internalBx += ((ty*viewz-tz*viewy)/distcube)*ds;
            internalBy += ((tz*viewx-tx*viewz)/distcube)*ds;
            internalBz += ((tx*viewy-ty*viewx)/distcube)*ds;
        }
    }

    Bx = internalBx;
    By = internalBy;
    Bz = internalBz;

    return;
}
void ComputeGradientAllPoints(const Link& Curve,vector<double>& Bx,vector<double>& By,vector<double>& Bz)
{

#pragma omp parallel default(none) shared(Bx,By,Bz,Curve,Nx,Ny,Nz,cout)
    {
        viewpoint Point;
        int progressbarcounter =0;

#pragma omp for
        for (int i=0; i<Nx; i++)
        {
            for (int j=0; j<Ny; j++)
            {
                for (int k=0; k<Nz; k++)
                {
                    Point.xcoord = x(i);
                    Point.ycoord = y(j);
                    Point.zcoord = z(k);

                    double bx,by,bz;
                    ComputeGradientOnePoint(Curve,Point,bx,by,bz);

                    int n = pt(i,j,k);
                    Bx[n]=bx;
                    By[n]=by;
                    Bz[n]=bz;
                }
            }
            // progress bar!
#ifdef _OPENMP
            if(omp_get_thread_num() == 0)
            {
                if ((  ( ( ((double)i)*omp_get_num_threads() )/Nx ) * 100 ) > progressbarcounter)
                {
                    cout<< "gradient " << progressbarcounter << "%" << endl;
                    progressbarcounter += 10;
                }
            }
#endif
#ifndef _OPENMP
            if (( ((double)(i)/Nx)*100 )> progressbarcounter )
            {
                cout<<"gradient " << progressbarcounter << "%" << endl;
                progressbarcounter += 10;
            }
#endif
        }
    }

}

int incp(int i, int p, int N)
{
    if(i+p<0) return (N+i+p);
    else return ((i+p)%N);
}

double x(int i)
{
    return (i+0.5-Nx/2.0)*h;
}

double y(int j)
{
    return (j+0.5-Ny/2.0)*h;
}

double z(int k)
{
    return (k+0.5-Nz/2.0)*h;
}

int pt( int i,  int j,  int k)
{
    return (i*Ny*Nz+j*Nz+k);
}
