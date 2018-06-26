#include <math.h>
#include <cmath>
#include <string.h>
#include "Wavefield.h"
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
        // all data read in or computed about the curve itself is stored in this data structure.
        Link Curve;
        if(0==InitialiseFromFile(Curve))
        {
            cout << "Filling in the Geometry of the Input Curve \n";
            ComputeGeometry(Curve);

            cout << "Beginning the Solid Angle Calculations \n";
            ComputeWavefieldAllPoints(Curve,omega);

            cout << "Printing the solid angle and the scaled knot \n";
            OutputWavefield(omega,knot_filename+"_Wavefield_offset");
            OutputScaledKnot(Curve);

            return 0;
        }
    }
    return 1;
}

double ComputeWavefieldOnePoint(const Link& Curve, const viewpoint& View)
{
    double totalomega = 0;
    // the shortest distance between our point and the link - initial value is greater than any possible interbox distance
    double WavefieldMinDistandOffset = 5*sqrt(Nx*h*Nx*h +Ny*h*Ny*h +Nz*h*Nz*h);

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

            // an offset function to modify the basic distance function
            //double offset = 20*sin(3*2*M_PI*Curve.Components[i].knotcurve[s].integratedlength/Curve.length); 
            double offset = 0; 
            double distandoffset = dist+offset; 

            if(distandoffset<WavefieldMinDistandOffset){WavefieldMinDistandOffset=distandoffset;}

            double ndotninfty = viewx*ninftyx + viewy*ninftyy + viewz*ninftyz;
            double tx = Curve.Components[i].knotcurve[s].tx;
            double ty = Curve.Components[i].knotcurve[s].ty;
            double tz = Curve.Components[i].knotcurve[s].tz;
            double ds = 0.5*(Curve.Components[i].knotcurve[s].length+Curve.Components[i].knotcurve[incp(s,-1,NP)].length);
            Integral += (ds/dist)*(ninftyz*(ty*viewx-tx*viewy)+ninftyx*(tz*viewy-ty*viewz)+ninftyy*(tx*viewz-tz*viewx))/(dist + ndotninfty);
        }

        totalomega += Integral;
    }

    double ScrollWavelength = 21.3;
    double ScrollPhase = 0.5*totalomega +(2*M_PI/ScrollWavelength)*WavefieldMinDistandOffset;

    while(ScrollPhase>M_PI) ScrollPhase -= 2*M_PI;
    while(ScrollPhase<-M_PI) ScrollPhase += 2*M_PI;
    return ScrollPhase;
}

void ComputeWavefieldAllPoints(const Link& Curve,vector<double>& omega)
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

                    double phi = ComputeWavefieldOnePoint(Curve,Point);

                    int n = pt(i,j,k);
                    omega[n]=phi;
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
