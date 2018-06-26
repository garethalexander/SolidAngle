#include "Geometry.h"
#include "Constants.h"
#include "Wavefield.h"

void ComputeGeometry(Link& Curve)
{
    ScaleKnot(Curve);
    ComputeLengths(Curve);
    ComputeTangent(Curve);
    ComputeKappaN(Curve);
    cout << "curve has size " << Curve.NumPoints << endl;
    cout << "Number of refinements requested is " << NumRefinements << endl;
    for(int counter=0;counter<NumRefinements;counter++)
    {
        RefineCurve(Curve);
        ComputeLengths(Curve);
        ComputeTangent(Curve);
        ComputeKappaN(Curve);
        cout << "curve has size " << Curve.NumPoints << endl;
    }
}

void ComputeLengths(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double dx = (Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - Curve.Components[i].knotcurve[s].xcoord);
            double dy = (Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - Curve.Components[i].knotcurve[s].ycoord);
            double dz = (Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - Curve.Components[i].knotcurve[s].zcoord);
            double deltas = sqrt(dx*dx+dy*dy+dz*dz);
            Curve.Components[i].knotcurve[s].length = deltas;
        }
    }
}

void ComputeTangent(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double dsp = Curve.Components[i].knotcurve[s].length;
            double dsm = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            Curve.Components[i].knotcurve[s].tx = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].xcoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].xcoord;
            Curve.Components[i].knotcurve[s].ty = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].ycoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].ycoord;
            Curve.Components[i].knotcurve[s].tz = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].zcoord - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].zcoord;
        }
    }
}

void ComputeKappaN(Link& Curve)
{
    for(int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for(int s=0; s<NP; s++)
        {
            double dsp = Curve.Components[i].knotcurve[s].length;
            double dsm = Curve.Components[i].knotcurve[incp(s,-1,NP)].length;
            double kappaNx = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].tx + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].tx - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].tx;
            double kappaNy = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].ty + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].ty - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].ty;
            double kappaNz = (dsm/(dsp*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,1,NP)].tz + ((dsp-dsm)/(dsp*dsm))*Curve.Components[i].knotcurve[s].tz - (dsp/(dsm*(dsp+dsm)))*Curve.Components[i].knotcurve[incp(s,-1,NP)].tz;
            Curve.Components[i].knotcurve[s].kappaNx = kappaNx;
            Curve.Components[i].knotcurve[s].kappaNy = kappaNy;
            Curve.Components[i].knotcurve[s].kappaNz = kappaNz;
        }
    }
}

void RefineCurve(Link& Curve)
{
    Link NewCurve;
    NewCurve.NumPoints= 0;
    NewCurve.NumComponents = NumComponents;
    NewCurve.Components.resize(NumComponents);

    for (int i=0; i<Curve.NumComponents; i++)
    {
        int NP = Curve.Components[i].knotcurve.size();
        for (int s=0; s<NP; s++)
        {
            knotpoint Point;
            // keep old point
            Point.xcoord = Curve.Components[i].knotcurve[s].xcoord;
            Point.ycoord = Curve.Components[i].knotcurve[s].ycoord;
            Point.zcoord = Curve.Components[i].knotcurve[s].zcoord;
            NewCurve.Components[i].knotcurve.push_back(Point);
            // create new point
            double ds = 0.5*Curve.Components[i].knotcurve[s].length;
            double x1 = Curve.Components[i].knotcurve[s].xcoord + ds*Curve.Components[i].knotcurve[s].tx + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNx;
            double x2 = Curve.Components[i].knotcurve[incp(s,1,NP)].xcoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].tx + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNx;
            Point.xcoord = 0.5*(x1+x2);
            double y1 = Curve.Components[i].knotcurve[s].ycoord + ds*Curve.Components[i].knotcurve[s].ty + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNy;
            double y2 = Curve.Components[i].knotcurve[incp(s,1,NP)].ycoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].ty + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNy;
            Point.ycoord = 0.5*(y1+y2);
            double z1 = Curve.Components[i].knotcurve[s].zcoord + ds*Curve.Components[i].knotcurve[s].tz + 0.5*ds*ds*Curve.Components[i].knotcurve[s].kappaNz;
            double z2 = Curve.Components[i].knotcurve[incp(s,1,NP)].zcoord - ds*Curve.Components[i].knotcurve[incp(s,1,NP)].tz + 0.5*ds*ds*Curve.Components[i].knotcurve[incp(s,1,NP)].kappaNz;
            Point.zcoord = 0.5*(z1+z2);
            NewCurve.Components[i].knotcurve.push_back(Point);
        }
        NewCurve.NumPoints += NewCurve.Components[i].knotcurve.size();
    }
    Curve = NewCurve;
}

void ScaleKnot(Link& Curve)
{
    if(Scaling)
    {
        double minxin = Curve.minx;
        double maxxin = Curve.maxx;
        double minyin = Curve.miny;
        double maxyin = Curve.maxy;
        double minzin = Curve.minz;
        double maxzin = Curve.maxz;
        double midpoint[3];
        midpoint[0] = 0.5*(maxxin+minxin);
        midpoint[1] = 0.5*(maxyin+minyin);
        midpoint[2] = 0.5*(maxzin+minzin);
        double scale[3];
        ScaleFunction(scale,maxxin,minxin,maxyin,minyin,maxzin,minzin);
        for(int i=0; i<Curve.NumComponents; i++)
        {
            for(int s=0; s<Curve.Components[i].knotcurve.size(); s++)
            {
                Curve.Components[i].knotcurve[s].xcoord = scale[0]*(Curve.Components[i].knotcurve[s].xcoord - midpoint[0]);
                Curve.Components[i].knotcurve[s].ycoord = scale[1]*(Curve.Components[i].knotcurve[s].ycoord - midpoint[1]);
                Curve.Components[i].knotcurve[s].zcoord = scale[2]*(Curve.Components[i].knotcurve[s].zcoord - midpoint[2]);
            }
        }
    }
}

void ScaleFunction(double *scale, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin)
{
     scale[0] = BoxFractionx/( (maxxin-minxin)/(Nx*h) );
     scale[1] = BoxFractiony/( (maxyin-minyin)/(Ny*h) );
     scale[2] = BoxFractionz/( (maxzin-minzin)/(Nz*h) );
    if(ScaleProportionally)
    {
        double minscale = scale[0];
        int imin=0;
        for(int i=0; i<3; i++)
        {
            if(scale[i] < minscale)
            {
                imin = i;
                minscale = scale[i];
            }
        }
        for(int i=0; i<3; i++){ scale[i] = scale[imin];}
    }
}
