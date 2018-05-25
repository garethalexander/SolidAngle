/*
    Copyright 2018
    Jack Binysh <j.binysh@warwick.ac.uk>,
    Gareth Alexander <g.p.alexadner@warwick.ack.uk>
    This software is provided under a BSD license. See LICENSE.txt.

    This file contains:
    -the datastuctures used to encode the link, and the viewpoint to look at it from.
    - a list of hardcoded vectors we succesively pick n_infty from if our initial choice was no good
    - the main functions to actually compute the solid angle function.
    - a few little helper math/logic functions

*/

#ifndef SOLIDANGLE_H
#define SOLIDANGLE_H

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
#include <vector>
#include <time.h>

using namespace std;

// A single point in the discretised embedding of the knot, with associated geometry.
struct knotpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
    double tx;
    double ty;
    double tz;
    double length;
    double kappaNx;
    double kappaNy;
    double kappaNz;
};

// A single component of the link, made up of knotpoints.
struct knotcurve
{
    std::vector<knotpoint> knotcurve;
};

// The main data struture - a single instance of this in the code contains all data associated to the knot itself.
// Made up of knotcurves, which in turn are made of knotpoints.
struct Link
{
    std::vector<knotcurve> Components;
    int NumComponents;
    int NumPoints;
    // bounding box
    double minx, maxx;
    double miny, maxy;
    double minz, maxz;
};

// A list of hardcoded directions to succesively try of the user inputted one gives n.n_infty exceeding the threshold.
// This is faster than actually selecting at random and sufficient for our purposes. And we don't have to worry about random numbers + threads = state that differs run to run
// read in threes:
const double hardcodedvectors [] = {1,0,0  ,  0,1,0  ,  0,0,1   ,  0,0.70710678118,0.70710678118  ,  0.70710678118,0,0.70710678118  ,  0, 0.70710678118,0.70710678118};
const int numhardcodedvectors = 6;

struct viewpoint
{
    double xcoord;
    double ycoord;
    double zcoord;
};

// Runs through our 3D grid omega, and computes the solid angle at each point.
void ComputeSolidAngleAllPoints(const Link& Curve,vector<double>& omega);

// All the real work of the code happens in this function - given a viewpoint, give me the solid angle subtended by the knot at that viewpoint.
// this function implements the formula for this solid angle, \omega({\bf x}) = \int \frac{{\bf n}_{\infty} \times {\bf n} }{1+{\bf n}_\infty \cdot {\bf n}} \cdot \mathrm{d}{\bf n} ,
// with a user defined threshold checking whether we are too close to a surface of discontinuity. If we are, we first swap from n_infty -> -n_infty and try again. If that doesn't work
// we pick a new one at random, and keep trying. Its output is standardised to the range 0 to 4 pi
double ComputeSolidAngleOnePoint(const Link& Curve, const viewpoint& View);

// periodic incrementer - increase (or decrease if -ve) i by p, wrapping around at N
int incp(int i, int p, int N);

// the functions x,y,z are used in the code to convert between grid locations and physical space locations.
// the physical space origin is in the middle of the grid, and then theres scaling by h the gridspacing.
double x(int i);
double y(int j);
double z(int k);

// getting the row major index for a 3D grid
int pt(int i,  int j,  int k);

#endif
