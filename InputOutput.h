/*
    Copyright 2018
    Jack Binysh <j.binysh@warwick.ac.uk>,
    Gareth Alexander <g.p.alexadner@warwick.ack.uk>
    This software is provided under a BSD license. See LICENSE.txt.

    all functions involving reading input from files, or writing output to files, are here.
*/

#ifndef INPUTOUTPUT_H
#define INPUTOUTPUT_H
#include <string>
#include <vector>

struct Link;

// reads in the user defined values set in "input.txt".
int InitialiseSystemParameters();
// Initialises the point data read in from the input Link embedding - performs no geometric calculations whatsoever
int InitialiseFromFile(struct Link& Curve);

// A word on output format: All files are output as human (or grad student) readable .vtk files, to be read into e.g. Paraview, https://www.paraview.org/ for visualisation.
// The format is the "Simple Legacy Format" described at https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf

// Output the same curve as read in, but scaled if scaling was asked for, and in .vtk format, so it can be imported to e.g. Paraview
void OutputScaledKnot(Link& Curve);
// Output the grid of omega values computed to a file with user specified filename, in legacy .vtk format
void OutputSolidAngle(const Link& Curve,const std::vector<double>& omega,const std::string filename);

#endif
