/*
    Copyright 2018
    Jack Binysh <j.binysh@warwick.ac.uk>,
    Gareth Alexander <g.p.alexadner@warwick.ack.uk>
    This software is provided under a BSD license. See LICENSE.txt.

    This file contains all functions associated with the geometry of the curve itself:
    scaling and refining the knot if the appropriate flags are set, computing the geometry from the point data - tangents, normals curvature etc.
*/


#ifndef GEOMETRY_H
#define GEOMETRY_H

struct Link;

// Runs scaling and refinement functions if the user has asked for them, then runs all geometric funtions to fill in the links geometry.
void ComputeGeometry(Link& Curve);
// ScaleKnot centers and scales the link to the box size if the user has requested it. ScaleFunction is a helper function for this purpose.
void ScaleKnot(Link& Curve);
void ScaleFunction(double *scale, double maxxin, double minxin, double maxyin, double minyin, double maxzin, double minzin);
// Subdivides each component of the link a user-defined number of times, placing a new point between each two exising points, with simple linear interpolation
void RefineCurve(struct Link& Curve);
void ComputeLengths(struct Link& Curve);
void ComputeTangent(struct Link& Curve);
void ComputeKappaN(struct Link& Curve);

#endif
