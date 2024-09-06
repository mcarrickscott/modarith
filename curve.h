
// API header file for Edwards and Weierstrass curves 

#ifndef CURVE_H
#define CURVE_H

// point.h generated from the curve.py script, and required for definition of point structure 
/*** Insert automatically generated point definition point.h here ***/

@point@

/*** End of automatically generated code ***/

// api functions. char* parameters are big-endian integers of fixed length
extern int ecnget(point *P,char *x, char *y);  // extract from point
extern void ecnset(int s,const char *x,const char *y,point *P); // set point
extern void ecninf(point *P);  // set point-at-infinity 
extern int ecnisinf(point *P); // check for point-at-infinity
extern void ecnneg(point *P);  // negate point
extern void ecnadd(point *Q,point *P); // add Q to P 
extern void ecnsub(point *Q,point *P); // subtract Q from P
extern void ecndbl(point *P);          // double P
extern void ecngen(point *P);          // create generator point
extern void ecnmul(const char *e,point *P); // multiply P by e
extern void ecnmul2(const char *e,point *P,const char *f,point *Q,point *R); // R=eP+fQ
extern int ecncmp(point *P,point *Q);  // compare points for equality
extern void ecnaffine(point *P);       // convert from projective (x,y,z) to (x,y,1)
extern void ecncpy(point *Q,point *P); // copy Q to P
extern void ecncof(point *P);  // multiply point by small curve co-factor

#endif