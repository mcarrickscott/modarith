
// API header file for Edwards and Weierstrass curves 

#ifndef CURVE_H
#define CURVE_H

// point.h generated from the curve.py script, and required for definition of point structure 
/*** Insert automatically generated point definition point.h here ***/

@point@

/*** End of automatically generated code ***/

// api functions. char* parameters are big-endian integers of fixed length
extern int ecnXXXget(point *P,char *x, char *y);  // extract from point
extern void ecnXXXset(int s,const char *x,const char *y,point *P); // set point
extern void ecnXXXinf(point *P);  // set point-at-infinity 
extern int ecnXXXisinf(point *P); // check for point-at-infinity
extern void ecnXXXneg(point *P);  // negate point
extern void ecnXXXadd(point *Q,point *P); // add Q to P 
extern void ecnXXXsub(point *Q,point *P); // subtract Q from P
extern void ecnXXXdbl(point *P);          // double P
extern void ecnXXXgen(point *P);          // create generator point
extern void ecnXXXran(int r,point *P);    // randomize projective point (side channel noise)
extern void ecnXXXmul(const char *e,point *P); // multiply P by e
extern void ecnXXXmul2(const char *e,point *P,const char *f,point *Q,point *R); // R=eP+fQ
extern int ecnXXXcmp(point *P,point *Q);  // compare points for equality
extern void ecnXXXaffine(point *P);       // convert from projective (x,y,z) to (x,y,1)
extern void ecnXXXcpy(point *Q,point *P); // copy Q to P
extern void ecnXXXcof(point *P);  // multiply point by small curve co-factor

#endif