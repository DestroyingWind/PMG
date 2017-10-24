
void exactinit();

double orient2d(double *pa, double *pb, double *pc);
double incircle(double *pa, double *pb, double *pc, double *pd);
double orient3d(double *pa, double *pb, double *pc, double *pd);

// Return 1 if p3 is CCW to the line p1->p2
//       -1 if p3 is CW ...
//        0 if p3 is collinear with p1->p2.
int Orient2d(double *p1, double *p2, double *p3);

// Return 1 if d lies inside the circumcircle of [a,b,c]
//       -1 if d lies outside
//        0 if d lies on.
// Assumption is: Orient2d(a,b,c) > 0
int InCircle(double *pa, double *pb, double *pc, double *pd);

int Orient3d(double *pa, double *pb, double *pc, double *pd);