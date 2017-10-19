#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "detri2.h"
#include "pred2d.h"

#pragma warning (disable:4996)
using  namespace detri2;

const char *vertextypename[] = {
  "UNUSEDVERTEX",     //  0  (default)
  "DUPLICATEDVERTEX", //  1
  "RIDGEVERTEX",      //  2
  "ACUTEVERTEX",      //  3
  "FACETVERTEX",      //  4
  "VOLVERTEX",        //  5
  "FREESEGVERTEX",    //  6
  "FREEFACETVERTEX",  //  7
  "FREEVOLVERTEX",    //  8
  "HIDDENVERTEX",     //  9
  "DEADVERTEX"        // 10
};

void Triangulation::pvert(Vertex *v)
{
  printf("  Vertex: t(0x%lx)", (unsigned long) v);
  if (v->idx == -1) {
    printf(" deleted.\n");
    return;
  }
  printf("  [%g, %g] [%g]", v->crd[0], v->crd[1], v->crd[2]);
  if (weighted && !weight_is_height) { // -w but not -wh
    printf(" h[%g]", v->crd[0]*v->crd[0]+v->crd[1]*v->crd[1]-v->crd[2]);
  }
  printf(" idx(%d)", v->idx);
  int t = v->get_vrttype();
  printf(" (%s)", vertextypename[t]);
  if (v->is_infected()) {
    printf(" (marked)");
  }
  printf("\n");
  if (v->adj != NULL) {
    printf("  adj: (0x%lx) [%d,%d,%d]\n", (unsigned long) v->adj,
           v->adj->vrt[0]->idx, v->adj->vrt[1]->idx, v->adj->vrt[2]->idx);
  } else {
    printf("  NO ADJCENT TRIANGLE.\n");
  }
  if (v->on != NULL) {
    printf("  on: (0x%lx) [%d,%d]\n", (unsigned long) v->on,
           v->on->vrt[0]->idx, v->on->vrt[1]->idx);
  }
}

void Triangulation::pverti(int idx)
{
  // Search the vertex with the given idx.
  Vertex *vrt = NULL;
  int i, c = 0;
  for (i = 0; c < vrts->objects; i++) {
    vrt = (Vertex *) (* vrts)[i];
    if (!vrts->is_deleted(vrt)) {
      if (vrt->idx == idx) break;
      c++;
    }
  }
  if (c == vrts->objects) {
    printf("  NOT EXIST!\n");
    return;
  }
  pvert(vrt);
}

void Triangulation::ptriv(TriEdge &E)
{
  if (E.tri->vrt[0] != NULL) {
    if (E.tri->vrt[2] != NULL) {
      printf("TriEdge: t(0x%lx) v(%d) [%d, %d, %d] tag(%d)",
             (unsigned long) E.tri, E.ver,
             (E.org())->idx, (E.dest())->idx, (E.apex())->idx, E.tri->tag);
      if (E.is_segment()) {
        printf(" (seg)");
      }
      if (E.tri->is_infected()) {
        printf(" (marked)");
      }
      if (E.tri->is_hulltri()) {
        printf(" (hulltri)");
      }
      if (E.tri->is_exterior()) {
        printf(" (extri)");
      }
      printf("\n");
    } else {
      printf("Segment: s(0x%lx) [%d, %d]\n",
             (unsigned long) E.tri, E.tri->vrt[0]->idx, E.tri->vrt[1]->idx);
    }
  } else {
    printf("TriEdge: t(0x%lx), deleted!\n", (unsigned long) E.tri);
  }
}

void Triangulation::ptriv(Triangle *t)
{
  TriEdge E = TriEdge(t, 0);
  ptriv(E);
}

void Triangulation::ptri(TriEdge &E)
{
  if (E.tri->vrt[0] == NULL) {
    printf("TriEdge: t(0x%lx), deleted!\n", (unsigned long) E.tri);
    return;
  }
  if (E.tri->vrt[2] != NULL) {
    printf("TriEdge: t(0x%lx) v(%d) [%d, %d, %d] tag(%d)",
           (unsigned long) E.tri, E.ver,
           (E.org())->idx, (E.dest())->idx, (E.apex())->idx, E.tri->tag);
    if (E.is_segment()) {
      printf(" (seg)");
    }
    if (E.tri->is_infected()) {
      printf(" (marked)");
    }
    if (E.tri->is_hulltri()) {
      printf(" (hulltri)");
    }
    if (E.tri->is_exterior()) {
      printf(" (extri)");
    }
    printf("\n");
  } else {
    printf("Segment: s(0x%lx) [%d, %d]\n",
           (unsigned long) E.tri, E.tri->vrt[0]->idx, E.tri->vrt[1]->idx);
  }
  for (int i = 0; i < 3; i++) {
    if (E.tri->nei[i]) {
      TriEdge N = TriEdge(E.tri->nei[i], E.tri->getnver(i));
      if (N.tri->vrt[0] != NULL) {
        printf("  N[%d]: t(0x%lx) v(%d) [%d, %d, %d]", i,
               (unsigned long) N.tri,  N.ver,
               (N.org())->idx, (N.dest())->idx, (N.apex())->idx);
        if (N.is_segment()) {
          printf(" (seg)");
        }
        if (N.tri->is_infected()) {
          printf(" (marked)");
        }
      } else {
        printf("  N[%d]: t(0x%lx), deleted!", i, (unsigned long) N.tri);
      }
    } else {
      printf("  N[%d]: NOT CONNECTED!", i);
    }
    if (i == E.ver) {
      printf(" (*)");
    }
    printf("\n");
  }
}

void Triangulation::ptri(Triangle *t)
{
  TriEdge E = TriEdge(t, 0);
  ptri(E);
}

TriEdge Triangulation::get_trii(int v1, int v2, int v3)
{
  // Search the triangle with vertices v1, v2, and v3;
  TriEdge E;
  Vertex *v;
  int i, c = 0;
  for (i = 0; c < tris->objects; i++) {
    E.tri = (Triangle *)(* tris)[i];
    if (!tris->is_deleted(E.tri)) {
      for (E.ver = 0; E.ver < 3; E.ver++) {
        v = E.org();
        if (v->idx == v1) {
          break;
        }
      }
      c++;
    }
  }
  if (c == tris->objects) {
    printf("  VERTEX %d NOT EXIST.\n", v1);
  }

  TriEdge N = E;
  do {
    v = N.dest();
    if (v->idx == v2) break;
    N.symself(); N.enextself();
  } while (N.tri != E.tri);
  if (v->idx != v2) {
    printf("  EDGE [%d,%d] NOT EXIST.\n", v1, v2);
  }

  v = N.apex();
  if (v->idx == v3) {
    ptri(N); // Found.
  } else {
    N.symself();
    v = N.apex();
    if (v->idx == v3) {
      ptri(N); // Found.
    } else {
      printf("  TRI [%d,%d,%d] NOT EXIST.\n", v1, v2, v3);
      N.tri = NULL;
    }
  }

  return N;
}

int Triangulation::ptri_orthocenter(int v1, int v2, int v3)
{
  // Make sure the index-to-vertex map is valid.
  build_index_to_vertex_map();

  return get_tri_orthocenter(idx_to_vert_list[v1], idx_to_vert_list[v2],
                             idx_to_vert_list[v3], NULL, NULL, NULL);
}

int Triangulation::pseg_orthocenter(int v1, int v2)
{
  // Make sure the index-to-vertex map is valid.
  build_index_to_vertex_map();

  return get_seg_orthocenter(idx_to_vert_list[v1], idx_to_vert_list[v2],
                             NULL, NULL, NULL);
}

// Check the validity of the mesh.
int Triangulation::check_mesh()
{
  TriEdge E, N, S;
  int horror = 0;
  int ori;
  int i, c = 0;

  for (i = 0; c < tris->objects; i++) {
    E.tri = (Triangle *)(* tris)[i];
    if (!tris->is_deleted(E.tri)) {
      if (!E.tri->is_hulltri()) {
        // Check if this triangle is inverted.
        ori = Orient2d(E.tri->vrt[0]->crd, E.tri->vrt[1]->crd, E.tri->vrt[2]->crd);
        if (ori <= 0) {
          if (!quiet) {
            printf("  !!!! Inverted Triangle [%d,%d,%d]\n", E.tri->vrt[0]->idx,
                   E.tri->vrt[1]->idx, E.tri->vrt[2]->idx);
          }
          horror++;
        }
      }
      // Check if this tri is correctly connected.
      for (E.ver = 0; E.ver < 3; E.ver++) {
        N.tri = E.tri->nei[E.ver];
        if (N.tri != NULL) {
          N.ver = E.tri->getnver(E.ver);
          if (!((E.org() == N.dest()) && (E.dest() == N.org()))) {
            if (!quiet) {
              printf("  !!!! WRONG CONNECTION AT [%d,%d,%d] - [%d,%d,%d]\n",
                     E.org()->idx, E.dest()->idx, E.apex()->idx,
                     N.org()->idx, N.dest()->idx, N.apex()->idx);
            }
            horror++;
          }
        } else {
          if (!quiet) {
            printf("  !! NO CONNECTION at [%d,%d]-%d\n", E.org()->idx,
                   E.dest()->idx, E.apex()->idx);
          }
          horror++;
        }
        if (E.is_segment()) {
          if (!N.is_segment()) {
            if (!quiet) {
				ptri(E); ptri(N); ptri(E.sym());
				printf("\n\n");
				printf("%d,%d   %d,%d  %d %d\n\n\n", E.org()->idx, E.dest()->idx, N.org()->idx, N.dest()->idx, E.is_segment(), N.is_segment());
              printf("  !!!! WRONG TRI->SEG AT [%d,%d] - %d\n",
                       E.org()->idx, E.dest()->idx, E.apex()->idx);
            }
            horror++;
          }
        }
      }
      c++;
    } // if (!tris->is_deleted
  } // i

  if (!quiet) {
    if (horror == 0) {
      printf("The mesh is consistent.\n");
    }
  }

  return horror;
}

// If 'queflag' is set, save locally non-Delaunay edges.
// It asseums that the arraypool 'chktris' was allocated, with entry 'TriEdge'.
int Triangulation::check_delaunay(int queflag)
{
  TriEdge E, N;
  int horror = 0;
  int ori;
  int i, c = 0;

  for (i = 0, c = 0; c < tris->objects; i++) {
    E.tri = (Triangle *)(* tris)[i];
    if (!tris->is_deleted(E.tri)) {
      if (!E.tri->is_hulltri()) {
        for (E.ver = 0; E.ver < 3; E.ver++) {
          N = E.sym();
          if (!N.tri->is_hulltri() && !N.tri->is_infected()) {
            if (!weighted) {
              ori = InCircle( E.org()->crd, E.dest()->crd,
                             E.apex()->crd, N.apex()->crd) * delaunay;
            } else { // weighted
              if (weight_is_height) {
                ori = Orient3d( E.org()->crd, E.dest()->crd,
                               E.apex()->crd, N.apex()->crd) * delaunay;
              } else {
                double *pa = E.org ()->crd;
                double *pb = E.dest()->crd;
                double *pc = E.apex()->crd;
                double *pd = N.apex()->crd;
                double wa = pa[2];  // Bakup
                double wb = pb[2];
                double wc = pc[2];
                double wd = pd[2];
                pa[2] = pa[0]*pa[0]+pa[1]*pa[1]-wa;
                pb[2] = pb[0]*pb[0]+pb[1]*pb[1]-wb;
                pc[2] = pc[0]*pc[0]+pc[1]*pc[1]-wc;
                pd[2] = pd[0]*pd[0]+pd[1]*pd[1]-wd;
                ori = Orient3d(pa, pb, pc, pd) * delaunay;
                pa[2] = wa;  // Restore
                pb[2] = wb;
                pc[2] = wc;
                pd[2] = wd;
              }
            }
            if (ori > 0) {
              if (!quiet) {
                printf("  !!!! NON-LOCALLY DELAUNAY AT [%d,%d,%d]-[%d,%d,%d]\n",
                       E.org()->idx, E.dest()->idx, E.apex()->idx,
                       N.org()->idx, N.dest()->idx, N.apex()->idx);
              }
              if (queflag) {
                // Save this edge.
                TriEdge *paryE = (TriEdge *) chktris->alloc();
                *paryE = E;
              }
              horror++;
            }
          } // if (!N.tri->is_hulltri())
        } // E.ver
        E.tri->set_infect();
      } // if (!E.tri->is_hulltri())
      c++;
    }
  }

  // Clear infect-flag.
  for (i = 0, c = 0; c < tris->objects; i++) {
    E.tri = (Triangle *)(* tris)[i];
    if (!tris->is_deleted(E.tri)) {
      if (!E.tri->is_hulltri()) {
        assert(E.tri->is_infected());
        E.tri->clear_infect();
      } // if (!E.tri->is_hulltri())
      c++;
    }
  }

  if (!quiet) {
    if (horror > 0) {
      printf("!! Found %d edges which are not locally", horror);
    } else {
      printf("The mesh is");
    }
    if (delaunay < 0) {
      printf(" furthest");
    }
    if (weighted) {
      printf(" weighted");
    }
    printf(" Delaunay.\n");
  }

  return horror;
}

void Triangulation::print_tri_array(arraypool *triarray)
{
  for (int i = 0; i < triarray->objects; i++) {
    Triangle **tri = (Triangle **) (* triarray)[i];
    ptriv(*tri);
  }
}

void Triangulation::print_triedge_array(arraypool *tearray)
{
  for (int i = 0; i < tearray->objects; i++) {
    TriEdge *tri = (TriEdge *) (* tearray)[i];
    ptriv(*tri);
  }
}

// Print the list of edges (TriEdges) saved in "edgearray"
void Triangulation::print_edge_array(arraypool *edgearray, const char* fn)
{
  FILE *fout;

  if (fn != NULL) {
    fout = fopen(fn, "w");
    if (fout == NULL) return;
    printf("  Print %ld edges into file %s.\n", edgearray->objects, fn);
  } else {
    fout = stdout;
  }

  for (int i = 0; i < edgearray->objects; i++) {
    TriEdge *tri = (TriEdge *) (* edgearray)[i];
    fprintf(fout, "p:draw_subseg(%d, %d) -- %d  #%d\n", tri->org()->idx,
            tri->dest()->idx, tri->apex()->idx, i);
  }

  if (fn != NULL) {
    fclose(fout);
  }
}

void Triangulation::print_edge_list(TriEdge *elist, int n, const char *fn)
{
  FILE *fout;

  if (fn != NULL) {
    fout = fopen(fn, "w");
    if (fout == NULL) return;
    printf("  Print %d edges into file %s.\n", n, fn);
  } else {
    fout = stdout;
  }

  for (int i = 0; i < n; i++) {
    TriEdge *tri = &(elist[i]);
    fprintf(fout, "p:draw_subseg(%d, %d) -- %d  #%d\n", tri->org()->idx,
            tri->dest()->idx, tri->apex()->idx, i);
  }

  if (fn != NULL) {
    fclose(fout);
  }
}

// Mark exterior (and hole) triangles
int Triangulation::mark_exterior_elements()
{
  arraypool *deadtris = new arraypool(sizeof(Triangle **), 8);
  TriEdge E, N;
  int c, i;

  // Loop through triangles, Mark all hull triangles.
  // Since some triangles may be deleted, use c to count the live triangle.
  for (i = 0, c = 0; c < tris->objects; i++) {
    Triangle* tri = (Triangle *)(* tris)[i];
    if (!tris->is_deleted(tri)) {
      // Look for a hull triangle.
      if (tri->is_hulltri()) {
        tri->set_infect();
        Triangle **pparytri = (Triangle **) deadtris->alloc();
        *pparytri = tri;
      }
      c++;
    }
  }

  // Loop through regions, find and infect a tet in holes.
  Triangle **pparytri;
  for (i = 0; i < regions->objects; i++) {
    Vertex *pt = (Vertex *) (*regions)[i];
    if ((pt->tag == -9999) || (pt->tag == 0)) {
      // Locate the point.
      int loc = locate_point(pt, E);  
      if (loc == 1) { // Inside triangle
        E.tri->set_infect();
        pparytri = (Triangle **) deadtris->alloc();
        *pparytri = E.tri;
      } else if (loc == 2) { // On edge.
        N = E.sym();
        if ((!E.tri->is_hulltri()) && (!N.tri->is_hulltri())) {
          E.tri->set_infect();
          pparytri = (Triangle **) deadtris->alloc();
          *pparytri = E.tri;
        }
      } else if (loc == 3) { // On vertex
        E.tri->set_infect();
        pparytri = (Triangle **) deadtris->alloc();
        *pparytri = E.tri;
      }
    }
  }

  // Find all exterior triangles by depth searching.
  for (i = 0; i < deadtris->objects; i++) {
    Triangle** pparytri = (Triangle **)(* deadtris)[i];
    E.tri = *pparytri;
    for (E.ver = 0; E.ver < 3; E.ver++) {
      if (!E.is_segment()) {
        N = E.sym();
        if (!N.tri->is_infected()) {
          N.tri->set_infect();
          Triangle **pnexttri = (Triangle **) deadtris->alloc();
          *pnexttri = N.tri;
        }
      }
    }
  }

  // Mark all exterior triangles.
  for (i = 0; i < deadtris->objects; i++) {
    Triangle** pparytri = (Triangle **)(* deadtris)[i];
    (*pparytri)->clear_infect();
    if (!(*pparytri)->is_hulltri()) {
      (*pparytri)->set_exterior();
      exterior_triangles++;
    }
  }
  assert((hullsize + exterior_triangles) == deadtris->objects);

  delete deadtris;
  return 1;
}

// The plane is polar to the pole (orthocenter).
int Triangulation::get_lifted_plane(REAL Ux, REAL Uy, REAL U_weight,
                                    REAL Vx, REAL Vy, REAL V_weight,
                                    REAL Wx, REAL Wy, REAL W_weight,
                                    REAL* a, REAL* b, REAL* c, REAL* d)
{
  // Translate weights to heights.
  REAL Uz = Ux*Ux + Uy*Uy - U_weight;
  REAL Vz = Vx*Vx + Vy*Vy - V_weight;
  REAL Wz = Wx*Wx + Wy*Wy - W_weight;

  // Calculate the normal (a,b,c) = A cross B of the plane.
  REAL Ax = Vx - Ux;
  REAL Ay = Vy - Uy;
  REAL Az = Vz - Uz;
  REAL Bx = Wx - Ux;
  REAL By = Wy - Uy;
  REAL Bz = Wz - Uz;

  if (a != NULL) {
    *a =   Ay*Bz - Az*By;
    *b = -(Ax*Bz - Az*Bx);
    *c =   Ax*By - Ay*Bx;
    // Calculate coff. d = -(ax + by + cz)
    *d = -(*a * Ux + *b * Uy + *c * Uz);
    // For debugging.
    REAL t1 = *a * Vx + *b * Vy + *c * Vz + *d;
    REAL t2 = *a * Wx + *b * Wy + *c * Wz + *d;
    assert(fabs(t1 - t2) < 1e-8);
  } else {
    // Only print the results.
    printf(" Plane (a,b,c,d): %g, %g, %g, %g\n", 
           Ay*Bz - Az*By, -(Ax*Bz - Az*Bx), Ax*By - Ay*Bx,
           -((Ay*Bz - Az*By) * Ux + (-(Ax*Bz - Az*Bx)) * Uy + (Ax*By - Ay*Bx) * Uz));
  }

  return 1;
}

int Triangulation::get_tri_lifted_plane(Triangle *tri, REAL* a, REAL* b, REAL* c, REAL* d)
{
  REAL Ux, Uy, U_weight;
  REAL Vx, Vy, V_weight;
  REAL Wx, Wy, W_weight;

  Ux = tri->vrt[0]->crd[0];
  Uy = tri->vrt[0]->crd[1];
  Vx = tri->vrt[1]->crd[0];
  Vy = tri->vrt[1]->crd[1];
  Wx = tri->vrt[2]->crd[0];
  Wy = tri->vrt[2]->crd[1];

  if (weighted) { // -w
    if (weight_is_height) { // -wh
      // Translate heights to weights.
      U_weight = Ux*Ux + Uy*Uy - tri->vrt[0]->crd[2];
      V_weight = Vx*Vx + Vy*Vy - tri->vrt[1]->crd[2];
      W_weight = Wx*Wx + Wy*Wy - tri->vrt[2]->crd[2];
    } else {
      U_weight = tri->vrt[0]->crd[2];
      V_weight = tri->vrt[1]->crd[2];
      W_weight = tri->vrt[2]->crd[2];
    }
  } else {
    U_weight = V_weight = W_weight = 0;
  }

  return get_lifted_plane(Ux, Uy, U_weight,
                          Vx, Vy, V_weight,
                          Wx, Wy, W_weight,
                          a, b, c, d);
}

int Triangulation::get_circumcenter(double* v0, double *v1, double *v2, 
                                    double *cx, double *cy, double *radius)
{
  REAL xdo = v1[0] - v0[0];
  REAL ydo = v1[1] - v0[1];
  REAL xao = v2[0] - v0[0];
  REAL yao = v2[1] - v0[1];
  REAL dodist = xdo * xdo + ydo * ydo;
  REAL aodist = xao * xao + yao * yao;

  REAL denominator = 0.5 / (xdo * yao - xao * ydo); // area.
  REAL dx = (yao * dodist - ydo * aodist) * denominator;
  REAL dy = (xdo * aodist - xao * dodist) * denominator;

  *cx = v0[0] + dx;
  *cy = v0[1] + dy;

  *radius = sqrt(dx * dx + dy * dy);

  return 1;
}

//=============================================================================
/*
void Triangulation::get_dual_diagram()
{
  if (dual_verts != NULL) {
    delete [] dual_verts; 
  }

  int nv = (int) tris->objects;
  dual_verts = new Vertex[nv];

  Vertex *e1, *e2;
  double cx, cy, r;
  int i, idx = 1;
  int count = 0;

  // Loop through all triangles (including hull triangles).
  for (i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      if (!tri->is_hulltri()) {
        get_circumcenter(tri->vrt[0]->crd, tri->vrt[1]->crd, tri->vrt[2]->crd,
                         &cx, &cy, &r);
      } else {
        if (tri->vrt[0] == infvrt) {
          e1 = tri->vrt[1]; e2 = tri->vrt[2];
        } else if (tri->vrt[1] == infvrt) {
          e1 = tri->vrt[2]; e2 = tri->vrt[0];
        } else { // 
          e1 = tri->vrt[0]; e2 = tri->vrt[1];
        }
        // Calucate the intersection point at the bounding circle.
        cx = 0.5 * (e1->crd[0] + e2->crd[0]);
        cy = 0.5 * (e1->crd[1] + e2->crd[1]);
      }
      dual_verts[count].crd[0] = cx;
      dual_verts[count].crd[1] = cy;
      dual_verts[count].crd[2] = cx * cx + cy * cy;
      dual_verts[count].idx = idx;
      dual_verts[count].adj = tri;
      tri->idx = count;
      count++; idx++;
    }
  }

  TriEdge tt[4];
  int t_count = 0;

  // Loop through all Delaunay vertices, calculate dual cells.
  count = 0;
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = vrts->get(i);
    if (!vrts->is_deleted(vrt)) {
      tt[0] = vrt->get_triorg();
      // Calculate the edge lengths and area of the cell.
      t_count = 1;
      tt[1] = tt[0];
      while (1) {
        tt[2] = tt[1].enext2().sym();
        if (!tt[1].tri->is_hulltri() || !tt[2].tri->is_hulltri()) {
          // length
          double x1 = dual_verts[tt[0].tri->idx].crd[0];
          double y1 = dual_verts[tt[0].tri->idx].crd[1];
          double x2 = dual_verts[tt[1].tri->idx].crd[0];
          double y2 = dual_verts[tt[1].tri->idx].crd[1];
          double len = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)); 
          // area of triangle
          double x3 = vrt->crd[0];
          double y3 = vrt->crd[1];
          //double area = ;
        } else {
          // Calucalte the area outside the convex hull.
        }
        if (tt[2].tri == tt[0].tri) break; // End of loop.
        tt[1] = tt[2];
        t_count++;
      }

      count++;
    }
  }

}
*/

// Save the orthocircles and lifted planes
void Triangulation::save_pole_and_polar(int meshidx)
{
  // Computer the set of dual (Voronoi) vertices.
  // A (non hull) triangle is dual to a vertex (its circumcenter).
  // A hull triangle is dual to a vertex on its hull edge.
  Vertex *tmp_dual_verts = new Vertex[tris->objects];
  REAL *lifted_planes = new REAL[tris->objects * 4];

  Vertex *e1, *e2;
  double Cx, Cy, radius2;
  double *P, A, B, C, D;
  int i, idx = 1; // UCD assumes vertex index starts from 1. 
  int count = 0;

  // Calculate dual vertices (orthocenters)
  for (i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      if (!tri->is_hulltri()) {
        get_tri_orthocenter(tri->vrt[0], tri->vrt[1], tri->vrt[2],
                            &Cx, &Cy, &radius2);
        get_tri_lifted_plane(tri, &A, &B, &C, &D);
      } else {
        if (tri->vrt[0] == infvrt) {
          e1 = tri->vrt[1]; e2 = tri->vrt[2];
        } else if (tri->vrt[1] == infvrt) {
          e1 = tri->vrt[2]; e2 = tri->vrt[0];
        } else { // 
          e1 = tri->vrt[0]; e2 = tri->vrt[1];
        }
        get_seg_orthocenter(e1, e2, &Cx, &Cy, &radius2);
        A = B = C = D = 0;
      }
      tmp_dual_verts[count].crd[0] = Cx;
      tmp_dual_verts[count].crd[1] = Cy;
      tmp_dual_verts[count].crd[2] = radius2; // w (weight, NOT height)
      tmp_dual_verts[count].idx = idx;
      P = &(lifted_planes[count]);
      P[0] = A;
      P[1] = B;
      P[2] = C;
      P[3] = D;
      tri->idx = count;
      count++; idx++;
    }
  }

  char filename[256];
  sprintf(filename, "%s_%d.txt", outfilename, meshidx);
  FILE *outfile = fopen(filename, "w");
  int nt = (int) vrts->objects;
  int nv = (int) tris->objects;
  int ne = (int) tris->objects * 3 / 2; // 2e = 3f
  printf("Saving poles and polars to file %s.\n", filename);

  fprintf(outfile, "# Size of the dual diagram: #cells=%d, #edges=%d #verts=%d\n", nv, ne, nt);
  fprintf(outfile, "#   Pole  (orthocircle) : Cx, Cy, Z(height), weight, radius\n");
  fprintf(outfile, "#   Polar (lifted plane): Ax+By+Cz+D: A, B, C, D\n");

  REAL x, y, h, w, r;
  int j;
  count = 0; idx = 1;
  for (i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      if (!tri->is_hulltri()) {
        fprintf(outfile, "Tri [%d]: %d, %d, %d\n", idx, 
                tri->vrt[0]->idx, tri->vrt[1]->idx, tri->vrt[2]->idx);
        for (j = 0; j < 3; j++) {
          x = tri->vrt[j]->crd[0];
          y = tri->vrt[j]->crd[1];
          if (weighted) { // -w
            if (weight_is_height) { // -wh
              h = tri->vrt[j]->crd[2];
              w = x*x + y*y - h;
            } else {
              h = x*x + y*y - tri->vrt[j]->crd[2];
              w = tri->vrt[j]->crd[2];
            }
            r = sqrt(fabs(w));
          } else {
            h = x*x + y*y;
            w = 0;
            r = 0;
          }
          fprintf(outfile, "  V[%d]: %g %g %g w(%g) r(%g)\n", j, x, y, h, w, r);
        }
        // Orthocircle (pole).
        x = tmp_dual_verts[count].crd[0];
        y = tmp_dual_verts[count].crd[1];
        w = tmp_dual_verts[count].crd[2]; // radius2
        h = x*x + y*y - w;
        r = sqrt(fabs(w));
        fprintf(outfile, "  Pole:  %g, %g, %g, w(%g) r(%g)\n", x, y, h, w, r);
        // Lifted plane (polar).
        P = &(lifted_planes[count]);
        fprintf(outfile, "  Polar: %g, %g, %g, %g\n", P[0], P[1], P[2], P[3]);
        fprintf(outfile, "  Evaluation: %g\n", P[0]*x + P[1]*y + P[2]*h + P[3]);
      } else {
        if (tri->vrt[0] == infvrt) {
          e1 = tri->vrt[1]; e2 = tri->vrt[2];
        } else if (tri->vrt[1] == infvrt) {
          e1 = tri->vrt[2]; e2 = tri->vrt[0];
        } else { // 
          e1 = tri->vrt[0]; e2 = tri->vrt[1];
        }
        fprintf(outfile, "Tri [%d]: %d, %d, -1\n", idx, e1->idx, e2->idx);
        // Orthocircle (pole).
        x = tmp_dual_verts[count].crd[0];
        y = tmp_dual_verts[count].crd[1];
        w = tmp_dual_verts[count].crd[2];
        h = x*x + y*y - w;
        r = sqrt(fabs(w));
        fprintf(outfile, "  Bisector:  %g, %g, %g, w(%g) r(%g)\n", x, y, h, w, r);
      }
      count++; idx++;
    }
  }

  fclose(outfile);

  delete [] tmp_dual_verts;
  delete [] lifted_planes;
}

//=============================================================================

// http://mathworld.wolfram.com/Circle-LineIntersection.html
// http://stackoverflow.com/questions/23016676/line-segment-and-circle-intersection

// If return 2, the two intersecting points are returned in pt1, and pt2,
//   where pt1 is more close to e1.
int line_circle_inter(double cx, double cy, double radius, 
                      double* e1, double *e2,
                      double* pt1, double* pt2)
{
  double dx, dy, A, B, C, det, t;

  dx = e2[0] - e1[0];
  dy = e2[1] - e1[1];

  A = dx * dx + dy * dy;
  B = 2 * (dx * (e1[0] - cx) + dy * (e1[1] - cy));
  C = (e1[0] - cx) * (e1[0] - cx) + 
      (e1[1] - cy) * (e1[1] - cy) - radius * radius;

  det = B * B - 4 * A * C;
  if ((A <= 0.0000001) || (det < 0)) {
    // No real solutions.
    return 0;
  } else if (det == 0) {
    // One solution.
    t = -B / (2 * A);
    pt1[0] = e1[0] + t * dx;
    pt1[1] = e1[1] + t * dy;
    return 1;
  } else {
    // Two solutions.
    t = (double)((-B + sqrt(det)) / (2 * A));
    pt1[0] = e1[0] + t * dx;
    pt1[1] = e1[1] + t * dy;
    t = (double)((-B - sqrt(det)) / (2 * A));
    pt2[0] = e1[0] + t * dx;
    pt2[1] = e1[1] + t * dy;
    return 2;
  }
        
  return 0;
}

int Triangulation::recover_edge_conforming(Vertex *e1, Vertex *e2, TriEdge &E) 
{
  TriEdge tt[2];
  tt[0].tri = e1->adj;
  for (tt[0].ver = 0; tt[0].ver < 3; tt[0].ver++) {
    if (tt[0].org() == e1) break;
  }
  assert(tt[0].ver < 3);
  assert(lawsonflag > 0); // Use Lawsonflip after inserting vertex.

  if (verbose > 2) {
    printf("    Recovering edge [%d, %d]\n", e1->idx, e2->idx);
  }

  int loc = find_direction(tt[0], e2);

  if (loc == 1) {
    // Flip the edge [b,c]
    //tt[0].enextself();
    //tt[1] = tt[0].sym();
    //flip22(tt);
    // Split the edge [e1,e2]
    double cx, cy, radius, pt1[3], pt2[3];
    cx = cy = radius = 0.0;
    get_circumcenter(tt[0].tri->vrt[0]->crd, tt[0].tri->vrt[1]->crd,
                     tt[0].tri->vrt[2]->crd, &cx, &cy, &radius);
    int chk = line_circle_inter(cx, cy, radius, e1->crd, e2->crd, pt1, pt2);
    assert(chk == 2);

    Vertex *newvrt = vrts->alloc_vrt();
    newvrt->crd[0] = pt1[0];
    newvrt->crd[1] = pt1[1];
    newvrt->idx = firstindex + (vrts->objects - 1);
    unused_vertices++;
  
    loc = insert_point(newvrt, tt, 0); // loc = 0, Locate point first.
    assert(loc > 0);
    //st_segref_count++;
  
    // Continue to recover edge (starting from pt2).
    return recover_edge_conforming(newvrt, e2, E);
  } else if (loc == 2) {
    if (tt[0].dest() == e2) {
      // Found the edge.
      E = tt[0];
      return 1;
    } else {
      assert(0); // Not handled yet.
    }
  }

  return 0;
}

int Triangulation::recover_segments_conforming()
{
  TriEdge E;

  if (verbose) {
    printf("Incrementally recovering segments.\n");
  }
  lawsonflag = 1;

  for (int i = 0; i < segs->objects; i++) {
    Triangle *paryseg = segs->get(i);
    if (verbose > 1) {
      printf("  Recovering segment [%d, %d]\n", paryseg->vrt[0]->idx,
             paryseg->vrt[1]->idx);
    }
    recover_edge_conforming(paryseg->vrt[0], paryseg->vrt[1], E);
  }

  lawsonflag = 0;

  return 1;
}

///////////////////////////////////////////////////////////////////////////////

// The furthest Delaunay triangulation T of a point set S is the triangulation
//   of the sites of the convex hull of S, such that the circumcircle of any
//   triangle in T has no vertex of S in its exterior.
// NOTE: different to the Delaunay triangulation of S, the furthest DT is
//   only defined for the convex hull points of S (a subset of S).
//
// Useful links:
//   LEDA: http://www.algorithmic-solutions.info/leda_guide/geo_algs/delaunay_triangulations.html
//
//   Book: LEDA: A Platform for Combinatorial and Geometric Computing,
//     by Kurt Mehlhorn,Stefan Näher.
//
//     (page 675) ... For points in convex hull position, there is also a
//      so-called the furthest Delaunay triangulation...
//      ... The (Lawson) flipping algorithm can also be used to construct the
//      furthest site Delaunay triangulation. We start with an arbitrary
//      triangulation of a set of points in convex position and flip as long as
//      the triangulation is not furthest site Delaunay. Of course, this time,
//      we flip the diagonal when the last vertex is outside the circumcircle
//      of the first three vertices.
//
//      (page 680) ... There are O(n \log n) algorothm for constructing the
//      (nearest site) Delaunay triangulation, i.e., incremental and
//      divide-and-conquer, etc.
//      ... For the furthest site DT, we only have the flip algorithm.
//      (MY COMMENT) Another algorithm is to use the convex hull construction
//      (see qhull below).
//
//   qhull: http://www.qhull.org/html/qdelau_f.htm
//     As with the Delaunay triangulation, Qhull computes the furthest-site
//     Delaunay triangulation by lifting the input sites to a paraboloid.
//     The lower facets correspond to the Delaunay triangulation while the upper
//     facets correspond to the furthest-site triangulation. Neither triangulation
//     includes "vertical" facets (i.e., facets whose last hyperplane coefficient
//     is nearly zero). Vertical facets correspond to input sites that are
//     coplanar to the convex hull of the input.
//     An example is points on the boundary of a lattice.

// Release the memory.
void Triangulation::cyclist_free()
{
  if (cyclist != NULL) {
    for (int i = 0; i < cyclist->objects; i++) {
      CycleEdges *parycyc = (CycleEdges *) (*cyclist)[i];
      delete [] parycyc->cycedges;
    }
    delete cyclist;
    cyclist = NULL;
  }
}

void Triangulation::get_cycles()
{
  if (cyclist) {
    cyclist_free();
  }

  cyclist = new arraypool(sizeof(CycleEdges), 2);

  TriEdge nextlist[100], workcycles[1000];
  TriEdge tt[4], startE, spinE, nextE;
  int nlen, clen, cidx;
  int flipflag;
  int i, j;

  int bak_lawsonflag = lawsonflag;
  lawsonflag = 1;

  // Used to save the list of locally non-Delaunay edges.
  assert(chktris == NULL);
  chktris = new arraypool(sizeof(TriEdge), 8);
  check_delaunay(1);

  printf("-- Get cycles from %ld locally non-Delaunay edges.\n", chktris->objects);
  if (chktris->objects > 0) {
    print_edge_array(chktris, "dump_nl_edges.lua");
  }

  for (i = 0; i < chktris->objects; i++) {
    TriEdge *paryE = (TriEdge *) (*chktris)[i];
    tt[0] = *paryE;
    if (flip(tt, flipflag, 0)) {
      //assert(0); // flippable?
      continue; // Skip a flippable face.
    }
    // It must be NOT locally (weighted) Delaunay.
    assert((flipflag & 2) == 0);

    // tt[0].dest() is the locally non-convex vertex.
    if (tt[0].is_edge_infected()) {
      // tt[0] already exists in a cycle.
      continue; // skip it.
    }

    if (verbose > 1) {
      printf("-- Start searching cycle from edge [%d,%d]\n",
             tt[0].org()->idx, tt[0].dest()->idx);
    }

    // Put it into workcycles.
    workcycles[0] = tt[0];
    clen = 1; cidx = 0;

    while (1) {
      if (verbose > 1) {
        printf("  p:draw_edge(%d,%d) -- %d\n", workcycles[clen-1].org()->idx,
               workcycles[clen-1].dest()->idx, clen - 1);
      }
      // workcycles[clen-1].dest() is the locally non-convex vertex.
      // Get all edges at this vertex.
      Vertex *pb = workcycles[clen-1].dest(); // The locally non-convex vertex.
      startE = workcycles[clen-1].enext();
      assert(startE.org() == pb);
      // Rotate CCW direction.
      nlen = 0;
      spinE = startE;
      while (1) {
        tt[0] = spinE;
        if (!flip(tt, flipflag, 0)) {
          // It is not flippable.
          if (((flipflag & 2) == 0) &&  // It is not locally (weighted) Delaunay
              (tt[0].dest() != pb)) { // and the non-convex vertex is not pb.
            nextlist[nlen] = tt[0];
            nlen++;
          }
        }
        // Go to the next edge.
        spinE.symself();
        spinE.enextself();
        assert(spinE.org() == pb);
        if (spinE.tri == startE.tri) {
          break; // The inner while (1) loop.
        }
      } // while (1)
      assert(nlen >= 1);

      // Randomly pick an edge to continue.
      if (nlen > 1) {
        if (verbose > 1) {
          for (j = 0; j < nlen; j++) {
            printf("  -- [%d] [%d,%d]\n", j+1, nextlist[j].org()->idx,
                   nextlist[j].dest()->idx);
          }
        }
        int rndix = (rand() % nlen); // tndidx is in [0,..,nlen-1]
        nextE = nextlist[rndix];
      } else {
        nextE = nextlist[0];
      }

      if (verbose > 1) {
        printf("  -- next [%d,%d]\n", nextE.org()->idx, nextE.dest()->idx);
      }

      if (nextE.is_edge_infected()) {
        // Found an edge in cycle. No need to conitnue.
        break;
      }

      // Check if this edge forms a cycle in workcycles,
      // i.e., check if nextE.dest() appears as any workcycles[j].org()
      Vertex *pc = nextE.dest();
      for (j = 0; j < clen; j++) {
        if (workcycles[j].org() == pc) {
          workcycles[clen] = nextE;
          clen++;
          cidx = j;
          break;
        }
      }
      if (j < clen) {
        // Found a new cycle, from workcycles[cidx] to workcycles[clen-1].
        // Save it.
        CycleEdges *parycyc = (CycleEdges *) cyclist->alloc();
        parycyc->cyclen = clen - cidx;
        parycyc->cycedges = new TriEdge[parycyc->cyclen];
        for (j = 0; j < parycyc->cyclen; j++) {
          parycyc->cycedges[j] = workcycles[cidx + j];
          parycyc->cycedges[j].set_edge_infect(); // Mark this edge in cycle.
        }
        break;
      } else {
        // Add an edge to the workcycles.
        workcycles[clen] = nextE;
        clen++;
        assert(clen < 100); // Overflow or in an endless loop.
      }
    } // while (1)

    if (verbose > 1) {
      printf("-- End searching cycle from edge [%d,%d] clen = %d\n",
             workcycles[0].org()->idx,
             workcycles[0].dest()->idx, clen);
    }
  } // i

  // Unmark all cycle edges.
  for (i = 0; i < cyclist->objects; i++) {
    CycleEdges *parycyc = (CycleEdges *)(* cyclist)[i];
    for (j = 0; j < parycyc->cyclen; j++) {
      assert(parycyc->cycedges[j].is_edge_infected());
      parycyc->cycedges[j].clear_edge_infect();
    }
  }
// \iffalse
  // Debug: make sure no edge is infected.
  {
    int c = 0;
    for (i = 0; c < tris->objects; i++) {
      tt[0].tri = tris->get(i);
      if (!tris->is_deleted(tt[0].tri)) {
        for (tt[0].ver = 0; tt[0].ver < 3; tt[0].ver++) {
          assert(!tt[0].is_edge_infected());
        }
        c++;
      }
    }
  }
// \fi

  printf("Found %ld cycles.\n", cyclist->objects);
  for (i = 0; i < cyclist->objects; i++) {
    CycleEdges *parycyc = (CycleEdges *) (*cyclist)[i];
    printf("--  C[%d] (%d edges):\n", i, parycyc->cyclen);
    for (j = 0; j < parycyc->cyclen; j++) {
      printf("p:draw_subseg(%d,%d) -- %d\n",
             parycyc->cycedges[j].org()->idx, parycyc->cycedges[j].dest()->idx,
             parycyc->cycedges[j].apex()->idx);
    }
  }

  lawsonflag = bak_lawsonflag;

  delete chktris;
  chktris = NULL;
}

// Constructing the furthest Delaunay triangulation using Lawson flips
int Triangulation::lawson_flip_furthest_delaunay()
{
  assert(flip_queue->objects == 0l);
  int bak_delaunay = delaunay;
  int bak_quiet = quiet;
  int bak_vertcies = unused_vertices;
  int nlcount = 0, nflips = 0, nvert = 0;
  int it = 0;

  quiet = 1;

  if (dump_to_ucd) { // -Gu
    save_to_ucd(totalflipcount, 0); // no_reindex
  }

  while (1) {
    printf("iteration: %d\n", it);
    nvert = unused_vertices;

    save_polyhedron(it);

    nflips = lawson_flip(NULL); it++;

    printf("  number of flips  = %d\n", nflips);
    printf("  removed vertices = %d\n", unused_vertices - nvert);

    if (dump_to_ucd) { // -Gu
      save_to_ucd(totalflipcount, 0); // no_reindex
    }

    nlcount = check_delaunay(0);

    if (nlcount > 0) {
      printf("  unflippable edges = %d\n", nlcount);
      save_polyhedron(it);

      delaunay = -delaunay;
      nvert = unused_vertices;

      nflips = lawson_flip(NULL); it++;

      printf("  number of reversed flips  = %d\n", nflips);
      printf("  removed vertices = %d\n", unused_vertices - nvert);
    
      if (dump_to_ucd) { // -Gu
        save_to_ucd(totalflipcount, 0); // no_reindex
      }

      delaunay = -delaunay;
    } else {
      break;
    }
  }

  printf("Finished after %d iterations\n", it);
  printf("  Removed %d vertices.\n", unused_vertices - bak_vertcies);

  delaunay = bak_delaunay;
  quiet = bak_quiet;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// [2016-09-21, 21:49]
int Triangulation::lawson_flip_furthest_delaunay_local()
{
  assert(flip_queue->objects == 0l);
  //int bak_quiet = quiet;
  int bak_vertcies = unused_vertices;
  int nlcount = 0, nflips = 0, nvert = 0;
  int it = 0;

  quiet = 1;

  //save_polyhedron(it);
  //save_triangulation(0, 1, it); // reindex(0), notag(1).

  while (1) {
    printf("iteration: %d\n", it);
    nvert = unused_vertices;

    nflips = lawson_flip(NULL); it++;
    printf("  number of flips  = %d\n", nflips);
    printf("  removed vertices = %d\n", unused_vertices - nvert);

    nlcount = check_delaunay(0);

    if (nlcount > 0) {
      printf("  %d locally non-Delaunay edges %s.\n", nlcount,
             delaunay < 0 ? "(furthest)" : "");

      //save_polyhedron(it);
      //save_triangulation(0, 1, it); // reindex(0), notag(1).

      get_cycles();

      // Randomly select a cycle.
      int rndcycle = rand() % (int) cyclist->objects;
      // Randomly select a vertex.
      CycleEdges *parycyc = (CycleEdges *) (*cyclist)[rndcycle];
      int rndvert = rand() % parycyc->cyclen;
      Vertex *delpt = parycyc->cycedges[rndvert].dest();
      TriEdge tt[4];

      printf("  Remove vertex %d\n", delpt->idx);
      int bak_lawsonflag = lawsonflag;
      lawsonflag = 0;
      remove_point(delpt, tt);
      lawsonflag = bak_lawsonflag;
    } // if (nlcount > 0)
    else {
      break; // Done!
    }
  } // while (1)

  printf("Finished after %d iterations\n", it);
  printf("  Removed %d vertices.\n", unused_vertices - bak_vertcies);

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// [2016-11-04, 10:39, Moscow]

int Triangulation::lawson_flip_furthest_delaunay_new()
{
  assert(flip_queue->objects == 0l);
  TriEdge tt[4], *qedge;
  int flipflag, flippable, flipped;
  int i;

  //int bak_quiet = quiet;
  //int bak_vertcies = unused_vertices;
  int bak_flip31count;
  int bak_lawsonflag = lawsonflag;
  lawsonflag = 1; // Use Delaunay check to determine flips.

  chktris = new arraypool(sizeof(TriEdge), 8);

  // Get all locally non-Delaunay edges.  
  check_delaunay(1);

  int success = 0;

  while (chktris->objects > 0) {
    // Push all queued edges, and try to flip them.
    for (i = 0; i < chktris->objects; i++) {
      tt[0] = * (TriEdge *) (* chktris)[i];
      qedge = (TriEdge *) flip_queue->alloc();
      *qedge = tt[0];
    }
    //print_edge_array(chktris, "dump_nl_edges.lua");
    chktris->restart();

    if (verbose > 2) {
      printf("    Lawson flip %ld edges.\n", flip_queue->objects);
    }
    bak_flip31count = flip31count;

    if (0) {
      save_polyhedron(1);
      save_triangulation(0, 1, 1); // reindex(0), notag(1).
    }

    for (i = 0; i < flip_queue->objects; i++) {
      tt[0] = * (TriEdge *) (* flip_queue)[i];
      if (tris->is_deleted(tt[0].tri)) continue;
      // Check whether this edge is locally Delaunay or not.
      flippable = flip(tt, flipflag, 1);
      if (((flipflag & 2) == 0) &&
          ((flipflag & 4) == 0)) { // It is not a constraint edge.
        // It is not locally Delaunay.
        flipped = 0;
        if (flippable) {
          // It is also flippable. 
          // Only flip it if it is a 3-to-1 flip.
          if (flipflag & 2048) {
            Vertex *pb = tt[0].apex();
            flip31(tt);
            pb->set_vrttype(0); // UNUSEDVERTEX
            unused_vertices++;
            flipped = 1;
          }
        }
        if (!flipped) {
          // Queue this edge for later flipping.
          qedge = (TriEdge *) chktris->alloc();
          *qedge = tt[0];
        }
      } // if ((flipflag & 2) == 0)
    } // i
    flip_queue->restart();

    if (chktris->objects == 0) {
      // The mesh is Delaunay!
      success = 1;
      break; // Done!
    }

    // Try 3-to-1 flip first.
    if (flip31count > bak_flip31count) {
      if (verbose > 2) {
        printf("    Removed %d locally non-Delaunay vertices %s.\n", 
               (int) flip31count - bak_flip31count,
               delaunay < 0 ? "(furthest)" : "");
      }
      continue;
    }

    // Push all queued edges, and try to flip them.
    for (i = 0; i < chktris->objects; i++) {
      tt[0] = * (TriEdge *) (* chktris)[i];
      qedge = (TriEdge *) flip_queue->alloc();
      *qedge = tt[0];
    }
    chktris->restart();

    // Try a 2-to-2 flip (an edge).
    for (i = 0; i < flip_queue->objects; i++) {
      tt[0] = * (TriEdge *) (* flip_queue)[i];
      if (tris->is_deleted(tt[0].tri)) continue;
      // Check whether this edge is locally Delaunay or not.
      flippable = flip(tt, flipflag, 1);
      if ((flipflag & 2) == 0) {
        // It is not locally Delaunay.
        flipped = 0;
        if (flippable) {
          // It is also flippable. 
          // It must be a 2-to-2 flip.
          if ((flipflag & 512) != 0) {
            flip22(tt);
            flipped = 1;
            i++; // Skip the current edge.
            break;
          } else {
            assert(0);
          }
        }
        if (!flipped) {
          // Queue this edge for later flipping.
          qedge = (TriEdge *) chktris->alloc();
          *qedge = tt[0];
        }
      } // if ((flipflag & 2) == 0)
    } // i
    if (i == flip_queue->objects) {
      // No edge is flippable.
      printf("!!DEBUG ME!!\n");
      //save_polyhedron(1);
      //save_triangulation(0, 1, 1); // reindex(0), notag(1).
      break; // Failed.
    }
    // Put the rest of the edges into queue.
    for (; i < flip_queue->objects; i++) {
      tt[0] = * (TriEdge *) (* flip_queue)[i];
      if (tris->is_deleted(tt[0].tri)) continue;
      qedge = (TriEdge *) chktris->alloc();
      *qedge = tt[0];
    }
    flip_queue->restart();

  } // while (flip_queue->objects > 0)

  if (!success) {
    printf("  Failed!!.\n");
  }

  lawsonflag = bak_lawsonflag;

  delete chktris;
  chktris = NULL;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////
// [2016-11-13, 21:28, Berlin]
// Removing a vertex by dynamically changing its weight.
// - Only work for weighted DT.

REAL Triangulation::get_flip_time(Vertex *a, Vertex *b, Vertex *c, Vertex *d)
{
  // pa is the deleting vertex, its weight is changing as
  //   new_a_weight = a_x^2 + a_y^2 - a_weight + t,
  // where "a_weight" is the original weight of pa, and "t" is the changing
  //   amount of the weight (the flip time). We use the '+' sign means we
  //   want to move the lifted pa' upwards. Note that 't' may still be negative.
  // The flip time (t) is when the four lifted points coplanar, which is
  //   t = orient3d(pa', pb', pc', pd') / orient2d(pb, pc, pd);
  REAL o3d, o2d;

  if (!weighted) {
    o3d = incircle(a->crd, b->crd, c->crd, d->crd);
  } else { // -w
    if (weight_is_height) { // -wh
      o3d = orient3d(a->crd, b->crd, c->crd, d->crd);
    } else {
      double *pa = a->crd;
      double *pb = b->crd;
      double *pc = c->crd;
      double *pd = d->crd;
      double wa = pa[2];  // Bakup
      double wb = pb[2];
      double wc = pc[2];
      double wd = pd[2];
      pa[2] = pa[0]*pa[0]+pa[1]*pa[1]-wa;
      pb[2] = pb[0]*pb[0]+pb[1]*pb[1]-wb;
      pc[2] = pc[0]*pc[0]+pc[1]*pc[1]-wc;
      pd[2] = pd[0]*pd[0]+pd[1]*pd[1]-wd;
      o3d = orient3d(pa, pb, pc, pd);
      pa[2] = wa;  // Restore
      pb[2] = wb;
      pc[2] = wc;
      pd[2] = wd;
    }
  }

  o2d = orient2d(b->crd, c->crd, d->crd);

  if (fabs(o2d) == 0) {
    // This is a possible 4-to-4 flip.
    assert(0); // Debug
  }
  if (fabs(o2d) < 1e-6) {
    if (fabs(o3d) > 1) {
      // A very strange time.
      assert(0); // Debug
    }
  }

  REAL t = -o3d / o2d;

  return t;
}

// "elist" is an array of TriEdges which are in the star of "pt". Moreover,
//   these edges are arranged in a cycle with CCW orientation. "n" is the
//   the current length of "elist".
// By moving the lifted vertex of "pt" upwards (or downwards), the set of
//   edges are ordered by their "flip-times", which are the time when four
//   lifted vertices (one is lifted "pt") become co-planar. The vertex removal
//   algorithm at each time flip the edge in "elist" which has the shortest
//   (or largest) flip time.
// The returned "n" value (integer) indicates the current length of "elist".
//   Note if the vertex is removed, "elist" only contains one triangle which
//   is the last triangle resuled by a 3-to-1 flip, and the value is '1'.
//   If the return value is larger than '3', it means this algorithm fails.
// Return 1 if the vertex is removed, otherwise, return 0.
// [NOTE]:
//   - This function assumes that "lawsonflag" is 0.
int Triangulation::remove_point_regular(Vertex *pt, TriEdge *elist, int &n)
{
  assert(lawsonflag == 0);
  TriEdge tt[4];
  int i;

  if (elist == NULL) {
    // Collect all edges at this vertex.
    n = 0;
    // Count the number of triagles (edges) at pt.
    tt[0] = pt->get_triorg();
    tt[1] = tt[0];
    do {
      n++;
      tt[1].enext2self();
      tt[1].symself();
    } while (tt[1].tri != tt[0].tri);

    elist = new TriEdge[n];
    tt[1] = tt[0]; i = 0;
    do {
      elist[i] = tt[1]; i++;
      tt[1].enext2self();
      tt[1].symself();
    } while (tt[1].tri != tt[0].tri);
    assert(i == n);
  } // elist

  if (verbose > 2) {
    printf("  Removing vertex %d, n=%d\n", pt->idx, n);
  }

  // Select edge.
  REAL t_min = (double) (1 << 30); // A very big value
  REAL t_max = 0.0;
  //t_min *= delaunay;
  int t_min_idx = 0; // default
  int t_max_idx = 0; // default

  for (i = 0; i < n; i++) {
    tt[0] = elist[i];
    tt[1] = tt[0].sym();
    // Skip a hull edge.
    if (tt[0].tri->is_hulltri() || tt[1].tri->is_hulltri()) continue;
    assert(tt[0].org() == pt);
    Vertex *pb = tt[0].dest();
    Vertex *pc = tt[0].apex();
    Vertex *pd = tt[1].apex();
    REAL t = get_flip_time(pt, pb, pc, pd);
    if (verbose > 2) {
      printf("    [e%d] [%d,%d]-[%d,%d] t=%g\n", i,
             pt->idx, pb->idx, pc->idx, pd->idx, t);
    }
    // Ignore this edge if t is negative.
    if (t < 0) continue;
    if (t < t_min) {
      t_min = t;
      t_min_idx = i;
    }
    if (t > t_max) {
      t_max = t;
      t_max_idx = i;
    }
  } // i

  // Flip the selected edge.
  int flipflag = 0;
  if (delaunay > 0) {
    tt[0] = elist[t_max_idx];
    tt[1] = tt[0].sym();
  } else { // -u furtheset Delaunay
    tt[0] = elist[t_min_idx];
    tt[1] = tt[0].sym();
  }
  if (verbose > 2) {
    printf("    Flip [%d,%d]-[%d,%d] t=%g\n",
           tt[0].org()->idx, tt[0].dest()->idx, tt[0].apex()->idx,
           tt[1].apex()->idx, t_min);
  }
  if (flip(tt, flipflag, 0)) {
    // The edge is flipped.
    if ((flipflag & 2048) || (flipflag & 4096)) {
      // The vertex has been removed.
      return 1;
    } else {
      // Re-fill the reduced array "elist".
      tt[0] = pt->get_triorg();
      tt[1] = tt[0]; i = 0;
      do {
        elist[i] = tt[1]; i++;
        tt[1].enext2self();
        tt[1].symself();
      } while (tt[1].tri != tt[0].tri);
      assert(i == (n-1));
      n -= 1;
      return remove_point_regular(pt, elist, n);
    }
  }

  assert(n > 3);
  return 0; // Fail to remove this vertex.
}

///////////////////////////////////////////////////////////////////////////////

// [2016-11-14, 19:28, Berlin]
int Triangulation::monotone_flip_furthest_delaunay_test()
{
  // First remove all interior vertices
  int bak_lawsonflag = lawsonflag;
  lawsonflag = 0;
  TriEdge tt[4];
  int fail_count = 0;
  int i, count = 0;
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = vrts->get(i);
    if (!vrts->is_deleted(vrt)) {
      if (vrt->get_vrttype() != 0) { // NOT UNUSEDVERTEX
        // Check if it is an interior vertex.
        tt[0] = vrt->get_triorg();
        tt[1] = tt[0];
        do {
          if (tt[1].tri->is_hulltri()) break;
          tt[1].enext2self();
          tt[1].symself();
        } while (tt[1].tri != tt[0].tri);
        if (!tt[1].tri->is_hulltri()) {
          // Remove an interior vertex.
          TriEdge *elist = NULL;
          int n = 0;
          int success = remove_point_regular(vrt, elist, n);
          delete elist;
          if (!success) {
            printf("  !! Failed to remove vertex %d.\n", vrt->idx);
            fail_count++;
          }
        }
      }
      count++;
    }
  } // i

  lawsonflag = bak_lawsonflag;
  if (fail_count > 0) {
    printf(" !! %d interior vertices are not removed.\n", fail_count);
  }

  assert(flip_queue->objects == 0);

  // Perform Lawson's flip
  return lawson_flip(NULL);
}

///////////////////////////////////////////////////////////////////////////////
// [2016-11-15, 22:15, Berlin] // -GW
// Based on an input triangulation, randomly generates weights at each vertex.
//   each weight is calculated by: weight = r * L_ave * scale_weight,
//   where, L_ave is the average edge length at this vertex, scale is a value
//   between 0 and 1 (to scale L_ave), r is a random value between 0 and 1.
//   default scale_weight is 0.5, it can be changed by -GW#.

REAL Triangulation::get_distance(Vertex *e1, Vertex *e2)
{
  REAL dx = e2->crd[0] - e1->crd[0];
  REAL dy = e2->crd[1] - e1->crd[1];
  return sqrt(dx * dx + dy * dy);
}

int Triangulation::generate_vertex_weights_random()
{
  REAL *vert_weights = new REAL[vrts->objects];

  TriEdge tt[4];
  Vertex *endpt;
  int ecount = 0; // Count the number of edges.
  REAL L_tot, L_ave, r;

  int i, count = 0;
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = vrts->get(i);
    if (!vrts->is_deleted(vrt)) {
      if (vrt->get_vrttype() != 0) { // NOT UNUSEDVERTEX
        ecount = 0;
        L_tot = 0.;
        tt[0] = vrt->get_triorg();
        tt[1] = tt[0];
        do {
          endpt = tt[1].dest();
          if (endpt != infvrt) {
            L_tot = get_distance(vrt, endpt);
            ecount++;
          }
          tt[1].enext2self();
          tt[1].symself();
        } while (tt[1].tri != tt[0].tri);
        assert(ecount > 0);
        L_ave = L_tot / (REAL) ecount;
        // Random value.
        int R = rand() % 100;
        r = (double) R / 100.0;
        vert_weights[count] = r * L_ave;
      } else {
        vert_weights[count] = 0.;
      }
      count++;
    }
  } // i

  char filename[256];
  strcpy(filename, outfilename);
  strcat(filename, ".weight");

  FILE *outfile = fopen(filename, "w");
  int nv = vrts->objects;
  printf("Writing %d weights to file %s.\n", nv, filename);

  fprintf(outfile, "%d\n", nv);

  int idx = firstindex; // Reindex the vertices.
  for (i = 0; i < nv; i++) {
    fprintf(outfile, "%d  %g\n", idx, vert_weights[i]);
    idx++;
  }

  fclose(outfile);

  delete  [] vert_weights;
  return 0;
}