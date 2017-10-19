#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "detri2.h"

using  namespace detri2;

///////////////////////////////////////////////////////////////////////////////

void Triangulation::randomly_permute(Vertex **permutarray, int arraysize)
{
  Vertex *swap;
  int randidx, i;

  for (i = 0; i < arraysize; i++) {
    randidx = rand() % (i+1);
    swap = permutarray[i];
    permutarray[i] = permutarray[randidx];
    permutarray[randidx] = swap;
  }
}

///////////////////////////////////////////////////////////////////////////////
// Print the binary representation of the integer 'I'.
// 'm' is the number of bits to represent 'I'.
void print_b(int I, int m)
{
  int b, i;

  printf("  %d: ", I);
  for (i = m - 1; i >= 0; i--) {
    b = I & (1 << i);
    printf("%d", b ? 1 : 0);
  }
  printf("\n");
}

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// hilbert_init()    Initialize the Gray code permutation table.             //
//                                                                           //
// The table 'transgc' has 8 x 3 x 8 entries. It contains all possible Gray  //
// code sequences traveled by the 1st order Hilbert curve in 3 dimensions.   //
// The first column is the Gray code of the entry point of the curve, and    //
// the second column is the direction (0, 1, or 2, 0 means the x-axis) where //
// the exit point of curve lies.                                             //
//                                                                           //
// The table 'tsb1mod3' contains the numbers of trailing set '1' bits of the //
// indices from 0 to 7, modulo by '3'. The code for generating this table is //
// from: http://graphics.stanford.edu/~seander/bithacks.html.                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

int transgc[8][3][8], tsb1mod3[8];

void Triangulation::hilbert_init(int n)
{
  int gc[8], N, mask, travel_bit;
  int e, d, f, k, g;
  int v, c;
  int i;

  N = (n == 2) ? 4 : 8;
  mask = (n == 2) ? 3 : 7;

  // Generate the Gray code sequence.
  for (i = 0; i < N; i++) {
    gc[i] = i ^ (i >> 1);
  }

  for (e = 0; e < N; e++) {
    for (d = 0; d < n; d++) {
      // Calculate the end point (f).
      f = e ^ (1 << d);  // Toggle the d-th bit of 'e'.
      // travel_bit = 2**p, the bit we want to travel. 
      travel_bit = e ^ f;
      for (i = 0; i < N; i++) {
        // // Rotate gc[i] left by (p + 1) % n bits.
        k = gc[i] * (travel_bit * 2);
        g = ((k | (k / N)) & mask);
        // Calculate the permuted Gray code by xor with the start point (e).
        transgc[e][d][i] = (g ^ e);
      }
// \iffalse
      assert(transgc[e][d][0] == e);
      assert(transgc[e][d][N - 1] == f);
// \fi
    } // d
  } // e

  // Count the consecutive '1' bits (trailing) on the right.
  tsb1mod3[0] = 0;
  for (i = 1; i < N; i++) {
    v = ~i; // Count the 0s.
    v = (v ^ (v - 1)) >> 1; // Set v's trailing 0s to 1s and zero rest
    for (c = 0; v; c++) {
      v >>= 1;
    }
    tsb1mod3[i] = c % n;
  }
}

int Triangulation::
  hilbert_split(Vertex **vertexarray, int arraysize, int gc0, int gc1,
                double bxmin, double bxmax, double bymin, double bymax)
{
  Vertex *swapvert;
  int axis, d;
  REAL split;
  int i, j;

  // \iffalse
  //if (b->debug) {
    //if (b->debug_hsort) { // -Gh
      //printf("    gc0 :");
      //print_b(gc0, 2);
      //printf("    gc1 :");
      //print_b(gc1, 2);
    //}
  //}
// \fi

  // Find the current splitting axis. 'axis' is a value 0, or 1, or 2, which 
  //   correspoding to x-, or y- or z-axis.
  axis = (gc0 ^ gc1) >> 1; 

  // Calulate the split position along the axis.
  if (axis == 0) {
    split = 0.5 * (bxmin + bxmax);
  } else if (axis == 1) {
    split = 0.5 * (bymin + bymax);
  } else { // == 2
    assert(0); //split = 0.5 * (bzmin + bzmax);
  }

  // Find the direction (+1 or -1) of the axis. If 'd' is +1, the direction
  //   of the axis is to the positive of the axis, otherwise, it is -1.
  d = ((gc0 & (1<<axis)) == 0) ? 1 : -1;

  // Partition the vertices into left- and right-arrays such that left points
  //   have Hilbert indices lower than the right points.
  i = 0;
  j = arraysize - 1;

  // Partition the vertices into left- and right-arrays.
  if (d > 0) {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i]->crd[axis] >= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j]->crd[axis] < split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  } else {
    do {
      for (; i < arraysize; i++) {      
        if (vertexarray[i]->crd[axis] <= split) break;
      }
      for (; j >= 0; j--) {
        if (vertexarray[j]->crd[axis] > split) break;
      }
      // Is the partition finished?
      if (i == (j + 1)) break;
      // Swap i-th and j-th vertices.
      swapvert = vertexarray[i];
      vertexarray[i] = vertexarray[j];
      vertexarray[j] = swapvert;
      // Continue patitioning the array;
    } while (true);
  }

// \iffalse
  //if (b->debug) {
    //if (b->debug_hsort) { // -Gh
      //printf("    leftsize = %d, rightsize = %d\n", i, arraysize - i);
    //}
  //}
// \fi
  return i;
}

// 'hvertices' returns the list of Hilbert points on the path.
// It is only for debugging.
void Triangulation::
  hilbert_sort(Vertex **vertexarray, int arraysize, int e, int d,
               double bxmin, double bxmax, double bymin, double bymax,
               int depth, int order, arraypool *hvertices)
{
  REAL x1, x2, y1, y2;
  int p[5], w, e_w, d_w, k, ei, di;
  int n = 2, mask = 3;

  if (vertexarray != NULL) {
    p[0] = 0;
    p[4] = arraysize;

    // Sort the points according to the 1st order Hilbert curve in 2d.
    p[2] = hilbert_split(vertexarray, p[4], transgc[e][d][1], transgc[e][d][2],
                         bxmin, bxmax, bymin, bymax);
    p[1] = hilbert_split(vertexarray, p[2], transgc[e][d][0], transgc[e][d][1],
                         bxmin, bxmax, bymin, bymax);
    p[3] = hilbert_split(&(vertexarray[p[2]]), p[4] - p[2], transgc[e][d][2], transgc[e][d][3],
                         bxmin, bxmax, bymin, bymax);
  }

  if (order > 0) {
    // A maximum order is prescribed. 
    if ((depth + 1) == order) {
      // The maximum prescribed order is reached.
// \iffalse
      if (hvertices != NULL) {
        // Output the centers of sub-boxes.
        for (w = 0; w < 4; w++) {
          if (transgc[e][d][w] & 1) { // x-axis
            x1 = 0.5 * (bxmin + bxmax);
            x2 = bxmax;
          } else {
            x1 = bxmin;
            x2 = 0.5 * (bxmin + bxmax);
          }
          if (transgc[e][d][w] & 2) { // y-axis
            y1 = 0.5 * (bymin + bymax);
            y2 = bymax;
          } else {
            y1 = bymin;
            y2 = 0.5 * (bymin + bymax);
          }
          Vertex *parypt = (Vertex *) hvertices->alloc();
          parypt->crd[0] = 0.5 * (x1 + x2);
          parypt->crd[1] = 0.5 * (y1 + y2);
        }
      } // if (hvertices != NULL)
// \fi
      return;
    }
  }

  // Recursively sort the points in sub-boxes.
  for (w = 0; w < 4; w++) {
    // w is the local Hilbert index (NOT Gray code).
    // Sort into the sub-box either there are more than 2 points in it, or
    //   the prescribed order of the curve is not reached yet.
    if ((p[w+1] - p[w] > 2) || (order > 0)) {
      // Calculcate the start point (ei) of the curve in this sub-box.
      //   update e = e ^ (e(w) left_rotate (d+1)).
      if (w == 0) {
        e_w = 0;
      } else {
        //   calculate e(w) = gc(2 * floor((w - 1) / 2)).
        k = 2 * ((w - 1) / 2); 
        e_w = k ^ (k >> 1); // = gc(k).
      }
      k = e_w;
      e_w = ((k << (d+1)) & mask) | ((k >> (n-d-1)) & mask);
      ei = e ^ e_w;
      // Calulcate the direction (di) of the curve in this sub-box.
      //   update d = (d + d(w) + 1) % n
      if (w == 0) {
        d_w = 0;
      } else {
// \iffalse
        /*if ((w % 2) == 0) {
          d_w = tsb1(w - 1) % n;
        } else {
          d_w = tsb1(w) % n;
	} */
// \fi
        d_w = ((w % 2) == 0) ? tsb1mod3[w - 1] : tsb1mod3[w];
      }
      di = (d + d_w + 1) % n;
      // Calculate the bounding box of the sub-box.
      if (transgc[e][d][w] & 1) { // x-axis
        x1 = 0.5 * (bxmin + bxmax);
        x2 = bxmax;
      } else {
        x1 = bxmin;
        x2 = 0.5 * (bxmin + bxmax);
      }
      if (transgc[e][d][w] & 2) { // y-axis
        y1 = 0.5 * (bymin + bymax);
        y2 = bymax;
      } else {
        y1 = bymin;
        y2 = 0.5 * (bymin + bymax);
      }
      //if (transgc[e][d][w] & 4) { // z-axis
      //  z1 = 0.5 * (bzmin + bzmax);
      //  z2 = bzmax;
      //} else {
      //  z1 = bzmin;
      //  z2 = 0.5 * (bzmin + bzmax);
      //}
// \iffalse
      //if (b->debug) {
        //if (b->debug_hsort) { // -Gh
          //printf("  array length = %d.\n", p[w+1] - p[w]);
          //printf("  e[%d]: ", w);
          //print_b(ei, 2);
          //printf("  d[%d]: ", w);
          //print_b(di, 2);
        //}
      //}
// \fi
      hilbert_sort(&(vertexarray[p[w]]), p[w+1] - p[w], ei, di,
                   x1, x2, y1, y2, depth+1, order, hvertices);
    } // if (p[w+1] - p[w] > 1)
  } // w
}

/*
  int threshold = 64;
  double ratio = 0.5;
  int depth = 0;
  int order = 0;
  randomly_permute(ptlist, ptnum);
  brio_multiscale_sort(ptlist, ptnum, threshold, ratio, &depth, order);
*/

void Triangulation::
  brio_multiscale_sort(Vertex **vertexarray, int arraysize,
                       int threshold, double ratio, int *depth, int order)
{
  int middle;

  middle = 0;
  if (arraysize >= threshold) {
    (*depth)++;
    middle = arraysize * ratio;
    brio_multiscale_sort(vertexarray, middle, threshold, ratio, depth, order);
  }
  // Sort the right-array (rnd-th round) using the Hilbert curve.
  hilbert_sort(&(vertexarray[middle]), arraysize - middle, 0, 0, // e, d
               xmin, xmax, ymin, ymax, 0, order, NULL); // depth.
}

// Sort the set of points along a given vector (vx, vy).
void Triangulation::sweepline_sort(Vertex **ptlist, int ptnum, double vx, double vy)
{
  int n = 0, i;
  while (n < ptnum) {
    double min = ptlist[n]->crd[0] * vx + ptlist[n]->crd[1] * vy;
    int min_idx = n;
    for (i = n + 1; i < ptnum; i++) {
      double value = ptlist[i]->crd[0] * vx + ptlist[i]->crd[1] * vy;
      if (value < min) {
        min = value;
        min_idx = i;
      }
    }
    // Swap i and n;
    Vertex *swappt = ptlist[min_idx];
    ptlist[min_idx] = ptlist[n];
    ptlist[n] = swappt;
    n++;
  }
  /*
  if (ptnum == 2) {
    double min = ptlist[0]->crd[0] * vx + ptlist[0]->crd[1] * vy;
    double val = ptlist[1]->crd[0] * vx + ptlist[1]->crd[1] * vy;
    if (val < min) {
      // Swap i and n;
      Vertex *swappt = ptlist[0];
      ptlist[0] = ptlist[1];
      ptlist[1] = swappt;
    }
    return;
  }

  int pivot = (rand() % ptnum);
  double pivot_val = ptlist[pivot]->crd[0] * vx + ptlist[pivot]->crd[1] * vy;
  //double val;
  int left = -1;
  int right = ptnum;
  while (left < right) {
    do {
      left++;
      //val = ptlist[left]->crd[0] * vx + ptlist[left]->crd[1] * vy;
    } while ((left <= right) &&
             ((ptlist[left]->crd[0] * vx + ptlist[left]->crd[1] * vy) <= pivot_val));
    do {
      right--;
      //val = ptlist[right]->crd[0] * vx + ptlist[right]->crd[1] * vy;
    } while ((left <= right) &&
             ((ptlist[right]->crd[0] * vx + ptlist[right]->crd[1] * vy) >= pivot_val));
    if (left < right) {
      // Swap left and right;
      Vertex *swappt = ptlist[left];
      ptlist[left] = ptlist[right];
      ptlist[right] = swappt;
    }
  } // while (left < right)

  if (left > 1) {
    sweepline_qsort(ptlist, left, vx, vy);
  }
  if (right < (ptnum - 2)) {
    sweepline_qsort(&(ptlist[right+1]), ptnum-right-1, vx, vy);
  }
  */
}