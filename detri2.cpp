
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include "detri2.h"
#include "pred2d.h"

using  namespace detri2;

#define UNDEFIND -99999

//=============================================================================
// TriEdge implementation

// It uses the 1st bit.
bool Vertex::is_infected() {return flags & 1;}
void Vertex::set_infect() {flags |= 1;}
void Vertex::clear_infect() {flags &= ~1;}

// Vertex types are stored from 9th-12th bits (4 bits)
// See table ``const char *vertextypename[]'' in dbg.cpp
int Vertex::get_vrttype() {return ((flags >> 8) & 15);}
void Vertex::set_vrttype(int val) {
  assert(val <= 15);
  flags &=  ~(15 << 8); // Clear the old bits.
  flags |= (val << 8); // Set the val.
}

TriEdge Vertex::get_triorg()
{
  int ver = 0;
  if (adj->vrt[0] == this) {
    ver = 2;
  } else if (adj->vrt[1] == this) {
    ver = 0;
  } else {
    assert(adj->vrt[2] == this);
    ver = 1;
  }
  return TriEdge(adj, ver);
}

// Traingle::flag (integer) is actually a collection of bits.
//   (In C/C++, an integer has at least 32 bits).
//
// The currently used bits are summaried below:
//   0th -- 11th, encode the three neighbor edges's ver,
//      used by Triangle::getnver(), setnver().
//   13th, 14th, 15th, store the segment-flag for its three edges,
//      used by TriEdge::is_segment(), ...
//   16th, stroes the hull (faked) flag for this triangle,
//      used by Triangle::is_hulltri(), ...
//   17th, stroes the infection flag for this triangle,
//      used by Triangle::is_infected(), ...
//   20th, 21th, 22th, store the infection flags for its three edges,
//      used by TriEdge::is_edge_infected(), ...

// Use 0th-11th bits
int Triangle::getnver(int ver) {return (flags >> (ver * 4)) & 15;}
void Triangle::setnver(int ver, int val) {
  flags &=  ~(15 << (ver * 4)); // Clear the old bits.
  flags |= (val << (ver * 4)); // Set the val.
}

// Use the 16th bit.
bool Triangle::is_hulltri() {return (flags & (1 << 15));}
void Triangle::set_hullflag() {
  assert(!is_hulltri());
  flags |= (1 << 15); // Mark the 16th bit.
}
void Triangle::clear_hullflag() {
  assert(is_hulltri());
  flags &= (~(1 << 15)); // Unmark the 16th bit.
}

// Use the 17th bit.
void Triangle::set_infect() {flags |= (1 << 16);}
bool Triangle::is_infected() {return (flags & (1 << 16));}
void Triangle::clear_infect() {flags &= (~(1 << 16));}

// Exterior triangles are marked as tag = -9999
void Triangle::set_exterior() {tag = -9999;}
bool Triangle::is_exterior() {return tag == -9999;}
void Triangle::clear_exterior() {tag = 0;}

// Initialize tables
static int _vo[3] = {1, 2, 0};
static int _vd[3] = {2, 0, 1};
static int _va[3] = {0, 1, 2};

Vertex* TriEdge::org() {return tri->vrt[_vo[ver]];}
Vertex* TriEdge::dest() {return tri->vrt[_vd[ver]];}
Vertex* TriEdge::apex() {return tri->vrt[_va[ver]];}

void TriEdge::set_vertices(Vertex *pa, Vertex *pb, Vertex *pc)
{
  tri->vrt[0] = pa;
  tri->vrt[1] = pb;
  tri->vrt[2] = pc;
  ver = 2; // Edge [a,b]
}

TriEdge TriEdge::enext() {return TriEdge(tri, _vo[ver]);}
void TriEdge::enextself() {ver = _vo[ver];}
TriEdge TriEdge::enext2() {return TriEdge(tri, _vd[ver]);}
void TriEdge::enext2self() {ver = _vd[ver];}

TriEdge TriEdge::sym() {return TriEdge(tri->nei[ver], tri->getnver(ver));}
void TriEdge::symself() {
  int nver = tri->getnver(ver);
  tri = tri->nei[ver];
  ver = nver;
}

// Connecting t1->t2 (one direction)
void TriEdge::bond1(const TriEdge &ntri)
{
  tri->nei[ver] = ntri.tri;
  tri->setnver(ver, ntri.ver);
}

// Connecting t1<->t2 (Both direction)
void TriEdge::bond2(const TriEdge &ntri)
{
  tri->nei[ver] = ntri.tri;
  tri->setnver(ver, ntri.ver);
  ntri.tri->nei[ntri.ver] = tri;
  ntri.tri->setnver(ntri.ver, ver);
}

void TriEdge::dissolve() {tri->nei[ver] = NULL;}

// check the 13th+ver bit.
bool TriEdge::is_segment() { return (tri->flags & (1 << (12+ver)));}
void TriEdge::set_segment() {tri->flags |= (1 << (12+ver));}
void TriEdge::clear_segment() {tri->flags &= (~(1 << (12+ver)));}

// check the 20th+ver bit.
bool TriEdge::is_edge_infected() { return (tri->flags & (1 << (20+ver)));}
void TriEdge::set_edge_infect() {tri->flags |= (1 << (20+ver));}
void TriEdge::clear_edge_infect() {tri->flags &= (~(1 << (20+ver)));}

bool EncEle::is_changed()
{
  if (Ele.tri != NULL) {
    if ((end[0] == Ele.tri->vrt[0]) &&
        (end[1] == Ele.tri->vrt[1]) &&
        (end[2] == Ele.tri->vrt[2])) return false;
  }
  return true;
}

//=============================================================================
// arraypool implementation

void arraypool::restart()
{
  objects = 0;
  deaditemstack = NULL;
}

void arraypool::poolinit(int sizeofobject, int log2objperblk)
{
  // Each object must be at least one byte long.
  objectbytes = sizeofobject > 1 ? sizeofobject : 1;

  log2objectsperblock = log2objperblk;
  // Compute the number of objects in each block.
  objectsperblock = ((int) 1) << log2objectsperblock;
  objectsperblockmark = objectsperblock - 1;

  // No memory has been allocated.
  totalmemory = 0l;
  // The top array has not been allocated yet.
  toparray = (char **) NULL;
  toparraylen = 0;

  // Ready all indices to be allocated.
  restart();
}

arraypool::arraypool() {}

arraypool::arraypool(int sizeofobject, int log2objperblk)
{
  poolinit(sizeofobject, log2objperblk);
}

arraypool::~arraypool()
{
  int i;

  // Has anything been allocated at all?
  if (toparray != (char **) NULL) {
    // Walk through the top array.
    for (i = 0; i < toparraylen; i++) {
      // Check every pointer; NULLs may be scattered randomly.
      if (toparray[i] != (char *) NULL) {
        // Free an allocated block.
        free((void *) toparray[i]);
      }
    }
    // Free the top array.
    free((void *) toparray);
  }

  // The top array is no longer allocated.
  toparray = (char **) NULL;
  toparraylen = 0;
  objects = 0;
  totalmemory = 0;
  deaditemstack = NULL;
}

// Used by newindex() (alloc)
char* arraypool::getblock(int objectindex)
{
  char **newarray;
  char *block;
  int newsize;
  int topindex;
  int i;

  // Compute the index in the top array (upper bits).
  topindex = objectindex >> log2objectsperblock;
  // Does the top array need to be allocated or resized?
  if (toparray == (char **) NULL) {
    // Allocate the top array big enough to hold 'topindex', and NULL out
    //   its contents.
    newsize = topindex + 128;
    toparray = (char **) malloc((size_t) (newsize * sizeof(char *)));
    toparraylen = newsize;
    for (i = 0; i < newsize; i++) {
      toparray[i] = (char *) NULL;
    }
    // Account for the memory.
    totalmemory = newsize * (uintptr_t) sizeof(char *);
  } else if (topindex >= toparraylen) {
    // Resize the top array, making sure it holds 'topindex'.
    newsize = 3 * toparraylen;
    if (topindex >= newsize) {
      newsize = topindex + 128;
    }
    // Allocate the new array, copy the contents, NULL out the rest, and
    //   free the old array.
    newarray = (char **) malloc((size_t) (newsize * sizeof(char *)));
    for (i = 0; i < toparraylen; i++) {
      newarray[i] = toparray[i];
    }
    for (i = toparraylen; i < newsize; i++) {
      newarray[i] = (char *) NULL;
    }
    free(toparray);
    // Account for the memory.
    totalmemory += (newsize - toparraylen) * sizeof(char *);
    toparray = newarray;
    toparraylen = newsize;
  }

  // Find the block, or learn that it hasn't been allocated yet.
  block = toparray[topindex];
  if (block == (char *) NULL) {
    // Allocate a block at this index.
    block = (char *) malloc((size_t) (objectsperblock * objectbytes));
    toparray[topindex] = block;
    // Account for the memory.
    totalmemory += objectsperblock * objectbytes;
  }

  // Return a pointer to the block.
  return block;
}

// Return the pointer to the object with a given index, or NULL
//   if the object's block doesn't exist yet.
void* arraypool::lookup(int objectindex)
{
  char *block;
  int topindex;

  // Has the top array been allocated yet?
  if (toparray == (char **) NULL) {
    return (void *) NULL;
  }

  // Compute the index in the top array (upper bits).
  topindex = objectindex >> log2objectsperblock;
  // Does the top index fit in the top array?
  if (topindex >= toparraylen) {
    return (void *) NULL;
  }

  // Find the block, or learn that it hasn't been allocated yet.
  block = toparray[topindex];
  if (block == (char *) NULL) {
    return (void *) NULL;
  }

  // Compute a pointer to the object with the given index.  Note that
  //   'objectsperblock' is a power of two, so the & operation is a bit mask
  //   that preserves the lower bits.
  return (void *)(block + (objectindex & (objectsperblock - 1)) * objectbytes);
}

void* arraypool::alloc()
{
  void *newptr;
  if (deaditemstack != (void *) NULL) {
    newptr = deaditemstack;                     // Take first item in list.
    deaditemstack = * (void **) deaditemstack;
  } else {
    // Allocate an object at index 'firstvirgin'.
    newptr = (void *) (getblock(objects) + 
                       (objects & (objectsperblock - 1)) * objectbytes);
  }
  objects++;
  return newptr;
}

void arraypool::dealloc(void *dyingitem)
{
  // Push freshly killed item onto stack.
  *((void **) dyingitem) = deaditemstack;
  deaditemstack = dyingitem;
  objects--;
}

void* arraypool::operator[](int index) 
{
  return (void *) (toparray[index >> log2objectsperblock] + 
            ((index) & objectsperblockmark) * objectbytes);
}

//=============================================================================

Triangle* Triangles::alloc_tri()
{
  Triangle *newtri = (Triangle *) alloc();
  newtri->init(); // Initialize.
  return newtri;
}

void Triangles::dealloc_tri(Triangle* tri)
{
  tri->vrt[0] = NULL;
  dealloc((void *) tri);
}

bool Triangles::is_deleted(Triangle* tri)
{
  return (tri->vrt[0] == NULL);
}

Triangle* Triangles::get(int index)
{
  return (Triangle *) (toparray[index >> log2objectsperblock] +
            ((index) & objectsperblockmark) * objectbytes);
}

Vertex* Vertices::alloc_vrt()
{
  Vertex *newvrt = (Vertex *) alloc();
  newvrt->init();
  return newvrt;
}

void Vertices::dealloc_vrt(Vertex* vrt)
{
  vrt->idx = UNDEFIND;
  dealloc((void *) vrt);
}

bool Vertices::is_deleted(Vertex* vrt)
{
  return (vrt->idx == UNDEFIND);
}

Vertex* Vertices::get(int index)
{
  return (Vertex *) (toparray[index >> log2objectsperblock] +
            ((index) & objectsperblockmark) * objectbytes);
}

//=============================================================================
// Triangulation implementation

Triangulation::Triangulation()
{
  tris = new Triangles(16);
  segs = new Triangles(8);
  vrts = new Vertices(8);
  regions = new Vertices(4);
  infvrt = new Vertex();
  infvrt->idx = -1;
  vertex_data = NULL;
  vertex_data_size = 0;
  triangle_data = NULL;
  triangle_data_size = 0;
  flip_queue = new arraypool(sizeof(TriEdge), 8);
  idx2seglist = NULL;
  idx_to_vert_list = NULL;
  idx_to_vert_list_len = 0;
  segperverlist = NULL;
  dual_verts = NULL;
  tet_verts = NULL;
  moving_vrts = new arraypool(sizeof(MovVertex), 4);//ADD

  lawsonflag = lawson_hullflag = 0;
  chktris = chktris2 = NULL;
  encsegs = encsegs2 = NULL;
  encsegflag = enqtriflag = 0;
  input_vertices = input_triangles = 0;
  hullsize = 0;
  unused_vertices = input_triangles = 0;
  flip13count = flip31count = flip22count = totalflipcount = 0;
  st_segref_count = st_triref_count = 0;
  infilename[0] = '\0';
  outfilename[0] = '\0';
  firstindex = 0;

  quiet = 0; // -Q
  verbose = 1; // -V
  norandomize = 0; // -R
  convex = 0;
  delaunay = 1; // -u
  weighted = 0; // -w
  weight_is_height = 0; // -wh
  transfer_weight_to_height = 0; // -wt
  gabriel = 0; // -G
  quality = 0; // -q
  minangle = 20.0;
  maxratio2 = 0;
  maxarea = 0.0;
  minedgelen = 0.0;
  no_reindex = 0; // -I
  no_node_tag = 0; // -B
  output_edges = 0;

  // Debug options
  checkmesh = 0;
  checkdelaunay = 0;
  test_hilbert_sort = 0; // -Gh
  sweep_sort = 0; // -Gs
  sweep_sort_vx = 0.0; // Default x.
  sweep_sort_vy = 1.0;
  do_lawson_flip = 0; // -Gl
  do_lawson_flip_origin = 0; // -GF
  do_lawson_flip_test = 0; // -Gt
  do_generate_vertex_weights = 0; // -GW
  dump_intermediate_meshes = 0;
  dump_to_mesh = 0;
  dump_to_ucd = 0; // -Gu
  dump_dual_to_ucd = 0; // -GU
  dump_lifting_map = 0; // -L
  lift_map_offset_z = 0.0; // -L# Default (for visualization)
  dump_dual_cell_areas = 0; // -GA
  dump_lifting_poly = 0; // -GP
  dump_lifting_tets = 0; // -GT
  dump_cycles = 0; // =GS
  exterior_triangles = 0;
  debug_conforming_delaunay = 0;

  for (int i = 0; i < 4; i++)//ADD for hi_di function parameter
	  para[i] = 0;

  PI = 3.14159265358979323846264338327950288419716939937510582;

  // Initialize exact arthimetics.
  exactinit();
}

Triangulation::~Triangulation()
{
  delete tris;
  delete segs;
  delete vrts;
  delete regions;
  delete infvrt;
  delete flip_queue;
  delete moving_vrts;//ADD
  if (vertex_data) delete vertex_data;
  if (triangle_data) delete triangle_data;
  if (dual_verts) delete [] dual_verts;
  if (tet_verts) delete tet_verts;
  if (idx_to_vert_list) delete [] idx_to_vert_list;
}

bool Triangulation::is_hulledge(TriEdge &E)
{
  return (E.apex() == infvrt);
}

void Triangulation::flip13(Vertex *pt, TriEdge *tt)
{
  // On input
  //   tt[0] = [a,b,c]
  // On ouput
  //   tt[0] = [a,b,p]
  //   tt[1] = [b,c,p]
  //   tt[2[ = [c,a,p]
  int hullflag = tt[0].tri->is_hulltri();
  int ttag = tt[0].tri->tag;
  int i;

  Vertex *pa = tt[0].org();
  Vertex *pb = tt[0].dest();
  Vertex *pc = tt[0].apex();

  if (verbose > 2) {
    printf("    flip13 [%d]-[%d,%d,%d].\n", pt->idx, pa->idx, pb->idx, pc->idx);
  }
  if (dump_lifting_tets) { // -GT
    Vertex **pTet = (Vertex **) tet_verts->alloc();
    pTet[0] = pa; pTet[1] = pb; pTet[2] = pc; pTet[3] = pt;
  }

  TriEdge nn[3];
  for (i = 0; i < 3; i++) {
    nn[i] = tt[0].sym();
    tt[0].enextself();
  }

  tris->dealloc_tri(tt[0].tri);

  tt[0].tri = tris->alloc_tri();
  tt[1].tri = tris->alloc_tri();
  tt[2].tri = tris->alloc_tri();

  tt[0].set_vertices(pa, pb, pt);
  tt[1].set_vertices(pb, pc, pt);
  tt[2].set_vertices(pc, pa, pt);

  if (hullflag) {
    // pa, pb, or pc may be infvert
    if (pa == infvrt) {
      tt[0].tri->set_hullflag();
      tt[2].tri->set_hullflag();
      tt[1].tri->tag = ttag;
    } else if (pb == infvrt) {
      tt[0].tri->set_hullflag();
      tt[1].tri->set_hullflag();
      tt[2].tri->tag = ttag;
    } else {
      assert(pc == infvrt);
      tt[1].tri->set_hullflag();
      tt[2].tri->set_hullflag();
      tt[0].tri->tag = ttag;
    }
    hullsize += 1; // Removed 1 added 2. (Convex Assumption)
  } else {
    for (i = 0; i < 3; i++) {
      tt[i].tri->tag = ttag;
    }
  }

  // Connect the three exterior neighbors.
  for (i = 0; i < 3; i++) {
    tt[i].bond2(nn[i]);
    if (nn[i].is_segment()) {
      tt[i].set_segment();
    }
  }

  // Connect the three interior neighbors. Re-use nn
  for (i = 0; i < 3; i++) {
    nn[0] = tt[i].enext();
    nn[1] = tt[(i+1)%3].enext2();
    nn[0].bond2(nn[1]);
  }

  // Update the vrt-to-tri map.
  pt->adj = tt[0].tri;
  pa->adj = tt[0].tri;
  pb->adj = tt[0].tri;
  pc->adj = tt[1].tri;

  if (lawsonflag) {
    // Push three exterior edges into flip queue.
    TriEdge *qedge;
    for (i = 0; i < 3; i++) {
      qedge = (TriEdge *) flip_queue->alloc();
      *qedge = tt[i];
    }
    if (weighted || (delaunay < 0)) { // weighted (or furthest) Delaunay
      // Push three interior edge as well.
      for (i = 0; i < 3; i++) {
        qedge = (TriEdge *) flip_queue->alloc();
        *qedge = tt[i].enext();
      }
    }
  }

  flip13count++;
  totalflipcount++;
  if (dump_intermediate_meshes) {
    save_to_ucd(totalflipcount);
  }
}

void Triangulation::flip31(TriEdge *tt)
{
  // On input
  //   tt[0] = [a,b,p]
  //   tt[1] = [b,c,p]
  //   tt[2[ = [c,a,p]
  // On ouput
  //   tt[0] = [a,b,c]
  int hullflag = (tt[0].tri->is_hulltri() ||
                  tt[1].tri->is_hulltri() ||
                  tt[2].tri->is_hulltri());
  int ttag = 0;
  if (!tt[0].tri->is_hulltri()) {
    ttag = tt[0].tri->tag;
  } else if (!tt[1].tri->is_hulltri()) {
    ttag = tt[1].tri->tag;
  } else if (!tt[2].tri->is_hulltri()) {
    ttag = tt[2].tri->tag;
  }
  int i;

  Vertex *pa = tt[0].org();
  Vertex *pb = tt[0].dest();
  Vertex *pc = tt[1].dest();
  Vertex *pt = tt[1].apex();

  if (verbose > 2) {
    printf("    flip31 [%d,%d,%d]-[%d].\n", pa->idx, pb->idx, pc->idx, pt->idx);
  }
  if (dump_lifting_tets) { // -GT
    Vertex **pTet = (Vertex **) tet_verts->alloc();
    pTet[0] = pa; pTet[1] = pb; pTet[2] = pc; pTet[3] = pt;
  }

  TriEdge nn[3];
  for (i = 0; i < 3; i++) {
    nn[i] = tt[i].sym();
  }

  tris->dealloc_tri(tt[0].tri);
  tris->dealloc_tri(tt[1].tri);
  tris->dealloc_tri(tt[2].tri);

  tt[0].tri = tris->alloc_tri();
  tt[0].set_vertices(pa, pb, pc);

  if (hullflag) {
    // One of pa, pb, and pc must be infvrt.
    assert((pa == infvrt) || (pb == infvrt) || (pc == infvrt));
    tt[0].tri->set_hullflag();
    hullsize -= 1; // Removed 2 hull edges, added one.
  }

  tt[0].tri->tag = ttag;

  // Connect the three exterior neighbors.
  for (i = 0; i < 3; i++) {
    tt[0].bond2(nn[i]);
    if (nn[i].is_segment()) {
      tt[0].set_segment();
    }
    tt[0].enextself();
  }

  // Update the vrt-to-tri map.
  pa->adj = tt[0].tri;
  pb->adj = tt[0].tri;
  pc->adj = tt[0].tri;

  if (lawsonflag) {
    // Push three edges into flip queue.
    TriEdge *qedge;
    for (i = 0; i < 3; i++) {
      qedge = (TriEdge *) flip_queue->alloc();
      *qedge = tt[0];
      tt[0].enextself();
    }
  }

  flip31count++;
  totalflipcount++;
  if (dump_intermediate_meshes) {
    save_to_ucd(totalflipcount);
  }
}

void Triangulation::flip22(TriEdge *tt)
{
  // On input,
  //   tt[0] is [a,b,c], where [a,b] to be flipped.
  //   tt[1] is [b,a,d]
  // On output:
  //   tt[0] is [c,d,b],
  //   tt[1] is [d,c,a].
  int hullflag = (tt[0].tri->is_hulltri() ||
                  tt[1].tri->is_hulltri());
  int ttag = 0;
  if (!tt[0].tri->is_hulltri()) {
    ttag = tt[0].tri->tag;
  } else if (!tt[1].tri->is_hulltri()) {
    ttag = tt[1].tri->tag;
  }

  Vertex *pa = tt[0].org();
  Vertex *pb = tt[0].dest();
  Vertex *pc = tt[0].apex();
  Vertex *pd = tt[1].apex();

  if (verbose > 2) {
    printf("    flip22 [%d,%d]-[%d,%d].\n", pa->idx, pb->idx, pc->idx, pd->idx);
  }
  if (dump_lifting_tets) { // -GT
    Vertex **pTet = (Vertex **) tet_verts->alloc();
    pTet[0] = pa; pTet[1] = pb; pTet[2] = pc; pTet[3] = pd;
  }

  TriEdge nn[4];
  nn[0] = (tt[0].enext()).sym();  // [b,c]
  nn[1] = (tt[0].enext2()).sym(); // [c,a]
  nn[2] = (tt[1].enext()).sym();  // [a,d]
  nn[3] = (tt[1].enext2()).sym(); // [d,b]

  tris->dealloc_tri(tt[0].tri);
  tris->dealloc_tri(tt[1].tri);

  tt[0].tri = tris->alloc_tri();
  tt[1].tri = tris->alloc_tri();

  tt[0].set_vertices(pc, pd, pb);
  tt[1].set_vertices(pd, pc, pa);
  if (hullflag) {
    if (pa == infvrt) {
      tt[1].tri->set_hullflag();
      hullsize -= 1; // Remove 2 hull edges, add 1. (Convex Assumption)
    } else if (pb == infvrt) {
      tt[0].tri->set_hullflag();
      hullsize -= 1;
    } else if (pc == infvrt) {
      tt[0].tri->set_hullflag();
      tt[1].tri->set_hullflag();
      hullsize += 1; // Remove 1 hull edges, add 2. (Convex Assumption)
    } else if (pd == infvrt) {
      tt[0].tri->set_hullflag();
      tt[1].tri->set_hullflag();
      hullsize += 1;
    } else {
      assert(0); // Not possible.
    }
  }

  tt[0].tri->tag = ttag;
  tt[1].tri->tag = ttag;

  tt[0].bond2(tt[1]);
  nn[0].bond2(tt[0].enext2()); // [b,c]
  nn[1].bond2(tt[1].enext());  // [c,a]
  nn[2].bond2(tt[1].enext2()); // [a,d]
  nn[3].bond2(tt[0].enext());  // [d,b]

  // Set the connections to segments (if any).
  for (int i = 0; i < 4; i++) {
    if (nn[i].is_segment()) {
      (nn[i].sym()).set_segment();
    }
  }

  // Update the vrt-to-tri map.
  pa->adj = tt[1].tri;
  pb->adj = tt[0].tri;
  pc->adj = tt[0].tri;
  pd->adj = tt[0].tri;

  if (lawsonflag) {
    // Push four edges into flip queue.
    TriEdge *qedge;
    for (int i = 0; i < 4; i++) {
      qedge = (TriEdge *) flip_queue->alloc();
      *qedge = nn[i];
    }
  }

  totalflipcount++;
  if (dump_intermediate_meshes) {
    save_to_ucd(totalflipcount);
  }
}

// This function does both the check and flip an edge.
// On input,
//   tt[0]: [a,b,c], where [a,b] is the edge to be checked/flipped.
//   tt[1]: unknown
//   tt[2]: unknown
//   tt[3]: unknown
//   flipflag, return the type of flip.
//   check_only, only do check, but not do flip.
// Return 1 if this edge is flippable, and is flipped (if check_only = 0).
int Triangulation::flip(TriEdge *tt, int &flipflag, int check_only)
{
  flipflag = 0; // Clear iit.

  if (tt[0].is_segment()) {
    flipflag |= 4; // Constraint flag.
    if (encsegflag) {
      check_segment(NULL, tt[0]);
    }
    return 0; // Do not flip a segment.
  }

  tt[1] = tt[0].sym();
  int doflip = 0;
  double ori;

  if (!tt[0].tri->is_hulltri()) {
    if (!tt[1].tri->is_hulltri()) {
      // An interior flip.
      if (lawsonflag) {
        // Only do flip if this edge is not locally (weighted) Delaunay.
        if (!lawson_hullflag) { // Not only exterior
          if (!weighted) {
            ori = InCircle( tt[0].org()->crd, tt[0].dest()->crd,
                           tt[0].apex()->crd, tt[1].apex()->crd) * delaunay;
          } else { // weighted
            if (weight_is_height) {
              ori = Orient3d( tt[0].org()->crd, tt[0].dest()->crd,
                             tt[0].apex()->crd, tt[1].apex()->crd) * delaunay;  
            } else {
              double *pa = tt[0].org ()->crd;
              double *pb = tt[0].dest()->crd;
              double *pc = tt[0].apex()->crd;
              double *pd = tt[1].apex()->crd; 
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
          //doflip = (ori > 0);
          if (ori > 0) {
            doflip = 1;
          } else {
            flipflag |= 2; // Locally Delaunay.
          }
        } // if (!lawson_hullflag)
      } else {
        doflip = 1;
      }
    }
  } else {
    if (tt[1].tri->is_hulltri()) {
      // An exterior flip.
      flipflag |= (1 << 24); // hullflag
      if (lawsonflag) {
        // Only do flip is this hull edge is not locally convex.
        if (tt[0].org() == infvrt) {
          TriEdge swap = tt[0];
          tt[0] = tt[1];
          tt[1] = swap;
        }
        assert(tt[0].dest() == infvrt);
        assert(tt[1].org() == infvrt);
        // Only do flip if both hull edges are not segments.
        if (!(tt[0].enext2()).is_segment()) {
          if (Orient2d((tt[0].apex())->crd, (tt[0].org())->crd,
                       (tt[1].apex())->crd) > 0) {
            doflip = 1;
          } else {
            flipflag |= 2; // Locally Delaunay.
          }
        }
      } // if (lawsonflag)
    }
  }

  if (!doflip) return 0; // Return if it is not flippable.

  if (lawsonflag) {
    // It is not locally Delaunay.
    if ((!weighted && (delaunay > 0)) || // not weighted Delaunay or
        (flipflag & (1 << 24))) { // it's an hull edge.
      // In this case, this edge must be 2-to-2 flippable.
      flipflag |= 512; // A 2-to-2 flip is found.
      if (!check_only) {
        flip22(tt);
      }
      return 1;
    }
  }

// \iffalse
  // Only for debugging.
  assert(!tt[0].is_segment());
  Vertex *pa = tt[0].org();
  Vertex *pb = tt[0].dest();
  Vertex *pc = tt[0].apex();
  Vertex *pd = tt[1].apex();
  assert((flipflag & (1 << 24)) == 0); // not a hull flip
// \fi

  int ori1 = Orient2d(tt[0].org()->crd, tt[0].apex()->crd, tt[1].apex()->crd);

  if (ori1 >= 0) {
    // Swap tt[0] and tt[1];
    flipflag |= (1 << 25); // Spivot = 1
    tt[0] = tt[1];
    tt[1] = tt[0].sym();
// \iffalse
    pa = tt[0].org();
    pb = tt[0].dest();
    pc = tt[0].apex();
    pd = tt[1].apex();
// \fi
    ori1 = -1;
  }

  // assert(ori1 < 0);
  int ori2 = Orient2d(tt[0].dest()->crd, tt[0].apex()->crd, tt[1].apex()->crd);

  if (ori2 > 0) {
    flipflag |= 512; // A 2-to-2 flip is found (set the 10th bit).
    if (!check_only) {
      flip22(tt);
    }
    return 1;
  } else {
    if (pb->get_vrttype() != 2) { // RIDGEVERTEX
      if ((ori2 != 0) && (pb->get_vrttype() == 6)) { // FREESEGVERTEX
        // [a,b] must not be segment, then [b,c] and [b,d] must be.
        assert(tt[0].enext().is_segment()); // [b,c]
        assert(tt[1].enext2().is_segment()); // [b,d]
        ori2 = 0;
      }
      if (ori2 < 0) {
        // Check if a 3-to-1 at [b] flip is possible.
        tt[2] = tt[1].enext2().sym(); // [b,d,#]
        if (tt[2].apex() == tt[0].apex()) { // pc
          // Found a 3-to-1 flip.
          flipflag |= 2048; // set the 12th bit.
          tt[0].enext2self(); // [c,a,b]
          tt[1].enextself();  // [a,d,b]
          tt[2].enextself();  // [d,c,b]
          if (!check_only) {
            flip31(tt);
            pb->set_vrttype(0); // UNUSEDVERTEX
            unused_vertices++;
          }
          return 1;
        }
      } else { // ori2 == 0
        // Check if a 4-to-2 flip is possible.
        tt[2] = tt[1].enext2().sym();
        tt[3] = tt[0].enext().sym();
        if (tt[2].apex() == tt[3].apex()) {
          // Found a 4-to-2 flip.
          flipflag |= 4096; // set the 13th bit.
          // tt[0] is [a,b,c]
          // tt[1] is [b,a,d]
          // tt[2] is [b,d,f] <== use it
          if (!check_only) {
            //assert(0); // Debug this case.
            flip22(tt);
            // Handle segments.
            tt[0] = tt[2]; // [b,d,f]
            tt[1] = tt[0].enext2().sym(); // [b,f,c]
            tt[2] = tt[1].enext2().sym(); // [b,c,d] (degenerate)
// \iffalse
            assert(tt[0].org() == pb);
            assert(tt[1].org() == pb);
            assert(tt[2].org() == pb);
            assert(tt[0].dest() == pd);
            assert(tt[1].apex() == pc);
            assert(tt[2].apex() == pd);
            // Handle segments.
            assert(!tt[1].is_segment());
// \fi
            if (tt[0].is_segment()) {
              assert(tt[2].is_segment());
              // tt[0] and tt[2] will be deleted.
              // Do not need to clear/set their flags.
              tt[3] = tt[2].enext().sym(); // [d,c,a]
// \iffalse
              assert(tt[3].apex() == pa);
              assert(!tt[3].is_segment());
// \fi
              tt[3].set_segment();
            }
            for (int i = 0; i < 3; i++) {
              tt[i].enextself();
            }
            flip31(tt);
            if (pb->get_vrttype() == 6) { // FREESEGVERTEX
              st_segref_count--;
            }
            pb->set_vrttype(0); // UNUSEDVERTEX
            unused_vertices++;
          } // if (!check_only)
          return 1;
        }
      }
    } else {
      flipflag |= 4; // Constraint
    }
  } // ori2 <= 0

  return 0;
}

// If "flip_queue" is empty, then do lawson_flip on all edges of the current
//   triangulation. 
// If "pt != NULL", it is the newly inserted vertex.
int Triangulation::lawson_flip(Vertex *pt)
{
  TriEdge tt[4], *qedge;
  int flipcount = 0;
  int flipflag;
  int i;

  if (flip_queue->objects == 0l) {
    // Push all edges of the mesh into queue.
    int c = 0;
    for (i = 0; c < tris->objects; i++) {
      tt[0].tri = tris->get(i);
      if (!tris->is_deleted(tt[0].tri)) {
        tt[0].tri->set_infect();
        for (tt[0].ver = 0; tt[0].ver < 3; tt[0].ver++) {
          if (!(tt[0].sym().tri->is_infected())) {
            qedge = (TriEdge *) flip_queue->alloc();
            *qedge = tt[0];
          }
        }
        c++;
      }
    }
    // Uninfect all triangles.
    c = 0;
    for (i = 0; c < tris->objects; i++) {
      tt[0].tri = tris->get(i);
      if (!tris->is_deleted(tt[0].tri)) {
        tt[0].tri->clear_infect();
        c++;
      }
    }
  }

  if (verbose > 2) {
    printf("    Lawson flip %ld edges.\n", flip_queue->objects);
  }
  int bak_flag = lawsonflag;
  lawsonflag = 1;

  for (i = 0; i < flip_queue->objects; i++) {
    tt[0] = * (TriEdge *) (* flip_queue)[i];
    if (tris->is_deleted(tt[0].tri)) continue;
    if (flip(tt, flipflag, 0)) {
      flipcount++;
    }
  } // i

  if (verbose > 2) {
    printf("    Performed %d flips.\n", flipcount);
  }
  lawsonflag = bak_flag;
  flip_queue->restart();
  return flipcount;
}

// "tt" is an array of 4 TriEdges, used for 1-3, 2-2, 2-4 flips.
// If "loc != 0", tt[0] contains the triangle containing "pt".
int Triangulation::insert_point(Vertex *pt, TriEdge *tt, int loc)
{
  assert(pt->get_vrttype() == 0); // UNUSEDVERTE

  if (loc == 0) {
    loc = locate_point(pt, tt[0]); // Locate point.
  }

  if (loc == 1) {
    flip13(pt, tt);
    pt->set_vrttype(4); // FACETVERTEX
  } else if (loc == 2) {
    // Perform a 2-to-4 flip.
    tt[3] = tt[0].sym(); // The adjacent edge to be flipped.
    flip13(pt, tt);
    if (tt[3].is_segment()) {
      // Split a segment [a,b].
      assert(tt[0].is_segment());
      pt->set_vrttype(6); // FREESEGVERTEX
      pt->on = get_segment(tt[3]);
      tt[0].clear_segment();
      tt[3].clear_segment();
      tt[1] = tt[0].enext();  // Goto subsegment [b,p]
      tt[1].set_segment();
      (tt[1].sym()).set_segment();
      tt[2] = tt[0].enext2(); // Goto subsegment [p,a]
      tt[2].set_segment();
      (tt[2].sym()).set_segment();
    } else {
      pt->set_vrttype(4); // FACETVERTEX
    }
    tt[0] = tt[3];
    tt[1] = tt[0].sym();
    flip22(tt);
  } else if (loc == 3) {
    printf("!! Warning: Vertex %d is coincident with %d. Skipped.\n",
           pt->idx, (tt[0].org())->idx);
    //pt->set_vrttype(1);  // DUPLICATEDVERTEX
    return loc;
  } else if (loc == 4) {
    // Save an encroached segment.
    if (encsegflag) {
      EncEle *paryEle = (EncEle *) encsegs->alloc();
      paryEle->init();
      paryEle->end[0] = tt[0].org();
      paryEle->end[1] = tt[0].dest();
      paryEle->encpt = pt; // Save the encroached point.  
    }
    return loc;
  } else {
    assert(0); // Unknown case.
  }

  assert(pt->get_vrttype() != 0); // UNUSEDVERTEX
  unused_vertices--;

  if (lawsonflag) {
    lawson_flip(pt);
  }

  return loc;
}

// We assume that this vertex is either a free vertex or a Steiner point
//   that was inserted on a segment. In the latter case, the original
//   segment is recovered.
// On input, "tt" is an array of 4 TriEdges, no need to initialize them.
// On output, "tt[0]" is the new triangle resulted by flip31().
int Triangulation::remove_point(Vertex *pt, TriEdge *tt)
{
  TriEdge tlist[100];
  int tnum = 0, i;
  int ori1, ori2;

  if (verbose > 1) {
    printf("  Removing point [%d]\n", pt->idx);
  }
  while (1) {
    // Collect all triangles at this vertex.
    tlist[0] = pt->get_triorg();
    tnum = 0;
    do {
      tnum++; assert(tnum < 100);
      tlist[tnum] = tlist[tnum - 1];
      tlist[tnum].symself();
      tlist[tnum].enextself(); // CW direction
      assert(tlist[tnum].org() == pt);
    } while (tlist[tnum].tri != tlist[0].tri);

    if (tnum == 3) {
      // A 3-to-1 flip.
      tt[0] = tlist[0].enext();
      tt[1] = tlist[2].enext(); // For CCW orientation
      tt[2] = tlist[1].enext();
      flip31(tt);
	  break;
    }

    for (i = 0; i < tnum; i++) {
      if (!tlist[i].is_segment()) {
        tt[0] = tlist[i];
        tt[1] = tt[0].sym();
		if (pt->get_vrttype() != 6)
		{
			ori1 = Orient2d(tt[1].apex()->crd, tt[0].apex()->crd, tt[0].org()->crd);
			ori2 = Orient2d(tt[1].apex()->crd, tt[0].apex()->crd, tt[0].dest()->crd);
			if (ori1 * ori2 < 0) {
				flip22(tt);
				break;
			}
		}
		else
		{
			double OR1, OR2;
			OR1 = orient2d(tt[1].apex()->crd, tt[0].apex()->crd, tt[0].org()->crd);
			OR2 = orient2d(tt[1].apex()->crd, tt[0].apex()->crd, tt[0].dest()->crd);
			if (OR1*OR2 < -1e-12)
			{
				flip22(tt);
				break;
			}
		}
      }
    }

	if (i == tnum) {
		//system("pause");
      // No 2-to-2 flip is done.
	  if (tnum != 4)
	  {
		  pvert(pt);
		  printf("\n");
	      for (int k = 0; k < tnum; k++)
		  {
			  pvert(tlist[k].dest());
		  }
		  assert(0);
	  }

	  //printf("the removing point's idx is: %d\n", pt->idx);
      //assert(0); // To be debugged.  ADD "//" for remove point on straight line.
      // It is a 4-to-2 flip. And the vertex lies on the cross of two
      //   straight lines. We could select an arbitrary edge and do a
      //   2-to-2 flip. However, avoid to flip a subsegment.
      //   Let the link vertices (in CW order) be p0, p1, p2, p3;
      if (!tlist[0].is_segment()) { // ADD the condition whether tlist[0] or tlist[1] should be flip
        tt[0] = tlist[0];  // [pt, p0, p3]
        tt[3] = tlist[3];  // [pt, p3, p2]
      } else {
        assert(!tlist[1].is_segment());
        tt[0] = tlist[1];  // [pt, p1, p2]
        tt[3] = tlist[0];  // [pt, p0, p3] (is_segment)
      }
      tt[1] = tt[0].sym();
      flip22(tt);
      // Handle segments before do 3-to-1 flip.
      tt[0] = tt[3];
      tt[1] = tt[0].enext2().sym();
      tt[2] = tt[1].enext2().sym();
	  tlist[0] = tt[2].enext().sym();
      if (tt[0].is_segment()) {
        assert(tt[2].is_segment());
        if (verbose > 1) {
          printf("    From subsegment [%d,%d]\n",
                 tt[0].dest()->idx, tt[2].dest()->idx);
        }
        tt[3] = tt[2].enext();
		tt[3].symself();
        tt[0].clear_segment();
        tt[2].clear_segment();
        tt[3].set_segment();
      }
      for (i = 0; i < 3; i++) {
        tt[i].enextself();
      }
      flip31(tt);
	  tlist[1] = tlist[0].sym().enext();
	  tlist[2] = tlist[0].sym().enext2();
	  if (Orient2d(tlist[1].org()->crd, tlist[1].dest()->crd, tlist[2].dest()->crd) <= 0.0 && !tlist[0].is_segment())
	  {
		  tlist[1] = tlist[0].sym();
		  flip22(tlist);
	  }
	  
      break;
    }
  } // while (1)

  //vrts->dealloc_vrt(pt); // Delete the vertex.
  pt->set_vrttype(0); // UNUSEDVERTEX
  unused_vertices++;

  if (lawsonflag) {
    lawson_flip(NULL);
  }

  return 1;
}

//=============================================================================
// Delaunay triangulation

int Triangulation::first_tri(Vertex **ptlist, int ptnum, TriEdge &E)
{
  int i, j, it = 0;
  int ori = 0;
  Vertex *swappt;

  if (!norandomize && !sweep_sort) { // Default
    // Randomly select 3 points.
    do {
      // Randomly select three vertices from the iuput list.
      for (i = 0; i < 3; i++) {
        // Swap ith and jth element.
        j = rand() % (ptnum - i);
        swappt = ptlist[i];
        ptlist[i] = ptlist[j];
        ptlist[j] = swappt;
      }
      ori = Orient2d(ptlist[0]->crd, ptlist[1]->crd, ptlist[2]->crd);
      if (ori != 0) break;
      it++;
    } while (it < 10);
    if (it >= 10) return UNDEFIND; // Failed.
  } else { // with -R or -s
    // Assume pt[0] and pt[1] are distinct.
    // Search the third non-collinear vertex.
    for (i = 2; i < ptnum; i++) {
      ori = Orient2d(ptlist[0]->crd, ptlist[1]->crd, ptlist[i]->crd);
      if (ori != 0) break;
    }
    if (i == ptnum) return UNDEFIND; // Failed
    // If i > 2, swap it to 2.
    if (i > 2) {
      swappt = ptlist[2];
      ptlist[2] = ptlist[i];
      ptlist[i] = swappt;
    }
  }

  if (ori < 0) {
    // Swap the first two vertices.
    swappt = ptlist[0];
    ptlist[0] = ptlist[1];
    ptlist[1] = swappt;
  }

  if (verbose > 1) {
    printf("  First triangle [%d,%d,%d]\n", ptlist[0]->idx, ptlist[1]->idx,
           ptlist[2]->idx);
  }

  TriEdge tt[4];

  // Create 4 new triangles.
  for (i = 0; i < 4; i++) {
    tt[i].tri = tris->alloc_tri();
  }

  tt[0].set_vertices(ptlist[0], ptlist[1], ptlist[2]);
  tt[1].set_vertices(ptlist[1], ptlist[0], infvrt);
  tt[2].set_vertices(ptlist[2], ptlist[1], infvrt);
  tt[3].set_vertices(ptlist[0], ptlist[2], infvrt);

  for (i = 1; i < 4; i++) {
    tt[i].tri->set_hullflag();
  }
  hullsize += 3;

  for (i = 0; i < 3; i++) {
    ptlist[i]->set_vrttype(4); // FACETVERTEX
  }
  unused_vertices -= 3;

  // Connect the 4 triangles.
  tt[1].bond2(tt[0]);
  tt[2].bond2(tt[0].enext());
  tt[3].bond2(tt[0].enext2());

  (tt[1].enext2()).bond2(tt[2].enext());
  (tt[2].enext2()).bond2(tt[3].enext());
  (tt[3].enext2()).bond2(tt[1].enext());

  // Setup the vertex-to-triangle map.
  for (i = 0; i < 3; i++) {
    ptlist[i]->adj = tt[0].tri;
  }

  if (dump_intermediate_meshes) {
    save_to_ucd(0);
  }

  E = tt[0];
  return 1;
}

int Triangulation::locate_point(Vertex *pt, TriEdge &E)
{
  if (verbose > 1) {
    printf("  Locating point [%d]\n", pt->idx);
  }

  int i;

  if ((E.tri == NULL) || tris->is_deleted(E.tri)) {
    // Find a libe triangle to start.
    int c = 0;
    for (i = 0; c < tris->objects; i++) {
      E.tri = tris->get(i);
      if (!tris->is_deleted(E.tri)) {
        E.ver = 0; c++;
        break;
      }
    }
  }

  if (E.tri->is_hulltri()) {
    // Get a non-hull triangle.
    for (i = 0; i < 3; i++) {
      if (is_hulledge(E)) break;
      E.enextself();
    }
    assert(i < 3);
    E.symself();
    assert(!E.tri->is_hulltri());
  }

  // Select an edge such that pt lies to CCW of it.
  for (i = 0; i < 3; i++) {
    if (Orient2d((E.org())->crd, (E.dest())->crd, pt->crd) > 0) break;
    E.enextself();
  }

  // Let E = [a,b,c] and p lies to the CCW of [a->b].
  int ori1, ori2;

  do {

    ori1 = Orient2d((E.dest())->crd, (E.apex())->crd, pt->crd);
    ori2 = Orient2d((E.apex())->crd, ( E.org())->crd, pt->crd);

    if (ori1 > 0) {
      if (ori2 > 0) {
        break; // Found.
      } else if (ori2 < 0) {
        E.enext2self();
      } else { // ori2 == 0
        E.enext2self(); // p lies on edge [c,a]
        return 2; // ONEDGE
      }
    } else if (ori1 < 0) {
      if (ori2 > 0) {
        E.enextself();
      } else if (ori2 < 0) {
        // Randomly choose one.
        if (rand() % 2) { // flipping a coin.
          E.enextself();
        } else {
          E.enext2self();
        }
      } else { // ori2 == 0
        E.enextself();
      }
    } else { // ori1 == 0
      if (ori2 > 0) {
        E.enextself(); // p lies on edge [b,c].
        return 2; // ONEDGE
      } else if (ori2 < 0) {
        E.enext2self();
      } else { // ori2 == 0
        E.enext2self(); // p is coincident with apex.
        return 3; // ONVERTEX Org(E)
      }
    }

    if (encsegflag) {
      if (E.is_segment()) {
        return 4;
      }
    }

    E.symself();
  } while (!E.tri->is_hulltri());

  return 1; // P lies inside E = [a,b,c]
}

int Triangulation::incremental_Delaunay()
{
  TriEdge tt[4]; // 1-3, 2-2, 2-4 flips.
  int i;

  int ptnum = vrts->objects;
  Vertex **ptlist = new Vertex*[ptnum];
  for (i = 0; i < ptnum; i++) {
    ptlist[i] = vrts->get(i);
  }

  // Sort points;
  if (sweep_sort) { // -Gs
    sweepline_sort(ptlist, ptnum, sweep_sort_vx, sweep_sort_vy);
    lawson_hullflag = 1; // Only do hull flips.
  } else {
    if (!norandomize) { // -R
      randomly_permute(ptlist, ptnum);
    }
  }
  lawsonflag = 1;

  if (verbose) {
    printf("Incrementally inserting vertices.\n");
  }

  if (!first_tri(ptlist, ptnum, tt[0])) {
    return 0;
  }

  for (i = 3; i < ptnum; i++) {
    insert_point(ptlist[i], tt, 0);
  }

  lawsonflag = 0;
  if (sweep_sort) {
    lawson_hullflag = 0;
  }

  delete [] ptlist;
  return (int) tris->objects;
}

//=============================================================================
// Constrained Delaunay triangulation

int Triangulation::find_direction(TriEdge& E, Vertex *pt)
{
  int loc = 0;
  Vertex *pa = E.org();
  int ori1 = Orient2d(pa->crd, (E.dest())->crd, pt->crd);
  if (ori1 == 0) {
    // Randomly select an edge to start.
    if (rand() % 2) {
      E.enext2self();
      E.symself();
    } else {
      E.symself();
      E.enextself();
    }
    assert(E.org() == pa);
    ori1 = Orient2d(pa->crd, (E.dest())->crd, pt->crd);
    assert(ori1 != 0);
  }

  if (ori1 > 0) {
    while (!loc) {
      int ori2 = Orient2d(pa->crd, (E.apex())->crd, pt->crd);
      if (ori2 > 0) {
        E.enext2self();
        E.symself();
        assert(E.org() == pa);
      } else if (ori2 < 0) {
        loc = 1; // Cross edge [b,c]
      } else { // ori2 == 0
        E.enext2self();
        E.symself();
        loc = 2; // Cross vertex [c].
      }
    }
  } else { //if (ori1 < 0) {
    while(!loc) {
      E.symself();
      E.enextself();
      assert(E.org() == pa);
      int ori2 = Orient2d(pa->crd, (E.dest())->crd, pt->crd);
      if (ori2 > 0) {
        loc = 1; // Cross edge [b,c]
      } else if (ori2 < 0) {
        // continue.
      } else { // ori2 == 0
        loc = 2; // Cross vertex [b].
      }
    }
  }

  return loc;
}

int Triangulation::recover_edge_flips(Vertex *e1, Vertex *e2, TriEdge &E)
{
  TriEdge tt[2];
  tt[0].tri = e1->adj;
  for (tt[0].ver = 0; tt[0].ver < 3; tt[0].ver++) {
    if (tt[0].org() == e1) break;
  }
  assert(tt[0].ver < 3);

  if (verbose > 2) {
    printf("    Recovering edge [%d, %d]\n", e1->idx, e2->idx);
  }

  int loc = find_direction(tt[0], e2);

  if (loc == 1) {
    // Flip the edge [b,c]
    tt[0].enextself();
    if (tt[0].is_segment()) {
      printf("!! Two segments are intersecting:\n");
      printf("  1st [%d, %d]\n", e1->idx, e2->idx);
      printf("  2nd [%d, %d]\n", tt[0].org()->idx, tt[0].dest()->idx);
      return 0;
    }
    tt[1] = tt[0].sym();
    flip22(tt);
    // Continue to recover edge (switch e1<->e2).
    return recover_edge_flips(e2, e1, E);
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

int Triangulation::recover_segments()
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
    if (recover_edge_flips(paryseg->vrt[0], paryseg->vrt[1], E)) {
      // DEBUG BEGIN
      if (E.org() != paryseg->vrt[0]) E.symself();
      assert(E.org() == paryseg->vrt[0]);
      assert(E.dest() == paryseg->vrt[1]);
      // DEBUG END
      E.set_segment();
      (E.sym()).set_segment();
      E.org ()->set_vrttype(2);  // RIDGEVERTEX
      E.dest()->set_vrttype(2);  // RIDGEVERTEX
      lawson_flip(NULL);
    } else {
      printf("!! Failed to recover segment [%d,%d].\n", paryseg->vrt[0]->idx,
             paryseg->vrt[1]->idx);
      // assert(0);
    }
  }

  lawsonflag = 0;

  return 1;
}

// Create a map from vertex to segments incident at it.
// The map is returned in two arrays 'idx2seglist' and 'segpervrtlist'.  All
// segments incident at i-th vertex (i is counted from 0) are found in the
// array segperverlist[j], where idx2seglist[i] <= j < idx2seglist[i + 1].
// Each entry in segperverlist[j] is a segment containing the vertex.
int Triangulation::make_vertex_to_segment_map()
{
  if (verbose > 1) {
    printf("  Build the vertex-to-segments map.\n");
  }
  int i, j, k;

  assert(idx2seglist == NULL);
  assert(segperverlist == NULL);

  idx2seglist = new int[vrts->objects + 1];
  for (i = 0; i < vrts->objects + 1; i++) idx2seglist[i] = 0;

  // Loop all segments, counter the number of segments incident at a vertex.
  for (i = 0; i < segs->objects; i++) {
    Triangle *seg = segs->get(i);
    j = seg->vrt[0]->idx - firstindex;
    idx2seglist[j]++;
    j = seg->vrt[1]->idx - firstindex;
    idx2seglist[j]++;
    assert(seg->vrt[2] == NULL); // SELF_CHECK
  }

  // Calculate the total length of array 'facperverlist'.
  j = idx2seglist[0];
  idx2seglist[0] = 0;  // Array starts from 0 element.
  for (i = 0; i < vrts->objects; i++) {
    k = idx2seglist[i + 1];
    idx2seglist[i + 1] = idx2seglist[i] + j;
    j = k;
  }

  // The total length is in the last unit of idx2seglist.
  segperverlist = new Triangle*[idx2seglist[i]];

  // Loop all segments again, remember the subfaces at each vertex.
  for (i = 0; i < segs->objects; i++) {
    Triangle *seg = segs->get(i);
    j = seg->vrt[0]->idx - firstindex;
    segperverlist[idx2seglist[j]] = seg;
    idx2seglist[j]++;
    j = seg->vrt[1]->idx - firstindex;
    segperverlist[idx2seglist[j]] = seg;
    idx2seglist[j]++;
  }

  // Contents in 'idx2faclist' are shifted, now shift them back.
  for (i = vrts->objects - 1; i >= 0; i--) {
    idx2seglist[i + 1] = idx2seglist[i];
  }
  idx2seglist[0] = 0;

  return idx2seglist[i];
}

// Return the segment that this TriEdge lies on.
Triangle* Triangulation::get_segment(TriEdge &E)
{
  assert(E.is_segment()); // It must be marked as a segment.
  Vertex *e1 = E.org();
  Vertex *e2 = E.dest();
  if (e1->get_vrttype() == 6) {
    return e1->on; // FREESEGVERTEX
  } else if (e2->get_vrttype() == 6) {
    return e2->on; // FREESEGVERTEX
  } else {
    // Both endpoints must be RIDGEVERTEX
    assert(e1->get_vrttype() == 2);
    assert(e2->get_vrttype() == 2);
    // Search the segment from the vertex-to-segment map.
    int i = e1->idx - firstindex, j;
    for (j = idx2seglist[i]; j < idx2seglist[i+1]; j++) {
      Triangle *seg = segperverlist[j];
      if (((seg->vrt[0] == e1) && (seg->vrt[1] == e2)) ||
          ((seg->vrt[0] == e2) && (seg->vrt[1] == e1))) break;
    }
    assert(j < idx2seglist[i+1]);
    return segperverlist[j];
  }
}

//=============================================================================

int Triangulation::build_index_to_vertex_map()
{
  int build_map = 0;

  if (idx_to_vert_list == NULL) {
    // Build the index-to-vertex map.
    build_map = 1;
  } else {
    if (idx_to_vert_list_len != vrts->objects) {
      // The current list is invalid, re-build it.
      delete [] idx_to_vert_list;
      build_map = 1;
    }
  }

  if (build_map) {
    // Make a map from idx to vertices.
    idx_to_vert_list = new Vertex*[vrts->objects + firstindex];
    for (int i = 0; i < vrts->objects; i++) {
      assert(vrts->get(i)->idx == (i + firstindex));
      idx_to_vert_list[i + firstindex] = vrts->get(i);
    }
    idx_to_vert_list_len = vrts->objects;
  }
  return 0;
}

//=============================================================================

// Removing exterior (and hole) elements, and recreating the hull connections. 
int Triangulation::remove_exterior_elements()
{
  arraypool *deadtris = new arraypool(sizeof(Triangle **), 8);
  TriEdge E, N;
  int c, i;

  // Loop through triangles, Mark all hull triangles.
  // Since some triangles may be deleted, use c to count the live triangle.
  for (i = 0, c = 0; c < tris->objects; i++) {
    Triangle* tri = tris->get(i);
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
    // Locate the point.
    int loc = locate_point(pt, E);
    if (loc == 1) { // Inside triangle
      if (!E.tri->is_hulltri()) {
        if ((pt->tag == -9999) || (pt->tag == 0) || 1) {//ADD "|| 1" for set_infect always
          E.tri->set_infect();
          pparytri = (Triangle **) deadtris->alloc();
          *pparytri = E.tri;
        } else {
          pt->adj = E.tri;
        }
      }
    } else if (loc == 2) { // On edge.
      N = E.sym();
      if ((!E.tri->is_hulltri()) && (!N.tri->is_hulltri())) {
        // It is an interior point.
        if ((pt->tag == -9999) || (pt->tag == 0) || 1) {//ADD "|| 1" for set_infect always
          E.tri->set_infect();
          pparytri = (Triangle **) deadtris->alloc();
          *pparytri = E.tri;
        } else {
          pt->adj = E.tri;
        }
      }
    } else if (loc == 3) { // On vertex
      if ((pt->tag == -9999) || (pt->tag == 0) || 1) {//ADD "|| 1" for set_infect always
        E.tri->set_infect();
        pparytri = (Triangle **) deadtris->alloc();
        *pparytri = E.tri;
      } else {
        pt->adj = E.tri;
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

  if (!convex) {  // Not -c
    // Delete all exterior triangles (including hull triangles)
    for (i = 0; i < deadtris->objects; i++) {
      Triangle** pparytri = (Triangle **)(* deadtris)[i];
      tris->dealloc_tri(*pparytri);
    }

    // Recreate the hull triangles.
    arraypool *hulltris = new arraypool(sizeof(TriEdge), 8);
    TriEdge htri;

    for (i = 0, c = 0; c < tris->objects; i++) {
      E.tri = tris->get(i);
      if (!tris->is_deleted(E.tri)) {
        for (E.ver = 0; E.ver < 3; E.ver++) {
          if (E.is_segment()) {
            N = E.sym();
            if (tris->is_deleted(N.tri)) {
              TriEdge *parytri = (TriEdge *) hulltris->alloc();
              *parytri = E;
            }
          }
        }
        c++;
      }
    }

    hullsize = (int) hulltris->objects; // Update the hullsize.

    // Create and connect hull triangles to neighbor (non-hull) triangles.
    for (i = 0; i < hulltris->objects; i++) {
      TriEdge *parytri = (TriEdge *) (* hulltris)[i];
      assert(parytri->is_segment()); // Must be a segment.
      htri.tri = tris->alloc_tri();
      htri.tri->set_hullflag();
      htri.set_vertices(parytri->dest(), parytri->org(), infvrt);
      htri.bond2(*parytri);
      htri.set_segment();
      // Update the vertex-to-triangle map.
      (htri.org()) ->adj = htri.tri;
      (htri.dest())->adj = htri.tri;
      *parytri = htri; // Save the new hull triangle.
    }

    // Connect hull triangles to each other.
    for (i = 0; i < hulltris->objects; i++) {
      TriEdge *parytri = (TriEdge *) (* hulltris)[i];
      E = parytri->enext();
      if (E.tri->nei[E.ver] == NULL) {
        N = parytri->sym();
        N.enext2self();
        while (N.tri->nei[N.ver]) {
          N.symself();
          N.enext2self();
        }
        E.bond2(N);
      }
      E = parytri->enext2();
      if (E.tri->nei[E.ver] == NULL) {
        N = parytri->sym();
        N.enextself();
        while (N.tri->nei[N.ver]) {
          N.symself();
          N.enextself();
        }
        E.bond2(N);
      }
    }

    // Some exterior vertices may be removed as well.
    for (i = 0, c = 0; c < vrts->objects; i++) {
      Vertex *vrt = (Vertex *)(* vrts)[i];
      if (!vrts->is_deleted(vrt)) {
        if (vrt->get_vrttype() != 0) { // NOT UNUSEDVERTEX
          if (tris->is_deleted(vrt->adj) ||
              ((vrt->adj->vrt[0] != vrt) &&
               (vrt->adj->vrt[1] != vrt) &&
               (vrt->adj->vrt[2] != vrt))) {
            vrt->set_vrttype(0); // UNUSEDVERTEX
            unused_vertices++;
          }
        }
        c++;
      }
    }

    delete hulltris;
  } else { // -c convex
    // Mark all exterior triangles.
    int count_exterior_triangles = 0;
    for (i = 0; i < deadtris->objects; i++) {
      Triangle** pparytri = (Triangle **)(* deadtris)[i];
      (*pparytri)->clear_infect();
      if (!(*pparytri)->is_hulltri()) {
        (*pparytri)->set_exterior();
        count_exterior_triangles++;
      }
    }
    assert(count_exterior_triangles == deadtris->objects);
  }

  delete deadtris;

  // Create regions 
  int regionattrib = 1; int j;
  arraypool *hulltris = new arraypool(sizeof(TriEdge), 8);

  for (i = 0; i < regions->objects; i++) {
    Vertex *pt = (Vertex *)(* regions)[i];
    if ((pt->tag != -9999) && (pt->tag != 0) && 0) {//ADD "&& 0" for always do not creat region
      E.tri = pt->adj; E.ver = 0;
      E.tri->tag = pt->tag;
      E.tri->set_infect();
      TriEdge *parytri = (TriEdge *) hulltris->alloc();
      *parytri = E;
      for (j = 0; j < hulltris->objects; j++) {
        parytri = (TriEdge *) (* hulltris)[j];
        for (parytri->ver = 0; parytri->ver < 3; parytri->ver++) {
          if (!parytri->is_segment()) {
            E = parytri->sym();
            if (!E.tri->is_infected()) {
              E.tri->tag = pt->tag;
              E.tri->set_infect();
              TriEdge *parytri2 = (TriEdge *) hulltris->alloc();
              *parytri2 = E;
            }
          }
        }
      }
      for (j = 0; j < hulltris->objects; j++) {
        parytri = (TriEdge *) (* hulltris)[j];
        parytri->tri->clear_infect();
      }
      hulltris->restart();
      // Record the largest positive integer.
      if (pt->tag > regionattrib) {
        regionattrib = pt->tag + 1;
      }
    }
  }

  delete hulltris;
  return 1;
}

//=============================================================================

EncEle* Triangulation::enq_triangle(arraypool* Eleque, TriEdge& E)
{
  EncEle *paryEle = (EncEle *) Eleque->alloc();
  paryEle->init(); // Initialize
  paryEle->Ele = E;
  for (int j = 0; j < 3; j++) {
    paryEle->end[j] = E.tri->vrt[j];
  }
  return paryEle;
}

void Triangulation::enq_vertex_star(Vertex *pt)
{
  TriEdge E = pt->get_triorg();
  TriEdge N = E;
  do {
    enq_triangle(chktris, N);
    N.symself(); N.enextself();
  } while (N.tri != E.tri);
}

// If 'pt' is given, only check if S is encroached by pt.
// Otherwise, check the apex of S.
// In Delaunay refinement case, 'pt' must be given.
int Triangulation::check_segment(Vertex *pt, TriEdge &S)
{
  assert(S.is_segment());

  Vertex *e1 = S.org();
  Vertex *e2 = S.dest();
  assert(e1 != infvrt);
  assert(e2 != infvrt);

  Vertex *chkpt;
  if (pt != NULL) {
    chkpt = pt;
  } else {
    chkpt = S.apex();
  }

  if (chkpt != infvrt) {
    double v1[2], v2[2];
    double costh = 0;
    v1[0] = e1->crd[0] - chkpt->crd[0];
    v1[1] = e1->crd[1] - chkpt->crd[1];
    v2[0] = e2->crd[0] - chkpt->crd[0];
    v2[1] = e2->crd[1] - chkpt->crd[1];
    costh = v1[0] * v2[0] + v1[1] * v2[1];
    if (costh < 0) {
      if (encsegflag) {
        EncEle *paryEle = (EncEle *) encsegs->alloc();
        paryEle->init();
        paryEle->end[0] = e1;
        paryEle->end[1] = e2;
        paryEle->encpt = chkpt; // Save the encroached point.
      }
      return 1;
    }
  }

  return 0;
}

int Triangulation::get_edge(Vertex *e1, Vertex *e2, TriEdge &N)
{
  TriEdge E = e1->get_triorg();
  N = E;
  do {
    if (N.dest() == e2) return 1;
    N.symself(); N.enextself();
  } while (N.tri != E.tri);
  return 0;
}

REAL distance(Vertex *e1, Vertex *e2)
{
  REAL vx = e2->crd[0] - e1->crd[0];
  REAL vy = e2->crd[1] - e1->crd[1];
  return sqrt(vx * vx + vy * vy);
}

Vertex* Triangulation::get_split_segment_vertex(EncEle *S)
{
  assert(S->Ele.is_segment());

  if (S->encpt != NULL) {
    if (!vrts->is_deleted(S->encpt) &&
        (S->encpt->get_vrttype() != 0)) { // NOT UNUSEDVERTEX
      if (S->encpt->on != NULL) {
        assert(S->encpt->get_vrttype() == 6); // FREESEGVERTEX
        // Check are these two segments incident.
        Triangle *seg = get_segment(S->Ele);
        if (seg != S->encpt->on) {
          // Check if these two segments shara a common vertex.
          if ((seg->vrt[0] == S->encpt->on->vrt[0]) ||
              (seg->vrt[0] == S->encpt->on->vrt[1])) {
            // Shared at seg->vrt[0]
            REAL split = distance(seg->vrt[0], S->encpt) / 
                         distance(seg->vrt[0], seg->vrt[1]);
            Vertex *newvrt = vrts->alloc_vrt();
            newvrt->crd[0] = seg->vrt[0]->crd[0] + 
              split * (seg->vrt[1]->crd[0] - seg->vrt[0]->crd[0]);
            newvrt->crd[1] = seg->vrt[0]->crd[1] + 
              split * (seg->vrt[1]->crd[1] - seg->vrt[0]->crd[1]);
            newvrt->idx = firstindex + (vrts->objects - 1);
            unused_vertices++;
            return newvrt;
          } else if ((seg->vrt[1] == S->encpt->on->vrt[0]) ||
                     (seg->vrt[1] == S->encpt->on->vrt[1])) {
            // Shared at seg->vrt[1]
            REAL split = distance(seg->vrt[1], S->encpt) / 
                         distance(seg->vrt[0], seg->vrt[1]);
            Vertex *newvrt = vrts->alloc_vrt();
            newvrt->crd[0] = seg->vrt[1]->crd[0] + 
              split * (seg->vrt[0]->crd[0] - seg->vrt[1]->crd[0]);
            newvrt->crd[1] = seg->vrt[1]->crd[1] + 
              split * (seg->vrt[0]->crd[1] - seg->vrt[1]->crd[1]);
            newvrt->idx = firstindex + (vrts->objects - 1);
            unused_vertices++;
            return newvrt;
          }
        }
      }
    }
  }

  // Split the segment at the middle.
  Vertex *e1 = S->Ele.org();
  Vertex *e2 = S->Ele.dest();
  Vertex *newvrt = vrts->alloc_vrt();
  newvrt->crd[0] = 0.5 * (e1->crd[0] + e2->crd[0]);
  newvrt->crd[1] = 0.5 * (e1->crd[1] + e2->crd[1]);
  newvrt->idx = firstindex + (vrts->objects - 1);
  unused_vertices++;
  return newvrt;
}

int Triangulation::repair_encsegments()
{
  TriEdge tt[4];
  int loc;

  while (encsegs->objects > 0l) {
    assert(encsegs2->objects == 0l);
    // Swap these two arrays.
    arraypool *swapary = encsegs;
    encsegs = encsegs2;
    encsegs2 = swapary;
    // Process chktris2, new tris are queued in chktris.
    for (int i = 0; i < encsegs2->objects; i++) {
      EncEle *paryEle = (EncEle *) (* encsegs2)[i];
      if (get_edge(paryEle->end[0], paryEle->end[1], paryEle->Ele)) {
        Vertex *newvrt = get_split_segment_vertex(paryEle);
        tt[0] = paryEle->Ele;
        loc = insert_point(newvrt, tt, 2); // loc = 2, ONEDGE
        assert(loc == 2);
        st_segref_count++;
        if (enqtriflag) {
          enq_vertex_star(newvrt);
        }
      }
    }
    encsegs2->restart();
  } // while

  return 1;
}

///////////////////////////////////////////////////////////////////////////////

// Get the circumcenter (or orthocenter) of 3 (weighted) vertices, U, V, W.
int Triangulation::get_orthocenter(REAL Ux, REAL Uy, REAL U_weight,
                                   REAL Vx, REAL Vy, REAL V_weight,
                                   REAL Wx, REAL Wy, REAL W_weight,
                                   REAL* Cx, REAL* Cy, REAL* r2)
{
  REAL Px = Vx - Ux;
  REAL Py = Vy - Uy;
  REAL Qx = Wx - Ux;
  REAL Qy = Wy - Uy;
  REAL Px2_plus_Py2 = Px * Px + Py * Py + U_weight - V_weight; // + w_0 - w_p
  REAL Qx2_plus_Qy2 = Qx * Qx + Qy * Qy + U_weight - W_weight; // + w_0 - w_q

  REAL det = Px * Qy - Qx * Py;
  if ((fabs(det) / (Px2_plus_Py2 + Qx2_plus_Qy2)) < 1e-8) {
    // The triangle is nearly degenerated.
    if (verbose) {
      printf("!! Warning: triangle is degenerated.\n");
    }
    return 0;
  }

  REAL denominator = 0.5 / det; //  0.5 / (Px * Qy - Qx * Py);
  REAL dx = (Qy * Px2_plus_Py2 - Py * Qx2_plus_Qy2) * denominator;
  REAL dy = (Px * Qx2_plus_Qy2 - Qx * Px2_plus_Py2) * denominator;

  // The power (of radius).
  REAL rr2 = dx * dx + dy * dy - U_weight;
  if (rr2 < 0) {
    //if (verbose) {
    //  printf("!! Warning: negative weights r^2 = %g.\n", rr2);
    //}
  }

  if (Cx != NULL) {
    *Cx = Ux + dx;
    *Cy = Uy + dy;
    *r2 = rr2;
  } else {
    // Only print the result (for debug).
    REAL cx, cy, r, h;
    cx = Ux + dx;
    cy = Uy + dy;
    r = sqrt(fabs(rr2));
    h = cx*cx + cy*cy - rr2;
    printf("  (Cx, Cy): %g %g r(%g) r2(%g) h(%g)\n", cx, cy, r, rr2, h);
  }

  return rr2 > 0;
}

int Triangulation::get_tri_orthocenter(Vertex *v1, Vertex *v2, Vertex *v3,
                                       REAL* Cx, REAL* Cy, REAL* r2)
{
  REAL Ux, Uy, U_weight;
  REAL Vx, Vy, V_weight;
  REAL Wx, Wy, W_weight;

  Ux = v1->crd[0];
  Uy = v1->crd[1];
  Vx = v2->crd[0];
  Vy = v2->crd[1];
  Wx = v3->crd[0];
  Wy = v3->crd[1];

  if (weighted) { // -w
    if (weight_is_height) { // -wh
      // Translate heights to weights.
      U_weight = Ux*Ux + Uy*Uy - v1->crd[2];
      V_weight = Vx*Vx + Vy*Vy - v2->crd[2];
      W_weight = Wx*Wx + Wy*Wy - v3->crd[2];
    } else {
      U_weight = v1->crd[2];
      V_weight = v2->crd[2];
      W_weight = v3->crd[2];
    }
  } else {
    U_weight = V_weight = W_weight = 0;
  }

  if (verbose > 2) {
    printf(" [%d]: %g,%g,%g, r(%g)\n", v1->idx, Ux, Uy, U_weight, sqrt(fabs(U_weight)));
    printf(" [%d]: %g,%g,%g, r(%g)\n", v2->idx, Vx, Vy, V_weight, sqrt(fabs(V_weight)));
    printf(" [%d]: %g,%g,%g, r(%g)\n", v3->idx, Wx, Wy, W_weight, sqrt(fabs(W_weight)));
  }

  return get_orthocenter(Ux, Uy, U_weight,
                         Vx, Vy, V_weight,
                         Wx, Wy, W_weight,
                         Cx, Cy, r2);
}

int Triangulation::get_seg_orthocenter(Vertex *e1, Vertex *e2,
                                       REAL* Cx, REAL* Cy, REAL* r2)
{
  REAL Ux, Uy, U_weight;
  REAL Vx, Vy, V_weight;

  Ux = e1->crd[0];
  Uy = e1->crd[1];
  Vx = e2->crd[0];
  Vy = e2->crd[1];

  REAL Qx = Vx - Ux;
  REAL Qy = Vy - Uy;
  REAL t = 0;

  if (weighted) { // -w
    if (weight_is_height) { // -wh
      // Translate heights to weights.
      U_weight = Ux*Ux + Uy*Uy - e1->crd[2];
      V_weight = Vx*Vx + Vy*Vy - e2->crd[2];
    } else {
      U_weight = e1->crd[2];
      V_weight = e2->crd[2];
    }
    REAL Qx2_plus_Qy2 = Qx*Qx + Qy*Qy;
    t = 0.5 * (1 + (U_weight - V_weight) / Qx2_plus_Qy2);
  } else {
    U_weight = V_weight = 0;
    t = 0.5; // Just middle point.
  }

  REAL rr2 = t*t * (Qx*Qx + Qy*Qy) - U_weight;
  if (rr2 < 0) {
    //if (verbose) {
    //  printf("!! Warning: negative weights r^2 = %g.\n", rr2);
    //}
  }

  if (Cx != NULL) {
    *Cx = Ux + t * Qx;
    *Cy = Uy + t * Qy;
    *r2 = rr2;
  } else {
    // Only print the result.
    REAL cx, cy, r, h;
    cx = Ux + t * Qx;
    cy = Uy + t * Qy;
    r = sqrt(fabs(rr2));
    h = cx*cx + cy*cy - rr2;
    printf("  (Cx, Cy): %g %g r(%g) r2(%g) h(%g)\n", cx, cy, r, rr2, h);
  }

  return 1;
}

REAL Triangulation::get_tri_area(Vertex* pa, Vertex* pb, Vertex* pc)
{
  REAL detleft = (pa->crd[0] - pc->crd[0]) * (pb->crd[1] - pc->crd[1]);
  REAL detright = (pa->crd[1] - pc->crd[1]) * (pb->crd[0] - pc->crd[0]);
  REAL det = detleft - detright;
  return 0.5 * fabs(det);
}

// Return the area of the Voronoi cell dual to the vertex pt.
REAL Triangulation::get_voronoi_cell_area(Vertex* pt)
{
  if (dual_verts == NULL) {
    // Computer the set of dual (Voronoi) vertices (only once).
    // A (non hull) triangle is dual to a vertex (its circumcenter).
    // A hull triangle is dual to a vertex on its hull edge.
    dual_verts = new Vertex[tris->objects];

    Vertex *e1, *e2;
    double Cx, Cy, radius2;
    int i, idx = 1; // UCD assumes vertex index starts from 1.
    int count = 0;

    // Calculate dual vertices (orthocenters)
    for (i = 0; count < tris->objects; i++) {
      Triangle* tri = tris->get(i);
      if (!tris->is_deleted(tri)) {
        if (!tri->is_hulltri()) {
          get_tri_orthocenter(tri->vrt[0], tri->vrt[1], tri->vrt[2],
                              &Cx, &Cy, &radius2);
        } else {
          if (tri->vrt[0] == infvrt) {
            e1 = tri->vrt[1]; e2 = tri->vrt[2];
          } else if (tri->vrt[1] == infvrt) {
            e1 = tri->vrt[2]; e2 = tri->vrt[0];
          } else { //
            e1 = tri->vrt[0]; e2 = tri->vrt[1];
          }
          get_seg_orthocenter(e1, e2, &Cx, &Cy, &radius2);
        }
        dual_verts[count].crd[0] = Cx;
        dual_verts[count].crd[1] = Cy;
        dual_verts[count].crd[2] = radius2;
        dual_verts[count].idx = idx;
        tri->idx = count;
        count++; idx++;
      }
    }
  } // if (dual_verts == NULL)

  TriEdge tt[4];
  REAL area = 0.0;

  tt[0] = pt->get_triorg();
  tt[1] = tt[0];
  do {
    tt[2] = tt[1]; // Bakup the current triangle.
    // Go to the next triangle (at its counter-clockwise direction).
    tt[1].enext2self();
    tt[1].symself();
    // Only skip the case when both triangles are hull triangles.
    if (!tt[1].tri->is_hulltri() || !tt[2].tri->is_hulltri()) {
      // Incremental the area by adding a triangle.
      Vertex *v1 = &dual_verts[tt[1].tri->idx];
      Vertex *v2 = &dual_verts[tt[2].tri->idx];
      area += get_tri_area(pt, v1, v2);
    }
  } while (tt[1].tri != tt[0].tri);

  return area;
}

///////////////////////////////////////////////////////////////////////////////

int Triangulation::check_triangle(Triangle *tri, double ccent[3])
{
  // Skip a hull and exterior triangle.
  if (tri->is_hulltri() || 
      tri->is_exterior()) return 0; 

  // torg  -> tri->vrt[0]->crd
  // tdest -> tri->vrt[1]->crd
  // tapex -> tri->vrt[2]->crd
  REAL xdo = tri->vrt[1]->crd[0] - tri->vrt[0]->crd[0];
  REAL ydo = tri->vrt[1]->crd[1] - tri->vrt[0]->crd[1];
  REAL xao = tri->vrt[2]->crd[0] - tri->vrt[0]->crd[0];
  REAL yao = tri->vrt[2]->crd[1] - tri->vrt[0]->crd[1];
  REAL dodist = xdo * xdo + ydo * ydo;
  REAL aodist = xao * xao + yao * yao;
  REAL dadist = (tri->vrt[1]->crd[0] - tri->vrt[2]->crd[0]) *
                (tri->vrt[1]->crd[0] - tri->vrt[2]->crd[0]) +
                (tri->vrt[1]->crd[1] - tri->vrt[2]->crd[1]) *
                (tri->vrt[1]->crd[1] - tri->vrt[2]->crd[1]);
  REAL denominator = 0.5 / (xdo * yao - xao * ydo); // area.
  REAL dx = (yao * dodist - ydo * aodist) * denominator;
  REAL dy = (xdo * aodist - xao * dodist) * denominator;

  ccent[0] = tri->vrt[0]->crd[0] + dx;
  ccent[1] = tri->vrt[0]->crd[1] + dy;

  if (maxarea > 0) {
    REAL area = 0.5 * (xdo * yao - xao * ydo);
    if (area > maxarea) {
      return 1; // Split triangle.
    }
	if (area < maxarea*0.01) //ADD a lower bound
		return 0;
  }

  TriEdge E; E.tri = tri;
  REAL minlen2;
  if (dodist < aodist) {
    if (dodist < dadist) { 
      E.ver = 2; // org->dest
      minlen2 = dodist;
    } else {
      E.ver = 0; // dest->apex
      minlen2 = dadist;
    }
  } else {
    if (aodist < dadist) {
      E.ver = 1; // apex->org
      minlen2 = aodist;
    } else {
      E.ver = 0; // dest->apex
      minlen2 = dadist;
    }
  }
  if (minedgelen > 0) {
    // TO DO: Check if we need to delete a point.
  }
  REAL radius2 = dx * dx + dy * dy;
  REAL ratio2 = radius2 / minlen2;
  if (ratio2 > maxratio2) {
    // Check if there are two segments.
    TriEdge E1 = E.enext();
    TriEdge E2 = E.enext2();
    if (E1.is_segment() && E2.is_segment()) {
      // Do not split a triangle containing a small angle.
      // Miller, Pav, and Walkington's adaptive version.
      return 0;
    }
    return 1;
  }

  return 0;
}

// Assume all triangles to be checked are in chktris, and are infected.
int Triangulation::repair_triangles()
{
  TriEdge tt[4];
  double ccent[3];
  int loc, i;
  while (chktris->objects > 0l) {
    assert(chktris2->objects == 0);
    // Swap the two arrays.
    arraypool *swapary = chktris;
    chktris = chktris2;
    chktris2 = swapary;
    // Process chktris2, new tris are queued in chktris.
    for (i = 0; i < chktris2->objects; i++) {
      EncEle *paryEle = (EncEle *)(* chktris2)[i];
      if (!paryEle->is_changed()) {
        if (check_triangle(paryEle->Ele.tri, ccent)) {
          Vertex *newvrt = vrts->alloc_vrt();
          newvrt->crd[0] = ccent[0];
          newvrt->crd[1] = ccent[1];
          newvrt->idx = firstindex + (vrts->objects - 1);
          unused_vertices++;
          assert(encsegs->objects == 0l);
          tt[0].tri = paryEle->Ele.tri;
          loc = insert_point(newvrt, tt, 0);
          if (encsegs->objects > 0l) {
            if (loc != 4) {
              remove_point(newvrt, tt);
            } else {
              unused_vertices--;
              vrts->dealloc_vrt(newvrt); // Delete the vertex.
            }
            repair_encsegments();
          } else {
            st_triref_count++;
            if (enqtriflag) {
              enq_vertex_star(newvrt);
            }
          }
        } // if (check_triangle
      } // if (!paryEle->is_changed())
    }
    // clear the array.
    chktris2->restart();
  } // while

  return 1;
}

int Triangulation::delaunay_refinement()
{
  if (verbose) {
    printf("Adding Steiner points to enforce quality.\n");
  }
  make_vertex_to_segment_map();
  encsegs  = new arraypool(sizeof(EncEle), 4);
  encsegs2 = new arraypool(sizeof(EncEle), 4);

  lawsonflag = 1;
  encsegflag = 1; 

  if (gabriel) { // -G
    // Put all encroached segments into queue.
    TriEdge E;
    int count = 0;
    for (int i = 0; count < tris->objects; i++) {
      E.tri = tris->get(i);
      if (!tris->is_deleted(E.tri)) {
        if (!E.tri->is_hulltri()) { 
          for (E.ver = 0; E.ver < 3; E.ver++) {
            if (E.is_segment()) {
              check_segment(NULL, E);
            }
          }
        }
        count++;
      }
    }

    repair_encsegments();
  }

  if (quality) { // -q
    chktris  = new arraypool(sizeof(EncEle), 8);
    chktris2 = new arraypool(sizeof(EncEle), 8);

    // Calculate the squre of the maximum radius-edge ratio.
    // r/d = 1/(2*sin(theta_min));
    maxratio2 = 0.5 / sin(PI * minangle / 180.0);
    maxratio2 *= maxratio2;

    enqtriflag = 1;

    // Put all triangles into queue.
    TriEdge E;
    int count = 0;
    for (int i = 0; count < tris->objects; i++) {
      E.tri = tris->get(i);
      if (!tris->is_deleted(E.tri)) {
        enq_triangle(chktris, E);
        count++;
      }
    }

    repair_triangles();

    enqtriflag = 0;

    delete chktris;
    delete chktris2;
    chktris  = NULL;
    chktris2 = NULL;
  }

  lawsonflag = 0;
  encsegflag = 0; 

  delete encsegs;
  delete encsegs2;
  encsegs  = NULL;
  encsegs2 = NULL;

  delete [] idx2seglist;
  delete [] segperverlist;
  idx2seglist = NULL;
  segperverlist = NULL;
  return 1;
}

//=============================================================================

// It assumes that the vertcies and triangles have been read.
//   Segments may be read (and will be created automatically).
//   Triangles are not connected to each other yet.
int Triangulation::reconstructmesh()
{
  // Allocate an array for creating the vertex-to-triangles map.
  // For each triangle, we allocate additional 3+1 pointers to store the
  //   three link lists (stacks) of triangle sharing at its three vertices.
  //   The extra pointer is to make the calulation faster (use bit operation)?.
  //   The total size of this array is: number of triangles times 3.
  Triangle **tristacks = new Triangle*[tris->objects * 4];
  // We will use the tags as index pointing to tristacks.
  int *bak_tags = new int[tris->objects];
  TriEdge E, N;
  int i;

  for (i = 0; i < tris->objects; i++) {
    E.tri = tris->get(i);
    bak_tags[i] = E.tri->tag;
    E.tri->tag = i; // Set the index-to-tristacks.
    for (E.ver = 0; E.ver < 3; E.ver++) {
      N.tri = E.org()->adj;
      if (N.tri == NULL) {
        // This vertex is used.
        E.org()->set_vrttype(4); // FACETVERTEX
        unused_vertices--;
      }
      // Link the current triangle to the next one in the stack.
      tristacks[(i << 2) + E.ver] = N.tri;
      // Push the current triangle onto the stack.
      E.org()->adj = E.tri;
      // Find a neighbor at E.org() and E.dest().
      while (N.tri != NULL) {
        for (N.ver = 0; N.ver < 3; N.ver++) {
          if (N.org() == E.org()) break; // Find the origin in N.
        }
        assert(N.ver < 3);
        // E and N have the same origin.
        if (E.dest() == N.apex()) {
          if (E.sym().tri == NULL) {
            assert(N.enext2().sym().tri == NULL);
            E.bond2(N.enext2());
          }
        } else if (E.apex() == N.dest()) {
          if (N.sym().tri == NULL) {
            assert(E.enext2().sym().tri == NULL);
            E.enext2().bond2(N);
          }
        }
        // Find the next triangle in the stack.
        N.tri = tristacks[(N.tri->tag << 2) + N.ver];
      }
    } // E.ver
  } // i

  delete [] tristacks;

  // Create hull triangles, and connect them together.
  arraypool *hulltris = new arraypool(sizeof(TriEdge), 8);

  for (i = 0; i < tris->objects; i++) {
    E.tri = tris->get(i);
    E.tri->tag = bak_tags[i]; // Restore the tags.
    for (E.ver = 0; E.ver < 3; E.ver++) {
      if (E.sym().tri == NULL) {
        // Save the hull triangle.
        TriEdge *parytri = (TriEdge *) hulltris->alloc();
        *parytri = E;
      }
    }
  }

  delete [] bak_tags;

  assert(hullsize == 0);
  hullsize = hulltris->objects;

  // Create the hull triangles.
  for (i = 0; i < hulltris->objects; i++) {
    TriEdge *parytri = (TriEdge *) (* hulltris)[i];
    E = *parytri;
    N.tri = tris->alloc_tri();
    N.tri->set_hullflag();
    N.set_vertices(E.dest(), E.org(), infvrt);
    N.bond2(E);
    *parytri = N; // Save the hull triangle.
  }

  // Connect hull triangles to each other.
  for (i = 0; i < hulltris->objects; i++) {
    TriEdge *parytri = (TriEdge *) (* hulltris)[i];
    E = parytri->enext();
    if (E.tri->nei[E.ver] == NULL) {
      N = parytri->sym();
      N.enext2self();
      while (N.tri->nei[N.ver]) {
        N.symself();
        N.enext2self();
      }
      E.bond2(N);
    }
    E = parytri->enext2();
    if (E.tri->nei[E.ver] == NULL) {
      N = parytri->sym();
      N.enextself();
      while (N.tri->nei[N.ver]) {
        N.symself();
        N.enextself();
      }
      E.bond2(N);
    }
  }

  delete hulltris;

  // Segments will be introduced.

  // First mark all segments which are already given.
  for (i = 0; i < segs->objects; i++) {
    Vertex *e1 = segs->get(i)->vrt[0];
    Vertex *e2 = segs->get(i)->vrt[1];
    if (get_edge(e1, e2, E)) {
      E.set_segment();
      E.sym().set_segment();
      e1->set_vrttype(2); // RIDGEVERTEX
      e2->set_vrttype(2); // RIDGEVERTEX
    } else {
      printf("!! Input segment [%d,%d] does not match mesh edge. Skipped.\n",
             e1->idx, e2->idx);
    }
  }

  // Create segments at exterior and interior boundary.
  for (i = 0; i < tris->objects; i++) {
    if (!tris->get(i)->is_hulltri()) {
      E.tri = tris->get(i);
      for (E.ver = 0; E.ver < 3; E.ver++) {
        if (!E.is_segment()) {
          N = E.sym();
          if (N.tri->is_hulltri() || (E.tri->tag != N.tri->tag)) {
            assert(!N.is_segment());
            E.set_segment();
            N.set_segment();
            // Create a new segment.
            Triangle *seg = segs->alloc_tri();
            seg->vrt[0] = E.org();
            seg->vrt[1] = E.dest();
            seg->vrt[0]->set_vrttype(2); // RIDGEVERTEX
            seg->vrt[1]->set_vrttype(2); // RIDGEVERTEX
          }
        }
      }
    } // not hulltri and not exterior tri
  }

  return 1;
}