#ifndef __DETRI2_H__  // Include this file only once!
#define __DETRI2_H__

// In case NULL is not defined (e.g., on Windows).
#ifndef NULL
  #define NULL 0
#endif

namespace detri2 {

// Mesh Data Strutcure 

#define REAL double

class Triangle;
class TriEdge;

class Vertex
{
 public:
  REAL  crd[3]; // x, y, w (weight)
  int   idx;    // Its index (0 or 1-based).
  int   tag;    // Boundary marker.
  int   flags;  // Vertex type, infect, etc.
  Triangle* adj;// Adjacent triangle.
  Triangle* on; // Segment containing this vertex.
  REAL  *data;  // An array of vertex data (length is vertex_data_size).

  double highdi[3];//ADD for additional dimension crd.

  void init() {
    crd[0] = crd[1] = crd[2] = 0.0;
    idx = 0; tag = 0; flags = 0;
    adj = on = (Triangle *) 0;
    data = NULL;
  }

  Vertex() {init();}

  bool is_infected();
  void set_infect();
  void clear_infect();
  int  get_vrttype();
  void set_vrttype(int); 
  TriEdge get_triorg();
};

class Triangle
{
 public:
  Triangle* nei[3];
  Vertex*   vrt[3];
  int       flags;  // Encode the orientation of each neighbor triangle.
                    // Each neigbor uses 2 bits.
  int       tag;    // Boundary marker.
  int       idx; 
  REAL      *data;  // An array of triangle data (length is triangle_data_size).

  void init() {
    nei[0] = nei[1] = nei[2] = (Triangle *) 0;
    vrt[0] = vrt[1] = vrt[2] = (Vertex *) 0;
    flags = 0; tag = 0; idx = 0;
    data = NULL;
  }

  Triangle() { init();}
  
  int  getnver(int ver);
  void setnver(int ver, int val);
  bool is_hulltri();
  void set_hullflag();
  void clear_hullflag();
  void set_infect();
  bool is_infected();
  void clear_infect();
  void set_exterior();
  bool is_exterior();
  void clear_exterior();
};

class TriEdge
{
 public:
  Triangle* tri;
  int       ver; // = 0,1,2

  TriEdge() {tri = (Triangle *) 0; ver = 0;}
  TriEdge(Triangle* _t, int _v) {tri = _t; ver = _v;}

  // Primitives
  Vertex* org();
  Vertex* dest();
  Vertex* apex();
  void set_vertices(Vertex *pa, Vertex *pb, Vertex *pc);

  TriEdge enext();
  TriEdge enext2();
  TriEdge sym();
  void enextself();
  void enext2self();
  void symself();
  void bond1(const TriEdge &Enei); // Only one direction
  void bond2(const TriEdge &Enei); // Both direction
  void dissolve();

  // Triangle(t)->segment(s)
  bool is_segment();
  void set_segment();
  void clear_segment();

  void set_edge_infect();
  bool is_edge_infected();
  void clear_edge_infect();
};

// An encroached segment or a bad-qaulity triangle.
class EncEle
{
 public:
  TriEdge Ele;
  Vertex *end[3];
  Vertex *encpt; // the encroaching vertex (for splitting a segment)
  EncEle *next;  // for making a link list, stack, or queue.
  double key;    // for makeing a priority queue.

  void init() {
    Ele.tri = NULL, Ele.ver = 0;
    end[0] = end[1] = end[2] = (Vertex *) NULL;
    encpt = (Vertex *) NULL;
    next = (EncEle *) NULL;
    key = 0.0;
  }
  bool is_changed();
};

// A dynamic array
class arraypool {
 public:
  int objectbytes;
  int objectsperblock;
  int log2objectsperblock;
  int objectsperblockmark;
  int toparraylen;
  char **toparray;
  long objects;
  unsigned long totalmemory;
  void *deaditemstack;

  void restart();
  void poolinit(int sizeofobject, int log2objperblk);
  char* getblock(int objectindex);
  void* alloc();
  void dealloc(void *dyingitem);
  void* lookup(int objectindex);
  void* operator[](int index); // fast lookup

  arraypool();
  arraypool(int sizeofobject, int log2objperblk);
  ~arraypool();
};

// A dynamic array of vertices
class Vertices : public arraypool
{
 public:
  Vertices(int val) {poolinit(sizeof(Vertex), val);}
  Vertex *alloc_vrt();
  void dealloc_vrt(Vertex* vrt);
  bool is_deleted(Vertex* vrt);
  Vertex *get(int index); // fast lookup
};

// A dynamic array of triangles / segments
class Triangles : public arraypool
{
 public: 
  Triangles(int val) {poolinit(sizeof(Triangle), val);}
  Triangle* alloc_tri();
  void dealloc_tri(Triangle* tri);
  bool is_deleted(Triangle* tri);
  Triangle* get(int index); // fast lookup
};

// A Triangulation class.
class Triangulation
{
 public:
  // Arrays of vertices, segments, and triangles.
  Vertices*  vrts;
  Triangles* segs;
  Triangles* tris;
  Vertices*  regions; // A list of regions (holes).

  Vertex* infvrt; // The infinite vertex.

  arraypool *vertex_data; // The array storing all vertex data. 
  int vertex_data_size;   // The size (bytes) of a single vertex data. 

  arraypool *triangle_data; // The array storing all triangle data.
  int triangle_data_size;   // The size (bytes) of a single triangle data.

  arraypool *flip_queue; // The flip queue (a linear array).
  int lawsonflag, lawson_hullflag;

  Vertex *dual_verts; // Dual vertices (bary-, circum-, or orthocenters).
  arraypool *tet_verts;  // Tetrahedra created by flips.

  // Index-to-Vertex map.
  Vertex **idx_to_vert_list;
  int idx_to_vert_list_len; // Current length (used to validate the map).

  // Vertex-to-segment map.
  int *idx2seglist;
  Triangle **segperverlist;

  arraypool *chktris, *chktris2; // The triangles lists (for quality checking).
  arraypool *encsegs, *encsegs2; // The encroached segments lists (to be split).
  int encsegflag, enqtriflag;

  REAL xmin, xmax, ymin, ymax; // The bounding box.

  // Counters
  int input_vertices;
  int input_triangles;
  int hullsize;
  int unused_vertices;
  int flip13count, flip31count, flip22count, totalflipcount;
  int st_segref_count, st_triref_count;

  // Options
  int quiet; // -Q
  int verbose; // -V
  int norandomize; // -R
  int convex; // -c
  int delaunay; // -u +1 or -1;
  int weighted; // -w
  int weight_is_height; // -wh
  int transfer_weight_to_height; // -wt for visualization
  int gabriel; // -g (conforming Delaunay with Gabriel segments)
  int quality; // -q
  double minangle, maxratio2;  // Default 20 degree
  double maxarea; // -a#
  double minedgelen; // -l# The smallest edge length to be accepted.
  int hdflag;//ADDing for high dimensional embedding flag
  //double revtol; // -T# A relative tolerance. (Default 1e-6).

  // Output options
  int no_reindex; // -I
  int no_node_tag; // -B (for Tetview visualization)
  int output_edges; // -e

  double PI; // = 3.14159265358979323846264338327950288419716939937510582;

  // Input / Output
  char infilename[1024];
  char outfilename[1024];
  int firstindex;

  // Construction / Destruction
  Triangulation();
  ~Triangulation();

  // Ptimitives (which are Triangulation-dependent).
  bool is_hulledge(TriEdge &E);

  // Flips
  void flip13(Vertex *pt, TriEdge *tt);
  void flip31(TriEdge *tt);
  void flip22(TriEdge *tt);
  int flip(TriEdge *tt, int &flipflag, int check_only);
  int lawson_flip(Vertex *pt);
  int insert_point(Vertex *pt, TriEdge *tt, int loc);
  int remove_point(Vertex *pt, TriEdge *tt);

  // Delaunay triangulation.
  int first_tri(Vertex **ptlist, int ptnum, TriEdge &E);
  int locate_point(Vertex *pt, TriEdge &E);
  int incremental_Delaunay();

  // Constrained Delaunay triangulation.
  int find_direction(TriEdge& E, Vertex *pt);
  int recover_edge_flips(Vertex *e1, Vertex *e2, TriEdge &E);
  int recover_segments();

  // Create the index-to-vertex map.
  int build_index_to_vertex_map();
  // Mark/remove exterior (and hole) triangles.
  int remove_exterior_elements();
  // Mesh reconstruction.
  int reconstructmesh();

  // Improve mesh quality by adding points.
  int make_vertex_to_segment_map();
  Triangle *get_segment(TriEdge &E);
  EncEle* enq_triangle(arraypool*, TriEdge&);
  void enq_vertex_star(Vertex *pt);
  int check_segment(Vertex *pt, TriEdge &S);
  int get_edge(Vertex *e1, Vertex *e2, TriEdge &E);
  Vertex* get_split_segment_vertex(EncEle *S);
  int repair_encsegments();

  // Get the circumcenter (or orthocenter) of 3 (weighted) vertices, U, V, W.
  int get_orthocenter(REAL, REAL, REAL,     // Ux, Uy, U_weight
                      REAL, REAL, REAL,     // Vx, Vy, V_weight
                      REAL, REAL, REAL,     // Wx, Wy, W_weight
                      REAL*, REAL*, REAL*); // Cx, Cy, radius^2
  int get_tri_orthocenter(Vertex *v1, Vertex *v2, Vertex *v3, REAL*, REAL*, REAL*);
  int get_seg_orthocenter(Vertex *e1, Vertex *e2, REAL*, REAL*, REAL*);
  // Return the area of the Voronoi cell dual to pt.
  REAL get_tri_area(Vertex* pa, Vertex* pb, Vertex* pc);
  REAL get_voronoi_cell_area(Vertex* pt);

  int check_triangle(Triangle *tri, double ccent[3]);
  int repair_triangles();
  int delaunay_refinement();

  // Maximal Poisson Disk Sampling [Ebeida, Mitchell, et al 2011]

  // Mesh smoothing, CVT, ODT, MMS, ...

  // Mesh adaptation, Metrics, HDE, ...

  // Mesh statistics.
  void meshstatic();

  // Input / Output
  int parse_commands(int argc, char* argv[]);
  int read_nodes();
  int read_poly();
  int read_mesh();
  void save_triangulation(int reindex, int notag, int meshidx);
  void save_polyhedron(int meshidx);
  void save_tetrahedra(int meshidx);
  void save_to_ucd(int meshidx, int reindex = 0);
  void save_dual_to_ucd(int meshidx);
  void save_voronoi_cell_area(int meshidx);

// \iffalse
  // Debug functions
  int check_mesh();
  int check_delaunay(int queflag);
  void pvert(Vertex *v);
  void pverti(int idx);
  void ptriv(TriEdge &E);
  void ptriv(Triangle *t);
  void ptri(TriEdge &E);
  void ptri(Triangle *t);
  TriEdge get_trii(int v1, int v2, int v3);
  int  ptri_orthocenter(int v1, int v2, int v3);
  int  pseg_orthocenter(int v1, int v2);
  void print_tri_array(arraypool *triarray);
  void print_triedge_array(arraypool *tearray);
  void print_edge_array(arraypool *edgearray, const char *fn);
  void print_edge_list(TriEdge *elist, int n, const char *fn);

  int checkmesh; // -GC
  int checkdelaunay; // -GD
  int test_hilbert_sort; // -Gh
  int sweep_sort; // -Gs option
  double sweep_sort_vx, sweep_sort_vy;
  int do_lawson_flip; // -Gl option (only testing)
  int do_lawson_flip_origin; // -GF (the original Lawson's flip algo).
  int do_lawson_flip_test; // -Gt (testing different ideas).
  int do_generate_vertex_weights; // -GW
  int dump_intermediate_meshes; // -Gd
  int dump_to_mesh; // -Gm (medit, gmsh) default
  int dump_to_ucd; // -Gu (vtk, paraview)
  int dump_dual_to_ucd; // -GU (vtk, paraview)
  int dump_lifting_map; // -GL
  double lift_map_offset_z; // -GL# Vertically (z) shifted distance (for visualization)
  int dump_dual_cell_areas; // -GA // cell areas of the Power diagram.
  int dump_lifting_poly; // -GP
  int dump_lifting_tets; // -GT
  int dump_cycles; // -GS
  int exterior_triangles;
  int debug_conforming_delaunay;

  // Vertex permutation and sorting.
  void randomly_permute(Vertex **permutarray, int arraysize);
  void hilbert_init (int n);
  int  hilbert_split(Vertex **sortarray, int arraysize, int gc0, int gc1,
                     double, double, double, double);
  void hilbert_sort (Vertex **sortarray, int arraysize, int e, int d,
                     double bxmin, double bxmax, double bymin, double bymax,
                     int depth, int order, arraypool *hvertices);
  void brio_multiscale_sort(Vertex **sortarray, int arraysize,
                            int threshold, double ratio, int *depth, int order);
  void sweepline_sort(Vertex **sortarray, int arraysize, double vx, double vy);

  int mark_exterior_elements();

  // The plane is polar to the pole (orthocenter).
  int get_lifted_plane(REAL, REAL, REAL,     // Ux, Uy, U_weight
                       REAL, REAL, REAL,     // Vx, Vy, V_weight
                       REAL, REAL, REAL,     // Wx, Wy, W_weight
                       REAL*, REAL*, REAL*, REAL*); // Coeffs. a, b, c, d
  int get_tri_lifted_plane(Triangle *tri, REAL*, REAL*, REAL*, REAL*);
  int get_circumcenter(double*, double*, double*, double*, double*, double*);

  // Calculate dual diagram (lengths, cell areas, etc).
  //void get_dual_diagram();
  // Save the orthocircles and lifted planes.
  void save_pole_and_polar(int meshidx); 

  int recover_edge_conforming(Vertex *e1, Vertex *e2, TriEdge &E); // Test
  int recover_segments_conforming();

  // [2016-10-12, 00:46] Data structure for saving cycles.
  //   filled by get_cycles().
  class CycleEdges {
    public:
      int cyclen;
      TriEdge *cycedges;
      void init() {cyclen = 0; cycedges = NULL;}
  };
  arraypool *cyclist; 
  void cyclist_free();

  void get_cycles();

  // Constructing the furthest Delaunay triangulation using Lawson flips
  int lawson_flip_furthest_delaunay();

  // [2016-09-21, 21:49]
  int lawson_flip_furthest_delaunay_local(); // -Gl

  // [2016-11-03, 23:14, Moscow]
  int lawson_flip_furthest_delaunay_new(); // -Gt

  // [2016-11-13, 21:28, Berlin]
  // Removing a vertex by dynamically changing its weight.
  // - Only work for weighted DT.
  REAL get_flip_time(Vertex *pa, Vertex *pb, Vertex *pc, Vertex *pd);
  int remove_point_regular(Vertex *pt, TriEdge *elist, int &n);

  // [2016-11-14, 19:28, Berlin]
  int monotone_flip_furthest_delaunay_test();

  // [2016-11-15, 22:15, Berlin] // -GW
  // Based on an input triangulation, randomly generates weights at each vertex.
  //   each weight is calculated by: weight = r * L_ave * scale_weight,
  //   where, L_ave is the average edge length at this vertex, scale is a value
  //   between 0 and 1 (to scale L_ave), r is a random value between 0 and 1.
  REAL get_distance(Vertex *e1, Vertex *e2);
  int generate_vertex_weights_random();

  //ADD for high dimensional embedding
  void high_dimen_embedding();
  void lift_vrt();                                         //lift all the vertices to the higher space
  void edge_flip();                                        //flip edges in the high space
    bool check_spaceangle(EncEle *paryEle);                  //check the angle in high dimensional space if the edge should be flip
    void flip_space22(EncEle *paryEle);                      //flip the edge in high dimensional space
    double cal_angle(TriEdge &tt1,TriEdge &tt2);             //calaulate the degree space angle
  void edge_contraction();                                 //remove too short edges
    int chk_contra(EncEle *paryEle);                         //check the edge if it should be contracted
	void contract_edge(EncEle *paryEle,int contracase);      //contracte the edge paryEle,the contracase determine which point should be removed.0 for no contra ,
	                                                         //1 for dest ,2 for org,and 3 for free.4 & 5 for segment case.4 for org,5 for dest and 6 for free.
  void edge_splitting();                                   //cut the long edges
    bool chk_split(EncEle *paryEle);                         //check if the edge should be split
    int split_edge(TriEdge &tt);                            //split the edge and change the arraypool of triangals
  void node_smoothing(double alph);                        //move vertices to better place......to be debugged
  double interpolation_error_sum();                        //calculate the sum of the piecewise linear interpolation error
    double interpolation_error(Triangle *tri);               //calculate the linear interpolation error on one grid
	int formula_solver(double *solve, Vertex* *vrt);         //solver the 3-D formula for the interpolation
	double inline detmination(double mat[3][3]);             //fast cal the det of a 3-D mat
	double extra_func(Vertex* temp[3],double coef[3],double x,double y);//calculate the value of the point pt in rectangle region
	double* area_to_axis(double area[2],Triangle tr);        //trans the coordinate from area to (x,y)
  double embed_fun(double x,double y);                          //embedding function
  double para[4];                                          //paraments for func
  double s;                                                //parameter for embedding space
  double L;                                                //length of standar edge
  double hi_distance(TriEdge E);                           //calculate the hi_length of the TriEdge E
  double hi_distance(Vertex* pt1, Vertex* pt2);
  double fmax, fmin, gxmax, gxmin, gymax, gymin;           //for normalize the function value and the gradiance value

  //ADD for locally remesh
  arraypool *moving_vrts=NULL;                                  //queue for pointer to all the moving vrts
  int count = 0;                                                //forgot this para has what mean
  void remesh_Triangulation(double *x, double *y, int *c, int k);
  void load_move_vertices(double *x, double *y, int *c, int k);   //initial,read all the moving vrts to the queue
  void findandremove_segpoint();                                //find the moving segments and remove the points on it
  void get_flip_queue();                                        //collect all the edge that should be push to the flip queue
  void delete_collision_points();                               //remove all the points accross to the moving segment
  void change_vrts();
  void locally_improvement();                                   //Locally Delaunay refinement
  void angle_statistic();                                       //print the statistic of all the angles
  Vertex* jilu[2][100];
  int ji = 0;
  //ADD fi
// fi
};

class MovVertex{
public:
	Vertex *mov_vrt;                                       //point to the vertex that should be moved
	double vec[2];                                         //moving vector
	int cla = 0;                                           //movvertex type.cla==0 common moving vrt,cla==1 moving holes
	MovVertex()
	{
		mov_vrt = NULL;
		vec[0] = vec[1] = 0.0;
	};
	MovVertex(Vertex* pt, double a, double b)
	{
		mov_vrt = pt;
		vec[0] = a;
		vec[1] = b;
	}
	~MovVertex()
	{
		delete mov_vrt;
		delete (void*)vec;
	};
};


} // namespace detri2 {

#endif  // #ifndef __DETRI2_H__
