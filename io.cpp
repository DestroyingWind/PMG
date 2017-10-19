#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "detri2.h"
#include "pred2d.h"

#pragma warning (disable:4996)
using  namespace detri2;

int Triangulation::parse_commands(int argc, char* argv[])
{
  infilename[0] = '\0';
  outfilename[0] = '\0';
  firstindex = 0;
  
  char workstring[256];
  int j, k;

  for (int i = 1; i < argc; i++) {
    // Is this string a filename?
    if (argv[i][0] != '-') {
      strncpy(infilename, argv[i], 1024 - 1);
      infilename[1024 - 1] = '\0';
      continue;                     
    } else {
      // Parse the individual switch from the string.
      for (j = 1; argv[i][j] != '\0'; j++) {
        if (argv[i][j] == 'Q') {
          quiet = 1;
        } else if (argv[i][j] == 'V') {
          verbose++;
        } else if (argv[i][j] == 'R') {
          norandomize = 1;
        } else if (argv[i][j] == 'c') {
          convex = 1;
        } else if (argv[i][j] == 'u') {
          delaunay = -1;
        } else if (argv[i][j] == 'w') {
          weighted = 1;
          if (argv[i][j+1] == 'h') { // -wh
            weight_is_height = 1; j++;
          } else if (argv[i][j+1] == 't') { // -wt
            transfer_weight_to_height = 1; j++;
          }
        } else if (argv[i][j] == 'g') {
          gabriel = 1;
        } else if (argv[i][j] == 'q') {
          quality = 1;
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            minangle = (REAL) strtod(workstring, (char **) NULL);
          }
        } else if (argv[i][j] == 'a') {
          if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
              (argv[i][j + 1] == '.')) {
            k = 0;
            while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                   (argv[i][j + 1] == '.')) {
              j++;
              workstring[k] = argv[i][j];
              k++;
            }
            workstring[k] = '\0';
            maxarea = (REAL) strtod(workstring, (char **) NULL);
          }
        } else if (argv[i][j] == 'I') {
          no_reindex = 1;
        } else if (argv[i][j] == 'B') {
          no_node_tag = 1;
        } else if (argv[i][j] == 'e') {
          output_edges = 1;
		} else if (argv[i][j] == 'H') {//ADDing for high dimensional embedding flag
			hdflag = 1;
		}
// \iffalse  // Debug options
        else if (argv[i][j] == 'G') { // -G
          if (argv[i][j+1] == 'C') { // -GC
            checkmesh = 1; j++;
          } else if (argv[i][j+1] == 'D') { // -GD
            checkdelaunay = 1; j++;
          } else if (argv[i][j+1] == 'h') { // -Gh
            test_hilbert_sort = 1; j++;
          } else if (argv[i][j+1] == 's') { // -Gs
            sweep_sort = 1; j++;
            if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
              k = 0;
              while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.')) {
                j++;
                workstring[k] = argv[i][j];
                k++;
              }
              workstring[k] = '\0';
              sweep_sort_vx = (REAL) strtod(workstring, (char **) NULL);
            }
            if ((argv[i][j + 1] == '/') || (argv[i][j + 1] == ',')) {
              j++;
              if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                  (argv[i][j + 1] == '.')) {
                k = 0;
                while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.') || (argv[i][j + 1] == 'e') ||
                     (argv[i][j + 1] == '-') || (argv[i][j + 1] == '+')) {
                  j++;
                  workstring[k] = argv[i][j];
                  k++;
                }
                workstring[k] = '\0';
                sweep_sort_vy = (REAL) strtod(workstring, (char **) NULL);
              }
            }
          } else if (argv[i][j+1] == 't') { // -Gt
            do_lawson_flip_test = 1; j++;
          } else if (argv[i][j+1] == 'F') { // -GF
            do_lawson_flip_origin = 1; j++;
          } else if (argv[i][j+1] == 'l') { // -Gl
            do_lawson_flip = 1; j++;
          } else if (argv[i][j+1] == 'W') { // -GW
            do_generate_vertex_weights = 1; j++;
          } else if (argv[i][j+1] == 'd') { // -Gd
            dump_intermediate_meshes = 1; j++;
          } else if (argv[i][j+1] == 'm') { // -Gm
            dump_to_mesh = 1; j++; j++;
          } else if (argv[i][j+1] == 'u') { // -Gu
            dump_to_ucd = 1; j++;
          } else if (argv[i][j+1] == 'U') { // -GU
            dump_dual_to_ucd = 1; j++;  
          } else if (argv[i][j+1] == 'L') { // -GL
            dump_lifting_map = 1; j++;
            if (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                (argv[i][j + 1] == '.')) {
              k = 0;
              while (((argv[i][j + 1] >= '0') && (argv[i][j + 1] <= '9')) ||
                     (argv[i][j + 1] == '.')) {
                j++;
                workstring[k] = argv[i][j];
                k++;
              }
              workstring[k] = '\0';
              lift_map_offset_z = (REAL) strtod(workstring, (char **) NULL);
            }
          } else if (argv[i][j+1] == 'A') { // -GA
            dump_dual_cell_areas = 1; j++;
          } else if (argv[i][j+1] == 'S') { // -GS
            dump_cycles = 1; j++;
            //debug_conforming_delaunay = 1; j++;
          } else if (argv[i][j+1] == 'P') { // -GP
            dump_lifting_poly = 1; j++;
          } else if (argv[i][j+1] == 'T') { // -GT
            dump_lifting_tets = 1; j++;
            //assert(tet_verts == NULL);
            tet_verts = new arraypool(sizeof(Vertex*)*4, 8);
          }
        } // -G
// \fi
      } // for j
    }
  }

  if (infilename[0] == '\0') {
    // No input file name. Print the syntax and exit.
    return 0;
  }

  // Reconginze any file format.
  if (!strcmp(&infilename[strlen(infilename) - 5], ".node")) {
    infilename[strlen(infilename) - 5] = '\0';
  } else if (!strcmp(&infilename[strlen(infilename) - 5], ".poly")) {
    infilename[strlen(infilename) - 5] = '\0';
  } else if (!strcmp(&infilename[strlen(infilename) - 4], ".ele")) {
    infilename[strlen(infilename) - 4] = '\0';
  } else if (!strcmp(&infilename[strlen(infilename) - 5], ".edge")) {
    infilename[strlen(infilename) - 5] = '\0';
  } else if (!strcmp(&infilename[strlen(infilename) - 7], ".region")) {
    infilename[strlen(infilename) - 7] = '\0';
  }

  int increment = 0;
  strcpy(workstring, infilename);
  j = 1;
  while (workstring[j] != '\0') {
    if ((workstring[j] == '.') && (workstring[j + 1] != '\0')) {
      increment = j + 1;
    }
    j++;
  }
  int meshnumber = 0;
  if (increment > 0) {
    j = increment;
    do {
      if ((workstring[j] >= '0') && (workstring[j] <= '9')) {
        meshnumber = meshnumber * 10 + (int) (workstring[j] - '0');
      } else {
        increment = 0;
      }
      j++;
    } while (workstring[j] != '\0');
  }
  if (increment == 0) {
    meshnumber = 0;
  } else {
    workstring[increment-1] = '\0';
  }
  sprintf(outfilename, "%s.%d", workstring, meshnumber + 1);

  return 1;
}

// Read the file "infilename.node"
int Triangulation::read_nodes()
{
  FILE *infile;
  char line[1024], *pstr;
  char delim[] = " ,";
  int pnum = 0, i;

  char filename[256];
  strcpy(filename, infilename);
  strcat(filename, ".node");
  infile = fopen(filename, "r");
  if (infile == NULL) {
    printf("Unable to open file %s.\n", filename);
    return 0;
  }

  // Read the number of nodes.
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }

  pnum = atoi(pstr);
  if (pnum <= 0) {
    printf("No points in file %s.\n", filename);
    fclose(infile);
    return 0;
  }

  printf("Reading %d points from file %s\n", pnum, filename);
  REAL x, y;

  for (i = 0; i < pnum; i++) {
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '#')) break;
    }
    if (feof(infile)) break;
    Vertex *vrt = vrts->alloc_vrt();
    unused_vertices++;
    vrt->idx = atoi(pstr);
    if (i == 0) {
      // Set the first index (0 or 1).
      firstindex = vrt->idx;
    }
    pstr = strtok(NULL, delim);
    x = vrt->crd[0] = atof(pstr);
    pstr = strtok(NULL, delim);
    y = vrt->crd[1] = atof(pstr);
    pstr = strtok(NULL, delim);
    if (pstr != NULL) {
      vrt->tag = atoi(pstr);
    }
    // Determine the smallest and largest x, and y coordinates.
    if (i == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
    } else {
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
    }
  }

  if (i < pnum) {
    printf("Missing %d points from file %s.\n", pnum - i, filename);
    fclose(infile);
    return 0;
  }
  input_vertices = pnum;

  fclose(infile);
  return 1;
}

int Triangulation::read_poly()
{
  FILE *infile;
  char line[1024], *pstr;
  char delim[] = " ,";
  int i;

  char filename[256];
  strcpy(filename, infilename);
  strcat(filename, ".poly");
  infile = fopen(filename, "r");
  if (infile == NULL) {
    printf("Unable to open file %s.poly\n", infilename);
    return 0;
  }

  // Read the number of nodes.
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }
  if (feof(infile)) {
    printf("Failed to read the point list.\n");
    fclose(infile);
    return 0;
  }

  int pnum = atoi(pstr);
  if (pnum == 0) {
    if (!read_nodes()) {
      printf("Failed to read the point list.\n");
      fclose(infile);
      return 0;
    }
  } else {
    printf("Reading %d points from file %s\n", pnum, filename);
    REAL x, y;
    for (i = 0; i < pnum; i++) {
      while (fgets(line, 1024, infile)) {
        //printf("%s", line);
        pstr = strtok(line, delim);
        if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
            (pstr[0] != '#')) break;
      }
      if (feof(infile)) break;
      Vertex *vrt = vrts->alloc_vrt();
      unused_vertices++;
      vrt->idx = atoi(pstr);
      if (i == 0) {
        // Set the first index (0 or 1).
        firstindex = vrt->idx;
      }
      pstr = strtok(NULL, delim);
      x = vrt->crd[0] = atof(pstr);
      pstr = strtok(NULL, delim);
      y = vrt->crd[1] = atof(pstr);
      pstr = strtok(NULL, delim);
      if (pstr != NULL) {
        vrt->tag = atoi(pstr);
      }
      // Determine the smallest and largest x, and y coordinates.
      if (i == 0) {
        xmin = xmax = x;
        ymin = ymax = y;
      } else {
        xmin = (x < xmin) ? x : xmin;
        xmax = (x > xmax) ? x : xmax;
        ymin = (y < ymin) ? y : ymin;
        ymax = (y > ymax) ? y : ymax;
      }
    }
    if (i < pnum) {
      printf("Missing %d points from file %s.\n", pnum - i, filename);
      fclose(infile);
      return 0;
    }
    input_vertices = pnum;
  }

  // Read the number of segments.
  int snum = 0, e1, e2, tag;
  
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }
  if (pstr != NULL) {
    snum = atoi(pstr);
  }
  if (snum <= 0) {
    printf("No segments in file %s.\n", filename);
    fclose(infile);
    return 0;
  }

  // Make a map from idx to vertices.
  arraypool *idx2vertlist = new arraypool(sizeof(Vertex *), 8);
  for (i = 0; i < vrts->objects; i++) {
    Vertex* paryvrt = (Vertex *) (*vrts)[i];
    // Find the location of this vertex in the map.
    Vertex** destvrt = (Vertex **) idx2vertlist->lookup(paryvrt->idx);
    if (destvrt == NULL) {
      // The memory for this location does not exist.
      idx2vertlist->getblock(paryvrt->idx);
      destvrt = (Vertex **) (*idx2vertlist)[paryvrt->idx];
    }
    *destvrt = paryvrt;
  }

  printf("Reading %d segments from file %s\n", snum, filename);

  for (i = 0; i < snum; i++) {
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (feof(infile)) break;
    //seg->idx = atoi(pstr);
    pstr = strtok(NULL, delim);
    e1 = atof(pstr);
    pstr = strtok(NULL, delim);
    e2 = atof(pstr);
    pstr = strtok(NULL, delim);
    if (pstr != NULL) {
      tag = atof(pstr);
    } else {
      tag = 0;
    }
    if (e1 != e2) {
      Triangle *seg = segs->alloc_tri();
      seg->init();
      seg->vrt[0] = * ((Vertex **) (*idx2vertlist)[e1]);
      seg->vrt[1] = * ((Vertex **) (*idx2vertlist)[e2]);
      seg->tag = tag;
    }
  }
  if (i < snum) {
    printf("Missing %d segments from file %s.\n", snum - i, filename);
    fclose(infile);
    return 0;
  }

  delete idx2vertlist;

  // Read the number of regions (holes).
  int rnum = 0;

  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }

  if (pstr != NULL) {
    rnum = atoi(pstr);
  }
  if (rnum <= 0) {
    fclose(infile);
    return 1;
  }

  printf("Reading %d regions (and holes) from file %s\n", rnum, filename);

  for (i = 0; i < rnum; i++) {
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (feof(infile)) break;
    Vertex *vrt = regions->alloc_vrt();
    vrt->idx = atoi(pstr);
    pstr = strtok(NULL, delim);
    vrt->crd[0] = atof(pstr);
    pstr = strtok(NULL, delim);
    vrt->crd[1] = atof(pstr);
    pstr = strtok(NULL, delim);
    if (pstr != NULL) {
      vrt->tag = atoi(pstr);
      pstr = strtok(NULL, delim);
      if (pstr != NULL) {
        vrt->crd[2] = atof(pstr); // region maxarea.
      }
    } else {
      vrt->tag = -9999; // Hole
    }
  }

  fclose(infile);
  return 1;
}

// Read the file(s): *.node, *.edge, *.ele, and *.region.
int Triangulation::read_mesh()
{
  FILE *infile = NULL;
  char line[1024], *pstr;
  char delim[] = " ,";
  int i;

  // For supporting the reader of Triangle's poly format.
  if (read_poly()) {
    return 1;
  }

  char filename[256];
  strcpy(filename, infilename);
  strcat(filename, ".node");
  infile = fopen(filename, "r");
  if (infile == NULL) {
    printf("Unable to open file %s.\n", filename);
    return 0;
  }

  // Read the number of nodes.
  while (fgets(line, 1024, infile)) {
    //printf("%s", line);
    pstr = strtok(line, delim);
    if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
        (pstr[0] != '#')) break;
  }

  int pnum = atoi(pstr);
  if (pnum <= 0) {
    printf("No points in file %s.\n", filename);
    fclose(infile);
    return 0;
  }

  printf("Reading %d points from file %s\n", pnum, filename);
  REAL x, y;

  for (i = 0; i < pnum; i++) {
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '#')) break;
    }
    //if (feof(infile)) break;
    Vertex *vrt = vrts->alloc_vrt();
    unused_vertices++;
    vrt->idx = atoi(pstr);
    if (i == 0) {
      // Set the first index (0 or 1).
      firstindex = vrt->idx;
      if ((firstindex != 0) && (firstindex != 1)) {
        printf("Wrong first index %d from file %s.\n", firstindex, filename);
        break;
      }
    } else {
      if (vrt->idx != (i + firstindex)) {
        printf("Wrong index %d at %d-th vertex from file %s.\n", vrt->idx,
               i + firstindex, filename);
        break;
      }
    }
    pstr = strtok(NULL, delim);
    x = vrt->crd[0] = atof(pstr);
    pstr = strtok(NULL, delim);
    y = vrt->crd[1] = atof(pstr);
    pstr = strtok(NULL, delim);
    if (pstr != NULL) {
      vrt->tag = atoi(pstr);
    }
    // Determine the smallest and largest x, and y coordinates.
    if (i == 0) {
      xmin = xmax = x;
      ymin = ymax = y;
    } else {
      xmin = (x < xmin) ? x : xmin;
      xmax = (x > xmax) ? x : xmax;
      ymin = (y < ymin) ? y : ymin;
      ymax = (y > ymax) ? y : ymax;
    }
  }

  if (i < pnum) {
    printf("Failed to read vertices from file %s.\n", filename);
    fclose(infile);
    return 0;
  }
  input_vertices = pnum;

  // Make a map from idx to vertices.
  Vertex **idx2vertlist = new Vertex*[input_vertices + firstindex];
  for (i = 0; i < vrts->objects; i++) {
    assert(vrts->get(i)->idx == (i + firstindex));
    idx2vertlist[i + firstindex] = vrts->get(i);
  }

  // Try to read a .weight file (if it exists).
  strcpy(filename, infilename);
  strcat(filename, ".weight");
  infile = fopen(filename, "r");
  if (infile != NULL) {
    // Read the number of segments.
    int wnum = 0;
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      wnum = atoi(pstr);
    }
    if (wnum == pnum) {
      // The number of points and weights must be equal.
      printf("Reading point weights from file %s\n", filename);
      for (i = 0; i < wnum; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        //if (feof(infile)) break;
        //seg->idx = atoi(pstr); // i
        pstr = strtok(NULL, delim);
        vrts->get(i)->crd[2] = atof(pstr);
      }
      if (i < wnum) {
        printf("Missing %d point weights from file %s.\n", wnum - i, filename);
      }
    } else {
      printf("Wrong number %d (should be %d) of point weights\n", wnum, pnum);
    }
    fclose(infile);
  } // Read point weights

  // Try to read a .edge file (if it exists).
  strcpy(filename, infilename);
  strcat(filename, ".edge");
  infile = fopen(filename, "r");
  if (infile != NULL) {
    // Read the number of segments.
    int snum = 0, e1, e2, tag;
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      snum = atoi(pstr);
    }
    if (snum > 0) {
      printf("Reading segments from %d edges from file %s\n", snum, filename);
      for (i = 0; i < snum; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        //if (feof(infile)) break;
        //seg->idx = atoi(pstr);
        e1 = e2 = 0;
        pstr = strtok(NULL, delim);
        e1 = atof(pstr);
        pstr = strtok(NULL, delim);
        e2 = atof(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          tag = atof(pstr);
        } else {
          tag = 0;
        }
        if (tag != 0) {
          // A segment must have non-zero tag.
          if ((e1 != e2) &&
              (e1 >= firstindex) && (e1 < input_vertices + firstindex) &&
              (e2 >= firstindex) && (e2 < input_vertices + firstindex)) {
            Triangle *seg = segs->alloc_tri();
            seg->vrt[0] = idx2vertlist[e1];
            seg->vrt[1] = idx2vertlist[e2];
            seg->tag = tag;
          } else {
            printf("Segment %d has invalid vertices.\n", i + firstindex);
          }
        }
      }
      if (i < snum) {
        printf("Missing %d edges from file %s.\n", snum - i, filename);
      } else {
        printf("Read %d segments from file %s\n", snum, filename);
      }
    } // snum > 0
    fclose(infile);
  } // Read segments.

  // Try to read a .ele file (if it exists).
  strcpy(filename, infilename);
  strcat(filename, ".ele");
  infile = fopen(filename, "r");
  if (infile != NULL) {
    // Read the number of segments.
    int trinum = 0, v1, v2, v3, tag, ori;
    Vertex *p1, *p2, *p3;
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      trinum = atoi(pstr);
    }
    if (trinum > 0) {
      printf("Reading %d triangles from file %s\n", trinum, filename);
      for (i = 0; i < trinum; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        //if (feof(infile)) break;
        //seg->idx = atoi(pstr);
        v1 = v2 = v3 = 0;
        pstr = strtok(NULL, delim);
        v1 = atof(pstr);
        pstr = strtok(NULL, delim);
        v2 = atof(pstr);
        pstr = strtok(NULL, delim);
        v3 = atof(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          tag = atof(pstr);
        } else {
          tag = 0;
        }
        if ((v1 != v2) && (v2 != v3) && (v3 != v1) &&
            (v1 >= firstindex) && (v1 < input_vertices + firstindex) &&
            (v2 >= firstindex) && (v2 < input_vertices + firstindex) &&
            (v3 >= firstindex) && (v3 < input_vertices + firstindex)) {
          p1 = idx2vertlist[v1];
          p2 = idx2vertlist[v2];
          p3 = idx2vertlist[v3];
          // Make sure all triangles are all CCW oriented.
          ori = Orient2d(p1->crd, p2->crd, p3->crd);
          if (ori < 0) {
            // Swap the points.
            Vertex *swap = p1;
            p1 = p2;
            p2 = swap;
          }
          if (ori != 0) {
            Triangle *tri = tris->alloc_tri();
            tri->vrt[0] = p1;
            tri->vrt[1] = p2;
            tri->vrt[2] = p3;
            tri->tag = tag;
          } else {
            printf("!! Triangle #%d [%d,%d,%d] is degenerated.\n", i+firstindex,
                   p1->idx, p2->idx, p3->idx);
          }
        } else {
          printf("!! Triangle #%d [%d,%d,%d] has invalid vertices.\n",
                 i + firstindex, v1, v2, v3);
        }
      }
      if (i < trinum) {
        printf("!! Missing %d triangles from file %s.\n", trinum - i, filename);
      }
      input_triangles = trinum;
    } // snum > 0
    fclose(infile);
  } // Read triangles.

  delete [] idx2vertlist;

  // Try to read a .region file (if it exists).
  strcpy(filename, infilename);
  strcat(filename, ".region");
  infile = fopen(filename, "r");
  if (infile != NULL) {
    // Read the number of regions (holes).
    int rnum = 0;
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      rnum = atoi(pstr);
    }
    if (rnum > 0) {
      printf("Reading %d regions (and holes) from file %s\n", rnum, filename);
      for (i = 0; i < rnum; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        if (feof(infile)) break;
        Vertex *vrt = regions->alloc_vrt();
        vrt->idx = atoi(pstr);
        pstr = strtok(NULL, delim);
        vrt->crd[0] = atof(pstr);
        pstr = strtok(NULL, delim);
        vrt->crd[1] = atof(pstr);
        pstr = strtok(NULL, delim);
        if (pstr != NULL) {
          vrt->tag = atoi(pstr);
          pstr = strtok(NULL, delim);
          if (pstr != NULL) {
            vrt->crd[2] = atof(pstr); // region maxarea.
          }
        } else {
          vrt->tag = -9999; // Hole
        }
      } // i
    } // rnum > 0
    fclose(infile);
  } // Read .region file.

  return 1;
}

// If the "reindex" is (1), then vertices with type UNUSEDVERTEX will be
//   removed from the output, and the vertices are re-indexed.
// If the "meshidx" is (> 0), then the output file name will contain this idex.
//   It is combined with the save_polyhedron() function.
// The "notag" option are only used for Tetview visualization.
void Triangulation::save_triangulation(int reindex, int notag, int meshidx)
{
  char filename[256];
  if (meshidx > 0) {
    sprintf(filename, "%s_%d.node", outfilename, meshidx);
  } else {
    strcpy(filename, outfilename);
    strcat(filename, ".node");
  }
  FILE *outfile = fopen(filename, "w");
  int nv = reindex ? (vrts->objects - unused_vertices) : vrts->objects;
  printf("Writing %d nodes to file %s.\n", nv, filename);

  fprintf(outfile, "%d 2 0 0\n", nv);

  int i, count = 0;
  int idx = firstindex; // Reindex the vertices.
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = (Vertex *)(* vrts)[i]; // = vrts->get(i)
    if (!vrts->is_deleted(vrt)) {
      if (reindex && (vrt->get_vrttype() == 0)) { // UNUSEDVERTEX
        if (verbose > 1) {
          printf("  Skip vertex %d\n", vrt->idx);
        }
        count++; continue; // Skip this vertex.
      }
      fprintf(outfile, "%d %g %g", idx, vrt->crd[0], vrt->crd[1]);
      if (!notag) {
        if (vrt->on) {
          fprintf(outfile, " %d", (vrt->on->tag != 0) ? vrt->on->tag : -1);
        } else {
          // An input vertex or facet vertex.
          fprintf(outfile, " %d", vrt->tag);
        }
      }
      fprintf(outfile, "\n");
      vrt->idx = idx;
      idx++; count++;
    }
  }
  fclose(outfile);

  if (weighted) {
    if (meshidx > 0) {
      sprintf(filename, "%s_%d.weight", outfilename, meshidx);
    } else {
      strcpy(filename, outfilename);
      strcat(filename, ".weight");
    }
    outfile = fopen(filename, "w");
    //int nv = reindex ? (vrts->objects - unused_vertices) : vrts->objects;
    printf("Writing %d weights to file %s.\n", nv, filename);
    fprintf(outfile, "%d\n", nv);
    count = 0;
    for (i = 0; count < vrts->objects; i++) {
      Vertex *vrt = vrts->get(i);
      if (!vrts->is_deleted(vrt)) {
        if (reindex && (vrt->get_vrttype() == 0)) { // UNUSEDVERTEX
          count++; continue; // Skip this vertex.
        }
        fprintf(outfile, " %d", vrt->idx);
        if (transfer_weight_to_height) { // -wt only for visulaiztion
          if (weight_is_height) {
            fprintf(outfile, " %g", vrt->crd[2]);
          } else {
            fprintf(outfile, " %g", 
              vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1] - vrt->crd[2]);
          }
        } else {
          fprintf(outfile, " %g", vrt->crd[2]);
        }
        fprintf(outfile, "\n");
        count++;
      }
    }
    fclose(outfile);
  }

  if ((segs->objects > 0l) || output_edges) { // -e
    // Output segments / edges.
    if (meshidx > 0) {
      sprintf(filename, "%s_%d.edge", outfilename, meshidx);
    } else {
      strcpy(filename, outfilename);
      strcat(filename, ".edge");
    }
    outfile = fopen(filename, "w");
    if (output_edges) { // -e option
      if (segs->objects > 0l) {
        // Since vertex may be reindexed, re-create this map.
        make_vertex_to_segment_map();
      }
      int nedg = (3 * ((int) tris->objects - hullsize) + hullsize) / 2;
      printf("Writing %d edges to file %s.\n", nedg, filename);
      fprintf(outfile, "%d 1\n", nedg);
      idx = firstindex;
      count = 0;
      TriEdge E;
      for (int i = 0; count < tris->objects; i++) {
        E.tri = tris->get(i);
        if (!tris->is_deleted(E.tri)) {
          if (!E.tri->is_hulltri()) { // Ignore a hull triangle.
            for (E.ver = 0; E.ver < 3; E.ver++) {
              if (!(E.sym().tri->is_infected())) {
                fprintf(outfile, "%d  %d %d", idx, E.org()->idx, E.dest()->idx);
                if (E.is_segment()) {
                  // Cannot use this function since the vertex may be reindexed.
                  Triangle *seg = get_segment(E);
                  // Make sure the tag is not non-zero (to distinguish it with
                  //   an interior edge).
                  int tag = (seg->tag != 0) ? seg->tag : -1;
                  fprintf(outfile, "  %d", tag);
                } else {
                  fprintf(outfile, "  0");
                }
                fprintf(outfile, "\n");
                idx++;
              }
            }
            E.tri->set_infect(); // Mark it.
          }
          count++;
        }
      }
      // Unmark all triangles.
      count = 0;
      for (int i = 0; count < tris->objects; i++) {
        E.tri = tris->get(i);
        if (!tris->is_deleted(E.tri)) {
          if (!E.tri->is_hulltri()) { // Ignore a hull triangle.
            E.tri->clear_infect();
          }
          count++;
        }
      }
      assert(idx == nedg + firstindex);
      if (segs->objects > 0l) {
        delete [] idx2seglist;
        delete [] segperverlist;
        idx2seglist = NULL;
        segperverlist = NULL;
      }
    } else {
      // Default output all segments
      printf("Writing %ld segments to file %s.\n", segs->objects, filename);
      fprintf(outfile, "%ld 1\n", segs->objects);
      idx = firstindex;
      count = 0;
      for (int i = 0; count < segs->objects; i++) {
        Triangle* seg = segs->get(i);
        if (!segs->is_deleted(seg)) {
          int tag = (seg->tag != 0) ? seg->tag : -1;
          fprintf(outfile, "%d  %d %d  %d\n", idx, seg->vrt[0]->idx,
                  seg->vrt[1]->idx, tag);
          idx++; count++;
        }
      }
    }
    fclose(outfile);
  } // Output segments / edges.

  int ntri = (int) tris->objects - hullsize;
  if (meshidx > 0) {
    sprintf(filename, "%s_%d.ele", outfilename, meshidx);
  } else {
    strcpy(filename, outfilename);
    strcat(filename, ".ele");
  }
  outfile = fopen(filename, "w");
  printf("Writing %d triangles to file %s.\n", ntri, filename);

  fprintf(outfile, "%d 3 0\n", ntri);
  idx = firstindex;
  count = 0;
  for (int i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      // ignore a hull triangle.
      if (!tri->is_hulltri()) {
        fprintf(outfile, "%d  %d %d %d  %d\n", idx, tri->vrt[0]->idx,
                tri->vrt[1]->idx, tri->vrt[2]->idx, tri->tag);
        idx++;
      }
      count++;
    }
  }
  fclose(outfile);
}

// Save the current triangulation (lower faces) together with the furthest
//   Delaunay triangulation (upper faces), must be pre-computed and saved
//   in the file "filename.u.ele".
// [Important] The faces saved in filename.u.ele must have the same indices
//   as the current nodes. Must use the "-I" option all the time. 
void Triangulation::save_polyhedron(int meshidx)
{
  char filename[256];
  FILE *infile = NULL;
  char line[1024], *pstr;
  char delim[] = " ,";
  int i;

  // Read the upperfaces.ele
  int trinum = 0, triidx = 0;
  int *trilist = NULL; // new int[trinum * 3];
  int *ptri = NULL;

  sprintf(filename, "%s.u.ele", infilename);
  infile = fopen(filename, "r");

  if (infile != NULL) {
    // Read the number of triangles.
    while (fgets(line, 1024, infile)) {
      //printf("%s", line);
      pstr = strtok(line, delim);
      if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
          (pstr[0] != '#')) break;
    }
    if (pstr != NULL) {
      trinum = atoi(pstr);
    }
    if (trinum > 0) {
      // Read the triangles.
      trilist = new int[trinum * 3];

      // Some triangles may be coincident with current ones,
      // They need to be identified and removed from the output.
      TriEdge tt[2];
      int loc, fexist;

      // Make a map from idx to vertices.
      Vertex **idx2vertlist = new Vertex*[input_vertices + firstindex];
      for (i = 0; i < vrts->objects; i++) {
        assert(vrts->get(i)->idx == (i + firstindex));
        idx2vertlist[i + firstindex] = vrts->get(i);
      }

      for (i = 0; i < trinum; i++) {
        while (fgets(line, 1024, infile)) {
          //printf("%s", line);
          pstr = strtok(line, delim);
          if ((pstr != NULL) && (pstr[0] != '\r') && (pstr[0] != '\n') &&
              (pstr[0] != '#')) break;
        }
        //if (feof(infile)) break;
        ptri = &(trilist[triidx * 3]);
        pstr = strtok(NULL, delim);
        ptri[0] = atof(pstr);
        pstr = strtok(NULL, delim);
        ptri[1] = atof(pstr);
        pstr = strtok(NULL, delim);
        ptri[2] = atof(pstr);

        // Check if this triangle exists in the current triangulation.
        fexist = 0;
        tt[0].tri = idx2vertlist[ptri[0]]->adj;
        for (tt[0].ver = 0; tt[0].ver < 3; tt[0].ver++) {
          if (tt[0].org() == idx2vertlist[ptri[0]]) break;
        }
        assert(tt[0].ver < 3);
        loc = find_direction(tt[0], idx2vertlist[ptri[1]]);
        if (loc == 2) {
          // Found the edge ptri[0]->ptri[1].
          if (tt[0].apex() == idx2vertlist[ptri[2]]) {
            // Found the triangle. Exists.
            tt[0].tri->set_infect();
            fexist = 1;
          } else {
            tt[1] = tt[0].sym();
            if (tt[1].apex() == idx2vertlist[ptri[2]]) {
              // Found the triangle. Exists.
              tt[1].tri->set_infect();
              fexist = 1;
            }
          }
        }

        if (!fexist) {
          triidx++; // Add a new triangle.
        }
      }

      delete [] idx2vertlist;
    } else { // if (trinum <= 0) {
      printf("Wrong number of faces %d in %s\n", trinum, filename);
    }
    fclose(infile);
  } else { // if (infile != NULL)
    printf("Unable to open file %s\n", filename);
  }

  if (meshidx > 0) {
    sprintf(filename, "%s_%d.smesh", outfilename, meshidx);
  } else {
    sprintf(filename, "%s.smesh", outfilename);
  }
  FILE *outfile = fopen(filename, "w");
  int ntri = (int) tris->objects - hullsize;

  int totaltri = ntri - (trinum - triidx) + triidx;

  printf("Writing %d triangles to file %s.\n", totaltri, filename);

  int nv = vrts->objects;
  fprintf(outfile, "%d 3 0 0\n", nv);

  int idx=1; // UCD index starts from 1.
  int count = 0;
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = (Vertex *)(* vrts)[i];
    if (!vrts->is_deleted(vrt)) {
      double h = 0; // height
      if (weighted) { // -w
        if (weight_is_height) { // -wh
          h = vrt->crd[2];
        } else {
          h = vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1] - vrt->crd[2];
        }
      } else {
        h = vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1];
      }
      fprintf(outfile, "%d %g %g %g\n", idx, vrt->crd[0], vrt->crd[1], 
              h + lift_map_offset_z);
      vrt->idx = idx; // Re-index.
      idx++; count++;
    }
  }

  // Output list of triangles.
  fprintf(outfile, "%d 1\n", totaltri);

  // UCD assumes vertex index starts from 1.
  int shift = (firstindex == 1 ? 0 : 1);
  idx = 1;
  count = 0;
  for (i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      // ignore a hull triangle.
      if (!tri->is_hulltri()) {
        if (!tri->is_infected()) {
          fprintf(outfile, "3 %d %d %d  1\n",
                  tri->vrt[0]->idx + shift,
                  tri->vrt[1]->idx + shift,
                  tri->vrt[2]->idx + shift);
          idx++;
        } else {
          // Skip a duplicated triangle.
          tri->clear_infect();
        }
      }
      count++;
    }
  }

  //for (i = 0; i < trinum; i++) {
  for (i = 0; i < triidx; i++) {
    ptri = &(trilist[i * 3]);
    fprintf(outfile, "3 %d %d %d  2\n", ptri[0], ptri[1], ptri[2]);
    idx++;
  }

  assert((idx - 1) == totaltri);

  fprintf(outfile, "0\n");
  fprintf(outfile, "0\n");

  fclose(outfile);

  if (trilist) {
    delete [] trilist;
  }
}

void Triangulation::save_tetrahedra(int meshidx)
{
  char filename[256];

  sprintf(filename, "%s_%d.node", outfilename, meshidx);
  FILE *outfile = fopen(filename, "w");
  int nv = vrts->objects;
  printf("Writing %d nodes to file %s.\n", nv, filename);

  fprintf(outfile, "%d 3 0 0\n", nv);

  int idx=1;
  int count = 0;
  for (int i = 0; count < vrts->objects; i++) {
    Vertex *vrt = (Vertex *)(* vrts)[i];
    if (!vrts->is_deleted(vrt)) {
      double h = 0; // height
      if (weighted) { // -w
        if (weight_is_height) { // -wh
          h = vrt->crd[2];
        } else {
          h = vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1] - vrt->crd[2];
        }
      } else {
        h = vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1];
      }
      fprintf(outfile, "%d %g %g %g\n", idx, vrt->crd[0], vrt->crd[1], 
              h + lift_map_offset_z);
      vrt->idx = idx; // Re-index.
      idx++; count++;
    }
  }

  fclose(outfile);

  sprintf(filename, "%s_%d.ele", outfilename, meshidx);
  outfile = fopen(filename, "w");
  int ntet = tet_verts->objects;
  printf("Writing %d tetrahedra to file %s.\n", ntet, filename);

  fprintf(outfile, "%d 4 0\n", ntet);

  for (int i = 0; i < tet_verts->objects; i++) {
    Vertex **pTet = (Vertex **)(* tet_verts)[i];
    fprintf(outfile, "%d  %d %d %d %d\n", i+1, pTet[0]->idx, pTet[1]->idx,
            pTet[2]->idx, pTet[3]->idx);
  }

  fclose(outfile);
}

void Triangulation::save_to_ucd(int meshidx, int reindex)
{
  char filename[256];
  sprintf(filename, "%s_%d.inp", outfilename, meshidx);
  FILE *outfile = fopen(filename, "w");
  int ntri = (int) tris->objects - hullsize;
  printf("Writing %d triangles to file %s.\n", ntri, filename);

  int nv = reindex ? (vrts->objects - unused_vertices) : vrts->objects;
  fprintf(outfile, "%d %d 0 0 0\n", nv, ntri);

  int i, idx=1; // UCD index starts from 1.
  int count = 0;
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = (Vertex *)(* vrts)[i];
    if (!vrts->is_deleted(vrt)) {
      if (reindex && (vrt->get_vrttype() == 0)) { // UNUSEDVERTEX
        count++; continue; // Skip this vertex.
      }
      fprintf(outfile, "%d %g %g 0\n", idx, vrt->crd[0], vrt->crd[1]);
      vrt->idx = idx; // Re-index.
      idx++; count++;
    }
  }

  // UCD assumes vertex index starts from 1.
  int shift = (firstindex == 1 ? 0 : 1);
  idx = 1;
  count = 0;
  for (i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      // ignore a hull triangle.
      if (!tri->is_hulltri()) {
        fprintf(outfile, "%d 0 tri %d %d %d\n", idx,
          tri->vrt[0]->idx + shift, 
          tri->vrt[1]->idx + shift, 
          tri->vrt[2]->idx + shift);
        idx++;
      }
      count++;
    }
  }
  fprintf(outfile, "\n");
  
  fclose(outfile);

  if (dump_lifting_map) { // -GL option
    sprintf(filename, "lift_%s_%d.inp", outfilename, meshidx);
    outfile = fopen(filename, "w");
    printf("Writing %d triangles to file %s.\n", ntri, filename);
    fprintf(outfile, "%d %d 0 0 0\n", nv, ntri);
    idx=1; count = 0;
    for (i = 0; count < vrts->objects; i++) {
      Vertex *vrt = (Vertex *)(* vrts)[i];
      if (!vrts->is_deleted(vrt)) {
        if (reindex && (vrt->get_vrttype() == 0)) { // UNUSEDVERTEX
          count++; continue; // Skip this vertex.
        }
        double h = 0; // height
        if (weighted) { // -w
          if (weight_is_height) { // -wh
            h = vrt->crd[2];
          } else {
            h = vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1] - vrt->crd[2];
          }
        } else {
          h = vrt->crd[0]*vrt->crd[0] + vrt->crd[1]*vrt->crd[1];
        }
        fprintf(outfile, "%d %g %g %g\n", idx, vrt->crd[0], vrt->crd[1], 
                h + lift_map_offset_z);
        vrt->idx = idx; // Re-index.
        idx++; count++;
      }
    }
    // UCD assumes vertex index starts from 1.
    shift = (firstindex == 1 ? 0 : 1);
    idx = 1;
    int count = 0;
    for (i = 0; count < tris->objects; i++) {
      Triangle* tri = tris->get(i);
      if (!tris->is_deleted(tri)) {
        // ignore a hull triangle.
        if (!tri->is_hulltri()) {
          fprintf(outfile, "%d 0 tri %d %d %d\n", idx,
            tri->vrt[0]->idx + shift,
            tri->vrt[1]->idx + shift,
            tri->vrt[2]->idx + shift);
          idx++;
        }
        count++;
      }
    }
    fprintf(outfile, "\n");
    fclose(outfile);
  }
}

void Triangulation::save_dual_to_ucd(int meshidx)
{
  // Computer the set of dual (Voronoi) vertices.
  // A (non hull) triangle is dual to a vertex (its circumcenter).
  // A hull triangle is dual to a vertex on its hull edge.
  Vertex *tmp_dual_verts = new Vertex[tris->objects];

  Vertex *e1, *e2;
  double Cx, Cy, radius2;
  int i, idx = 1; // UCD assumes vertex index starts from 1. 
  int count = 0;

  // Calculate dual vertices (orthocenters)
  for (i = 0; count < tris->objects; i++) {
    Triangle* tri = tris->get(i);
    if (!tris->is_deleted(tri)) {
      if (!tri->is_hulltri()) {
        get_tri_orthocenter(tri->vrt[0], tri->vrt[1], tri->vrt[2], &Cx, &Cy, &radius2);
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
      tmp_dual_verts[count].crd[0] = Cx;
      tmp_dual_verts[count].crd[1] = Cy;
      tmp_dual_verts[count].crd[2] = radius2; // w (weight, NOT height)
      tmp_dual_verts[count].idx = idx;
      tri->idx = count;
      count++; idx++;
    }
  }

  char filename[256];
  sprintf(filename, "%s_dual_%d.inp", outfilename, meshidx);
  FILE *outfile = fopen(filename, "w");
  int nv = (int) tris->objects;
  int ne = (int) tris->objects * 3 / 2; // 2e = 3f
  printf("Writing %d vertices and %d edges to file %s.\n", nv, ne, filename);

  fprintf(outfile, "%d %d 0 0 0\n", nv, ne);

  idx = 1; // // UCD index starts from 1.
  for (i = 0; i < nv; i++) {
    fprintf(outfile, "%d %g %g 0\n", idx, tmp_dual_verts[i].crd[0], tmp_dual_verts[i].crd[1]);
    idx++;
  }

  // Loop through all triangles, output edges.
  TriEdge tt;
  count = 0; idx = 1;
  for (i = 0; count < tris->objects; i++) {
    tt.tri = tris->get(i);
    if (!tris->is_deleted(tt.tri)) {
      for (tt.ver = 0; tt.ver < 3; tt.ver++) {
        if (tt.tri->idx < tt.sym().tri->idx) {
          fprintf(outfile, "%d 0 line %d %d\n", idx, 
                  tmp_dual_verts[tt.tri->idx].idx, 
                  tmp_dual_verts[tt.sym().tri->idx].idx);
          idx++;
        }
      }
      count++;
    }
  }
  fclose(outfile);

  if (dump_lifting_map) { // -GL option
    sprintf(filename, "lift_%s_dual_%d.inp", outfilename, meshidx);
    outfile = fopen(filename, "w");
    printf("Writing %d vertices and %d edges to file %s.\n", nv, ne, filename);

    fprintf(outfile, "%d %d 0 0 0\n", nv, ne);

    idx = 1; // // UCD index starts from 1.
    for (i = 0; i < nv; i++) {
      // Each lifted vertex is: (cx, cy, cx^2+cy^2-r^2 (=height)).
      REAL h = tmp_dual_verts[i].crd[0] * tmp_dual_verts[i].crd[0] 
             + tmp_dual_verts[i].crd[1] * tmp_dual_verts[i].crd[1]               
             - tmp_dual_verts[i].crd[2];
      fprintf(outfile, "%d %g %g %g\n", idx, 
              tmp_dual_verts[i].crd[0], tmp_dual_verts[i].crd[1], 
              h + lift_map_offset_z);
      idx++;
    }

    // Loop through all triangles, output edges.
    count = 0; idx = 1;
    for (i = 0; count < tris->objects; i++) {
      tt.tri = tris->get(i);
      if (!tris->is_deleted(tt.tri)) {
        for (tt.ver = 0; tt.ver < 3; tt.ver++) {
          if (tt.tri->idx < tt.sym().tri->idx) {
            fprintf(outfile, "%d 0 line %d %d\n", idx, 
                    tmp_dual_verts[tt.tri->idx].idx, 
                    tmp_dual_verts[tt.sym().tri->idx].idx);
            idx++;
          }
        }
        count++;
      }
    }
    fclose(outfile);
  }

  delete [] tmp_dual_verts;
}

// Save the cell area of the dual graph (Voronoi diagram or power diagram)
//   into a .area file.
void Triangulation::save_voronoi_cell_area(int meshidx)
{
  char filename[256];
  if (meshidx > 0) {
    sprintf(filename, "%s_%d.v.area", outfilename, meshidx);
  } else {
    strcpy(filename, outfilename);
    strcat(filename, ".v.area");
  }

  FILE *outfile = fopen(filename, "w");
  int nv = vrts->objects;
  printf("Writing %d dual cell areas to file %s.\n", nv, filename);
  fprintf(outfile, "%d\n", nv);

  REAL total_area = 0.0, area;
  int i, count = 0;
  int idx = firstindex; // Reindex the vertices.
  for (i = 0; count < vrts->objects; i++) {
    Vertex *vrt = (Vertex *)(* vrts)[i]; // = vrts->get(i)
    if (!vrts->is_deleted(vrt)) {
      if (vrt->get_vrttype() == 0) { // UNUSEDVERTEX
        area = 0.0;
      } else {
        area = get_voronoi_cell_area(vrt);
      }
      total_area += area;
      fprintf(outfile, "%d %g\n", idx, area);
      vrt->idx = idx;
      idx++; count++;
    }
  }

  printf("# Total area = %g.\n", total_area);
  fclose(outfile);
}

///////////////////////////////////////////////////////////////////////////////

void Triangulation::meshstatic()
{
  printf("\nStatistics:\n\n");
  printf("  Input vertices: %d\n", input_vertices);
  if (segs->objects > 0l) {
    printf("  Input segments: %ld\n", segs->objects);
  }
  if (input_triangles > 0) {
    printf("  Input triangles: %d\n", input_triangles);
  }
  if (regions->objects > 0l) {
    printf("  Input region (holes): %ld\n", regions->objects);
  }

  int nvrt = vrts->objects - unused_vertices;
  int ntri = (int) tris->objects - hullsize;
  int nedg = (3 * ((int) tris->objects - hullsize) + hullsize) / 2;
  printf("\n  Mesh vertices: %d\n", nvrt);
  printf("  Mesh edges: %d\n", nedg);
  printf("  Mesh triangles: %d\n", ntri);
  if (segs->objects > 0l) {
    printf("  Mesh exterior boundary edges: %d\n", hullsize);
  } else {
    printf("  Convex hull edges: %d\n", hullsize);
  }
  if (gabriel || quality) {
    printf("  Boundary Steiner points: %d\n", st_segref_count);
  }
  if (quality) {
    printf("  Interior Steiner points: %d\n", st_triref_count);
  }
  if (unused_vertices > 0) {
    if (no_reindex) {
      printf("  Number of unused vertices: %d\n", unused_vertices);
    }
  }
  printf("\n");
}