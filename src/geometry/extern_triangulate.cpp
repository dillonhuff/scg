#include "geometry/extern_triangulate.h"

#include "geometry/vtk_utils.h"
#include "geometry/vtk_debug.h"

#define ANSI_DECLARATORS

#ifdef SINGLE

#define REAL float

#else

#define REAL double

#endif

#include <string>
#include <stdio.h>
#include <stdlib.h>

//extern "C" {

//}

#include "triangle_lib/triangle.h"

namespace gca {

  void report(struct triangulateio *io,
	      int markers,
	      int reporttriangles,
	      int reportneighbors,
	      int reportsegments,
	      int reportedges,
	      int reportnorms) {

    int i, j;

    for (i = 0; i < io->numberofpoints; i++) {
      printf("Point %4d:", i);
      for (j = 0; j < 2; j++) {
	printf("  %.6g", io->pointlist[i * 2 + j]);
      }
      if (io->numberofpointattributes > 0) {
	printf("   attributes");
      }
      for (j = 0; j < io->numberofpointattributes; j++) {
	printf("  %.6g",
	       io->pointattributelist[i * io->numberofpointattributes + j]);
      }
      if (markers) {
	printf("   marker %d\n", io->pointmarkerlist[i]);
      } else {
	printf("\n");
      }
    }
    printf("\n");

    if (reporttriangles || reportneighbors) {
      for (i = 0; i < io->numberoftriangles; i++) {
	if (reporttriangles) {
	  printf("Triangle %4d points:", i);
	  for (j = 0; j < io->numberofcorners; j++) {
	    printf("  %4d", io->trianglelist[i * io->numberofcorners + j]);
	  }
	  if (io->numberoftriangleattributes > 0) {
	    printf("   attributes");
	  }
	  for (j = 0; j < io->numberoftriangleattributes; j++) {
	    printf("  %.6g", io->triangleattributelist[i *
						       io->numberoftriangleattributes + j]);
	  }
	  printf("\n");
	}
	if (reportneighbors) {
	  printf("Triangle %4d neighbors:", i);
	  for (j = 0; j < 3; j++) {
	    printf("  %4d", io->neighborlist[i * 3 + j]);
	  }
	  printf("\n");
	}
      }
      printf("\n");
    }

    if (reportsegments) {
      for (i = 0; i < io->numberofsegments; i++) {
	printf("Segment %4d points:", i);
	for (j = 0; j < 2; j++) {
	  printf("  %4d", io->segmentlist[i * 2 + j]);
	}
	if (markers) {
	  printf("   marker %d\n", io->segmentmarkerlist[i]);
	} else {
	  printf("\n");
	}
      }
      printf("\n");
    }

    if (reportedges) {
      for (i = 0; i < io->numberofedges; i++) {
	printf("Edge %4d points:", i);
	for (j = 0; j < 2; j++) {
	  printf("  %4d", io->edgelist[i * 2 + j]);
	}
	if (reportnorms && (io->edgelist[i * 2 + 1] == -1)) {
	  for (j = 0; j < 2; j++) {
	    printf("  %.6g", io->normlist[i * 2 + j]);
	  }
	}
	if (markers) {
	  printf("   marker %d\n", io->edgemarkerlist[i]);
	} else {
	  printf("\n");
	}
      }
      printf("\n");
    }
  }

  std::vector<triangle> trilib_to_tris(const polygon_3& original,
				       struct triangulateio mid) {
    DBG_ASSERT(original.vertices().size() > 2);

    double z = original.vertices().front().z;

    vector<point> point_vec;
    for (unsigned i = 0; i < mid.numberofpoints; i++) {
      REAL px = mid.pointlist[2*i];
      REAL py = mid.pointlist[2*i + 1];
      point_vec.push_back(point(px, py, z));
    }

    // cout << "POINTS" << endl;
    // for (auto p : point_vec) {
    //   cout << p << endl;
    // }
  
    vector<triangle> tris;
    for (unsigned i = 0; i < mid.numberoftriangles; i++) {
      vector<point> pts;
      for (unsigned j = 0; j < mid.numberofcorners; j++) {
	auto ind = mid.trianglelist[i * mid.numberofcorners + j];
	pts.push_back(point_vec[ind]);
      }

      DBG_ASSERT(pts.size() == 3);

      point v0 = pts[0];
      point v1 = pts[1];
      point v2 = pts[2];
    
      point q1 = v1 - v0;
      point q2 = v2 - v0;
      point norm = cross(q2, q1).normalize();
      tris.push_back(triangle(norm, v0, v1, v2));
    }

    // cout << "Number of triangles = " << tris.size() << endl;
    // for (auto t : tris) {
    //   cout << t << endl;
    // }

    return tris;
  }

  point interior_point(const std::vector<point>& pts) {
    polygon_3 hole_poly = build_clean_polygon_3(pts);
    auto bp = to_boost_poly_2(hole_poly);

    for (unsigned i = 0; i < pts.size(); i++) {

      for (unsigned j = 2; j < pts.size() - 1; j++) {
	unsigned i2 = (i + j) % pts.size();

	point pi0 = pts[i];

	point pi2 = pts[i2];

	point mid = (1.0 / 2.0) * (pi2 + pi0);
      
	auto pt_2 = bg::model::d2::point_xy<double>(mid.x, mid.y);
	if (bg::within(pt_2, bp)) { return mid; }
      }

    }

    vtk_debug_ring(pts);

    DBG_ASSERT(false);
  }

  vector<std::pair<int, int> > ring_segments(const std::vector<point>& pts) {
    vector<std::pair<int, int> > segs;
    for (unsigned i = 0; i < pts.size(); i++) {
      segs.push_back(make_pair(i, (i + 1) % pts.size() ));
    }
    return segs;
  }
  
  void set_points_and_segments(struct triangulateio* in, const polygon_3& p) {
    in->numberofpoints = p.vertices().size();

    for (auto h : p.holes()) {
      in->numberofpoints += h.size();
    }

    in->pointlist = (REAL *) malloc(in->numberofpoints * 2 * sizeof(REAL));
    unsigned vert_ind = 0;

    for (unsigned i = 0; i < p.vertices().size(); i++) {
      point pi = p.vertices()[i];
      in->pointlist[2*vert_ind] = pi.x;
      in->pointlist[2*vert_ind + 1] = pi.y;

      vert_ind++;
    }

    for (auto h : p.holes()) {
      for (unsigned i = 0; i < h.size(); i++) {
	point pi = h[i];
	in->pointlist[2*vert_ind] = pi.x;
	in->pointlist[2*vert_ind + 1] = pi.y;

	vert_ind++;
      }
    }
    
    in->numberofpointattributes = 0;
    in->pointattributelist = NULL;

    in->pointmarkerlist = (int *) malloc(in->numberofpoints * sizeof(int));
    for (unsigned i = 0; i < in->numberofpoints; i++) {
      in->pointmarkerlist[i] = 0;
    }

    in->numberofsegments = p.vertices().size();
    for (auto h : p.holes()) {
      in->numberofsegments += h.size();
    }

    vector<unsigned> buffer_offsets;
    vector<vector<std::pair<int, int> > > segments;

    buffer_offsets.push_back(0);
    segments.push_back(ring_segments(p.vertices()));

    auto last_size = p.vertices().size();
    for (auto h : p.holes()) {
      buffer_offsets.push_back(buffer_offsets.back() + last_size);
      segments.push_back(ring_segments(h));
      last_size = h.size();
    }

    in->segmentlist =
      static_cast<int*>(malloc(in->numberofsegments*sizeof(int)*2));

    cout << "# of points = " << in->numberofpoints << endl;
    unsigned seg_num = 0;
    for (unsigned i = 0; i < segments.size(); i++) {
      auto seg_group = segments[i];
      auto offset = buffer_offsets[i];
      for (unsigned j = 0; j < seg_group.size(); j++) {
	in->segmentlist[2*seg_num] = offset + seg_group[j].first;
	in->segmentlist[2*seg_num + 1] = offset + seg_group[j].second;

	DBG_ASSERT(in->segmentlist[2*seg_num] >= 0);
	DBG_ASSERT(in->segmentlist[2*seg_num + 1] >= 0);

	// cout << "segment " << seg_num << " start : " << in->segmentlist[2*seg_num] << endl;
	// cout << "segment " << seg_num << " end: " << in->segmentlist[2*seg_num + 1] << endl;
	
	DBG_ASSERT(in->segmentlist[2*seg_num] < in->numberofpoints);
	DBG_ASSERT(in->segmentlist[2*seg_num + 1] < in->numberofpoints);
	
	seg_num++;
      }
    }

    DBG_ASSERT(seg_num == in->numberofsegments);

    in->numberofholes = p.holes().size();

    in->holelist =
      static_cast<REAL*>(malloc(in->numberofholes*sizeof(REAL)*2));

    // TODO: Actually use interior point computation
    for (unsigned i = 0; i < in->numberofholes; i++) {
      point pt = interior_point(p.holes()[i]);

      polygon_3 hole_poly = build_clean_polygon_3(p.holes()[i]);
      auto bp = to_boost_poly_2(hole_poly);
      auto pt_2 = bg::model::d2::point_xy<double>(pt.x, pt.y);

      DBG_ASSERT(bg::within(pt_2, bp));

      in->holelist[2*i] = pt.x;
      in->holelist[2*i + 1] = pt.y;
    }

    in->numberofregions = 0;
    in->regionlist = 0;

  }

  void set_output(struct triangulateio* mid) {
    /* Make necessary initializations so that Triangle can return a */
    /*   triangulation in `mid' */

    mid->pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
    /* Not needed if -N switch used or number of point attributes is zero: */
    mid->pointattributelist = (REAL *) NULL;
    mid->pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
    mid->trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
    /* Not needed if -E switch used or number of triangle attributes is zero: */
    mid->triangleattributelist = (REAL *) NULL;
    mid->neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
    /* Needed only if segments are output (-p or -c) and -P not used: */
    mid->segmentlist = (int *) NULL;
    /* Needed only if segments are output (-p or -c) and -P and -B not used: */
    mid->segmentmarkerlist = (int *) NULL;
    mid->edgelist = (int *) NULL;             /* Needed only if -e switch used. */
    mid->edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

    /* Triangulate the points.  Switches are chosen to read and write a  */
  }

  std::vector<triangle> triangulate_flat_3d(const polygon_3& p) {
    // DANGER: This was inserted to eliminate a link error in the GUI
    // version of our UIST paper. Restore once that is done!
    //DBG_ASSERT(false);

    struct triangulateio in, mid;

    /* Define input points. */

    set_points_and_segments(&in, p);

    printf("Input point set:\n\n");
    //report(&in, 1, 0, 0, 0, 0, 0);

    set_output(&mid);


    /*   PSLG (p), preserve the convex hull (c), number everything from  */
    /*   zero (z), assign a regional attribute to each element (A), and  */
    /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
    /*   neighbor list (n).                                              */

    char settings_2[] = "pczAen";
    triangulate(&(settings_2[0]), &in, &mid, NULL);

    printf("Initial triangulation:\n\n");
    //report(&mid, 1, 1, 1, 1, 1, 0);

    DBG_ASSERT(mid.numberoftriangles > 0);

    /* Free all allocated arrays, including those allocated by Triangle. */

    vector<triangle> tris = trilib_to_tris(p, mid);

    free(in.pointlist);
    free(in.pointattributelist);
    free(in.pointmarkerlist);
    free(in.regionlist);

    free(mid.pointlist);
    free(mid.pointattributelist);
    free(mid.pointmarkerlist);
    free(mid.trianglelist);
    free(mid.triangleattributelist);

    free(mid.neighborlist);
    free(mid.segmentlist);
    free(mid.segmentmarkerlist);
    free(mid.edgelist);
    free(mid.edgemarkerlist);

    return tris;

  }

}
