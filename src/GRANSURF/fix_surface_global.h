/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(surface/global,FixSurfaceGlobal)
// clang-format on
#else

#ifndef LMP_FIX_SURFACE_GLOBAL_H
#define LMP_FIX_SURFACE_GLOBAL_H

#include <stdio.h>
#include "fix.h"
#include "fix.h"

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranularModel;
}

class FixSurfaceGlobal : public Fix {
 public:

  // neighbor lists for spheres with surfs and shear history
  // accessed by fix shear/history

  class NeighList *list;
  class NeighList *listhistory;

  FixSurfaceGlobal(class LAMMPS *, int, char **);
  ~FixSurfaceGlobal();
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void setup_pre_neighbor() override;
  void initial_integrate(int) override;
  void pre_neighbor() override;
  void post_force(int) override;

  int modify_param(int, char **x) override;
  void reset_dt() override;
  double memory_usage() override;

  void *extract(const char *, int &) override;
  int image(int *&, double **&) override;

 private:
  int dimension,firsttime,use_history;
  double dt,skin;

  // for granular model choices

  class Granular_NS::GranularModel *model;
  int history, size_history, heat_flag;

  double Twall;
  int tvar;
  char *tstr;

  // per-surf properties

  double **xsurf,**vsurf,**omegasurf,*radsurf;

  double triggersq;

  // group settings

  int ngroup;             // # of defined groups
  char **gnames;          // name of each group
  int *bitmask;           // one-bit mask for each group
  int *gmask;             // mask bits for each surf

  // motion settings

  struct Motion {
    int igroup;
    int mstyle;
    int vxflag,vyflag,vzflag;
    int axflag,ayflag,azflag;
    double vx,vy,vz;
    double ax,ay,az;
    double period;
    double point[3],axis[3],unit[3];
    double omega;
    double time_origin;
  };

  struct Motion *motions;
  int nmotion;
  int maxmotion;

  double **points_original,**xsurf_original;
  double **points_lastneigh;

  // storage of granular history info

  class FixNeighHistory *fix_history;
  double *zeroes;

  // rigid body masses for use in granular interactions

  class Fix *fix_rigid;    // ptr to rigid body fix, NULL if none
  double *mass_rigid;      // rigid mass for owned+ghost atoms
  int nmax;                // allocated size of mass_rigid

  // data structs for extracting surfs from molecule files

  struct Point {
    double x[3];
  };

  struct Line {
    int mol,type;           // molID and type of the line
    int p1,p2;              // indices of points in line segment
    double norm[3];         // unit normal to line = Z x (p2-p1)
  };

  struct Tri {
    int mol,type;           // modID and type of the triangle
    int p1,p2,p3;           // indices of points in triangle
    double norm[3];         // unit normal to tri plane = (p2-p1) x (p3-p1)
  };

  Point *points;              // global list of unique points
  Line *lines;                // global list of lines
  Tri *tris;                  // global list of tris
  int npoints,nlines,ntris;   // count of each
  int nsurf;                  // count of lines or tris for 2d/3d

                              // ragged 2d arrays for 2d connectivity
  int **plines;               // indices of lines which contain each point
  int **neigh_p1;             // indices of other lines connected to endpt1
  int **pwhich_p1;            // which point (0/1) on other line is endpt1
  int **nside_p1;             // consistency of other line normal
                              //   SAME_SIDE or OPPOSITE_SIDE
  int **aflag_p1;             // is this line + other line a FLAT,CONCAVE,CONVEX surf
                              //   surf = on normal side of this line
  int **neigh_p2;             // ditto for connections to endpt2
  int **pwhich_p2;            // ditto for endpt2
  int **nside_p2;             // ditto for endpt2
  int **aflag_p2;             // ditto for endpt2

  int **elist;                // ragged 2d array for global tri edge lists
  int **clist;                // ragged 2d array for global tri corner pt lists

  // 2d/3d connectivity

  struct Connect2d {      // line connectivity
    int np1,np2;          // # of lines connected to endpt1/2 (NOT including self)
    int *neigh_p1;        // indices of lines connected to emdpt1
    int *neigh_p2;        // ditto for connections to endpt2
    int *pwhich_p1;       // which point (0,1) on other line is endpt1
    int *pwhich_p2;       // ditto for endpt2
    int *nside_p1;        // consitency of other line normal
    int *nside_p2;        // ditto for endpt2
                          //   SAME_SIDE = 2 normals are on same side of surf
                          //   OPPOSITE_SIDE = opposite sides of surf
    int *aflag_p1;        // is this line + other line a FLAT,CONCAVE,CONVEX surf
    int *aflag_p2;        // ditto for endpt2
                          //   surf = on normal side of this line
                          //   aflag = FLAT, CONCAVE, CONVEX
  };

  struct Connect3d {      // tri connectivity
    int ne1,ne2,ne3;      // # of tris connected to edges 1,2,3 (including self)
    int nc1,nc2,nc3;      // # of tris connected to corner pts 1,2,3 (including self)
    int *neigh_e1;        // indices of all tris connected to 1-2 edge (if ne1 > 1)
    int *neigh_e2;        // ditto for 2-3 edge
    int *neigh_e3;        // ditto for 3-1 edge
    int *neigh_c1;        // indices of all tris connected to corner pt 1 (if nc1 > 1)
    int *neigh_c2;        // ditto for corner pt 2
    int *neigh_c3;        // ditto for corner pt 3
    int flags;            // future flags for edge and corner pt coupling
  };

  Connect2d *connect2d;             // 2d connection info
  Connect3d *connect3d;             // 3d connection info

  // struct for storing contact data

  struct ContactingSurf {
    int index, type;
    double r[3];
    double overlap;
  };

  ContactingSurf *contacting_surfs;
  int nmaxcontacts;

  // data for DumpImage

  int *imflag;
  double **imdata;
  int imax;

  // private methods

  void extract_from_molecules(char *);
  void extract_from_stlfile(char *);
  void connectivity2d_global();
  void connectivity3d_global();
  void surface_attributes();
  void move_init();
  void move_clear();

  int modify_params_group(int, int, int, char **);
  int modify_params_move(Motion *, int, char **);
  int find_group(const char *);
  int add_group(const char *);
};

}    // namespace LAMMPS_NS

#endif
#endif
