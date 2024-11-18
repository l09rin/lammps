// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// for both lines and tris:
// NOTE: allow for multiple instances of this fix or not ?
// NOTE: set flatthresh based on angle (degrees) set by user
// NOTE: need to order connections with FLAT first?
// NOTE: error checks for zero-length lines, duplicate lines
// NOTE: warn for too-small lines
// NOTE: more efficient neighbor lists
//       see Joel's 18 Nov email for idea
// NOTE: alter connection info if 2 lines are different types ?
// NOTE: allow setting a granular model for each type
// MOTE: how to set molecule IDs and types for STL and molecule input
//       STL files have no per-tri flags
//       molecule files can have types - how to apply offset?
//       could allow reading of multiple STL or molecule files in one fix
//         with an arg for mol ID or type
//         or could use different instances of FSG to do this
// NOTE: assign motion via fix modify move keyword with a type (range)
//       this means each type can have a Motion class assigned to it
//       integrate() should access the matching Motion class for each surf
// NOTE: should this fix produce any output
//         global array with force on each surf
//         or global array of forces per molecule ID (consecutive) ?
//         ditto for particle contact counts ?
// NOTE: should all fix move options be supported by this fix?
// NOTE: allow use of the scale keyword in the molecule command for lines/tris ?

#include "fix_surface_global.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "my_page.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "region.h"
#include "stl_reader.h"
#include "surf_extra.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <algorithm>
#include <map>
#include <tuple>
#include <unordered_set>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace Granular_NS;
using namespace MathConst;
using namespace MathExtra;
using namespace SurfExtra;

enum{TYPE,MOLECULE,ID};
enum{LT,LE,GT,GE,EQ,NEQ,BETWEEN};

enum{SPHERE,LINE,TRI};           // also in DumpImage
enum{LINEAR,WIGGLE,ROTATE,TRANSROT,VARIABLE};

enum{FLAT,CONCAVE,CONVEX};
enum{SAME_SIDE,OPPOSITE_SIDE};

#define FLATTHRESH 0.01
#define DELTA 128
#define DELTACONTACTS 4

/* ---------------------------------------------------------------------- */

FixSurfaceGlobal::FixSurfaceGlobal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), tstr(nullptr)
{
  if (narg < 11) error->all(FLERR,"Illegal fix surface/global command");

  if (!atom->omega_flag) error->all(FLERR,"Fix surface/global requires atom attribute omega");
  if (!atom->radius_flag) error->all(FLERR,"Fix surface/global requires atom attribute radius");

  // set interaction style
  // disable bonded/history option for now

  model = new GranularModel(lmp);
  model->contact_type = SURFACE;

  heat_flag = 0;
  int classic_flag = 1;
  if (strcmp(arg[4], "granular") == 0) classic_flag = 0;

  // wall/particle coefficients

  int iarg;
  if (classic_flag) {
    iarg = model->define_classic_model(arg, 4, narg);

    while (iarg < narg) {
      if (strcmp(arg[iarg],"limit_damping") == 0) {
        model->limit_damping = 1;
        iarg += 1;
      }
    }

  } else {
    iarg = 5;
    iarg = model->add_sub_model(arg, iarg, narg, NORMAL);

    while (iarg < narg) {
      if (strcmp(arg[iarg], "damping") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, DAMPING);
      } else if (strcmp(arg[iarg], "tangential") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, TANGENTIAL);
      } else if (strcmp(arg[iarg], "rolling") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, ROLLING);
      } else if (strcmp(arg[iarg], "twisting") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, TWISTING);
      } else if (strcmp(arg[iarg], "heat") == 0) {
        iarg = model->add_sub_model(arg, iarg + 1, narg, HEAT);
        heat_flag = 1;
      } else if (strcmp(arg[iarg],"limit_damping") == 0) {
        model->limit_damping = 1;
        iarg += 1;
      } else {
        break;
      }
    }
  }

  // define default damping sub model if unspecified, takes no args

  if (!model->damping_model) model->construct_sub_model("viscoelastic", DAMPING);
  model->init();

  size_history = model->size_history;
  if (model->beyond_contact) size_history += 1; //Need to track if particle is touching
  if (size_history == 0) use_history = 0;
  else use_history = 1;

  // optional args

  int scaleflag = 0;
  int Twall_defined = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix surface/global command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix surface/global command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"temperature") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix surface/global command");
      if (utils::strmatch(arg[iarg+1], "^v_")) {
        tstr = utils::strdup(arg[iarg+1] + 2);
      } else {
        Twall = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }
      Twall_defined = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix surface/global command");
  }

  if (heat_flag != Twall_defined)
    error->all(FLERR, "Must define wall temperature with heat model");

  // initializations

  dimension = domain->dimension;

  points = nullptr;
  lines = nullptr;
  tris = nullptr;

  nmotion = maxmotion = 0;
  motions = NULL;

  points_lastneigh = nullptr;
  points_original = nullptr;
  xsurf_original = nullptr;

  neigh_p1 = neigh_p2 = nullptr;
  pwhich_p1 = pwhich_p2 = nullptr;
  nside_p1 = nside_p2 = nullptr;
  aflag_p1 = aflag_p2 = nullptr;

  neigh_e1 = neigh_e2 = neigh_e3 = nullptr;
  ewhich_e1 = ewhich_e2 = ewhich_e3 = nullptr;
  nside_e1 = nside_e2 = nside_e3 =nullptr;
  aflag_e1 = aflag_e2 = aflag_e3 = nullptr;
  neigh_c1 = neigh_c2 = neigh_c3 = nullptr;
  cwhich_c1 = cwhich_c2 = cwhich_c3 = nullptr;

  connect2d = nullptr;
  connect3d = nullptr;

  xsurf = vsurf = omegasurf = nullptr;
  radsurf = nullptr;

  contact_surfs = nullptr;
  contact_forces = nullptr;
  nmax_contact_surfs = 0;
  nmax_contact_forces = 0;

  nmax = 0;
  mass_rigid = nullptr;

  fix_rigid = nullptr;
  fix_history = nullptr;

  list = new NeighList(lmp);
  if (use_history) {
    listhistory = new NeighList(lmp);
    zeroes = new double[size_history];
    for (int i = 0; i < size_history; i++) zeroes[i] = 0.0;
  } else {
    listhistory = nullptr;
    zeroes = nullptr;
  }

  imax = 0;
  imflag = nullptr;
  imdata = nullptr;

  firsttime = 1;

  // define points/lines/tris and their connectivity
  // done via molecule template ID or STL file
  // first check if arg is valid molecute template ID, else treat as STL file

  int imol = atom->find_molecule(arg[3]);
  if (imol >= 0) extract_from_molecules(arg[3]);
  else extract_from_stlfile(arg[3]);

  if (dimension == 2) connectivity2d_global();
  else connectivity3d_global();

  nsurf = nlines;
  if (dimension == 3) nsurf = ntris;

  surface_attributes();
}

/* ---------------------------------------------------------------------- */

FixSurfaceGlobal::~FixSurfaceGlobal()
{
  memory->sfree(points);
  memory->sfree(lines);
  memory->sfree(tris);

  // NOTE: need to free motions and their contents

  memory->destroy(points_lastneigh);
  memory->destroy(points_original);
  memory->destroy(xsurf_original);

  memory->destroy(neigh_p1);
  memory->destroy(neigh_p2);
  memory->destroy(pwhich_p1);
  memory->destroy(pwhich_p2);
  memory->destroy(nside_p1);
  memory->destroy(nside_p2);
  memory->destroy(aflag_p1);
  memory->destroy(aflag_p2);

  memory->destroy(neigh_e1);
  memory->destroy(neigh_e2);
  memory->destroy(neigh_e3);
  memory->destroy(ewhich_e1);
  memory->destroy(ewhich_e2);
  memory->destroy(ewhich_e3);
  memory->destroy(nside_e1);
  memory->destroy(nside_e2);
  memory->destroy(nside_e3);
  memory->destroy(aflag_e1);
  memory->destroy(aflag_e2);
  memory->destroy(aflag_e3);

  memory->destroy(neigh_c1);
  memory->destroy(neigh_c2);
  memory->destroy(neigh_c3);
  memory->destroy(cwhich_c1);
  memory->destroy(cwhich_c2);
  memory->destroy(cwhich_c3);

  memory->sfree(connect2d);
  memory->sfree(connect3d);

  memory->sfree(contact_surfs);
  memory->destroy(contact_surfs);

  memory->sfree(contact_forces);
  memory->destroy(contact_forces);

  memory->destroy(xsurf);
  memory->destroy(vsurf);
  memory->destroy(omegasurf);
  memory->destroy(radsurf);

  memory->destroy(mass_rigid);

  delete list;
  delete listhistory;
  delete [] zeroes;
  delete[] tstr;

  if (use_history)
    modify->delete_fix("NEIGH_HISTORY_SURFACE_GLOBAL_" + std::to_string(instance_me));

  memory->destroy(imflag);
  memory->destroy(imdata);
}

/* ----------------------------------------------------------------------
   create Fix needed for storing shear history if needed
   must be done in post_constructor()
------------------------------------------------------------------------- */

void FixSurfaceGlobal::post_constructor()
{
  if (use_history) {
    auto cmd = fmt::format("NEIGH_HISTORY_SURFACE_GLOBAL_" + std::to_string(instance_me) + " all NEIGH_HISTORY {}",size_history);
    fix_history = dynamic_cast<FixNeighHistory *>(modify->add_fix(cmd));
  } else
    fix_history = nullptr;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceGlobal::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::init()
{
  dt = update->dt;
  triggersq = 0.25 * neighbor->skin * neighbor->skin;

  // check for compatible heat conduction atom style

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR, "Heat conduction in fix surface/global requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR, "Heat conduction in fix surface/global requires atom style with heatflow property");
  }

  // define history indices

  int next_index = 0;
  if (model->beyond_contact) //next_index = 1;
    error->all(FLERR, "Beyond contact models not currenty supported");

  for (int i = 0; i < NSUBMODELS; i++) {
    model->sub_models[i]->history_index = next_index;
    next_index += model->sub_models[i]->size_history;
  }

  // one-time setup and allocation of neighbor list
  // wait until now, so neighbor settings have been made

  if (firsttime) {
    firsttime = 0;
    int pgsize = neighbor->pgsize;
    int oneatom = neighbor->oneatom;
    list->setup_pages(pgsize,oneatom);
    list->grow(atom->nmax,atom->nmax);

    if (use_history) {
      listhistory->setup_pages(pgsize,oneatom);
      listhistory->grow(atom->nmax,atom->nmax);
    }
  }

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) error->all(FLERR, "Variable {} for fix surface/global does not exist", tstr);
    if (! input->variable->equalstyle(tvar))
      error->all(FLERR, "Variable {} for fix surface/global must be an equal style variable", tstr);
  }
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   move surfaces via fix_modify setting
   similar to fix move operations
------------------------------------------------------------------------- */

void FixSurfaceGlobal::initial_integrate(int vflag)
{
  /*
  double ddotr;
  double a[3],b[3],c[3],d[3],disp[3],p12[3],p13[3];
  double *pt,*p1,*p2,*p3;

  double delta = (update->ntimestep - time_origin) * dt;

  // for rotate by right-hand rule around omega:
  // P = point = vector = point of rotation
  // R = vector = axis of rotation
  // w = omega of rotation (from period)
  // X0 = xoriginal = initial coord of atom
  // R0 = runit = unit vector for R
  // D = X0 - P = vector from P to X0
  // C = (D dot R0) R0 = projection of atom coord onto R line
  // A = D - C = vector from R line to X0
  // B = R0 cross A = vector perp to A in plane of rotation
  // A,B define plane of circular rotation around R line
  // X = P + C + A cos(w*dt) + B sin(w*dt)
  // V = w R0 cross (A cos(w*dt) + B sin(w*dt))

  if (mstyle == ROTATE) {
    double arg = omega_rotate * delta;
    double cosine = cos(arg);
    double sine = sin(arg);

    for (int i = 0; i < npoints; i++) {
      d[0] = points_original[i][0] - rpoint[0];
      d[1] = points_original[i][1] - rpoint[1];
      d[2] = points_original[i][2] - rpoint[2];

      ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
      c[0] = ddotr*runit[0];
      c[1] = ddotr*runit[1];
      c[2] = ddotr*runit[2];
      a[0] = d[0] - c[0];
      a[1] = d[1] - c[1];
      a[2] = d[2] - c[2];
      b[0] = runit[1]*a[2] - runit[2]*a[1];
      b[1] = runit[2]*a[0] - runit[0]*a[2];
      b[2] = runit[0]*a[1] - runit[1]*a[0];
      disp[0] = a[0]*cosine  + b[0]*sine;
      disp[1] = a[1]*cosine  + b[1]*sine;
      disp[2] = a[2]*cosine  + b[2]*sine;

      pt = points[i].x;
      pt[0] = rpoint[0] + c[0] + disp[0];
      pt[1] = rpoint[1] + c[1] + disp[1];
      pt[2] = rpoint[2] + c[2] + disp[2];
    }

    for (int i = 0; i < nsurf; i++) {
      d[0] = xsurf_original[i][0] - rpoint[0];
      d[1] = xsurf_original[i][1] - rpoint[1];
      d[2] = xsurf_original[i][2] - rpoint[2];
      ddotr = d[0]*runit[0] + d[1]*runit[1] + d[2]*runit[2];
      c[0] = ddotr*runit[0];
      c[1] = ddotr*runit[1];
      c[2] = ddotr*runit[2];
      a[0] = d[0] - c[0];
      a[1] = d[1] - c[1];
      a[2] = d[2] - c[2];
      b[0] = runit[1]*a[2] - runit[2]*a[1];
      b[1] = runit[2]*a[0] - runit[0]*a[2];
      b[2] = runit[0]*a[1] - runit[1]*a[0];
      disp[0] = a[0]*cosine  + b[0]*sine;
      disp[1] = a[1]*cosine  + b[1]*sine;
      disp[2] = a[2]*cosine  + b[2]*sine;

      xsurf[i][0] = rpoint[0] + c[0] + disp[0];
      xsurf[i][1] = rpoint[1] + c[1] + disp[1];
      xsurf[i][2] = rpoint[2] + c[2] + disp[2];
      vsurf[i][0] = omega_rotate * (runit[1]*disp[2] - runit[2]*disp[1]);
      vsurf[i][1] = omega_rotate * (runit[2]*disp[0] - runit[0]*disp[2]);
      vsurf[i][2] = omega_rotate * (runit[0]*disp[1] - runit[1]*disp[0]);
    }

    if (dimension == 3) {
      for (int i = 0; i < nsurf; i++) {
        p1 = points[tris[i].p1].x;
        p2 = points[tris[i].p2].x;
        p3 = points[tris[i].p3].x;
        MathExtra::sub3(p1,p2,p12);
        MathExtra::sub3(p1,p3,p13);
        MathExtra::cross3(p12,p13,tris[i].norm);
        MathExtra::norm3(tris[i].norm);
      }
    }
  }

  // trigger reneighbor if any point has moved skin/2 distance

  double dx,dy,dz,rsq;

  int triggerflag = 0;

  if (mstyle != NONE) {
    for (int i = 0; i < npoints; i++) {
      pt = points[i].x;
      dx = pt[0] - points_lastneigh[i][0];
      dy = pt[1] - points_lastneigh[i][1];
      dz = pt[2] - points_lastneigh[i][2];
      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq > triggersq) {
        triggerflag = 1;
        break;
      }
    }
  }

  if (triggerflag) next_reneighbor = update->ntimestep;

  */
}

/* ----------------------------------------------------------------------
   build neighbor list for sphere/surf interactions
   I = sphere, J = surf
   similar to methods in neigh_gran.cpp
------------------------------------------------------------------------- */

void FixSurfaceGlobal::pre_neighbor()
{
  int i,j,m,n,nn,dnum,dnumbytes;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,rsq,radsum,cutsq;
  int *neighptr,*touchptr;
  double *shearptr;

  int *npartner;
  tagint **partner;
  double **shearpartner;
  int **firsttouch;
  double **firstshear;
  MyPage<int> *ipage_touch;
  MyPage<double> *dpage_shear;

  double **x = atom->x;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double skin = neighbor->skin;

  list->grow(nlocal,nall);
  if (use_history) listhistory->grow(nlocal,nall);

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  if (use_history) {
    fix_history->nlocal_neigh = nlocal;
    //npartner = fix_history->npartner;
    //partner = fix_history->partner;
    //shearpartner = fix_history->shearpartner;
    firsttouch = fix_history->firstflag;
    firstshear = fix_history->firstvalue;
    //ipage_touch = listhistory->ipage;
    //dpage_shear = listhistory->dpage;
    //dnum = listhistory->dnum;
    dnumbytes = dnum * sizeof(double);
  }

  // store current point positions for future neighbor trigger check
  // check is performed in intitial_integrate()

  /*
  if (mstyle != NONE) {
    for (i = 0; i < npoints; i++) {
      points_lastneigh[i][0] = points[i].x[0];
      points_lastneigh[i][1] = points[i].x[1];
      points_lastneigh[i][2] = points[i].x[2];
    }
  }
  */

  int inum = 0;
  ipage->reset();
  if (use_history) {
    ipage_touch->reset();
    dpage_shear->reset();
  }

  for (i = 0; i < nlocal; i++) {
    n = 0;
    neighptr = ipage->vget();
    if (use_history) {
      nn = 0;
      touchptr = ipage_touch->vget();
      shearptr = dpage_shear->vget();
    }

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    // for now, loop over all surfs
    // NOTE: use a more sophisticated neighbor check

    for (j = 0; j < nsurf; j++) {
      delx = xtmp - xsurf[j][0];
      dely = ytmp - xsurf[j][1];
      delz = ztmp - xsurf[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radsum = radi + radsurf[j] + skin;
      cutsq = radsum*radsum;
      if (rsq <= cutsq) {
        neighptr[n] = j;

        if (use_history) {
          if (rsq < radsum*radsum) {
            for (m = 0; m < npartner[i]; m++)
              if (partner[i][m] == j) break;
            if (m < npartner[i]) {
              touchptr[n] = 1;
              memcpy(&shearptr[nn],&shearpartner[i][dnum*m],dnumbytes);
              nn += dnum;
            } else {
              touchptr[n] = 0;
              memcpy(&shearptr[nn],zeroes,dnumbytes);
              nn += dnum;
            }
          } else {
            touchptr[n] = 0;
            memcpy(&shearptr[nn],zeroes,dnumbytes);
            nn += dnum;
          }
        }

        n++;
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Fix surface/global neighbor list overflow, "
                 "boost neigh_modify one");

    if (use_history) {
      firsttouch[i] = touchptr;
      firstshear[i] = shearptr;
      ipage_touch->vgot(n);
      dpage_shear->vgot(nn);
    }
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   compute particle/surface interactions
   impart force and torque to spherical particles
------------------------------------------------------------------------- */

void FixSurfaceGlobal::post_force(int vflag)
{
  int i,j,k,n,m,ii,jj,inum,jnum,jflag,otherflag;
  int n_contact_surfs, n_contact_forces;
  double xtmp,ytmp,ztmp,radi,delx,dely,delz;
  double meff;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch,touch_flag;
  double rsq,radsum;
  double dr[3],contact[3],ds[3],vs[3],*forces,*torquesi;
  double *history,*allhistory,**firsthistory;

  model->history_update = 1;
  if (update->setupflag) model->history_update = 0;

  // if just reneighbored:
  // update rigid body masses for owned atoms if using FixRigid
  //   body[i] = which body atom I is in, -1 if none
  //   mass_body = mass of each rigid body

  if (neighbor->ago == 0 && fix_rigid) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body",tmp);
    double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid,nmax,"surface/global:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++) {
      if (body[i] >= 0) mass_rigid[i] = mass_body[body[i]];
      else mass_rigid[i] = 0.0;
    }
  }

  // loop over neighbors of my atoms
  // I is always sphere, J is always line

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (heat_flag) {
    if (tstr)
      Twall = input->variable->compute_equal(tvar);
    model->Tj = Twall;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firsttouch = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    model->xi = x[i];
    model->radi = radius[i];
    model->vi = v[i];
    model->omegai = omega[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }

    n_contact_surfs = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - xsurf[j][0];
      dely = xtmp - xsurf[j][1];
      delz = xtmp - xsurf[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // skip contact check if particle/surf are too far apart

      radsum = radi + radsurf[j];
      if (rsq > radsum * radsum) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // check for contact between particle and line/tri

      if (dimension == 2) {

        // check for overlap of sphere and line segment
        // jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
        // if no overlap, just continue
        // for overlap, also return:
        //   contact = nearest point on line to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = SurfExtra::
          overlap_sphere_line(x[i],radius[i],
                              points[lines[j].p1].x,points[lines[j].p2].x,
                              contact,dr,rsq);
      } else {

        // check for overlap of sphere and triangle
        // jflag = 0 for no overlap, 1 for interior line pt,
        //   -1/-2/-3 for 3 edges, -4/-5/-6 for 3 corner pts
        // if no overlap, just continue
        // for overlap, also returns:
        //   contact = nearest point on tri to sphere center
        //   dr = vector from contact pt to sphere center
        //   rsq = squared length of dr

        jflag = SurfExtra::
          overlap_sphere_tri(x[i],radius[i],
                             points[tris[j].p1].x,points[tris[j].p2].x,
                             points[tris[j].p3].x,tris[j].norm,
                             contact,dr,rsq);
      }

      if (!jflag) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // append surf to list of contacts

      if (n_contact_surfs + 1 > nmax_contact_surfs) {
        nmax_contact_surfs += DELTACONTACTS;
        memory->grow(contact_surfs, nmax_contact_surfs * sizeof(ContactSurf),
                                      "surface/global:contact_surfs");
      }

      if (dimension == 2) {
        contact_surfs[n_contact_surfs].index = j;
        contact_surfs[n_contact_surfs].type = lines[j].type;
        contact_surfs[n_contact_surfs].overlap = sqrt(rsq) - radsum;
        contact_surfs[n_contact_surfs].r[0] = dr[0];
        contact_surfs[n_contact_surfs].r[1] = dr[1];
        contact_surfs[n_contact_surfs].r[2] = dr[2];
      }

      n_contact_surfs += 1;
    }

    // Reduce set of contacts
    std::sort(contact_surfs, contact_surfs + n_contact_surfs, [](ContactSurf a, ContactSurf b) {return a.overlap > b.overlap;});

    n_contact_forces = 0;
    processed_contacts.clear();
    connected_contacts.clear();

    for (n = 0; n < n_contact_surfs; n++) {
      j = contact_surfs[n].index;

      if (processed_contacts.find(j) != processed_contacts.end()) continue;
      processed_contacts.insert(j);

      n_contact_forces += 1;

      if (n_contact_forces + 1 > nmax_contact_forces) {
        nmax_contact_forces += DELTACONTACTS;
        memory->grow(contact_forces, nmax_contact_forces * sizeof(ContactForce),
                                      "surface/global:contact_forces");
      }

      if (dimension == 2) {
        for (m = 0; m < connect2d[j].np1; m++) {
          k = connect2d[j].neigh_p1[m];
          if (processed_contacts.find(k) != processed_contacts.end())
            connected_contacts.insert(k);
        }
        for (m = 0; m < connect2d[j].np2; m++) {
          k = connect2d[j].neigh_p2[m];
          if (processed_contacts.find(k) != processed_contacts.end())
            connected_contacts.insert(k);
        }

        while (connected_contacts.size() != 0) {
          auto itr = connected_contacts.begin();
          k = *itr;
          connected_contacts.erase(itr);
          processed_contacts.insert(k);

          // if flat, convex, concave...
          //  need to identify which face is in contact,
          //  do I need to save j & m when walking across flats?
        }
      } else {
        continue; // 3D
      }
    }


    /*
    For contact in reduced contacts:
      // reset model and copy initial geometric data

      factor_lj = special_lj[sbmask(j)]; // presumably not necessary
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      model->xj = xsurf[j];
      model->radj = radsurf[j];
      if (use_history) model->touch = touch[jj];

      // unset non-touching neighbors
      // NOTE: in pair_surf_granular, this unsetting occurs twice ?
      // NOTE: maybe it should only be below, after call to overlap() methods ?

      touch_flag = model->check_contact();

      if (!touch_flag) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      rsq = model->rsq;

      // meff = effective mass of sphere
      // if I is part of rigid body, use body mass

      meff = rmass[i];
      if (fix_rigid && mass_rigid[i] > 0.0) meff = mass_rigid[i];

      // copy additional information and prepare force calculations

      model->meff = meff;

      ds[0] = contact[0] - xsurf[j][0];
      ds[1] = contact[1] - xsurf[j][1];
      ds[2] = contact[2] - xsurf[j][2];

      vs[0] = vsurf[j][0] + (omegasurf[j][1] * ds[2] - omegasurf[j][2] * ds[1]);
      vs[1] = vsurf[j][1] + (omegasurf[j][2] * ds[0] - omegasurf[j][0] * ds[2]);
      vs[2] = vsurf[j][2] + (omegasurf[j][0] * ds[1] - omegasurf[j][1] * ds[0]);

      model->vj = vs;
      model->omegaj = omegasurf[j];

      if (heat_flag) model->Ti = temperature[i];

      // pairwise interaction between sphere and surface element

      if (use_history) {
        touch[jj] = 1;
        history = &allhistory[3*jj];
        model->history = history;
      }

      model->dx[0] = dr[0];
      model->dx[1] = dr[1];
      model->dx[2] = dr[2];

      // need to add support coupled contacts
      // is this just multiplying forces (+torques?) by factor_couple?

      model->calculate_forces();

      forces = model->forces;
      torquesi = model->torquesi;

      // apply forces & torques

      add3(f[i], forces, f[i]);
      add3(torque[i], torquesi, torque[i]);
      if (heat_flag) heatflow[i] += model->dq;
    */
  }
}

/* ----------------------------------------------------------------------
   process fix_modify commands specific to fix surface/global
   when motion is active, set INTITIAL_INTEGRATE mask flag
------------------------------------------------------------------------- */

int FixSurfaceGlobal::modify_param(int narg, char **arg)
{
  // move options

  if (strcmp(arg[0],"move") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");

    int ifix = modify->find_fix(id);

    if (strcmp(arg[1],"none") == 0) {
      nmotion = 0;
      modify->fmask[ifix] &= ~INITIAL_INTEGRATE;
      force_reneighbor = 0;
      next_reneighbor = -1;

      // NOTE: this data exists for every motion instance ?

      move_clear();
      return 2;
    }

    // add a new motion operation

    modify->fmask[ifix] |= INITIAL_INTEGRATE;

    if (nmotion == maxmotion) {
      maxmotion++;
      motions = (Motion *)
        memory->srealloc(motions,maxmotion*sizeof(Motion),"fix_surface_global::motion");
    }

    int argcount = modify_params_move(&motions[nmotion],narg-1,&arg[1]);
    nmotion++;

    force_reneighbor = 1;
    next_reneighbor = -1;

    return argcount + 2;
  }

  // keyword not recognized

  return 0;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceGlobal::modify_params_move(Motion *motion, int narg, char **arg)
{
  if (strcmp(arg[0],"linear") == 0) {
    if (narg < 4) error->all(FLERR,"Illegal fix_modify move command");
    motion->mstyle = LINEAR;

    if (strcmp(arg[1], "NULL") == 0) motion->vxflag = 0;
    else {
      motion->vxflag = 1;
      motion->vx = utils::numeric(FLERR, arg[4], false, lmp);
    }
    if (strcmp(arg[2], "NULL") == 0) motion->vyflag = 0;
    else {
      motion->vyflag = 1;
      motion->vy = utils::numeric(FLERR, arg[2], false, lmp);
    }
    if (strcmp(arg[3], "NULL") == 0) motion->vzflag = 0;
    else {
      motion->vzflag = 1;
      motion->vz = utils::numeric(FLERR, arg[3], false, lmp);
    }

    return 4;
  }

  if (strcmp(arg[0],"wiggle") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal fix_modify move command");
    motion->mstyle = WIGGLE;

    if (strcmp(arg[1], "NULL") == 0) motion->axflag = 0;
    else {
      motion->axflag = 1;
      motion->ax = utils::numeric(FLERR, arg[4], false, lmp);
    }
    if (strcmp(arg[2], "NULL") == 0) motion->ayflag = 0;
    else {
      motion->ayflag = 1;
      motion->ay = utils::numeric(FLERR, arg[2], false, lmp);
    }
    if (strcmp(arg[3], "NULL") == 0) motion->azflag = 0;
    else {
      motion->azflag = 1;
      motion->az = utils::numeric(FLERR, arg[3], false, lmp);
    }

    motion->period = utils::numeric(FLERR, arg[7], false, lmp);
    if (motion->period <= 0.0) error->all(FLERR, "Illegal fix_modify move command");

    return 5;
  }

  if (strcmp(arg[0],"rotate") == 0) {
    if (narg < 8) error->all(FLERR,"Illegal fix_modify move command");
    motion->mstyle = ROTATE;

    motion->point[0] = utils::numeric(FLERR,arg[1],false,lmp);
    motion->point[1] = utils::numeric(FLERR,arg[2],false,lmp);
    motion->point[2] = utils::numeric(FLERR,arg[3],false,lmp);

    motion->axis[0] = utils::numeric(FLERR,arg[4],false,lmp);
    motion->axis[1] = utils::numeric(FLERR,arg[5],false,lmp);
    motion->axis[2] = utils::numeric(FLERR,arg[6],false,lmp);
    if (dimension == 2)
      if (motion->mstyle == ROTATE && (motion->axis[0] != 0.0 || motion->axis[1] != 0.0))
        error->all(FLERR,"Fix_modify move cannot rotate around "
                   "non z-axis for 2d problem");

    motion->period = utils::numeric(FLERR,arg[7],false,lmp);
    if (motion->period <= 0.0) error->all(FLERR,"Illegal fix_modify move command");

    motion->time_origin = update->ntimestep;
    motion->omega = MY_2PI / motion->period;

    // runit = unit vector along rotation axis

    double len = MathExtra::len3(motion->axis);
    if (len == 0.0)
      error->all(FLERR,"Fix_modify move zero length rotation vector");
    MathExtra::normalize3(motion->axis,motion->unit);

    // NOTE: how do these 2 operations now work

    move_clear();
    move_init();

    return 8;
  }

  if (strcmp(arg[0],"transrot") == 0) {
    if (narg < 11) error->all(FLERR,"Illegal fix_modify move command");
    motion->mstyle = TRANSROT;

    motion->vxflag = motion->vyflag = motion->vzflag = 1;
    motion->vx = utils::numeric(FLERR, arg[1], false, lmp);
    motion->vy = utils::numeric(FLERR, arg[2], false, lmp);
    motion->vz = utils::numeric(FLERR, arg[3], false, lmp);

    motion->point[0] = utils::numeric(FLERR,arg[4],false,lmp);
    motion->point[1] = utils::numeric(FLERR,arg[5],false,lmp);
    motion->point[2] = utils::numeric(FLERR,arg[6],false,lmp);

    motion->axis[0] = utils::numeric(FLERR,arg[7],false,lmp);
    motion->axis[1] = utils::numeric(FLERR,arg[8],false,lmp);
    motion->axis[2] = utils::numeric(FLERR,arg[9],false,lmp);
    if (dimension == 2)
      if (motion->mstyle == ROTATE && (motion->axis[0] != 0.0 || motion->axis[1] != 0.0))
        error->all(FLERR,"Fix_modify move cannot rotate around "
                   "non z-axis for 2d problem");

    motion->period = utils::numeric(FLERR,arg[10],false,lmp);
    if (motion->period <= 0.0) error->all(FLERR,"Illegal fix_modify move command");

    motion->time_origin = update->ntimestep;
    motion->omega = MY_2PI / motion->period;

    // runit = unit vector along rotation axis

    double len = MathExtra::len3(motion->axis);
    if (len == 0.0)
      error->all(FLERR,"Fix_modify move zero length rotation vector");
    MathExtra::normalize3(motion->axis,motion->unit);

    // NOTE: how do these 2 operations now work

    move_clear();
    move_init();

    return 11;
  }

  if (strcmp(arg[0],"variable") == 0) {
  }

  //error return;

  //return iarg;
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixSurfaceGlobal::reset_dt()
{
  /*
  if (mstyle != NONE)
    error->all(FLERR,"Resetting timestep size is not allowed with "
               "fix surface/global motion");
  */
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSurfaceGlobal::memory_usage()
{
  // NOTE: need to include neigh lists

  double bytes = 0.0;
  bytes += npoints*sizeof(Point);
  if (dimension == 2) {
    bytes += nlines*sizeof(Line);
    bytes = nlines*sizeof(Connect2d);
  } else {
    bytes += ntris*sizeof(Tri);
    bytes = ntris*sizeof(Connect3d);
  }
  return bytes;
}

/* ----------------------------------------------------------------------
   extract neighbor lists
------------------------------------------------------------------------- */

void *FixSurfaceGlobal::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"list") == 0) return list;
  else if (strcmp(str,"listhistory") == 0) return listhistory;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

int FixSurfaceGlobal::image(int *&ivec, double **&darray)
{
  int n;
  double *p1,*p2,*p3;

  if (dimension == 2) {
    n = nlines;

    if (imax == 0) {
      imax = n;
      memory->create(imflag,imax,"surface/global:imflag");
      memory->create(imdata,imax,7,"surface/global:imflag");
    }

    for (int i = 0; i < n; i++) {
      p1 = points[lines[i].p1].x;
      p2 = points[lines[i].p2].x;

      imflag[i] = LINE;
      imdata[i][0] = lines[i].type;
      imdata[i][1] = p1[0];
      imdata[i][2] = p1[1];
      imdata[i][3] = p1[2];
      imdata[i][4] = p2[0];
      imdata[i][5] = p2[1];
      imdata[i][6] = p2[2];
    }

  } else {
    n = ntris;

    if (imax == 0) {
      imax = n;
      memory->create(imflag,imax,"surface/global:imflag");
      memory->create(imdata,imax,10,"surface/global:imflag");
    }

    for (int i = 0; i < n; i++) {
      p1 = points[tris[i].p1].x;
      p2 = points[tris[i].p2].x;
      p3 = points[tris[i].p3].x;

      imflag[i] = TRI;
      imdata[i][0] = tris[i].type;
      imdata[i][1] = p1[0];
      imdata[i][2] = p1[1];
      imdata[i][3] = p1[2];
      imdata[i][4] = p2[0];
      imdata[i][5] = p2[1];
      imdata[i][6] = p2[2];
      imdata[i][7] = p3[0];
      imdata[i][8] = p3[1];
      imdata[i][9] = p3[2];
    }
  }

  ivec = imflag;
  darray = imdata;
  return n;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// initializiation of surfs
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

/* ----------------------------------------------------------------------
   extract lines or surfs from molecule template ID for one or more molecules
   concatenate into single list of lines and tris
   identify unique points using hash
   this creates points and lines or tris data structs (2d or 3d)
------------------------------------------------------------------------- */

void FixSurfaceGlobal::extract_from_molecules(char *molID)
{
  // populate global point/line/tri data structs

  npoints = nlines = ntris = 0;
  int maxpoints = 0;

  int imol = atom->find_molecule(molID);
  if (imol == -1)
    error->all(FLERR,"Molecule template ID for fix surface/global does not exist");

  // loop over one or more molecules in molID

  Molecule **onemols = &atom->molecules[imol];
  int nmol = onemols[0]->nset;

  for (int m = 0; m < nmol; m++) {
    if (dimension == 2)
      if (onemols[m]->lineflag == 0)
        error->all(FLERR,"Fix surface/global molecule must have lines");
    if (dimension == 3)
      if (onemols[m]->triflag == 0)
        error->all(FLERR,"Fix surface/global molecule must have triangles");

    int nl = onemols[m]->nlines;
    int nt = onemols[m]->ntris;

    nlines += nl;
    ntris += nt;
    lines = (Line *) memory->srealloc(lines,nlines*sizeof(Line),
                                      "surface/global:lines");
    tris = (Tri *) memory->srealloc(tris,ntris*sizeof(Tri),
                                    "surface/global:tris");

    // create a map
    // key = xyz coords of a point
    // value = index into unique points vector

    std::map<std::tuple<double,double,double>,int> hash;

    // offset line/tri index lists by previous npoints
    // pi,p2,p3 are C-style indices into points vector

    if (dimension == 2) {
      int *molline = onemols[m]->molline;
      int *typeline = onemols[m]->typeline;
      double **epts = onemols[m]->lines;
      int iline = nlines - nl;

      for (int i = 0; i < nl; i++) {
        lines[iline].mol = molline[i];
        lines[iline].type = typeline[i];

        auto key = std::make_tuple(epts[i][0],epts[i][1],0.0);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = epts[i][0];
          points[npoints].x[1] = epts[i][1];
          points[npoints].x[2] = 0.0;
          lines[iline].p1 = npoints;
          npoints++;
        } else lines[iline].p1 = hash[key];

        key = std::make_tuple(epts[i][2],epts[i][3],0.0);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = epts[i][2];
          points[npoints].x[1] = epts[i][3];
          points[npoints].x[2] = 0.0;
          lines[iline].p2 = npoints;
          npoints++;
        } else lines[iline].p2 = hash[key];

        iline++;
      }
    }

    if (dimension == 3) {
      int *moltri = onemols[m]->moltri;
      int *typetri = onemols[m]->typetri;
      double **cpts = onemols[m]->tris;
      int itri = ntris - nt;

      for (int i = 0; i < nt; i++) {
        tris[itri].mol = moltri[i];
        tris[itri].type = typetri[i];

        auto key = std::make_tuple(cpts[i][0],cpts[i][1],cpts[i][2]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][0];
          points[npoints].x[1] = cpts[i][1];
          points[npoints].x[2] = cpts[i][2];
          tris[itri].p1 = npoints;
          npoints++;
        } else tris[itri].p1 = hash[key];

        key = std::make_tuple(cpts[i][3],cpts[i][4],cpts[i][5]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][3];
          points[npoints].x[1] = cpts[i][4];
          points[npoints].x[2] = cpts[i][5];
          tris[itri].p2 = npoints;
          npoints++;
        } else tris[itri].p2 = hash[key];

        key = std::make_tuple(cpts[i][6],cpts[i][7],cpts[i][8]);
        if (hash.find(key) == hash.end()) {
          if (npoints == maxpoints) {
            maxpoints += DELTA;
            points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                                "surface/global:points");
          }
          hash[key] = npoints;
          points[npoints].x[0] = cpts[i][6];
          points[npoints].x[1] = cpts[i][7];
          points[npoints].x[2] = cpts[i][8];
          tris[itri].p3 = npoints;
          npoints++;
        } else tris[itri].p3 = hash[key];

        itri++;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   extract triangles from an STL file, can be text or binary
   identify unique points using hash
   this creates points and tris data structs (3d only)
------------------------------------------------------------------------- */

void FixSurfaceGlobal::extract_from_stlfile(char *filename)
{
  if (dimension == 2)
    error->all(FLERR,"Fix surface/global cannot use an STL file for 2d simulations");

  // read tris from STL file
  // stltris = tri coords internal to STL reader

  STLReader *stl = new STLReader(lmp);
  double **stltris;
  ntris = stl->read_file(filename,stltris);

  // create points and tris data structs

  npoints = 0;
  int maxpoints = 0;

  tris = (Tri *) memory->smalloc(ntris*sizeof(Tri),"surface/global:tris");

  // create a map
  // key = xyz coords of a point
  // value = index into unique points vector

  std::map<std::tuple<double,double,double>,int> hash;

  // loop over STL tris
  // populate points and tris data structs
  // set molecule and type of tri = 1

  for (int itri = 0; itri < ntris; itri++) {
    tris[itri].mol = 1;
    tris[itri].type = 1;

    auto key = std::make_tuple(stltris[itri][0],stltris[itri][1],stltris[itri][2]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][0];
      points[npoints].x[1] = stltris[itri][1];
      points[npoints].x[2] = stltris[itri][2];
      tris[itri].p1 = npoints;
      npoints++;
    } else tris[itri].p1 = hash[key];

    key = std::make_tuple(stltris[itri][3],stltris[itri][4],stltris[itri][5]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][3];
      points[npoints].x[1] = stltris[itri][4];
      points[npoints].x[2] = stltris[itri][5];
      tris[itri].p2 = npoints;
      npoints++;
    } else tris[itri].p2 = hash[key];

    key = std::make_tuple(stltris[itri][6],stltris[itri][7],stltris[itri][8]);
    if (hash.find(key) == hash.end()) {
      if (npoints == maxpoints) {
        maxpoints += DELTA;
        points = (Point *) memory->srealloc(points,maxpoints*sizeof(Point),
                                            "surface/global:points");
      }
      hash[key] = npoints;
      points[npoints].x[0] = stltris[itri][6];
      points[npoints].x[1] = stltris[itri][7];
      points[npoints].x[2] = stltris[itri][8];
      tris[itri].p3 = npoints;
      npoints++;
    } else tris[itri].p3 = hash[key];
  }

  // delete STL reader

  delete stl;
}

/* ----------------------------------------------------------------------
   create and initialize Connect2d info for all lines
   this creates plines and connect2d data structs
------------------------------------------------------------------------- */

void FixSurfaceGlobal::connectivity2d_global()
{
  connect2d = (Connect2d *)
    memory->smalloc(nlines*sizeof(Connect2d),
                    "surface/global:connect2d");

  // setup line end point connectivity lists
  // count # of lines containing each end point (including self)
  // plines = ragged 2d array with indices of lines which contain each point

  int *counts;
  memory->create(counts,npoints,"surface/global:counts");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    counts[lines[i].p1]++;
    counts[lines[i].p2]++;
  }

  memory->create_ragged(plines,npoints,counts,"surface/global:plines");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < nlines; i++) {
    plines[lines[i].p1][counts[lines[i].p1]++] = i;
    plines[lines[i].p2][counts[lines[i].p2]++] = i;
  }

  // allocate all ragged arrays which Connect2d will point to
  // first subtract self from each point's count
  //   b/c connect3d vectors will NOT include self

  for (int i = 0; i < npoints; i++) counts[i]--;

  memory->create_ragged(neigh_p1,npoints,counts,"surface/global:neigh_p1");
  memory->create_ragged(neigh_p2,npoints,counts,"surface/global:neigh_p2");
  memory->create_ragged(pwhich_p1,npoints,counts,"surface/global:pwhich_p1");
  memory->create_ragged(pwhich_p2,npoints,counts,"surface/global:pwhich_p2");
  memory->create_ragged(nside_p1,npoints,counts,"surface/global:nside_p1");
  memory->create_ragged(nside_p2,npoints,counts,"surface/global:nside_p2");
  memory->create_ragged(aflag_p1,npoints,counts,"surface/global:aflag_p1");
  memory->create_ragged(aflag_p2,npoints,counts,"surface/global:aflag_p2");

  // set connect2d vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < nlines; i++) {
    connect2d[i].np1 = counts[lines[i].p1];
    if (connect2d[i].np1 == 0) {
      connect2d[i].neigh_p1 = nullptr;
      connect2d[i].pwhich_p1 = nullptr;
      connect2d[i].nside_p1 = nullptr;
      connect2d[i].aflag_p1 = nullptr;
    } else {
      connect2d[i].neigh_p1 = neigh_p1[lines[i].p1];
      connect2d[i].pwhich_p1 = pwhich_p1[lines[i].p1];
      connect2d[i].nside_p1 = nside_p1[lines[i].p1];
      connect2d[i].aflag_p1 = aflag_p1[lines[i].p1];
    }

    connect2d[i].np2 = counts[lines[i].p2];
    if (connect2d[i].np2 == 0) {
      connect2d[i].neigh_p2 = nullptr;
      connect2d[i].pwhich_p2 = nullptr;
      connect2d[i].nside_p2 = nullptr;
      connect2d[i].aflag_p2 = nullptr;
    } else {
      connect2d[i].neigh_p2 = neigh_p2[lines[i].p2];
      connect2d[i].pwhich_p2 = pwhich_p2[lines[i].p2];
      connect2d[i].nside_p2 = nside_p2[lines[i].p2];
      connect2d[i].aflag_p2 = aflag_p2[lines[i].p2];
    }
  }

  // initialize connect2d neigh vectors for each end point of each line
  // do NOT include self

  int j,m;

  for (int i = 0; i < nlines; i++) {
    if (connect2d[i].np1) {
      j = 0;
      for (m = 0; m <= connect2d[i].np1; m++) {
        if (plines[lines[i].p1][m] == i) continue;
        connect2d[i].neigh_p1[j] = plines[lines[i].p1][m];
        j++;
      }
    }
    if (connect2d[i].np2) {
      j = 0;
      for (m = 0; m <= connect2d[i].np2; m++) {
        if (plines[lines[i].p2][m] == i) continue;
        connect2d[i].neigh_p2[j] = plines[lines[i].p2][m];
        j++;
      }
    }
  }

  // deallocate counts and plines

  memory->destroy(counts);
  memory->destroy(plines);

  // set connect2d pwhich/nside/aflag for each end point of each line
  // see fsg.h file for an explanation of each vector in Connect2d
  // aflag is based on dot and cross product of 2 connected line normals
  //   cross product is either along +z or -z direction

  double dotline,dotnorm;
  double *inorm,*jnorm;
  double icrossj[3];

  for (int i = 0; i < nlines; i++) {
    for (m = 0; m < connect2d[i].np1; m++) {
      j = connect2d[i].neigh_p1[m];

      inorm = lines[i].norm;
      jnorm = lines[j].norm;
      dotnorm = MathExtra::dot3(inorm,jnorm);

      if (lines[i].p1 == lines[j].p1) {
        connect2d[i].pwhich_p1[j] = 0;
        connect2d[i].nside_p1[j] = OPPOSITE_SIDE;
        if (dotnorm < -1.0+FLATTHRESH) connect2d[i].aflag_p1[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
          if (icrossj[2] > 0.0) connect2d[i].aflag_p1[j] = CONCAVE;
          else connect2d[i].aflag_p1[j] = CONVEX;
        }
      } else if (lines[i].p1 == lines[j].p2) {
        connect2d[i].pwhich_p1[j] = 1;
        connect2d[i].nside_p1[j] = SAME_SIDE;
        if (dotnorm > 1.0-FLATTHRESH) connect2d[i].aflag_p1[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
          if (icrossj[2] < 0.0) connect2d[i].aflag_p1[j] = CONCAVE;
          else connect2d[i].aflag_p1[j] = CONVEX;
        }
      }
    }

    for (m = 0; m < connect2d[i].np2; m++) {
      j = connect2d[i].neigh_p2[m];

      inorm = lines[i].norm;
      jnorm = lines[j].norm;
      dotnorm = MathExtra::dot3(inorm,jnorm);

      if (lines[i].p2 == lines[j].p1) {
        connect2d[i].pwhich_p2[j] = 0;
        connect2d[i].nside_p2[j] = SAME_SIDE;
        if (dotnorm > 1.0-FLATTHRESH) connect2d[i].aflag_p2[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
          if (icrossj[2] > 0.0) connect2d[i].aflag_p2[j] = CONCAVE;
          else connect2d[i].aflag_p2[j] = CONVEX;
        }
      } else if (lines[i].p2 == lines[j].p2) {
        connect2d[i].pwhich_p2[j] = 1;
        connect2d[i].nside_p2[j] = OPPOSITE_SIDE;
        if (dotnorm < -1.0+FLATTHRESH) connect2d[i].aflag_p2[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
          if (icrossj[2] < 0.0) connect2d[i].aflag_p2[j] = CONCAVE;
          else connect2d[i].aflag_p2[j] = CONVEX;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   create and initialize Connect3d info for all triangles
   this creates etris/ctris and connect3d data structs
------------------------------------------------------------------------- */

void FixSurfaceGlobal::connectivity3d_global()
{
  int p1,p2,p3;

  connect3d = (Connect3d *)
    memory->smalloc(ntris*sizeof(Connect3d),
                    "surface/global:connect3d");

  // create hash = map of unique edges
  //   key = <p1,p2> indices of 2 points, in either order
  //   value = index of the unique edge (0 to Nedge-1)
  // tri2edges[i][j] = index of unique edge for tri I and edge J
  // nedges = total count of unique edges

  int **tri2edge;
  memory->create(tri2edge,ntris,3,"surfface/global::tri2edge");

  std::map<std::tuple<int,int>,int> hash;
  int nedges = 0;

  for (int i = 0; i < ntris; i++) {
    p1 = tris[i].p1;
    p2 = tris[i].p2;
    p3 = tris[i].p3;

    auto key1 = std::make_tuple(p1,p2);
    auto key2 = std::make_tuple(p2,p1);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][0] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][0] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][0] = hash[key2];

    key1 = std::make_tuple(p2,p3);
    key2 = std::make_tuple(p3,p2);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][1] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][1] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][1] = hash[key2];

    key1 = std::make_tuple(p3,p1);
    key2 = std::make_tuple(p1,p3);

    if (hash.find(key1) == hash.end() && hash.find(key2) == hash.end()) {
      hash[key1] = nedges;
      tri2edge[i][2] = nedges;
      nedges++;
    }
    else if (hash.find(key1) != hash.end()) tri2edge[i][2] = hash[key1];
    else if (hash.find(key2) != hash.end()) tri2edge[i][2] = hash[key2];
  }

  // setup tri edge connectivity lists
  // count # of tris containing each edge (including self)
  // etris = ragged 2d array with indices of tris which contain each edge

  int *counts;
  memory->create(counts,nedges,"surface/global:count");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tri2edge[i][0]]++;
    counts[tri2edge[i][1]]++;
    counts[tri2edge[i][2]]++;
  }

  memory->create_ragged(etris,nedges,counts,"surface/global:etris");

  for (int i = 0; i < nedges; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    etris[tri2edge[i][0]][counts[tri2edge[i][0]]++] = i;
    etris[tri2edge[i][1]][counts[tri2edge[i][1]]++] = i;
    etris[tri2edge[i][2]][counts[tri2edge[i][2]]++] = i;
  }

  // allocate all edge ragged arrays which Connect3d will point to
  // first subtract self from each edge's count
  //   b/c connect3d vectors will NOT include self

  for (int i = 0; i < nedges; i++) counts[i]--;

  memory->create_ragged(neigh_e1,nedges,counts,"surface/global:neigh_e1");
  memory->create_ragged(neigh_e2,nedges,counts,"surface/global:neigh_e2");
  memory->create_ragged(neigh_e3,nedges,counts,"surface/global:neigh_e3");
  memory->create_ragged(ewhich_e1,nedges,counts,"surface/global:ewhich_e1");
  memory->create_ragged(ewhich_e2,nedges,counts,"surface/global:ewhich_e2");
  memory->create_ragged(ewhich_e3,nedges,counts,"surface/global:ewhich_e3");
  memory->create_ragged(nside_e1,nedges,counts,"surface/global:nside_e1");
  memory->create_ragged(nside_e2,nedges,counts,"surface/global:nside_e2");
  memory->create_ragged(nside_e3,nedges,counts,"surface/global:nside_e3");
  memory->create_ragged(aflag_e1,nedges,counts,"surface/global:aflag_e1");
  memory->create_ragged(aflag_e2,nedges,counts,"surface/global:aflag_e2");
  memory->create_ragged(aflag_e3,nedges,counts,"surface/global:aflag_e3");

  // set connect3d edge vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < ntris; i++) {
    connect3d[i].ne1 = counts[tri2edge[i][0]];
    if (connect3d[i].ne1 == 0) {
      connect3d[i].neigh_e1 = nullptr;
      connect3d[i].ewhich_e1 = nullptr;
      connect3d[i].nside_e1 = nullptr;
      connect3d[i].aflag_e1 = nullptr;
    } else {
      connect3d[i].neigh_e1 = neigh_e1[tri2edge[i][0]];
      connect3d[i].ewhich_e1 = ewhich_e1[tri2edge[i][0]];
      connect3d[i].nside_e1 = nside_e1[tri2edge[i][0]];
      connect3d[i].aflag_e1 = aflag_e1[tri2edge[i][0]];
    }

    connect3d[i].ne2 = counts[tri2edge[i][1]];
    if (connect3d[i].ne2 == 0) {
      connect3d[i].neigh_e2 = nullptr;
      connect3d[i].ewhich_e2 = nullptr;
      connect3d[i].nside_e2 = nullptr;
      connect3d[i].aflag_e2 = nullptr;
    } else {
      connect3d[i].neigh_e2 = neigh_e2[tri2edge[i][1]];
      connect3d[i].ewhich_e2 = ewhich_e2[tri2edge[i][1]];
      connect3d[i].nside_e2 = nside_e2[tri2edge[i][1]];
      connect3d[i].aflag_e2 = aflag_e2[tri2edge[i][1]];
    }

    connect3d[i].ne3 = counts[tri2edge[i][2]];
    if (connect3d[i].ne3 == 0) {
      connect3d[i].neigh_e3 = nullptr;
      connect3d[i].ewhich_e3 = nullptr;
      connect3d[i].nside_e3 = nullptr;
      connect3d[i].aflag_e3 = nullptr;
    } else {
      connect3d[i].neigh_e3 = neigh_e3[tri2edge[i][2]];
      connect3d[i].ewhich_e3 = ewhich_e3[tri2edge[i][2]];
      connect3d[i].nside_e3 = nside_e3[tri2edge[i][2]];
      connect3d[i].aflag_e3 = aflag_e3[tri2edge[i][2]];
    }
  }

  // initialize connect3d edge neigh vectors for each edge of each tri
  // do NOT include self

  int j,m;

  for (int i = 0; i < ntris; i++) {
    if (connect3d[i].ne1) {
      j = 0;
      for (m = 0; m <= connect3d[i].ne1; m++) {
        if (etris[tri2edge[i][0]][m] == i) continue;
        connect3d[i].neigh_e1[j] = etris[tri2edge[i][0]][m];
        j++;
      }
    }
    if (connect3d[i].ne2) {
      j = 0;
      for (m = 0; m <= connect3d[i].ne2; m++) {
        if (etris[tri2edge[i][1]][m] == i) continue;
        connect3d[i].neigh_e2[j] = etris[tri2edge[i][1]][m];
        j++;
      }
    }
    if (connect3d[i].ne3) {
      j = 0;
      for (m = 0; m <= connect3d[i].ne3; m++) {
        if (etris[tri2edge[i][2]][m] == i) continue;
        connect3d[i].neigh_e3[j] = etris[tri2edge[i][2]][m];
        j++;
      }
    }
  }

  // deallocate counts, tri2edge, etris

  memory->destroy(counts);
  memory->destroy(tri2edge);
  memory->destroy(etris);

  // set connect3d edge ewhich/nside/aflag for each edge of each tri
  // see fsg.h file for an explanation of each edge vector in Connect3d
  // aflag is based on dot and cross product of 2 connected tri normals
  //   cross product is either along itri edge or in opposite dir

  int jpfirst,jpsecond;
  double dotline,dotnorm;
  double *inorm,*jnorm;
  double icrossj[3],iedge[3];

  for (int i = 0; i < ntris; i++) {
    for (m = 0; m < connect3d[i].ne1; m++) {
      j = connect3d[i].neigh_e1[m];

      if (tris[i].p1 == tris[j].p1) jpfirst = 1;
      else if (tris[i].p1 == tris[j].p2) jpfirst = 2;
      else if (tris[i].p1 == tris[j].p3) jpfirst = 3;

      if (tris[i].p2 == tris[j].p1) jpsecond = 1;
      else if (tris[i].p2 == tris[j].p2) jpsecond = 2;
      else if (tris[i].p2 == tris[j].p3) jpsecond = 3;

      inorm = tris[i].norm;
      jnorm = tris[j].norm;
      dotnorm = MathExtra::dot3(inorm,jnorm);
      MathExtra::sub3(points[tris[i].p2].x,points[tris[i].p1].x,iedge);

      if ((jpfirst == 1 && jpsecond == 2) ||
	  (jpfirst == 2 && jpsecond == 3) ||
	  (jpfirst == 3 && jpsecond == 1)) {
	connect3d[i].ewhich_e1[j] = jpfirst - 1;
	connect3d[i].nside_e1[j] = OPPOSITE_SIDE;
	if (dotnorm < -1.0+FLATTHRESH) connect3d[i].aflag_e1[j] = FLAT;
	else {
	  MathExtra::cross3(inorm,jnorm,icrossj);
	  if (MathExtra::dot3(icrossj,iedge) > 0.0)
	    connect3d[i].aflag_e1[j] = CONCAVE;
	  else
	    connect3d[i].aflag_e1[j] = CONVEX;
	}
      } else {
	if (jpfirst == 2) connect3d[i].ewhich_e1[j] = 0;
	else if (jpfirst == 3) connect3d[i].ewhich_e1[j] = 1;
	else if (jpfirst == 1) connect3d[i].ewhich_e1[j] = 2;
	connect3d[i].nside_e1[j] = SAME_SIDE;
        if (dotnorm > 1.0-FLATTHRESH) connect3d[i].aflag_e1[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
          if (MathExtra::dot3(icrossj,iedge) < 0.0)
	    connect3d[i].aflag_e1[j] = CONCAVE;
          else
	    connect3d[i].aflag_e1[j] = CONVEX;
        }
      }
    }

    for (m = 0; m < connect3d[i].ne2; m++) {
      j = connect3d[i].neigh_e2[m];

      if (tris[i].p2 == tris[j].p1) jpfirst = 1;
      else if (tris[i].p2 == tris[j].p2) jpfirst = 2;
      else if (tris[i].p2 == tris[j].p3) jpfirst = 3;

      if (tris[i].p3 == tris[j].p1) jpsecond = 1;
      else if (tris[i].p3 == tris[j].p2) jpsecond = 2;
      else if (tris[i].p3 == tris[j].p3) jpsecond = 3;

      inorm = tris[i].norm;
      jnorm = tris[j].norm;
      dotnorm = MathExtra::dot3(inorm,jnorm);
      MathExtra::sub3(points[tris[i].p3].x,points[tris[i].p2].x,iedge);

      if ((jpfirst == 1 && jpsecond == 2) ||
	  (jpfirst == 2 && jpsecond == 3) ||
	  (jpfirst == 3 && jpsecond == 1)) {
	connect3d[i].ewhich_e2[j] = jpfirst - 1;
	connect3d[i].nside_e2[j] = OPPOSITE_SIDE;
	if (dotnorm < -1.0+FLATTHRESH) connect3d[i].aflag_e2[j] = FLAT;
	else {
	  MathExtra::cross3(inorm,jnorm,icrossj);
	  if (MathExtra::dot3(icrossj,iedge) > 0.0)
	    connect3d[i].aflag_e2[j] = CONCAVE;
	  else
	    connect3d[i].aflag_e2[j] = CONVEX;
	}
      } else {
	if (jpfirst == 2) connect3d[i].ewhich_e2[j] = 0;
	else if (jpfirst == 3) connect3d[i].ewhich_e2[j] = 1;
	else if (jpfirst == 1) connect3d[i].ewhich_e2[j] = 2;
	connect3d[i].nside_e2[j] = SAME_SIDE;
        if (dotnorm > 1.0-FLATTHRESH) connect3d[i].aflag_e2[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
          if (MathExtra::dot3(icrossj,iedge) < 0.0)
	    connect3d[i].aflag_e2[j] = CONCAVE;
          else
	    connect3d[i].aflag_e2[j] = CONVEX;
        }
      }
    }

    for (m = 0; m < connect3d[i].ne3; m++) {
      j = connect3d[i].neigh_e3[m];

      if (tris[i].p3 == tris[j].p1) jpfirst = 1;
      else if (tris[i].p3 == tris[j].p2) jpfirst = 2;
      else if (tris[i].p3 == tris[j].p3) jpfirst = 3;

      if (tris[i].p1 == tris[j].p1) jpsecond = 1;
      else if (tris[i].p1 == tris[j].p2) jpsecond = 2;
      else if (tris[i].p1 == tris[j].p3) jpsecond = 3;

      inorm = tris[i].norm;
      jnorm = tris[j].norm;
      dotnorm = MathExtra::dot3(inorm,jnorm);
      MathExtra::sub3(points[tris[i].p1].x,points[tris[i].p3].x,iedge);

      if ((jpfirst == 1 && jpsecond == 2) ||
	  (jpfirst == 2 && jpsecond == 3) ||
	  (jpfirst == 3 && jpsecond == 1)) {
	connect3d[i].ewhich_e3[j] = jpfirst - 1;
	connect3d[i].nside_e3[j] = OPPOSITE_SIDE;
	if (dotnorm < -1.0+FLATTHRESH) connect3d[i].aflag_e3[j] = FLAT;
	else {
	  MathExtra::cross3(inorm,jnorm,icrossj);
	  if (MathExtra::dot3(icrossj,iedge) > 0.0)
	    connect3d[i].aflag_e3[j] = CONCAVE;
	  else
	    connect3d[i].aflag_e3[j] = CONVEX;
	}
      } else {
	if (jpfirst == 2) connect3d[i].ewhich_e3[j] = 0;
	else if (jpfirst == 3) connect3d[i].ewhich_e3[j] = 1;
	else if (jpfirst == 1) connect3d[i].ewhich_e3[j] = 2;
	connect3d[i].nside_e3[j] = SAME_SIDE;
        if (dotnorm > 1.0-FLATTHRESH) connect3d[i].aflag_e3[j] = FLAT;
        else {
          MathExtra::cross3(inorm,jnorm,icrossj);
	  if (MathExtra::dot3(icrossj,iedge) < 0.0)
	    connect3d[i].aflag_e3[j] = CONCAVE;
          else
	    connect3d[i].aflag_e3[j] = CONVEX;
        }
      }
    }
  }

  // setup tri corner point connectivity lists
  // count # of tris containing each corner point (including self)
  // ctris =  ragged 2d array with indices of tris which contain each point

  memory->create(counts,npoints,"surface/global:counts");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    counts[tris[i].p1]++;
    counts[tris[i].p2]++;
    counts[tris[i].p3]++;
  }

  memory->create_ragged(ctris,npoints,counts,"surface/global:ctris");

  for (int i = 0; i < npoints; i++) counts[i] = 0;

  for (int i = 0; i < ntris; i++) {
    ctris[tris[i].p1][counts[tris[i].p1]++] = i;
    ctris[tris[i].p2][counts[tris[i].p2]++] = i;
    ctris[tris[i].p3][counts[tris[i].p3]++] = i;
  }

  // allocate all ragged corner arrays which Connect2d will point to
  // first subtract self from each point's count
  //   b/c connect2d vectors will NOT include self

  memory->create_ragged(neigh_c1,npoints,counts,"surface/global:neigh_c1");
  memory->create_ragged(neigh_c2,npoints,counts,"surface/global:neigh_c2");
  memory->create_ragged(neigh_c3,npoints,counts,"surface/global:neigh_c3");
  memory->create_ragged(cwhich_c1,npoints,counts,"surface/global:cwhich_c1");
  memory->create_ragged(cwhich_c2,npoints,counts,"surface/global:cwhich_c2");
  memory->create_ragged(cwhich_c3,npoints,counts,"surface/global:cwhich_c3");

  // set connect3d corner vector ptrs to rows of corresponding ragged arrays

  for (int i = 0; i < ntris; i++) {
    connect3d[i].nc1 = counts[tris[i].p1];
    if (connect3d[i].nc1 == 0) {
      connect3d[i].neigh_c1 = nullptr;
      connect3d[i].cwhich_c1 = nullptr;
    } else {
      connect3d[i].neigh_c1 = neigh_c1[tris[i].p1];
      connect3d[i].cwhich_c1 = cwhich_c1[tris[i].p1];
    }

    connect3d[i].nc2 = counts[tris[i].p2];
    if (connect3d[i].nc2 == 0) {
      connect3d[i].neigh_c2 = nullptr;
      connect3d[i].cwhich_c2 = nullptr;
    } else {
      connect3d[i].neigh_c2 = neigh_c2[tris[i].p2];
      connect3d[i].cwhich_c2 = cwhich_c2[tris[i].p2];
    }

    connect3d[i].nc3 = counts[tris[i].p3];
    if (connect3d[i].nc3 == 0) {
      connect3d[i].neigh_c3 = nullptr;
      connect3d[i].cwhich_c3 = nullptr;
    } else {
      connect3d[i].neigh_c3 = neigh_c3[tris[i].p3];
      connect3d[i].cwhich_c3 = cwhich_c3[tris[i].p3];
    }
  }

  // initialize connect3d corner neigh vectors for each corner point of each tri
  // do NOT include self

  for (int i = 0; i < ntris; i++) {
    if (connect3d[i].nc1) {
      j = 0;
      for (m = 0; m <= connect3d[i].nc1; m++) {
        if (ctris[tris[i].p1][m] == i) continue;
        connect3d[i].neigh_c1[j] = ctris[tris[i].p1][m];
        j++;
      }
    }
    if (connect3d[i].nc2) {
      j = 0;
      for (m = 0; m <= connect3d[i].nc2; m++) {
        if (ctris[tris[i].p2][m] == i) continue;
        connect3d[i].neigh_c2[j] = ctris[tris[i].p2][m];
        j++;
      }
    }
    if (connect3d[i].nc3) {
      j = 0;
      for (m = 0; m <= connect3d[i].nc3; m++) {
        if (ctris[tris[i].p3][m] == i) continue;
        connect3d[i].neigh_c3[j] = ctris[tris[i].p3][m];
        j++;
      }
    }
  }

  // deallocate counts and ctris

  memory->destroy(counts);
  memory->destroy(ctris);

  // set connect3d cwhich for each end point of each line
  // see fsg.h file for an explanation of each corner vector in Connect3d
  // aflag is based on dot and cross product of 2 connected line normals

  for (int i = 0; i < ntris; i++) {
    for (m = 0; m < connect3d[i].nc1; m++) {
      j = connect3d[i].neigh_c1[m];
      if (tris[i].p1 == tris[j].p1) connect3d[i].cwhich_c1[j] = 0;
      else if (tris[i].p1 == tris[j].p2) connect3d[i].cwhich_c1[j] = 1;
      else if (tris[i].p1 == tris[j].p3) connect3d[i].cwhich_c1[j] = 2;
    }
    for (m = 0; m < connect3d[i].nc2; m++) {
      j = connect3d[i].neigh_c2[m];
      if (tris[i].p2 == tris[j].p1) connect3d[i].cwhich_c2[j] = 0;
      else if (tris[i].p2 == tris[j].p2) connect3d[i].cwhich_c2[j] = 1;
      else if (tris[i].p2 == tris[j].p3) connect3d[i].cwhich_c2[j] = 2;
    }
    for (m = 0; m < connect3d[i].nc3; m++) {
      j = connect3d[i].neigh_c3[m];
      if (tris[i].p3 == tris[j].p1) connect3d[i].cwhich_c3[j] = 0;
      else if (tris[i].p3 == tris[j].p2) connect3d[i].cwhich_c3[j] = 1;
      else if (tris[i].p3 == tris[j].p3) connect3d[i].cwhich_c3[j] = 2;
    }
  }
}

/* ----------------------------------------------------------------------
   set attributes of all lines or tris
   xsurf,vsurf,omegasurf,norm
------------------------------------------------------------------------- */

void FixSurfaceGlobal::surface_attributes()
{
  double delta[3],p12[3],p13[3];
  double *p1,*p2,*p3;
  double zunit[3] = {0.0,0.0,1.0};

  memory->create(xsurf,nsurf,3,"surface/global:xsurf");
  memory->create(vsurf,nsurf,3,"surface/global:vsurf");
  memory->create(omegasurf,nsurf,3,"surface/global:omegasurf");
  memory->create(radsurf,nsurf,"surface/global:radsurf");

  if (dimension == 2) {
    for (int i = 0; i < nsurf; i++) {
      p1 = points[lines[i].p1].x;
      p2 = points[lines[i].p2].x;
      xsurf[i][0] = 0.5 * (p1[0]+p2[0]);
      xsurf[i][1] = 0.5 * (p1[1]+p2[1]);
      xsurf[i][2] = 0.0;

      MathExtra::sub3(p2,p1,p12);
      radsurf[i] = 0.5 * MathExtra::len3(p12);

      MathExtra::cross3(zunit,p12,lines[i].norm);
      MathExtra::norm3(lines[i].norm);
    }

  } else {

    for (int i = 0; i < nsurf; i++) {
      p1 = points[tris[i].p1].x;
      p2 = points[tris[i].p2].x;
      p3 = points[tris[i].p3].x;
      xsurf[i][0] = (p1[0]+p2[0]+p3[0]) / 3.0;
      xsurf[i][1] = (p1[1]+p2[1]+p3[1]) / 3.0;
      xsurf[i][2] = (p1[2]+p2[2]+p3[2]) / 3.0;

      MathExtra::sub3(p1,xsurf[i],delta);
      radsurf[i] = MathExtra::lensq3(delta);
      MathExtra::sub3(p2,xsurf[i],delta);
      radsurf[i] = MAX(radsurf[i],MathExtra::lensq3(delta));
      MathExtra::sub3(p3,xsurf[i],delta);
      radsurf[i] = MAX(radsurf[i],MathExtra::lensq3(delta));
      radsurf[i] = sqrt(radsurf[i]);

      MathExtra::sub3(p2,p1,p12);
      MathExtra::sub3(p3,p1,p13);
      MathExtra::cross3(p12,p13,tris[i].norm);
      MathExtra::norm3(tris[i].norm);
    }
  }

  for (int i = 0; i < nsurf; i++) {
    vsurf[i][0] = vsurf[i][1] = vsurf[i][2] = 0.0;
    omegasurf[i][0] = omegasurf[i][1] = omegasurf[i][2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   setup for surface motion
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_init()
{
  /*
  memory->create(points_lastneigh,npoints,3,"surface/global:points_lastneigh");
  memory->create(points_original,npoints,3,"surface/global:points_original");
  memory->create(xsurf_original,nsurf,3,"surface/global:xsurf_original");

  for (int i = 0; i < npoints; i++) {
    points_lastneigh[i][0] = points_original[i][0] = points[i].x[0];
    points_lastneigh[i][1] = points_original[i][1] = points[i].x[1];
    points_lastneigh[i][2] = points_original[i][2] = points[i].x[2];
  }

  for (int i = 0; i < nsurf; i++) {
    xsurf_original[i][0] = xsurf[i][0];
    xsurf_original[i][1] = xsurf[i][1];
    xsurf_original[i][2] = xsurf[i][2];
    omegasurf[i][0] = omega_rotate*runit[0];
    omegasurf[i][1] = omega_rotate*runit[1];
    omegasurf[i][2] = omega_rotate*runit[2];
  }
  */
}

/* ----------------------------------------------------------------------
   turn off surface motion and free memory
------------------------------------------------------------------------- */

void FixSurfaceGlobal::move_clear()
{
  /*
  // reset v,omega to zero

  for (int i = 0; i < nsurf; i++) {
    vsurf[i][0] = vsurf[i][1] = vsurf[i][2] = 0.0;
    omegasurf[i][0] = omegasurf[i][1] = omegasurf[i][2] = 0.0;
  }

  // deallocate memory

  memory->destroy(points_lastneigh);
  memory->destroy(points_original);
  memory->destroy(xsurf_original);
  points_lastneigh = nullptr;
  points_original = nullptr;
  xsurf_original = nullptr;
  */
}
