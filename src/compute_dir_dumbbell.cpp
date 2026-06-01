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

#include "compute_dir_dumbbell.h"

#include <cstring>
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"
#include "neighbor.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeDirDumb::ComputeDirDumb(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), ke(nullptr)
{
  if (narg != 3) error->all(FLERR, "Illegal compute dir/dumb command");

  peratom_flag = 1;
  size_peratom_cols = 0;
  comm_reverse = 2;
  comm_forward = 2;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeDirDumb::~ComputeDirDumb()
{
  memory->destroy(ke);
}

/* ---------------------------------------------------------------------- */

void ComputeDirDumb::init()
{
  if (modify->get_compute_by_style(style).size() > 1)
    if (comm->me == 0) error->warning(FLERR, "More than one compute {}", style);

  //LINO - define external variable for memorizing thetaold and differences in theta
  int flag, ncols;

  int index1 = atom->find_custom("dumbdir",flag,ncols);
  int index2 = atom->find_custom("dumbpass",flag,ncols);
  if (index1 < 0)
        error->all(FLERR,"p property/atom floating point "
                   "vector does not exist");
  if (index2 < 0)
        error->all(FLERR,"p property/atom floating point "
                   "vector does not exist");

  index_dumbdir = index1;
  index_dumbpass = index2;
  double *dumbdir = atom->dvector[index1];
  double *dumbpass = atom->dvector[index2];
  int nlocal = atom->nlocal;

  for (int n = 0; n < nlocal; n++) {
    dumbdir[n] = 0.0;
    dumbpass[n] = 0;
  }
  //LINO - end
}

/* ---------------------------------------------------------------------- */

void ComputeDirDumb::compute_peratom()
{
  double *dumbdir = atom->dvector[index_dumbdir];
  double *dumbpass = atom->dvector[index_dumbpass];

  invoked_peratom = update->ntimestep;

  // grow ke array if necessary

  if (atom->nmax > nmax) {
    memory->destroy(ke);
    nmax = atom->nmax;
    memory->create(ke, nmax, "ke/atom:ke");
    vector_atom = ke;
  }

  // compute kinetic energy for each atom in group

  double theta, thetaold, thetadiff, itheta;
  int i1, i2;
  double delx, dely, delz;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int *tag = atom->tag;

  for (int i = nlocal; i < nlocal + nghost; i++) {
    dumbdir[i] = 0;
    dumbpass[i] = 0;
  }
 
  for (int n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    if (tag[i1] > tag[i2]) {
      i1 = i2;
      i2 = bondlist[n][0];
    }

    delx = x[i2][0] - x[i1][0];
    dely = x[i2][1] - x[i1][1];
    delz = x[i2][2] - x[i1][2];

    // Update dipole direction and check if it passes 2PI
    domain->minimum_image(delx,dely,delz);
    theta = atan2(dely,delx);
    itheta = 0;

    // if the atom is owned, not a ghost atom and belonging to the group
    if (mask[i1] & groupbit && mask[i2] & groupbit) {
      if (i1 < nlocal) {
	thetaold = dumbdir[i1];
	thetadiff = theta - thetaold;
	if (thetadiff<-MY_PI) {
	  itheta = 1; // qui  passata in senso antiorario
	  thetadiff = MY_PI + thetadiff;
	} else if (thetadiff>MY_PI) {
	  itheta = -1;  // qui  passata in senso orario
	  thetadiff = thetadiff - MY_PI;
	}
	dumbdir[i1] = theta;
	dumbpass[i1] += itheta;

	if (i2 < nlocal) {
	  dumbdir[i2] = dumbdir[i1];
	  dumbpass[i2] = dumbpass[i1];
	} else {
	  dumbdir[i2] = thetadiff;
	  dumbpass[i2] = itheta;
	}
      } else if (i2 < nlocal) {
	thetaold = dumbdir[i2];
	thetadiff = theta - thetaold;
	if (thetadiff<-MY_PI) {
	  itheta = 1; // qui  passata in senso antiorario
	  thetadiff = MY_PI + thetadiff;
	} else if (thetadiff>MY_PI) {
	  itheta = -1;  // qui  passata in senso orario
	  thetadiff = thetadiff - MY_PI;
	}
	dumbdir[i2] = theta;
	dumbpass[i2] += itheta;

	if (i1 < nlocal) {
	  dumbdir[i1] = dumbdir[i2];
	  dumbpass[i1] = dumbpass[i2];
	} else {
	  dumbdir[i1] = thetadiff;
	  dumbpass[i1] = itheta;
	}
      }
    }
  }
  // To safely update the values relative to ghost atoms
  comm->reverse_comm(this);
  //comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int ComputeDirDumb::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = atom->dvector[index_dumbdir][i];
    buf[m++] = static_cast<double>(atom->dvector[index_dumbpass][i]);
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeDirDumb::unpack_reverse_comm(int n, int *list, double *buf) {
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    atom->dvector[index_dumbdir][j] += buf[m++];
    atom->dvector[index_dumbpass][j] += static_cast<int>(buf[m++]);
  }
}

/* ---------------------------------------------------------------------- */

int ComputeDirDumb::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, m = 0, j;
  for (i = 0; i < n; ++i) {
    j = list[i];
    buf[m++] = atom->dvector[index_dumbdir][j];
    buf[m++] = static_cast<double>(atom->dvector[index_dumbpass][j]);
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeDirDumb::unpack_forward_comm(int n, int first, double *buf)
{
  int i, last, m = 0;
  last = first + n;
  for (i = first; i < last; ++i) {
    atom->dvector[index_dumbdir][i] += buf[m++];
    atom->dvector[index_dumbpass][i] += static_cast<int>(buf[m++]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeDirDumb::memory_usage()
{
  double bytes = (double) nmax * sizeof(double);
  return bytes;
}
