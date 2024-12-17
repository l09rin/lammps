/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(fene/alpha, BondFENEALPHA);
// clang-format on
#else

#ifndef LMP_BOND_FENEALPHA_H
#define LMP_BOND_FENEALPHA_H
#define MY_GAMMA 3.173072867831619

#include "bond.h"

namespace LAMMPS_NS {

class BondFENEALPHA : public Bond {
 public:
  BondFENEALPHA(class LAMMPS *lmp) : Bond(lmp) {}
  virtual ~BondFENEALPHA();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  virtual void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);
  virtual void *extract(const char *, int &);

 protected:
  double *k, *r0, *epsilon, *sigma, *alpha;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

W: FENEALPHA bond too long: %ld %d %d %g

A FENEALPHA bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

E: Bad FENEALPHA bond

Two atoms in a FENEALPHA bond have become so far apart that the bond cannot
be computed.

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

W: Use special bonds = 0,1,1 with bond style fene/alpha

Most FENEALPHA models need this setting for the special_bonds command.

W: FENEALPHA bond too long: %ld %g

A FENEALPHA bond has stretched dangerously far.  It's interaction strength
will be truncated to attempt to prevent the bond from blowing up.

*/
