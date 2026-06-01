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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(dir/dumbbell,ComputeDirDumb);
// clang-format on
#else

#ifndef LMP_COMPUTE_DIR_FUMB_H
#define LMP_COMPUTE_DIR_DUMB_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeDirDumb : public Compute {
 public:
  ComputeDirDumb(class LAMMPS *, int, char **);
  ~ComputeDirDumb() override;
  void init() override;
  void compute_peratom() override;
  double memory_usage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:
  int nmax;
  int index_dumbdir, index_dumbpass;
  double *ke;
};

}    // namespace LAMMPS_NS

#endif
#endif
