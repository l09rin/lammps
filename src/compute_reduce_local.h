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
ComputeStyle(reduce/local,ComputeReduceLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_REDUCE_LOCAL_H
#define LMP_COMPUTE_REDUCE_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeReduceLocal : public Compute {
 public:
  ComputeReduceLocal(class LAMMPS *, int, char **);
  ~ComputeReduceLocal() {}
  void init() override;
  double compute_scalar() override;

 private:
  char *id_compute;
  int index_compute;
  int which_component;
};

}    // namespace LAMMPS_NS

#endif
#endif
