#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "input.h"
#include "compute_reduce_local.h"
#include <cstring>
#include <cmath>

using namespace LAMMPS_NS;

ComputeReduceLocal::ComputeReduceLocal(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR, "Illegal compute reduce/local command");

  id_compute = utils::strdup(arg[3]);
  which_component = utils::inumeric(FLERR, arg[4], true, lmp);

  scalar_flag = 1;
  extscalar = 0;
}

void ComputeReduceLocal::init()
{
  index_compute = modify->find_compute(id_compute);
  if (index_compute < 0)
    error->all(FLERR, "Compute ID for reduce/local not found");

  if (!modify->compute[index_compute]->local_flag)
    error->all(FLERR, "Compute does not provide local data");

  //  if (!modify->compute[index_compute]->array_flag)
  //    error->all(FLERR, "Compute must provide a local array");
}

double ComputeReduceLocal::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  modify->compute[index_compute]->compute_local();   // update data

  double **ldata = modify->compute[index_compute]->array_local;
  int nrows = modify->compute[index_compute]->size_local_rows;

  double sum = 0.0;
  for (int i = 0; i < nrows; i++) {
    sum += ldata[which_component][i];
  }

  scalar = sum / nrows;
  return scalar;
}
