.. index:: fix surface/global

fix surface/global command
===============

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID surface/global input args input args ... model args
   model args ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* surface/global = style name of this fix command
* input = one or more input keywords can be specified

  .. parsed-literal::

       *input* args = source source-args
         *source* = *mol* or *stl*
	    *mol* arg = template-ID
               template-ID = ID of molecule template specified in a
	       separate :doc:`molecule <molecule>` command, which defines a set of triangles or lines
            *stl* args = stype stlfile
	       stype = numeric type assigned to all triangles in STL file
	       stlfile = STL filename which defines a set of triangles

* model = one or more model keywords can be specified

  .. parsed-literal::
     
        *model* args = ptype stype fstyle fstyle-args
           ptype = numeric particle type (see asterisk form below)
           stype = numeric surface type (see asterisk form below)
           *fstyle* =  *hooke* or *hooke/history* or *hertz/history* or *granular*
             *hooke* or *hooke/history* or *hertz/history* args = Kn Kt gamma_n gamma_t xmu dampflag (limit_damping)
                Kn = elastic constant for normal particle repulsion (force/distance units or pressure units - see discussion below)
                Kt = elastic constant for tangential contact (force/distance units or pressure units - see discussion below)
                gamma_n = damping coefficient for collisions in normal direction (1/time units or 1/time-distance units - see discussion below)
                gamma_t = damping coefficient for collisions in tangential direction (1/time units or 1/time-distance units - see discussion below)
                xmu = static yield criterion (unitless value between 0.0 and 1.0e4)
                dampflag = 0 or 1 if tangential damping force is excluded or included
                limit_damping = optional argument to prevent attractive interactions
             *granular* args = same syntax for args as *pair_coeff* command of :doc:`pair_style granular <pair_granular>`

* zero or more keyword/value pairs may be appended
* keyword = *smax* or *flat* or *temperature*

  .. parsed-literal::

       *smaxtype* value = maximum surface type allowed
       *flat* value = maximum angle (degrees) between a pair of connected triangles/lines for a flat connection
       *temperature* value = surface temperature (degrees Kelvin), required if model with heat is used
	     
Examples
""""""""

.. code-block:: LAMMPS

   molecule lines surf.line
   fix 1 all surface/global input mol linesurf model 1 1 hooke 4000.0 NULL 100.0 NULL 0.5 1

   fix 1 all surface/global input stl 1 object.stl model * 1 hooke/history 4000.0 NULL 100.0 NULL 0.5 1

Description
"""""""""""

Enable granular surfaces to be used as boundary conditions on
particles in a granular simulation.  Granular surfaces are defined as
a set of triangles (3d) or a set of line segments (2d).

The :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
gives an overview of granular surfaces of two types, *global* and
*local*, with guidelines for how to use them.

This command is used for models with *global* surfaces.  The :doc:`fix
surface/local <fix_surface_local>` command is used for models with
*local* surfaces.  As explained on the :doc:`Howto granular surfaces
<Howto_granular_surfaces>` doc page, *global* surfaces are most
appropriate when there is a modest number of them.  Each surface
(triangle/line) can be of any size, even as large as a dimension of
the simulation box.  A copy of the list of *global* surfaces is stored
by each processor.

*Global* surfacescan be defined in 2 ways, which correspond to the 2
options listed above for the *source* argument of the *input* keyword:

* via a molecule file(s), read by the :doc:`molecule <molecule>` command
* via an STL file, read by this command

If triangles or lines were previously read in by the :doc:`molecule
<molecule>` command, then the *source* keyword is *mol* and its
argument is the molecule template ID used with the :doc:`molecule
<molecule>` command.  Note that a doc:`molecule <molecule>` command
can read and assign serveral molecule files to the same template-ID.
Each molecule file must define triangles or lines, not atoms.  For
multiple molecule files, the set of surfaces used by this command will
be the union of the triangles and lines from all the molecule files.
Note that each line/triangle in a molecule file is assigned a type and
molecule-ID.

An STL (stereolithography) file defines a set of triangles.  For use
with this command, the *source* argument of the *input* keyword is
*stl*.  The *stype* argument is the numeric type assigned to all the
triangles from the file.  Note that STL files do not contain types or
other flags for each triangle.  The *stlfile* argument is the name of
the STL file.  It can be in text or binary format; this command
auto-detects the format.  Note that STL files cannot be used for 2d
simulations since they only define triangles.

This `Wikepedia page
<https://en.wikipedia.org/wiki/STL_(file_format)>`_ describes the
format of both text and binary STL files.  Binary STL files can be
converted to ASCII for editing with the stl_bin2txt tool in the
lammps/tools directory.  Examples of text-based STL files are included
in the examples/gransurf directory.

Note that this command allows for multiple uses of the *input*
keyword, each with a *source* argument as either *mol* or *stl*.  The
surfaces used by this command are the union of the triangles and lines
from all the input files.

Once surfaces are defined, this command calculates their connectivity.
Two triangles are "connected" if they have a single corner point in
common or an edge in common (2 corner points).  Two line segments are
"connected" if the they have an end point in common.  More technical
details on connectivity and its significance for granular surface
simulations is given on :doc:`Howto granular surfaces
<Howto_granular_surfaces>` doc page.  In brief, a pair of connected
surfaces interact with a particle which contacts both of them
simultaneously according to a set of rules which are designed to
generate physically sensible forces on the particle.

Note that there is no requirement that all the surfaces be connected
to one another.  The surfaces can represent the surface of one or more
independent objects.  Particles in the specified group-ID interact
with the surface when they are close enough to overlap (touch) one or
more individual triangles or lines.  Both sides of a triangle or line
interact with particles.  Thus a surface can be infinitely thin,
e.g. the blade of a mixer.  See the :doc:`Howto granular surfaces
<Howto_granular_surfaces>` doc page for restrictions on the geometry
of a collection of triangles or lines.

The nature of individual surface/particle interactions are determined
by the *model* keyword.  Each use of the model keyword is applied to
one or more particle types interacting with one or more surface types.
The *ptype* argument is the pariticle type, *stype* is the surface
type.

Either *ptype* and *stype* can be specified as a single numeric value.
Or a wildcard asterisk can be used to specify multiple particle or
surface types.  This takes the form "\*" or "\*n" or "n\*" or "m\*n".
If :math:`N` is the number of particle or surface types, then an
asterisk with no numeric values means all types from 1 to :math:`N`.
A leading asterisk means all types from 1 to n (inclusive).  A
trailing asterisk means all types from n to :math:`N` (inclusive).  A
middle asterisk means all types from m to n (inclusive).

The model keywords must specify an interactions for each particle type
interacting with each surface type, otherwise an error is flagged.
If use of the model keywords specifies an individual particle/surface
type pair more than once, then the final specification is used.

The number of particle types is the number of atom types in the
system.  The number of surface types is determined by the maximum
surface type in any of files read by the *input* keyword(s) or by the
optional *smaxtype* keyword.  The latter can be useful if the
:doc:`fix_modify type/region <fix_modify>` command (described below)
is used to assign new types to surfaces after they are read in.  As
for particles, there is no requirement that triangles/lines exist for
every surface type.

The *fstyle* argument (for force style) can be any of the styles
defined by the :doc:`pair_style gran/\* <pair_gran>` or the more
general :doc:`pair_style granular <pair_granular>` commands.
Currently the options are *hooke*, *hooke/history*, or *hertz/history*
for the former, and *granular* with all the possible options of the
associated *pair_coeff* command for the latter.  The equation for the
force between a triangle/line and a particle touching it is the same
as the corresponding equation on the :doc:`pair_style gran/\*
<pair_gran>` and :doc:`pair_style granular <pair_granular>` doc pages,
in the limit of one of the two particles going to infinite radius and
mass (flat surface).  Specifically, delta = radius - r = overlap of
particle with triangle/line, m_eff = mass of particle, and the
effective radius of contact = RiRj/Ri+Rj is set to the radius of the
particle.

The parameters *Kn*, *Kt*, *gamma_n*, *gamma_t*, *xmu*, *dampflag*,
and the optional keyword *limit_damping* have the same meaning and
units as those specified with the :doc:`pair_style gran/\*
<pair_gran>` commands.  This means a NULL can be used for either *Kt*
or *gamma_t* as described on that page.  If a NULL is used for *Kt*,
then a default value is used where *Kt* = 2/7 *Kn*\ .  If a NULL is
used for *gamma_t*, then a default value is used where *gamma_t* = 1/2
*gamma_n*.

Note that the fix surface/global command can be used multiple times
though it is not typically necessary to do so.  Note that if it is
used multiple times, the surfaces defined by the different commands
will NOT be "connected" to each other in the manner described above or
on the :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc
page.

----------

These are the optional keywords and values.

The *smaxtype* keyword sets the number of surface types which can be
used.  By default, this is the maximum type of any surface defined by
the *input* keyword(s).  If the :doc:`fix_modify type/region
<fix_modify>` command (described below) will be used to change a
surface type to a larger value than the default, then the *smax*
keyword can allow this.

The *flat* keyword can be used to set a threshold for the angle (in
degrees) between two connected surfaces (triangles or line segments)
which will be treated as "flat" by the particle/surface interactions
models.  A flat connection means a single force will be applied to the
particle even if it is contact with both surfaces simultaneously.  See
the :doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page
for more details.  The default for the *flat* keword is one degree.

The *temperature* keyword is required if any of the granular models
used includes a heat model which depends on the surface temperature.
Otherwise it is ignored.  The value is the temperature of the surface
in degrees Kelvin.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

This fix defines two new keywords for the doc:`fix_modify
<fix_modify>` command, *move* and *type/region*.  Because they are
specific to this command, they are only described here, not on the
doc:`fix_modify <fix_modify>` doc page.  Both keywords can be used
multiple times.  In the description that follows, a surface means a
triangle (3d) or line segment (2d).

The *move* keyword can be used to make all or a subset of the surfaces
move in a prescribed manner, similar to the :doc:`fix move <fix_move>`
command.  The *type/region* keyword can be used to change the types of
surfaces which are within a geometric region.  Their syntax is as follows:

.. code-block:: LAMMPS

   fix_modify fix-ID keyword values ...

* fix-ID = ID of the fix to modify
* keyword (specific to this fix) = *move* or *type/region*

  .. parsed-literal::

       *move* values = stype mstyle args
          stype =  numeric surface type(s) (see comma-separated, asterisk form below)
          mstyle = *none* or *linear* or *wiggle* or *rotate* or *transrot* or *variable*
           *none* args = none
           *linear* args = Vx Vy Vz
           *wiggle* args = Ax Ay Az period
           *rotate* args = Px Py Pz Rx Ry Rz period
           *transrot* args = Vx Vy Vz Px Py Pz Rx Ry Rz period
           *variable* args = v_dx v_dy v_dz v_vx v_vy v_vz
       *type/region* values = stype region-ID
         stype = numeric surface type
	 region-ID = ID of a region previously defined by the :doc:`region <region>` command

The *stype* argument can specify one or more surface types.  It must
specify all the surface types within a connecte object(s).  If an
object is composed of surfaces of 2 or more types, it is an error to
use the *move* keyword and not specify all those types, since this
would break the connections.  Note that LAMMPS does NOT check that
this requirement is met.  It is likewise an error to use the *move*
keyword multiple times to induce motion which overlaps surfaces in
ways that violate the surface geometry restrictions explained on the
:doc:`Howto granular surfaces <Howto_granular_surfaces>` doc page.
Again, LAMMPS does NOT check that this requirement is met.

The general format of *stype* is st,st,...,st where one or more *st*
sub-arguments are separated by commas.  A single *st* sub-argument is
either a single numeric value or contains a wildcard asterisk.  The
asterisk is used in place of or in conjunction with numeric arguments
to specify multiple type values.  This takes the form "\*" or "\*n" or
"n\*" or "m\*n".  If :math:`N` is the number of atom types, then an
asterisk with no numeric values means all types from 1 to :math:`N`.
A leading asterisk means all types from 1 to n (inclusive).  A
trailing asterisk means all types from n to :math:`N` (inclusive).  A
middle asterisk means all types from m to n (inclusive).

The *mstyle* argument is one of the listed styles above.  The *none*
style turns off motion which was previously enabled, e.g. stops the
rotation of an object.  Again, the list of surface types must include
all the surfaces in a connected object.  The other styles and their
effects on motion are the same as those defined by the :doc:`fix move
<fix_move>` command.  Their arguments are also the same as those
documented by the :doc:`fix move <fix_move>` command.

Note that for *local* surfaces the same motion operations can be
performed using the :doc:`fix move <fix_move>` command with a group-ID
defined by the :doc:`group <group>` which includes the appropriate
particle types for triangle and line-segment particles.

The *type/region* keyword can be used to re-assign surface types to
surfaces after they have been initialized by the *input* keyword.
This is most useful for STL triangles since STL files do not allow for
assignment of types to individual triangles.

The *stype* argument is a single numeric value, which must be between
1 and maxtype inclusive.  Maxtype is either the maximum type of all
surfaces read in by the *input* keyword or the setting of the optional
*smax* keyword.  The *region-ID* is the ID of a geometric region
defined by the :doc:`region <region>` command.  Note that regions can
be the inside or outside of a geometric object, such as a sphere or
block.  The geometric center point of a triangle or line segment is
used to determine where a surface is in the region or not.  If it is,
its type is reset to *stype*.

Examples for both keywords are as follows:

.. code-block:: LAMMPS

   fix_modify 1 move 2 rotate 0 0 0 0 0 1 25
   fix_modify 1 move 1,3*5,8* rotate 0 0 0 0 0 1 25
   fix_modify 1 type/region 3 myBlock
   
No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.  No parameter of this
fix can be used with the *start/stop* keywords of the :doc:`run <run>`
command.  This fix is not invoked during :doc:`energy minimization
<minimize>`.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""

:doc:`fix surface/local <fix_surface_local>`

Default
"""""""

The keyword defaults are smax = max type of all surfaces defined by
the input keyword(s), flat = one degree, temperaure = none.
