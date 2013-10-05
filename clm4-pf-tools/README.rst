
CLM-PFLOTRAN Development Strategy
=================================

CLM-PFLOTRAN is making invasive changes to CLM in what is essentially
a permanent fork. An organized development strategey is necessary to
ensure long term maintenance and sustainability of CLM-PFLOTRAN and
its long term viability as a platform to enable new science.

Our plan is to adapt a trunk and science branch strategy similar to
CESM. In addition to the current clm4-pf “trunk”, we will be creating
a shared clm4-pf-ngee-sci repo for changes related to the ngee-arctic
project.

To make the ongoing changes to the code as easy and maintainable as
possible, we have put together the following list of suggested
practices.

CLM-PFLOTRAN trunk
------------------

1.0 - The main clm-pflotran repository, clm4-pf, should be a "pure"
implementation of clm-pflotran against clm-trunk. This version will
have only the minimal changes necessary to implement the clm-pflotran
interface and will always be based on clm trunk tags. All
non-interface changes, i.e. science, should go on a science branch.

1.1 - For testing and validation purposes, the main clm-pflotran
interface should be runtime configurable behavior. CLM-PFLOTRAN should
always compile with and without pflotran, and should produce the same
numerical results as the trunk tag if:

    * pflotran is not compiled into clm

    * pflotran is compiled into clm but turned off.

1.2 - CLM is working on removing all ifdefs from the code base, and we
are adopting that strategy and not allowing any new ifdefs in the
interface repo. The only exception to this is the CLM_PFLOTRAN ifdefs
inside clm_pfotran_interfaceMod.F90. These are required to maintain
the ability to compile without pflotran and ensure runtime
functionality. All new uses of ifdef CLM_PFLOTRAN should follow this
convention (see below).

1.3 - The clm-pflotran interface is adapting a strategy of "surgical
strikes", inserting the smallest possible amount of code into the clm
modules. An example of the preferred method is:

.. code-block:: fortran

    ! For example: Within Hydrology2() subroutine in Hydrology2Mod.F90
    if (use_pflotran) then
       call clm_pf_step_th(...) ! [This subroutine should be added in
                                ! clm_pflotran_interfaceMod.F90]
    else
       call SoilWater(...) ! standard clm function call
    end if

1.4 - All the major code development for the interface should go into
the clm-pflotran interface module.

.. code-block:: fortran

    ! Add clm_pf_step_th() in clm_pflotran_interaceMod.F90
    ! This is a dummy subroutine to check if user requested PFLOTRAN
    ! via usr_clm file, but the code was compiled without PFLOTRAN
    subroutine clm_pf_step_th()
    
    #ifdef CLM_PFLOTRAN
     call step_th_clm_pf( … ) ! worker subroutine
    #else
     call pflotran_not_available(subname) ! Error checking
    #endif
    
    end subroutine clm_pf_step_th()
    
    ! Add step_th_clm_pf() in clm_pflotran_interaceMod.F90
    ! This subroutine does the work of interfacing with pflotran
    subroutine step_th_clm_pf()
     ! Pass data from CLM to PFLOTRAN; Step PFLOTRAN;
     ! Get update states/fluxes from PFLOTRAN, etc
    end subroutine step_th_clm_pf

1.5 - The primary clm-pflotran repo will continue to have the
gatekeepers doing code reviews. Code acceptance requirements for this
repo will get stricter.

1.6 - An automated build and testing platform will be used to ensure
that changes are not breaking the expected functionality.


CLM-PFLOTRAN sciences branches
------------------------------

2.0 - All science work based on clm-pflotran should be done on a
science branch. Depending on the nature of the changes, this should
either be:

    2.0.1 - An NCAR svn science branch that gets pulled into a clone
    of clm-pflotran. We will provide tools to make this easier.
    
    2.0.2 - A separate clone of clm4-pf dedicated to the science branch.

2.1 - Science branches will be more open. Code maintenance on the
science branches will be left up to the scientists who are working on
the branches.

2.2 - We encourage scientists to adapt good coding practices like
properly integrating your code instead of adding ifdefs, not making
lots of invasive changes that will be hard to merge, making changes
directly to the source rather than copying files into source mods,
etc. But this is left up to the discretion of the maintainers of the
science branch.

2.3 - Science branches can be shared among collaborators according to
the sharing and NCAR-CESM guidelines
(http://www.cesm.ucar.edu/working_groups/Software/secp/repo/).

2.4 - Science changes must go through peer review and be accepted into
the CLM trunk before being merged into clm-pflotran trunk.

2.5 - We will provide assistance in setting up an automated testing
system for long lived science branches.

PFLOTRAN
--------

3.0 - PFLOTRAN is considered an external library for clm-pflotran, the
same way BLAS provides linear algebra routines or petsc provides mesh
and solver functionality to pflotran.

3.1 - The goal to provide a generic interface into pflotran for
external drivers in pflotran-dev, doing away with the need for a
pflotran_coupled repository. To couple directly to pflotran-dev, clm
specific code for things like data mapping between meshes should go
into the clm side of the code as it is unlikely to be benefit any
other project using pflotran-dev.

3.2 - Once clm-pflotran trunk migrates to building on pflotran-dev,
all changes to pflotran will have to go through the standard
acceptance procedures for pflotran, including pflotran’s testing and
coding standards.

3.3 - We will document a “supported” version of pflotran-dev and
petsc-dev that will work with clm-pflotran.

3.4 - Science branches can fork from the supported version of
pflotran-dev to add experimental science to the pflotran side of the
interface.

