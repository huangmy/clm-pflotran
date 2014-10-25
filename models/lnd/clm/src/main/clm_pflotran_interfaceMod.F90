module clm_pflotran_interfaceMod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: clm_pflotran_interfaceMod
  !
  ! !DESCRIPTION:
  ! Performs
  !
  ! !USES:

#ifdef CLM_PFLOTRAN
  use clm_pflotran_interface_data
  use pflotran_model_module
#endif

  use shr_log_mod    , only : errMsg => shr_log_errMsg

  ! !PUBLIC TYPES:
  implicit none


  save

  private    ! By default everything is private

  !type(clm_pflotran_interface_data_type),pointer,public   :: clm_pf_idata
#ifdef CLM_PFLOTRAN
  type(pflotran_model_type), pointer, public              :: pflotran_m
#endif
  !
  character(len=256), private :: pflotran_prefix = ''
  character(len=32), private :: restart_stamp = ''

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_pf_readnl

  ! wrappers around ifdef statements to maintain sane runtime behavior
  ! when pflotran is not available.
  public :: clm_pf_interface_init, &
       clm_pf_set_sflow_forcing, &
       clm_pf_update_soil_moisture, &
       clm_pf_update_soil_temperature, &
       clm_pf_update_soil_therm_cond, &
       clm_pf_update_drainage, &
       clm_pf_update_h2osfc, &
       clm_pf_step_th, &
       clm_pf_write_restart, clm_pf_set_restart_stamp, &
       clm_pf_vecget_gflux, clm_pf_vecrestore_gflux

#ifdef CLM_PFLOTRAN
  ! private work functions that truely require ifdef CLM_PFLOTRAN
  private :: interface_init_clm_pf, & ! Phase one initialization
       set_sflow_forcing_clm_pf, &
       update_soil_moisture_clm_pf, &
       update_soil_temperature_clm_pf, &
       step_th_clm_pf, &
       write_restart_clm_pf, &
       vecget_gflux_subsurf_clm_pf, vecrestore_gflux_subsurf_clm_pf
#endif

contains

!-----------------------------------------------------------------------
!
! public interface functions allowing runtime behavior regardless of
! whether pflotran is compiled in.
!
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_readnl
  !
  ! !INTERFACE:
  subroutine clm_pf_readnl( NLFilename )
  !
  ! !DESCRIPTION:
  ! Read namelist for clm-pflotran interface
  !
  ! !USES:
    use clm_varctl    , only : iulog
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use clm_nlUtilsMod, only : find_nlgroup_name
    use shr_nl_mod    , only : shr_nl_find_group_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun

    implicit none

  ! !ARGUMENTS:
    character(len=*), intent(IN) :: NLFilename ! Namelist filename
  ! !LOCAL VARIABLES:
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'clm_pf_readnl'  ! subroutine name
  !EOP
  !-----------------------------------------------------------------------
    namelist / clm_pflotran_inparm / pflotran_prefix

    ! ----------------------------------------------------------------------
    ! Read namelist from standard namelist file.
    ! ----------------------------------------------------------------------

    if ( masterproc )then

       unitn = getavu()
       write(iulog,*) 'Read in clm-pflotran namelist'
       call opnfil (NLFilename, unitn, 'F')
       call shr_nl_find_group_name(unitn, 'clm_pflotran_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, clm_pflotran_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading clm_pflotran_inparm namelist')
          end if
       end if
       call relavu( unitn )
       write(iulog, '(/, A)') " clm-pflotran namelist parameters :"
       write(iulog, '(A, " : ", A,/)') "   pflotran_prefix", trim(pflotran_prefix)
    end if

    ! Broadcast namelist variables read in
    call shr_mpi_bcast(pflotran_prefix, mpicom)
  end subroutine clm_pf_readnl



  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: clm_pf_set_restart_stamp
  !
  ! !INTERFACE:
  subroutine clm_pf_set_restart_stamp(clm_restart_filename)
  !
  ! !DESCRIPTION: Set the pflotran restart date stamp. Note we do NOT
  ! restart here, that gets handled by pflotran's internal
  ! initialization during interface_init_clm_pf()
  !
  ! !USES:
  ! !ARGUMENTS:
    character(len=256), intent(in) :: clm_restart_filename
  ! !LOCAL VARIABLES:
    integer :: name_length, start_pos, end_pos
    character(len=32) :: clm_stamp
  !EOP
  !-----------------------------------------------------------------------

    ! clm restart file name is of the form:
    !     ${CASE_NAME}.clm2.r.YYYY-MM-DD-SSSSS.nc
    ! we need to extract the: YYYY-MM-DD-SSSSS
    write(*, '("clm-pf : clm restart file name : ", A/)') trim(clm_restart_filename)
    name_length = len(trim(clm_restart_filename))
    start_pos = name_length - 18
    end_pos = name_length - 3
    clm_stamp = clm_restart_filename(start_pos : end_pos)
    write(*, '("clm-pf : clm date stamp : ", A/)') trim(clm_stamp)
    restart_stamp = clm_stamp

  end subroutine clm_pf_set_restart_stamp


  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: pflotran_not_available
  !
  ! !INTERFACE:
  subroutine pflotran_not_available(subname)
  !
  ! !DESCRIPTION:
  ! Print an error message and abort.
  !
  ! !USES:
    use abortutils    , only : endrun
  ! !ARGUMENTS:
    implicit none
    character(len=*), intent(in) :: subname
  ! !LOCAL VARIABLES:
  !EOP
  !-----------------------------------------------------------------------
    call endrun(trim(subname) // ": ERROR: CLM-PFLOTRAN interface has not been compiled " // &
         "into this version of CLM.")
  end subroutine pflotran_not_available


!-----------------------------------------------------------------------
!
! public interface function wrappers
!
!-----------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  subroutine clm_pf_interface_init(bounds, col, lun, grc, &
       soilstate_vars, waterstate_vars, temperature_vars, &
       soilhydrology_vars)

    use decompMod, only : bounds_type
    use GridcellType    , only : gridcell_type
    use ColumnType      , only : column_type
    use LandUnitType    , only : landunit_type
    use SoilStateType   , only : soilstate_type
    use WaterStateType  , only : waterstate_type
    use TemperatureType , only : temperature_type
    use SoilHydrologyType , only : soilhydrology_type

    implicit none

    type(bounds_type), intent(in) :: bounds  ! bounds
    type(gridcell_type), intent(in) :: grc
    type(column_type), intent(in) :: col
    type(landunit_type), intent(in) :: lun
    type(soilstate_type), intent(in) :: soilstate_vars
    type(waterstate_type), intent(in) :: waterstate_vars
    type(temperature_type), intent(in) :: temperature_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars

#ifdef CLM_PFLOTRAN
    call interface_init_clm_pf(bounds, col, lun, grc, &
       soilstate_vars, waterstate_vars, temperature_vars, &
       soilhydrology_vars)
#else
    character(len=256) :: subname = "clm_pf_interface_init()"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_interface_init


  !-----------------------------------------------------------------------------
  subroutine clm_pf_set_sflow_forcing(num_hydrologyc, filter_hydrologyc, &
       atm2lnd_vars, waterflux_vars)

    use filterMod, only : clumpfilter
    use atm2lndType, only : atm2lnd_type
    use WaterfluxType, only : waterflux_type

    implicit none

    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(waterflux_type), intent(in) :: waterflux_vars

#ifdef CLM_PFLOTRAN
    call set_sflow_forcing_clm_pf(num_hydrologyc, filter_hydrologyc, &
         atm2lnd_vars, waterflux_vars)
#else
    character(len=256) :: subname = "clm_pf_set_sflow_forcing"
    call pflotran_not_available(subname)
#endif

  end subroutine clm_pf_set_sflow_forcing


  !-----------------------------------------------------------------------------
  subroutine clm_pf_update_soil_moisture(waterstate_vars, soilstate_vars, bounds, &
       num_hydrologyc, filter_hydrologyc, col, soilhydrology_vars)

    use WaterStateType, only : waterstate_type
    use SoilStateType, only : soilstate_type
    use filterMod, only : clumpfilter
    use decompMod, only : bounds_type
    use ColumnType, only : column_type
    use SoilHydrologyType, only : soilhydrology_type

    implicit none

    type(waterstate_type), intent(inout) :: waterstate_vars
    type(soilstate_type), intent(inout) :: soilstate_vars
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(column_type), intent(in) :: col
    type(soilhydrology_type), intent(out) :: soilhydrology_vars

#ifdef CLM_PFLOTRAN
    call update_soil_moisture_clm_pf(waterstate_vars, soilstate_vars, bounds, num_hydrologyc, filter_hydrologyc, col, soilhydrology_vars)
#else
    character(len=256) :: subname = "clm_pf_update_soil_moisture"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_update_soil_moisture


  !-----------------------------------------------------------------------------
  subroutine clm_pf_update_soil_temperature(bounds, &
       num_nolakec, filter_nolakec, temperature_vars, col)

    use decompMod, only : bounds_type
    use TemperatureType, only : temperature_type
    use ColumnType, only : column_type

    implicit none

    type(bounds_type), intent(in) :: bounds  ! bounds
    integer , intent(in)  :: num_nolakec         ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)   ! column filter for non-lake points
    type(temperature_type), intent(in) :: temperature_vars
    type(column_type), intent(in) :: col

#ifdef CLM_PFLOTRAN
    call update_soil_temperature_clm_pf(bounds, &
       num_nolakec, filter_nolakec, temperature_vars, col)
#else
    character(len=256) :: subname = "clm_pf_update_soil_temperature"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_update_soil_temperature

  !-----------------------------------------------------------------------------
  subroutine clm_pf_update_soil_therm_cond(bounds, &
       num_urbanl, filter_urbanl, &
       num_nolakec, filter_nolakec, tk)

    use shr_kind_mod  , only : r8 => shr_kind_r8
    use filterMod, only : clumpfilter
    use decompMod, only : bounds_type
    use clm_varpar, only : nlevsno, nlevgrnd

    implicit none

    type(bounds_type), intent(in) :: bounds  ! bounds
    integer , intent(in)  :: num_nolakec         ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)   ! column filter for non-lake points
    integer , intent(in)  :: num_urbanl          ! number of urban landunits in clump
    integer , intent(in)  :: filter_urbanl(:)    ! urban landunit filter
    real(r8), intent(out) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)          ! thermal conductivity [W/(m K)]

    character(len=256) :: subname = "clm_pf_update_soil_therm_cond"

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

#ifdef CLM_PFLOTRAN
    call update_soil_therm_cond_clm_pf(bounds, &
       num_urbanl, filter_urbanl, &
       num_nolakec, filter_nolakec, &
       tk)
#else
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_update_soil_therm_cond


  !-----------------------------------------------------------------------------
  subroutine clm_pf_update_drainage(num_hydrologyc, filter_hydrologyc, waterflux_vars)

    use decompMod, only : bounds_type
    use WaterFluxType, only : waterflux_type

    implicit none

    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterflux_type), intent(out) :: waterflux_vars

#ifdef CLM_PFLOTRAN
    call update_drainage_clm_pf(num_hydrologyc, filter_hydrologyc, waterflux_vars)
#else
    character(len=256) :: subname = "clm_pf_update_drainage"
    call pflotran_not_available(subname)
#endif

  end subroutine clm_pf_update_drainage

  !-----------------------------------------------------------------------------
  subroutine clm_pf_update_h2osfc(num_hydrologyc, filter_hydrologyc, &
       waterstate_vars)

    use WaterStateType, only : waterstate_type

    implicit none

    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterstate_type), intent(out) :: waterstate_vars

#ifdef CLM_PFLOTRAN
    call update_h2osfc_clm_pf(num_hydrologyc, filter_hydrologyc, waterstate_vars)
#else
    character(len=256) :: subname = "clm_pf_update_h2osfc"
    call pflotran_not_available(subname)
#endif

  end subroutine clm_pf_update_h2osfc


  !-----------------------------------------------------------------------------
  subroutine clm_pf_step_th(bounds, &
       num_hydrologyc, filter_hydrologyc, &
       waterflux_vars, soilstate_vars, patch_vars)

    use WaterFluxType, only : waterflux_type
    use SoilStateType, only : soilstate_type
    use PatchType, only : patch_type
    use filterMod,            only : clumpfilter
    use decompMod, only : bounds_type

    implicit none

    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points

    type(waterflux_type), intent(in) :: waterflux_vars
    type(soilstate_type), intent(in) :: soilstate_vars
    type(patch_type), intent(in) :: patch_vars

#ifdef CLM_PFLOTRAN
    call step_th_clm_pf(bounds, &
       num_hydrologyc, filter_hydrologyc, &
       waterflux_vars, soilstate_vars, patch_vars)
#else
    character(len=256) :: subname = "clm_pf_step_th"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_step_th


  !-----------------------------------------------------------------------------
  subroutine clm_pf_vecget_gflux(gflux_clm_loc)

    use shr_kind_mod  , only : r8 => shr_kind_r8

    implicit none
    real(r8), pointer :: gflux_clm_loc(:)

#ifdef CLM_PFLOTRAN
    call vecget_gflux_subsurf_clm_pf(gflux_clm_loc)
#else
    character(len=256) :: subname = "clm_pf_vecget_gflux()"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_vecget_gflux



  !-----------------------------------------------------------------------------
  subroutine clm_pf_vecrestore_gflux(gflux_clm_loc)

    use shr_kind_mod  , only : r8 => shr_kind_r8

    implicit none
    real(r8), pointer :: gflux_clm_loc(:)

#ifdef CLM_PFLOTRAN
    call vecrestore_gflux_subsurf_clm_pf(gflux_clm_loc)
#else
    character(len=256) :: subname = "clm_pf_vecrestore_gflux()"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_vecrestore_gflux



  !-----------------------------------------------------------------------------
  subroutine clm_pf_write_restart(date_stamp)

    implicit none
    character(len=*), intent(in) :: date_stamp

#ifdef CLM_PFLOTRAN
    call write_restart_clm_pf(date_stamp)
#else
    character(len=32) :: subname = "clm_pf_write_restart"
    call pflotran_not_available(subname)
#endif
  end subroutine clm_pf_write_restart



!-----------------------------------------------------------------------
!
! private work functions requiring pflotran
!
!-----------------------------------------------------------------------
#ifdef CLM_PFLOTRAN
  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: write_restart_clm_pf
  !
  ! !INTERFACE:
  subroutine write_restart_clm_pf(date_stamp)
  !
  ! !DESCRIPTION:
  ! Trigger a pflotran checkpoint file to be written
  !
  ! !USES:
  ! !ARGUMENTS:
    character(len=32), intent(in) :: date_stamp ! file name date stamp
  ! !LOCAL VARIABLES:


  !EOP
  !-----------------------------------------------------------------------

    call pflotranModelStepperCheckpoint(pflotran_m, date_stamp)

  end subroutine write_restart_clm_pf



  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: vecget_gflux_subsurf_clm_pf
  !
  ! !INTERFACE:
  subroutine vecget_gflux_subsurf_clm_pf(gflux_clm_loc)
  !
  ! !DESCRIPTION:
  ! Wrapper around pflotran init
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    real(r8), pointer :: gflux_clm_loc(:)
  ! !LOCAL VARIABLES:
    PetscErrorCode ierr
  !EOP
  !-----------------------------------------------------------------------
    call VecGetArrayF90(clm_pf_idata%gflux_subsurf_clm, gflux_clm_loc, ierr); CHKERRQ(ierr)
  end subroutine vecget_gflux_subsurf_clm_pf



  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: vecrestore_gflux_subsurf_clm_pf
  !
  ! !INTERFACE:
  subroutine vecrestore_gflux_subsurf_clm_pf(gflux_clm_loc)
  !
  ! !DESCRIPTION:
  ! Wrapper around pflotran init
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    real(r8), pointer :: gflux_clm_loc(:)
    PetscErrorCode ierr
  ! !LOCAL VARIABLES:
  !EOP
  !-----------------------------------------------------------------------

    ! TODO(bja, 2014-02) make gflux = 0.0 a debugging runtime flag...
    !gflux_clm_loc = 0.0_r8
    call VecRestoreArrayF90(clm_pf_idata%gflux_subsurf_clm, gflux_clm_loc, ierr); CHKERRQ(ierr)

  end subroutine vecrestore_gflux_subsurf_clm_pf



  !-----------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: interface_init_clm_pf
  !
  ! !INTERFACE:
  subroutine interface_init_clm_pf(bounds, col, lun, grc, &
       soilstate_vars, waterstate_vars, temperature_vars, &
       soilhydrology_vars)
    !
    ! !DESCRIPTION:
    ! initialize the pflotran iterface
    !
    ! !USES:
    use shr_assert_mod , only : shr_assert
    use shr_log_mod    , only : errMsg => shr_log_errMsg
    use clm_varctl      , only : iulog, fsurdat
    use clm_varctl      , only : pflotran_surfaceflow, pflotran_th_mode, pflotran_th_freezing
    use decompMod       , only : bounds_type, get_proc_total, ldecomp
    use clm_varpar      , only : nlevsoi, nlevgrnd
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use domainMod       , only : ldomain

    use GridcellType    , only : gridcell_type
    use ColumnType      , only : column_type
    use LandUnitType    , only : landunit_type
    use SoilStateType   , only : soilstate_type
    use WaterStateType  , only : waterstate_type
    use TemperatureType , only : temperature_type
    use SoilHydrologyType , only : soilhydrology_type

    use fileutils       , only : getfil
    use spmdMod         , only : mpicom, masterproc, iam
    use organicFileMod  , only : organicrd
    use landunit_varcon , only : istsoil, istice, istdlak, istwet, istice_mec
    use column_varcon   , only : icol_roof, icol_sunwall, icol_shadewall, icol_road_perv, icol_road_imperv
    use clm_varcon      , only : zsoi, denice, denh2o
    use clm_varcon      , only : grlnd
    use abortutils      , only : endrun

    use ncdio_pio

    ! pflotran
    use Option_module, only : printErrMsg
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Realization_class, only : realization_type
    use Surface_Realization_class, only : surface_realization_type
    use PFLOTRAN_Constants_module
    !
    ! !ARGUMENTS:

    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"

    type(bounds_type), intent(in) :: bounds
    type(column_type), intent(in) :: col
    type(landunit_type), intent(in) :: lun
    type(gridcell_type), intent(in) :: grc
    type(soilstate_type), intent(in) :: soilstate_vars
    type(waterstate_type), intent(in) :: waterstate_vars
    type(temperature_type), intent(in) :: temperature_vars
    type(soilhydrology_type), intent(in) :: soilhydrology_vars

    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    integer  :: local_num_g      ! local number of gridcells across this processor
    integer  :: local_num_l      ! local number of landunits across this processor
    integer  :: local_num_c      ! local number of columns across this processor
    integer  :: local_num_p      ! local number of pfts across this processor
    integer  :: local_num_cohorts ! local number of ed cohorts across this processor
    integer  :: g,l,c,lev,j  ! indices
    integer  :: gcount

    !
    ! From iniTimeConst.F90
    !
    type(file_desc_t)  :: ncid   ! netcdf id
    real(r8) :: clay,sand        ! temporaries

    real(r8),pointer :: sand3d(:,:)    ! read in - soil texture: percent sand
    real(r8),pointer :: clay3d(:,:)    ! read in - soil texture: percent clay
    real(r8),pointer :: organic3d(:,:) ! read in - organic matter: kg/m3

    real(r8) :: om_frac                ! organic matter fraction
    real(r8) :: om_watsat    = 0.9_r8  ! porosity of organic soil
    real(r8) :: om_hksat     = 0.1_r8  ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8) :: om_tkm       = 0.25_r8 ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8) :: om_sucsat    = 10.3_r8 ! saturated suction for organic matter (Letts, 2000)
    real(r8) :: om_b         = 2.7_r8  ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8) :: organic_max  = 130._r8 ! organic matter (kg/m3) where soil is assumed to act like peat
    real(r8) :: pc           = 0.5_r8   ! percolation threshold
    real(r8) :: pcbeta       = 0.139_r8 ! percolation exponent
    real(r8) :: perc_frac               ! "percolating" fraction of organic soil
    real(r8) :: perc_norm               ! normalize to 1 when 100% organic soil
    real(r8) :: uncon_hksat             ! series conductivity of mineral/organic soil
    real(r8) :: uncon_frac              ! fraction of "unconnected" soil

    real(r8) :: watsat_tmp, bsw_tmp, sucsat_tmp, press_tmp
    real(r8) :: bd, tkm, bsw2_tmp,psisat_tmp
    real(r8) :: vwcsat_tmp, xksat, hksat_tmp

    character(len=256) :: locfn                    ! local filEname
    character(len= 32) :: subname = 'clm_pf_interface_init' ! subroutine name
    integer :: mxsoil_color                        ! maximum number of soil color classes
    
    integer :: nlevmapped

    logical :: readvar

    integer :: clm_npts
    integer :: clm_surf_npts
    integer :: num_active_columns

    class(realization_type), pointer    :: realization
    class(surface_realization_type), pointer    :: surf_realization

    !PetscViewer :: viewer
    PetscScalar, pointer :: hksat_x_clm_loc(:) ! hydraulic conductivity in x-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_y_clm_loc(:) ! hydraulic conductivity in y-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: hksat_z_clm_loc(:) ! hydraulic conductivity in z-dir at saturation (mm H2O /s)
    PetscScalar, pointer :: watsat_clm_loc(:)  ! minimum soil suction (mm)
    PetscScalar, pointer :: sucsat_clm_loc(:)  ! volumetric soil water at saturation (porosity)
    PetscScalar, pointer :: bsw_clm_loc(:)     ! Clapp and Hornberger "b"
    PetscScalar, pointer :: press_clm_loc(:)   ! Pressure
    PetscScalar, pointer :: temp_clm_loc(:)    ! Temperature
    PetscScalar, pointer :: sat_clm_loc(:)     ! Saturation
    PetscErrorCode :: ierr
    integer, pointer :: clm_cell_ids_nindex(:)
    integer, pointer :: clm_surf_cell_ids_nindex(:)

    associate( &
         ! Assign local pointers to derived subtypes components (landunit-level)
         ltype      =>  lun%itype      , & !  [integer (:)]  landunit type index
         urbpoi     =>  lun%urbpoi     , & !  [logical (:)]  true => landunit is an urban point
         ! Assign local pointer to derived subtypes components (column-level)
         clandunit  =>  col%landunit   , & !  [integer (:)]  landunit index of column
         cgridcell  =>  col%gridcell   , & !  [integer (:)]  gridcell index of column
         cwtgcell   =>  col%wtgcell    , & !  [real(r8) (:)]  weight (relative to gridcell
         ctype      =>  col%itype      , & !  [integer (:)]  column type index
         hksat      =>  soilstate_vars%hksat_col      , & !  [real(r8) (:,:)]  hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
         sucsat     =>  soilstate_vars%sucsat_col     , & !  [real(r8) (:,:)]  minimum soil suction (mm) (nlevgrnd)
         watsat     =>  soilstate_vars%watsat_col     , & !  [real(r8) (:,:)]  volumetric soil water at saturation (porosity) (nlevgrnd)
         h2osoi_vol =>  waterstate_vars%h2osoi_vol_col , & !  [real(r8) (:,:)]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         h2osoi_liq =>  waterstate_vars%h2osoi_liq_col , & !  [real(r8) (:,:)]  liquid water (kg/m2)
         h2osoi_ice =>  waterstate_vars%h2osoi_ice_col , & !  [real(r8) (:,:)]  ice lens (kg/m2)
         t_soisno   =>  temperature_vars%t_soisno_col   , & !  [real(r8) (:,:)]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
         topo       =>  ldomain%topo   , & !  [real(r8) (:)]  topography
         zwt        =>  soilhydrology_vars%zwt_col , & !  [real(r8) (:)]  water table depth (m)
         latdeg     =>  grc%latdeg     , & !  [real(r8) (:)]  latitude (radians)
         londeg     =>  grc%londeg     , & !  [real(r8) (:)]  longitude (radians)
         lakpoi     =>  lun%lakpoi     , & !  [logical (:)]  true => landunit is a lake point
         dz         =>  col%dz           & !  [real(r8) (:,:)]  layer thickness (m)
         )

    ! Determine necessary indices
    call get_proc_total(iam, local_num_g, local_num_l, local_num_c, local_num_p, local_num_cohorts)

    !------------------------------------------------------------------------
    allocate(pflotran_m)

    ! Create PFLOTRAN model
    pflotran_m => pflotranModelCreate(mpicom, pflotran_prefix)

    call pflotranModelSetupRestart(pflotran_m, restart_stamp)

    ! Initialize PETSc vector for data transfer between CLM and PFLOTRAN
    call CLMPFLOTRANIDataInit()

    ! Allocate memory for CLM-PFLOTRAN interface data
    call CreateCLMPFLOTRANInterfaceDate(bounds, clm_npts, clm_surf_npts)

    ! Create CLM-PFLOTRAN mapping files
    call CreateCLMPFLOTRANMaps(bounds, clm_npts, clm_surf_npts)

    nlevmapped     = clm_pf_idata%nzclm_mapped

    call VecGetArrayF90(clm_pf_idata%hksat_x_clm, hksat_x_clm_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_y_clm, hksat_y_clm_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%hksat_z_clm, hksat_z_clm_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%sucsat_clm,  sucsat_clm_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clm,  watsat_clm_loc,  ierr)
    call VecGetArrayF90(clm_pf_idata%bsw_clm,     bsw_clm_loc,     ierr)
    call VecGetArrayF90(clm_pf_idata%press_clm,   press_clm_loc,   ierr)

    write(iulog,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%          Within clm_pf_interface_init               %%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(iulog,*) ' '

    allocate(sand3d(bounds%begg:bounds%endg, nlevsoi), clay3d(bounds%begg:bounds%endg, nlevsoi))
    allocate(organic3d(bounds%begg:bounds%endg, nlevsoi))

    ! --------------------------------------------------------------------
    ! Read soil color, sand and clay from surface dataset
    ! --------------------------------------------------------------------
    if (masterproc) then
       write(iulog,*) 'Attempting to read soil color, sand and clay boundary data .....'
    endif

    call getfil (fsurdat, locfn, 0)
    call ncd_pio_openfile(ncid, locfn, 0)

    ! Determine number of soil color classes - if number of soil color classes is not
    ! on input dataset set it to 8
    call ncd_io(ncid=ncid, varname='mxsoil_color', flag='read', data=mxsoil_color, &
                  readvar=readvar)
    if ( .not. readvar) mxsoil_color = 8

    call ncd_io(ncid=ncid, varname='PCT_SAND', flag='read', data=sand3d, dim1name=grlnd, readvar=readvar)
    if(.not. readvar) call endrun( trim(subname)//' ERROR: PCT_SAND NOT on surfadata file' )

    call ncd_io(ncid=ncid, varname='PCT_CLAY', flag='read', data=clay3d, dim1name=grlnd,readvar=readvar)
    if(.not. readvar) call endrun( trim(subname)//' ERROR: PCT_CLAY NOT on surfadata file' )


    ! Ensure that there is only one active column in each grid cell.
    !
    ! NOTE(bja, 2014-02) should check each grid cell to ensure that
    ! there is only one active column. Suggestion from from Bill
    ! Sacks: "See reweightMod.F90 : checkWeights for
    ! inspiration. Basically, you could have an array of size
    ! begg:endg, where you accumulate a count of the number of active
    ! columns in that grid cell. Then you could ensure that those
    ! counts are 1 at every point."
    num_active_columns = 0
    do c = bounds%begc, bounds%endc
       if (col%active(c)) then
          num_active_columns = num_active_columns + 1
       end if
    end do
    if (local_num_g /= num_active_columns) then
      write (iulog,*), 'ERROR: CLM-PFLOTRAN requires one active column per grid cell: '
      write (iulog, *) "decomMod - num grid cells = ", local_num_g
      write (iulog, *) "num active columns = ", num_active_columns
      call endrun( trim(subname)//' ERROR: More than 1 col per grid cell' )
    endif

    ! --------------------------------------------------------------------
    ! If a organic matter dataset has been specified, read it
    ! --------------------------------------------------------------------

    call organicrd(organic3d)

    gcount = 0 ! assumption that only 1 soil-column per grid cell
    do c = bounds%begc, bounds%endc
      if (.not. col%active(c)) then
         ! only operate on active columns.
         cycle
      end if

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)
      gcount = g - bounds%begg

      if (ltype(l)==istdlak .or. ltype(l)==istwet .or. ltype(l)==istice .or. ltype(l)==istice_mec) then
        write (iulog,*), 'WARNING: Land Unit type Lake/Wet/Ice/Ice_mec ... within the domain'
        write (iulog,*), 'CLM-CN -- PFLOTRAN does not support this land unit presently'
        call endrun( trim(subname)//' ERROR: Land Unit type not supported' )
      else if (urbpoi(l) .and. (ctype(c) /= icol_road_perv) .and. (ctype(c) /= icol_road_imperv) )then
        ! Urban Roof, sunwall, shadewall properties set to special value
        write (iulog,*), 'ERROR: Unsupported Urban Land Unit type: '
        write (iulog, *) "    urbpoi(l) = ", urbpoi(l) 
        write (iulog, *) "    ctype(c) = ", ctype(c)
        write (iulog,*), '  Supported Urban Land Unit types: '
        write (iulog, *) "    icol_road_perv = ", icol_road_perv
        write (iulog, *) "    icol_road_imperv = ", icol_road_imperv
        write (iulog,*), 'CLM-CN -- PFLOTRAN does not support this land unit presently'
        call endrun( trim(subname)//' ERROR: Land Unit type not supported' )
      else  ! soil columns of both urban and non-urban types

        do lev = 1,nlevgrnd

          ! duplicate clay and sand values from 10th soil layer
          if (lev .le. nlevsoi) then
            clay    = clay3d(g,lev)
              sand    = sand3d(g,lev)
            om_frac = (organic3d(g,lev)/organic_max)**2._r8
          else
            clay    = clay3d(g,nlevsoi)
            sand    = sand3d(g,nlevsoi)
            om_frac = 0._r8
          endif

          ! No organic matter for urban
          if (urbpoi(l)) then
            om_frac = 0._r8
          end if

          watsat_tmp = 0.489_r8 - 0.00126_r8*sand
          bsw_tmp    = 2.91 + 0.159*clay
          sucsat_tmp = 10._r8 * ( 10._r8**(1.88_r8-0.0131_r8*sand) )
          bd            = (1._r8-watsat_tmp)*2.7e3_r8
          watsat_tmp = (1._r8 - om_frac)*watsat_tmp + om_watsat*om_frac
          tkm           = (1._r8-om_frac)*(8.80_r8*sand+2.92_r8*clay)/(sand+clay)+om_tkm*om_frac ! W/(m K)
          bsw_tmp    = (1._r8-om_frac)*(2.91_r8 + 0.159_r8*clay) + om_frac*om_b
          bsw2_tmp   = -(3.10_r8 + 0.157_r8*clay - 0.003_r8*sand)
          psisat_tmp = -(exp((1.54_r8 - 0.0095_r8*sand + 0.0063_r8*(100.0_r8-sand-clay))*log(10.0_r8))*9.8e-5_r8)
          vwcsat_tmp = (50.5_r8 - 0.142_r8*sand - 0.037_r8*clay)/100.0_r8
          sucsat_tmp = (1._r8-om_frac)*sucsat_tmp + om_sucsat*om_frac
          xksat         = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s

          ! perc_frac is zero unless perf_frac greater than percolation threshold
          if (om_frac > pc) then
            perc_norm=(1._r8 - pc)**(-pcbeta)
            perc_frac=perc_norm*(om_frac - pc)**pcbeta
          else
            perc_frac=0._r8
          endif

          !
          !                                ||
          !                                ||
          !                               \||/
          !
          !
          !            ******************************************
          !            *                                        *
          !            *                                        *
          !            *           ORGANIC PERCOLATING          *
          !            *               f_om*f_pre               *
          !            *                                        *
          ! ---\       ******************************************
          ! ---/       *                    *                   *
          !            *                    *                   *
          !            *                    *                   *
          !            *     ORGANIC        *       MINERAL     *
          !            *  NON-PERCOLATING   *                   *
          !            *                    *                   *
          !            *   f_om*(1-f_pre)   *       1-f_om      *
          !            *                    *                   *
          !            ******************************************
          !

          ! ---------------------------------------------------------------
          ! Hydraulic conductivity in Z-direction
          ! ---------------------------------------------------------------

          ! uncon_frac is fraction of mineral soil plus fraction of "nonpercolating" organic soil
          uncon_frac=(1._r8-om_frac)+(1._r8-perc_frac)*om_frac

          ! uncon_hksat is series addition of mineral/organic conductivites
          if (om_frac .lt. 1._r8) then
            uncon_hksat = uncon_frac/((1._r8-om_frac)/xksat &
                         + ((1._r8-perc_frac)*om_frac)/om_hksat)
          else
            uncon_hksat = 0._r8
          end if
          hksat_tmp  = uncon_frac*uncon_hksat + (perc_frac*om_frac)*om_hksat

          ! ---------------------------------------------------------------
          ! Hydraulic conductivity in X/Y-direction
          ! ---------------------------------------------------------------

          !
          if (om_frac .lt. 1._r8) then
            uncon_hksat = ( (1._r8 - om_frac)*xksat + &
                            (1._r8 - perc_frac)*om_frac*om_hksat )/uncon_frac
            hksat_tmp   = uncon_hksat*om_hksat/(om_frac*perc_frac*uncon_hksat &
                          + (1._r8 - om_frac*perc_frac)*om_hksat)
          else
            uncon_hksat = 0._r8
            hksat_tmp   = om_hksat
          end if

          press_tmp = 101325.0_r8 - 998.2_r8*9.81_r8*(zwt(c) - zsoi(lev))
          press_tmp = 101325.0_r8 - 998.2_r8*9.81_r8*(2.0_r8 - zsoi(lev))

          if (lev <= nlevmapped) then

            hksat_x_clm_loc(gcount*nlevmapped + lev ) = hksat_x_clm_loc(gcount*nlevmapped + lev ) + hksat_tmp*cwtgcell(c)
            hksat_y_clm_loc(gcount*nlevmapped + lev ) = hksat_y_clm_loc(gcount*nlevmapped + lev ) + hksat_tmp*cwtgcell(c)
            hksat_z_clm_loc(gcount*nlevmapped + lev ) = hksat_z_clm_loc(gcount*nlevmapped + lev ) + hksat(c,lev)*cwtgcell(c)
            sucsat_clm_loc( gcount*nlevmapped + lev ) = sucsat_clm_loc( gcount*nlevmapped + lev ) + sucsat(c,lev)*cwtgcell(c)
            watsat_clm_loc( gcount*nlevmapped + lev ) = watsat_clm_loc( gcount*nlevmapped + lev ) + watsat(c,lev)*cwtgcell(c)
            bsw_clm_loc(    gcount*nlevmapped + lev ) = bsw_clm_loc(    gcount*nlevmapped + lev ) + bsw_tmp*cwtgcell(c)
            press_clm_loc(  gcount*nlevmapped + lev ) = press_clm_loc(  gcount*nlevmapped + lev ) + press_tmp*cwtgcell(c)
          endif

        enddo
      endif
    enddo ! do c = bounds%begc, bounds%endc

    call VecRestoreArrayF90(clm_pf_idata%hksat_x_clm, hksat_x_clm_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_y_clm, hksat_y_clm_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%hksat_z_clm, hksat_z_clm_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%sucsat_clm,  sucsat_clm_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_clm,  watsat_clm_loc,  ierr)
    call VecRestoreArrayF90(clm_pf_idata%bsw_clm,     bsw_clm_loc,     ierr)
    call VecRestoreArrayF90(clm_pf_idata%press_clm,   press_clm_loc,   ierr)

    ! Set CLM soil properties onto PFLOTRAN grid
    call pflotranModelSetSoilProp(pflotran_m)
    !call pflotranModelSetICs(pflotran_m)

    ! Initialize PFLOTRAN states
    call pflotranModelStepperRunInit(pflotran_m)

    ! Get top surface area
    call pflotranModelGetTopFaceArea(pflotran_m)

    ! Get PFLOTRAN states
    call pflotranModelGetUpdatedData(pflotran_m)

    ! Initialize soil temperature
    if(pflotran_m%option%iflowmode==TH_MODE) then
      call VecGetArrayF90(clm_pf_idata%temp_clm, temp_clm_loc, ierr)
      do c = bounds%begc, bounds%endc
        if (col%active(c)) then
          l = clandunit(c)
          if (.not. lakpoi(l)) then  !not lake
            g = cgridcell(c)
            gcount = g - bounds%begg
            do j = 1, nlevmapped
              t_soisno(c, j) = temp_clm_loc(gcount*nlevmapped + j) + 273.15_r8
            enddo
            if ( nlevmapped /= nlevgrnd) then
              t_soisno(c, nlevmapped+1:nlevgrnd) = t_soisno(c, nlevmapped)
            end if
          endif
        endif
      enddo
      call VecRestoreArrayF90(clm_pf_idata%temp_clm, temp_clm_loc, ierr)
    endif

    ! Initialize soil moisture
    call VecGetArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr)
    do c = bounds%begc, bounds%endc
      if (col%active(c)) then
        l = clandunit(c)
        if (ltype(l) == istsoil) then
          g = cgridcell(c)
          gcount = g - bounds%begg
          do j = 1, nlevsoi
            h2osoi_liq(c,j) = sat_clm_loc(gcount*nlevmapped + j)*dz(c,j)*1.e3_r8
            h2osoi_vol(c,j) = h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                 h2osoi_ice(c,j)/dz(c,j)/denice
            h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
          enddo
        endif
      endif
    enddo
    call VecRestoreArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr)

    deallocate(sand3d,clay3d,organic3d)
    end associate
  end subroutine interface_init_clm_pf

  !-----------------------------------------------------------------------
  subroutine CreateCLMPFLOTRANInterfaceDate(bounds, clm_npts, clm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Allocates memory for CLM-PFLOTRAN interface
    !
    ! !USES:
    use clm_varctl      , only : pflotran_surfaceflow, pflotran_th_mode, pflotran_th_freezing
    use decompMod       , only : bounds_type, ldecomp
    use clm_varpar      , only : nlevsoi, nlevgrnd
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use spmdMod         , only : mpicom, masterproc, iam
    use abortutils      , only : endrun
    ! pflotran
    use Option_module, only : printErrMsg
    use Simulation_Base_class, only : simulation_base_type
    use Subsurface_Simulation_class, only : subsurface_simulation_type
    use Surface_Simulation_class, only : surface_simulation_type
    use Surf_Subsurf_Simulation_class, only : surfsubsurface_simulation_type
    use Realization_class, only : realization_type
    use Surface_Realization_class, only : surface_realization_type
    use PFLOTRAN_Constants_module
    !
    implicit none
    !
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(inout) :: clm_npts
    integer, intent(inout)  :: clm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer  :: g,l,c,lev,j  ! indices
    integer  :: gcount
    integer :: nlevmapped
    character(len= 128) :: subname = 'CreateCLMPFLOTRANInterfaceDate' ! subroutine name
    class(realization_type), pointer    :: realization
    class(surface_realization_type), pointer    :: surf_realization

    ! Initialize PETSc vector for data transfer between CLM and PFLOTRAN
    call CLMPFLOTRANIDataInit()

    select type (simulation => pflotran_m%simulation)
      class is (subsurface_simulation_type)
         realization => simulation%realization
         nullify(surf_realization)
      class is (surface_simulation_type)
         nullify(realization)
         surf_realization => simulation%surf_realization
      class is (surfsubsurface_simulation_type)
         realization => simulation%realization
         surf_realization => simulation%surf_realization
      class default
         pflotran_m%option%io_buffer = "clm-pflotran only works with surface and subsurface simulations."
         write(*, '(/A/)') pflotran_m%option%io_buffer
         call printErrMsg(pflotran_m%option)
    end select

    ! Compute number of cells in CLM domain.
    ! Assumption-1: One column per CLM grid cell.

    ! Check if the number of CLM vertical soil layers defined in the mapping
    ! file read by PFLOTRAN matches either nlevsoi or nlevgrnd
    clm_pf_idata%nzclm_mapped = pflotran_m%map_clm_sub_to_pf_sub%clm_nlevsoi
    nlevmapped                = clm_pf_idata%nzclm_mapped
    if ( (nlevmapped /= nlevsoi) .and. (nlevmapped /= nlevgrnd) ) then
       call endrun(trim(subname)//' ERROR: Number of layers PFLOTRAN thinks CLM should '// &
                'have do not match either nlevsoi or nlevgrnd. Abortting' )
    end if

    clm_npts = (bounds%endg - bounds%begg + 1)*nlevmapped
    clm_surf_npts = (bounds%endg - bounds%begg + 1)

    ! CLM: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_sub = clm_npts
    clm_pf_idata%ngclm_sub = clm_npts

    ! CLM: Surface of subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_2dsub = (bounds%endg - bounds%begg + 1)
    clm_pf_idata%ngclm_2dsub = (bounds%endg - bounds%begg + 1)
    ! For CLM: Same as surface of subsurface domain
    clm_pf_idata%nlclm_srf = clm_surf_npts
    clm_pf_idata%ngclm_srf = clm_surf_npts

    ! PFLOTRAN: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlpf_sub = realization%patch%grid%nlmax
    clm_pf_idata%ngpf_sub = realization%patch%grid%ngmax

    ! PFLOTRAN: Surface of subsurface domain (local and ghosted cells)
    if (pflotran_m%option%iflowmode == TH_MODE) then
      clm_pf_idata%nlpf_2dsub = pflotranModelNSurfCells3DDomain(pflotran_m)
      clm_pf_idata%ngpf_2dsub = pflotranModelNSurfCells3DDomain(pflotran_m)
    else
      clm_pf_idata%nlpf_2dsub = 0
      clm_pf_idata%ngpf_2dsub = 0
    endif

    ! PFLOTRAN: Surface domain (local and ghosted cells)
    if (associated(surf_realization) .and. pflotran_m%option%nsurfflowdof > 0) then
      clm_pf_idata%nlpf_srf = surf_realization%patch%grid%nlmax
      clm_pf_idata%ngpf_srf = surf_realization%patch%grid%ngmax
    else
      clm_pf_idata%nlpf_srf = 0
      clm_pf_idata%ngpf_srf = 0
    endif

    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
    call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

  end subroutine CreateCLMPFLOTRANInterfaceDate

  !-----------------------------------------------------------------------
  subroutine CreateCLMPFLOTRANMaps(bounds, clm_npts, clm_surf_npts)
    !
    ! !DESCRIPTION:
    ! Creates maps to transfer data between CLM and PFLOTRAN
    !
    ! !USES:
    use clm_varctl      , only : pflotran_surfaceflow, pflotran_th_mode, pflotran_th_freezing
    use decompMod       , only : bounds_type, ldecomp
    use clm_varpar      , only : nlevsoi, nlevgrnd
    use PFLOTRAN_Constants_module
    !
    implicit none
    !
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    integer, intent(inout) :: clm_npts
    integer, intent(inout)  :: clm_surf_npts
    !
    ! LOCAL VARAIBLES:
    integer  :: g,l,c,lev,j  ! indices
    integer  :: gcount
    character(len= 128) :: subname = 'CreateCLMPFLOTRANMaps' ! subroutine name
    integer :: nlevmapped
    integer, pointer :: clm_cell_ids_nindex(:)
    integer, pointer :: clm_surf_cell_ids_nindex(:)

    ! Save cell IDs of CLM grid
    allocate(clm_cell_ids_nindex(     1:clm_npts     ))
    allocate(clm_surf_cell_ids_nindex(1:clm_surf_npts))

    nlevmapped     = clm_pf_idata%nzclm_mapped
    clm_npts       = 0
    clm_surf_npts  = 0
    do g = bounds%begg, bounds%endg
       do j = 1,nlevmapped
          clm_npts = clm_npts + 1
          clm_cell_ids_nindex(clm_npts) = (ldecomp%gdc2glo(g)-1)*nlevmapped + j - 1
       enddo
       clm_surf_npts = clm_surf_npts + 1
       clm_surf_cell_ids_nindex(clm_surf_npts) = (ldecomp%gdc2glo(g)-1)*nlevmapped
    enddo

    ! Initialize maps for transferring data between CLM and PFLOTRAN.
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, CLM_SUB_TO_PF_SUB)
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, CLM_SUB_TO_PF_EXTENDED_SUB)
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, PF_SUB_TO_CLM_SUB)

    if (pflotran_m%option%iflowmode == TH_MODE) then
      pflotran_th_mode = .true.
      if (pflotran_m%option%use_th_freezing) pflotran_th_freezing  = .true.
    endif

    if (pflotran_m%option%nsurfflowdof > 0) then
      pflotran_surfaceflow = .true.
      call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, PF_SRF_TO_CLM_SRF)
      call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, CLM_SRF_TO_PF_SRF)
    else
      if (pflotran_m%option%iflowmode == TH_MODE) then
        call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids_nindex, &
                                      clm_surf_npts, CLM_SRF_TO_PF_2DSUB)
      endif
    endif

    deallocate(clm_cell_ids_nindex)
    deallocate(clm_surf_cell_ids_nindex)

  end subroutine CreateCLMPFLOTRANMaps

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: set_sflow_forcing_clm_pf
  !
  ! !INTERFACE:
  subroutine set_sflow_forcing_clm_pf(num_hydrologyc, filter_hydrologyc, &
       atm2lnd_vars, waterflux_vars)
  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
    use shr_kind_mod    , only: r8 => shr_kind_r8
    use atm2lndType, only : atm2lnd_type
    use WaterfluxType, only : waterflux_type
    use clm_time_manager, only : get_step_size
    use filterMod,            only : clumpfilter
    use abortutils  , only : endrun

    integer, intent(in) :: num_hydrologyc        ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:)  ! column filter for soil points
    type(atm2lnd_type), intent(in) :: atm2lnd_vars
    type(waterflux_type), intent(in) :: waterflux_vars

  ! !LOCAL VARIABLES:
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    integer  :: c,fc                       !indices
    integer  :: count
    PetscScalar, pointer :: rain_clm_loc(:)
    PetscScalar, pointer :: rain_temp_clm_loc(:)
    PetscErrorCode :: ierr
    !-----------------------------------------------------------------------

    associate( &
    qflx_snow_h2osfc  =>    waterflux_vars%qflx_snow_h2osfc_col  , & ! Input:  [real(r8) (:)]  snow falling on surface water (mm/s)
    qflx_floodc       =>    waterflux_vars%qflx_floodc_col       , & ! Input:  [real(r8) (:)]  column flux of flood water from RTM
    qflx_top_soil     =>    waterflux_vars%qflx_top_soil_col     , & ! Input:  [real(r8) (:)]  net water input into soil from top (mm/s)
    qflx_infl         =>    waterflux_vars%qflx_infl_col         , & ! Output: [real(r8) (:)] infiltration (mm H2O /s)
    qflx_surf         =>    waterflux_vars%qflx_surf_col         , & ! Output: [real(r8) (:)]  surface runoff (mm H2O /s)
    forc_t            =>    atm2lnd_vars%forc_t_downscaled_col & ! Input:  [real(r8) (:)]  atmospheric temperature (Kelvin)
    )

    call VecGetArrayF90(clm_pf_idata%rain_clm,rain_clm_loc,ierr)
    call VecGetArrayF90(clm_pf_idata%rain_temp_clm,rain_temp_clm_loc,ierr)
    count = 0
    do fc = 1, num_hydrologyc

      c = filter_hydrologyc(fc)
      qflx_surf(c) = 0._r8
      qflx_infl(c) = 0._r8
      count = count + 1

      ! Convert mm/s to m/s
      rain_clm_loc(count) = (qflx_top_soil(c) + qflx_snow_h2osfc(c) + &
                             qflx_floodc(c))/1000._r8

      ! Assumption: Rain water is at air-temperature
      ! Convert Kelvin to degC
      rain_temp_clm_loc(count) = forc_t(c) - 273.15_r8
    enddo
    call VecRestoreArrayF90(clm_pf_idata%rain_clm,rain_clm_loc,ierr)
    call VecRestoreArrayF90(clm_pf_idata%rain_temp_clm,rain_temp_clm_loc,ierr)

    end associate
  end subroutine set_sflow_forcing_clm_pf

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_soil_moisture_clm_pf
  !
  ! !INTERFACE:
  subroutine update_soil_moisture_clm_pf(waterstate_vars, soilstate_vars, bounds, &
       num_hydrologyc, filter_hydrologyc, col, soilhydrology_vars)
  !
  ! !DESCRIPTION:
  ! 
  !
  ! !USES:
    use ColumnType, only : column_type
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use WaterStateType, only : waterstate_type
    use SoilStateType, only : soilstate_type
    use SoilHydrologyType, only : soilhydrology_type
    use decompMod           , only : bounds_type
    use clm_time_manager    , only : get_nstep
    use clm_varcon,           only : denh2o, denice
    use clm_varpar,           only : nlevsoi
    use filterMod,            only : clumpfilter
    use PFLOTRAN_Constants_module

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(waterstate_type), intent(inout) :: waterstate_vars
    type(soilstate_type), intent(inout) :: soilstate_vars
    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc        ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:)  ! column filter for soil points
    type(column_type), intent(in) :: col
    type(soilhydrology_type), intent(out) :: soilhydrology_vars

  ! !LOCAL VARIABLES:
    integer  :: nstep                    ! time step number
    integer  :: fc, c, g, gcount   ! indices

    ! col%dz(:,:) == layer thickness depth (m)
    ! waterstate_vars%h2osoi_liq(:,:)  == liquid water (kg/m2)
    ! h2osoi_ice(:,:)  == ice lens (kg/m2)
    ! h2osoi_vol(:,:) == volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    ! col%cgridcell(:) == column's gridcell
    PetscScalar, pointer :: sat_clm_loc(:)
    PetscScalar, pointer :: watsat_clm_loc(:)
    PetscScalar, pointer :: sat_ice_clm_loc(:)
    PetscErrorCode :: ierr
    integer :: j
    integer :: nlevmapped
  !EOP
  !-----------------------------------------------------------------------

  ! =======================================================================
  ! For NSTEP=0; update the soil moisture values that were initialized in
  !              PFLOTRAN. Variables modified
  !   h2osoi_liq [kg/m^2]  
  !   h2osoi_vol[m^3/m^3] (water + ice)
  ! =======================================================================

    nstep   = get_nstep()

    call VecGetArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr); CHKERRQ(ierr)

    nlevmapped = clm_pf_idata%nzclm_mapped

    if (pflotran_m%option%iflowmode == TH_MODE .and. &
        pflotran_m%option%use_th_freezing) then

      call VecGetArrayF90(clm_pf_idata%sat_clm, sat_ice_clm_loc, ierr); CHKERRQ(ierr)

      do fc = 1,num_hydrologyc
        c = filter_hydrologyc(fc)
        g = col%gridcell(c)
        gcount = g - bounds%begg
        do j = 1, nlevsoi
          waterstate_vars%h2osoi_liq_col(c,j) = sat_clm_loc(gcount*nlevmapped + j)*soilstate_vars%watsat_col(c,j)*col%dz(c,j)*denh2o
          waterstate_vars%h2osoi_ice_col(c,j) = sat_ice_clm_loc(gcount*nlevmapped + j)*soilstate_vars%watsat_col(c,j)*col%dz(c,j)*denice
          waterstate_vars%h2osoi_vol_col(c,j) = waterstate_vars%h2osoi_liq_col(c,j)/col%dz(c,j)/denh2o + &
                                waterstate_vars%h2osoi_ice_col(c,j)/col%dz(c,j)/denice
          waterstate_vars%h2osoi_vol_col(c,j) = min(waterstate_vars%h2osoi_vol_col(c,j), soilstate_vars%watsat_col(c,j))
        enddo
      enddo

      call VecRestoreArrayF90(clm_pf_idata%sat_clm, sat_ice_clm_loc, ierr); CHKERRQ(ierr)

    else

      do fc = 1,num_hydrologyc
        c = filter_hydrologyc(fc)
        g = col%gridcell(c)
        gcount = g - bounds%begg
        do j = 1, nlevsoi
          waterstate_vars%h2osoi_liq_col(c,j) = sat_clm_loc(gcount*nlevmapped + j)*soilstate_vars%watsat_col(c,j)*col%dz(c,j)*denh2o
          waterstate_vars%h2osoi_vol_col(c,j) = waterstate_vars%h2osoi_liq_col(c,j)/col%dz(c,j)/denh2o + &
                                waterstate_vars%h2osoi_ice_col(c,j)/col%dz(c,j)/denice
          waterstate_vars%h2osoi_vol_col(c,j) = min(waterstate_vars%h2osoi_vol_col(c,j), soilstate_vars%watsat_col(c,j))
        enddo
      enddo

    endif

    call VecRestoreArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr); CHKERRQ(ierr)

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
      soilhydrology_vars%qcharge_col(c) = 0.0_r8
    end do

  end subroutine update_soil_moisture_clm_pf



  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: step_th_clm_pf
  !
  ! !INTERFACE:
  subroutine step_th_clm_pf(bounds, &
       num_hydrologyc, filter_hydrologyc, &
       waterflux_vars, soilstate_vars, patch_vars)
  !
  ! !DESCRIPTION:
  ! 
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8

    use ColumnType, only : col
    use PatchType, only : pft
    use WaterFluxType, only : waterflux_type
    use SoilStateType, only : soilstate_type
    use PatchType, only : patch_type

    use pflotran_model_module, only :pflotranModelUpdateFlowConds, &
         pflotranModelStepperRunTillPauseTime, &
         pflotranModelGetUpdatedData

    use clm_pflotran_interface_data
    use decompMod                  , only : bounds_type
    use clm_varpar                 , only : max_patch_per_col
    use clm_varpar      , only : nlevsoi
    use clm_time_manager, only : get_step_size, get_nstep, is_perpetual
    use abortutils  , only : endrun

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type), intent(in) :: bounds
    integer, intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer, intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterflux_type), intent(in) :: waterflux_vars
    type(soilstate_type), intent(in) :: soilstate_vars
    type(patch_type), intent(in) :: patch_vars

  ! !LOCAL VARIABLES:
    integer  :: c, fc, g, gcount, j, p   ! do loop indices
    integer  :: pftindex                        ! pft index
    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: nstep                      ! time step number

    !real(r8) :: den
    real(r8) :: temp(bounds%begc:bounds%endc) ! accumulator for rootr weighting

    PetscScalar, pointer :: sat_clm_loc(:)    !
    PetscScalar, pointer :: qflx_clm_loc(:)   !
    PetscScalar, pointer :: area_clm_loc(:)   !
    PetscErrorCode :: ierr
    integer :: nlevmapped
    real(r8) :: area

  !EOP
  !-----------------------------------------------------------------------
    !den = 998.2_r8 ! [kg/m^3]
    !den = 1000._r8 ! [kg/m^3]

    associate( &
         qflx_tran_veg_col => waterflux_vars%qflx_tran_veg_col  , & !  [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)
         ! Assign local pointers to derived type members (pft-level)
         qflx_tran_veg_patch => waterflux_vars%qflx_tran_veg_patch    , & !  [real(r8) (:)]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_evap_soi_patch => waterflux_vars%qflx_evap_soi_patch    , & !  [real(r8) (:)]  soil evaporation (mm H2O/s) (+ = to atm)
         rootr_patch         => soilstate_vars%rootr_patch            , & !  [real(r8) (:,:)]  effective fraction of roots in each soil layer
         pwtgcell          => pft%wtgcell          , & !  [real(r8) (:)]  weight relative to gridcell for each pft
         pwtcol            => pft%wtcol            , & !  [real(r8) (:)]  weight relative to column for each pft
         pfti              => col%pfti             , & !  [integer (:)]  beginning pft index for each column
         rootr_col         => soilstate_vars%rootr_col , & !  [real(r8) (:,:)]  effective fraction of roots in each soil layer
         wtcol             => patch_vars%wtcol              & !  [real(r8) (:)]  pft weight relative to column
         )

    nstep = get_nstep()
    dtime = get_step_size()

    call VecGetArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%qflx_clm, qflx_clm_loc, ierr); CHKERRQ(ierr)
    call VecGetArrayF90(clm_pf_idata%area_top_face_clm, area_clm_loc, ierr); CHKERRQ(ierr)

    nlevmapped = clm_pf_idata%nzclm_mapped

    ! Initialize to ZERO
    do g = bounds%begg, bounds%endg
      do j = 1,nlevmapped
        gcount = g - bounds%begg
        qflx_clm_loc(gcount*nlevmapped + j ) = 0.0_r8
      end do
    end do

    ! Compute the Infiltration - Evaporation at each grid-level
    ! qflx_infl [mm/sec],
    ! qflx_clm_loc [m^3/sec] (assuming top surf-area = 1 m^2)
    !
    ! [m^3/s] = [mm/s]/1000
    !

    ! Note: When surface-flows are turned on in PFLOTRAN, qflx_infl(c) is set
    !       to 0.0_r8.
    do c = bounds%begc, bounds%endc
      if (col%active(c)) then
        ! Set gridcell indices
        g = col%gridcell(c)
        gcount = g - bounds%begg
        j = 1
        area = area_clm_loc(gcount*nlevsoi+j)
        qflx_clm_loc(gcount*nlevmapped + j) = qflx_clm_loc(gcount*nlevmapped + j) + &
             waterflux_vars%qflx_infl_col(c)*col%wtgcell(c)*area
      end if
    enddo

    ! Compute the Transpiration sink at grid-level for each soil layer
    ! qflx_tran_veg_patch [mm/sec], while
    ! qflx_clm_loc      [kg/sec] (assuming top surf-area = 1 m^2)

    ! (i) Initialize root faction at column level to be zero
    do j = 1, nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        rootr_col(c,j) = 0._r8
      end do
    end do

    temp(:) = 0._r8

    ! (ii) Compute the root fraction at column level
    do pftindex = 1,max_patch_per_col
      do j = 1,nlevsoi
        do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (pftindex <= col%npfts(c)) then
            p = pfti(c) + pftindex - 1
            if (pwtgcell(p)>0._r8) then
              rootr_col(c,j) = rootr_col(c,j) + &
              rootr_patch(p,j) * qflx_tran_veg_patch(p) * pwtcol(p)
            end if
          end if
        end do
      end do


      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        if (pftindex <= col%npfts(c)) then
          p = pfti(c) + pftindex - 1
          if (pwtgcell(p)>0._r8) then
            temp(c) = temp(c) + qflx_tran_veg_patch(p) * pwtcol(p)
          end if
        end if
      end do
    end do

    ! (iii) Compute the Transpiration sink
    do j = 1, nlevsoi
      do fc = 1, num_hydrologyc
        c = filter_hydrologyc(fc)
        g = col%gridcell(c)
        gcount = g - bounds%begg
        if (temp(c) /= 0._r8) then
          rootr_col(c,j) = rootr_col(c,j)/temp(c)
          area = area_clm_loc(gcount*nlevsoi+j)
          qflx_clm_loc(gcount*nlevmapped + j ) = &
                              qflx_clm_loc(gcount*nlevmapped + j ) - &
                              qflx_tran_veg_col(c)*rootr_col(c,j)*area
        end if
      end do
    end do

    call VecRestoreArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%qflx_clm, qflx_clm_loc, ierr); CHKERRQ(ierr)
    call VecRestoreArrayF90(clm_pf_idata%area_top_face_clm, area_clm_loc, ierr); CHKERRQ(ierr)

    call pflotranModelUpdateFlowConds( pflotran_m )
    call pflotranModelStepperRunTillPauseTime( pflotran_m, (nstep+1.0d0)*dtime )
    call pflotranModelGetUpdatedData( pflotran_m )

    end associate
  end subroutine step_th_clm_pf

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_soil_temperature_clm_pf
  !
  ! !INTERFACE:
  subroutine update_soil_temperature_clm_pf(bounds, &
       num_nolakec, filter_nolakec, &
       temperature_vars, col)

  !
  ! !DESCRIPTION:
  ! 
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use clm_time_manager  , only : get_step_size
    use clm_varpar    , only : nlevsno, nlevgrnd
    use decompMod     , only : bounds_type
    use ColumnType, only : column_type
    use TemperatureType, only : temperature_type

    use clm_pflotran_interface_data, only : clm_pf_idata

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type) , intent(in)  :: bounds
    integer , intent(in)  :: num_nolakec         ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)   ! column filter for non-lake points
    type(temperature_type), intent(in) :: temperature_vars
    type(column_type), intent(in) :: col

  ! !LOCAL VARIABLES:
    integer  :: j,c,g                     !  indices
    integer  :: fc                        ! lake filtered column indices
    integer  :: gcount
    integer  :: nlevmapped

    PetscScalar, pointer :: temp_clm_loc(:)  !
    PetscErrorCode :: ierr
  !EOP
  !-----------------------------------------------------------------------

    associate( &
         t_soisno   => temperature_vars%t_soisno_col  , & !  [real(r8) (:,:)]  soil temperature (Kelvin)
         cgridcell  => col%gridcell    & !  [integer (:)]  column's gridcell
         )

    call VecGetArrayF90(clm_pf_idata%temp_clm, temp_clm_loc, ierr); CHKERRQ(ierr)

    nlevmapped = clm_pf_idata%nzclm_mapped

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       g = cgridcell(c)
       gcount = g - bounds%begg
       do j = 1, nlevmapped
          t_soisno(c,j) = temp_clm_loc(gcount*nlevmapped+j) + 273.15_r8
       enddo
       if ( nlevmapped /= nlevgrnd) then
          t_soisno(c, nlevmapped+1:nlevgrnd) = t_soisno(c, nlevmapped)
       end if
    enddo

    call VecRestoreArrayF90(clm_pf_idata%temp_clm, temp_clm_loc, ierr); CHKERRQ(ierr)
    end associate
  end subroutine update_soil_temperature_clm_pf

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_soil_therm_cond_clm_pf
  !
  ! !INTERFACE:
  subroutine update_soil_therm_cond_clm_pf(bounds, &
       num_urbanl, filter_urbanl, &
       num_nolakec, filter_nolakec, &
       tk)

  !
  ! !DESCRIPTION:
  ! Overwrites the effective thermal conductivity values computed by CLM with
  ! values obtained from PFLOTRAN.
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use clm_time_manager  , only : get_step_size
    use clm_varctl    , only : iulog
    use clm_varpar    , only : nlevsno, nlevgrnd, nlevsoi
    use decompMod     , only : bounds_type
    use ColumnType      , only : col

    use clm_pflotran_interface_data, only : clm_pf_idata

  ! !ARGUMENTS:
    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    type(bounds_type) , intent(in)  :: bounds
    integer , intent(in)  :: num_nolakec         ! number of column non-lake points in column filter
    integer , intent(in)  :: filter_nolakec(:)   ! column filter for non-lake points
    integer , intent(in)  :: num_urbanl          ! number of urban landunits in clump
    integer , intent(in)  :: filter_urbanl(:)    ! urban landunit filter
    real(r8), intent(out) :: tk (bounds%begc:bounds%endc,-nlevsno+1:nlevgrnd)          ! thermal conductivity [W/(m K)]

  ! !LOCAL VARIABLES:
    character(len=256) :: subname = 'Update_soil_temperature_clm_pf'
    integer  :: j,c,g                     !  indices
    integer  :: fc                        ! lake filtered column indices
    integer  :: gcount
    integer  :: nlevmapped

    PetscScalar, pointer :: eff_therm_cond_clm_loc(:)  !
    PetscErrorCode :: ierr
  !EOP
  !-----------------------------------------------------------------------

    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(tk)  == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    associate( &
      cgridcell  =>  col%gridcell    & !  [integer (:)]  column's gridcell
    )

    call VecGetArrayF90(clm_pf_idata%eff_therm_cond_clm, eff_therm_cond_clm_loc, ierr); CHKERRQ(ierr)

    nlevmapped = clm_pf_idata%nzclm_mapped

    do fc = 1,num_nolakec
       c = filter_nolakec(fc)
       g = cgridcell(c)
       gcount = g - bounds%begg
       do j = 1, nlevmapped
          tk(c,j) = eff_therm_cond_clm_loc(gcount*nlevmapped+j)
       enddo
       if ( nlevmapped /= nlevgrnd) then
          tk(c, nlevmapped+1:nlevgrnd) = tk(c,j)
       end if
    enddo

    call VecRestoreArrayF90(clm_pf_idata%eff_therm_cond_clm, eff_therm_cond_clm_loc, ierr); CHKERRQ(ierr)

    end associate
  end subroutine update_soil_therm_cond_clm_pf

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_drainage_clm_pf
  !
  ! !INTERFACE:
  subroutine update_drainage_clm_pf(num_hydrologyc, filter_hydrologyc, &
       waterflux_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8

    use WaterFluxType, only : waterflux_type

    implicit none
    integer , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    type(waterflux_type), intent(out) :: waterflux_vars

!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
    integer  :: c,fc                               ! indices

    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       waterflux_vars%qflx_drain_perched_col(c) = 0._r8   ! perched wt sub-surface runoff (mm H2O /s)
       waterflux_vars%qflx_drain_col(c)         = 0._r8   ! sub-surface runoff (mm H2O /s)
       waterflux_vars%qflx_irrig_patch(c)         = 0._r8   ! irrigation flux (mm H2O /s)
       waterflux_vars%qflx_qrgwl_col(c)         = 0._r8   ! qflx_surf at glaciers, wetlands, lakes (mm H2O /s)
       waterflux_vars%qflx_rsub_sat_col(c)      = 0._r8   ! soil saturation excess [mm h2o/s]
    end do

    end subroutine update_drainage_clm_pf

  !-----------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: update_h2osfc_clm_pf
  !
  ! !INTERFACE:
  subroutine update_h2osfc_clm_pf(num_hydrologyc, filter_hydrologyc, &
       waterstate_vars)

  !
  ! !DESCRIPTION:
  !
  !
  ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8

    use WaterStateType, only : waterstate_type

    implicit none

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"

    integer , intent(in) :: num_hydrologyc       ! number of column soil points in column filter
    integer , intent(in) :: filter_hydrologyc(:) ! column filter for soil points
    ! FIXME(bja, 201409) just pass in h2osfc_col instead of the entire structure!
    type(waterstate_type), intent(out) :: waterstate_vars

!
!EOP
!
! !OTHER LOCAL VARIABLES:
!

    integer  :: c,fc                               ! indices
    integer  :: ccount
    PetscScalar, pointer :: h2osfc_clm_loc(:)
    PetscErrorCode :: ierr

    associate(&
         h2osfc => waterstate_vars%h2osfc_col & ! Input:  [real(r8) (:)]  surface water (mm)
    )

    call VecGetArrayF90(clm_pf_idata%h2osfc_clm, h2osfc_clm_loc, ierr); CHKERRQ(ierr)
    ccount = 0
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       ccount = ccount + 1
       h2osfc(c) = h2osfc_clm_loc(ccount)
    end do
    call VecRestoreArrayF90(clm_pf_idata%h2osfc_clm, h2osfc_clm_loc, ierr); CHKERRQ(ierr)

    end associate

    end subroutine update_h2osfc_clm_pf

#endif
end module clm_pflotran_interfaceMod

