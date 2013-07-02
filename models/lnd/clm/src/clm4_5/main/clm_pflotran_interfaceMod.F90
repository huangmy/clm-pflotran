#ifdef CLM_PFLOTRAN

module clm_pflotran_interfaceMod

  !-----------------------------------------------------------------------
  !BOP
  !
  ! !MODULE: clm_pflotran_interfaceMod
  !
  ! !DESCRIPTION:
  ! Performs
  !
  ! !USES:
  use clmtype
  use clm_varctl      , only : iulog, fsurdat, scmlon, scmlat, single_column
  use decompMod       , only : get_proc_bounds, get_proc_global
  use clm_varpar      , only : nlevsoi, nlevgrnd
  use shr_kind_mod    , only: r8 => shr_kind_r8
  use decompMod       , only : ldecomp
  use domainMod       , only : ldomain

  use fileutils       , only : getfil
  use spmdMod         , only : mpicom, masterproc
  use organicFileMod  , only : organicrd
  use clm_varcon      , only : istice, istdlak, istwet, isturb, istice_mec,  &
                               icol_roof, icol_sunwall, icol_shadewall, &
                               icol_road_perv, icol_road_imperv, zisoi, zsoi, &
                               istsoil, denice, denh2o
  use abortutils      , only : endrun

  !use clm_pflotran_interface_type
  use clm_pflotran_interface_data
  use pflotran_model_module
  use ncdio_pio


  ! !PUBLIC TYPES:
  implicit none


  save

  !type(clm_pflotran_interface_data_type),pointer,public   :: clm_pf_idata
  type(pflotran_model_type),pointer,public                :: pflotran_m

  !
  private    ! By default everything is private

  character(len=256), public :: pflotran_prefix = ''

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_pf_interface_init  ! Phase one initialization
  public :: clm_pf_readnl


contains

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
    use spmdMod       , only : masterproc, mpicom
    use fileutils     , only : getavu, relavu, opnfil
    use clm_nlUtilsMod, only : find_nlgroup_name
    use shr_mpi_mod   , only : shr_mpi_bcast
    use abortutils    , only : endrun
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
       call find_nlgroup_name(unitn, 'clm_pflotran_inparm', status=ierr)
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
! !IROUTINE: clm_pf_interface_init
!
! !INTERFACE:
  subroutine clm_pf_interface_init()
!
! !DESCRIPTION:
! initialize the pflotran iterface
!
! !USES:
    use clm_varctl, only : use_pflotran
    !
    ! !ARGUMENTS:

    implicit none

#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscviewer.h"
#include "definitions.h"

    !
    ! !REVISION HISTORY:
    ! Created by Gautam Bisht
    !
    !EOP
    !
    ! LOCAL VARAIBLES:

    integer  :: begp, endp       ! per-proc beginning and ending pft indices
    integer  :: begc, endc       ! per-proc beginning and ending column indices
    integer  :: begl, endl       ! per-proc beginning and ending landunit indices
    integer  :: begg, endg       ! per-proc gridcell ending gridcell indices
    integer  :: nbeg, nend
    integer  :: numg             ! total number of gridcells across all processors
    integer  :: numl             ! total number of landunits across all processors
    integer  :: numc             ! total number of columns across all processors
    integer  :: nump             ! total number of pfts across all processors
    integer  :: n,g,l,c,p,lev,j  ! indices
    integer  :: gcount

    integer , pointer :: ctype(:)                   ! column type index
    real(r8), pointer :: hksat(:,:)                 ! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: sucsat(:,:)                ! minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: watsat(:,:)                ! volumetric soil water at saturation (porosity) (nlevgrnd)
    integer , pointer :: cgridcell(:)               ! gridcell index of column
    integer , pointer :: clandunit(:)               ! landunit index of column
    real(r8), pointer :: wtgcell(:)                 ! weight (relative to gridcell)
    integer , pointer :: ltype(:)                   ! landunit type index
    real(r8), pointer :: h2osoi_vol(:,:)            ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: h2osoi_liq(:,:)            ! liquid water (kg/m2)
    real(r8), pointer :: h2osoi_ice(:,:)  ! ice lens (kg/m2)
    real(r8), pointer :: t_soisno(:,:)              ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: topo(:)                    ! topography

    real(r8), pointer :: zwt(:)                     ! water table depth (m)

    real(r8), pointer :: hksat_x_loc(:)         ! hydraulic conductivity in x-dir at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: hksat_y_loc(:)         ! hydraulic conductivity in y-dir at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: hksat_z_loc(:)         ! hydraulic conductivity in z-dir at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: sucsat_loc(:)          ! minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: watsat_loc(:)          ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: bsw_loc(:)             ! Clapp and Hornberger "b" (nlevgrnd)
    real(r8), pointer :: zwt_2d_loc(:)          ! water table depth (m)
    real(r8), pointer :: topo_2d_loc(:)         ! Topogrpahy
    real(r8), pointer :: dz(:,:)                ! layer thickness (m)
    integer :: index
    
    real(r8), pointer :: latdeg(:)             ! latitude (radians)
    real(r8), pointer :: londeg(:)             ! longitude (radians)

    !
    ! From iniTimeConst.F90
    !
    type(file_desc_t)  :: ncid   ! netcdf id
    real(r8) :: clay,sand        ! temporaries

    real(r8),pointer :: arrayl(:)      ! generic global array
    integer ,pointer :: irrayg(:)      ! generic global array
    integer ,pointer :: soic2d(:)      ! read in - soil color
    real(r8),pointer :: sand3d(:,:)    ! read in - soil texture: percent sand
    real(r8),pointer :: clay3d(:,:)    ! read in - soil texture: percent clay
    real(r8),pointer :: organic3d(:,:) ! read in - organic matter: kg/m3
    real(r8),pointer :: gti(:)         ! read in - fmax
    integer  :: varid                  ! netCDF id's
    integer  :: ret

    real(r8) :: om_frac                ! organic matter fraction
    real(r8) :: om_watsat    = 0.9_r8  ! porosity of organic soil
    real(r8) :: om_hksat     = 0.1_r8  ! saturated hydraulic conductivity of organic soil [mm/s]
    real(r8) :: om_tkm       = 0.25_r8 ! thermal conductivity of organic soil (Farouki, 1986) [W/m/K]
    real(r8) :: om_sucsat    = 10.3_r8 ! saturated suction for organic matter (Letts, 2000)
    real(r8) :: om_csol      = 2.5_r8  ! heat capacity of peat soil *10^6 (J/K m3) (Farouki, 1986)
    real(r8) :: om_tkd       = 0.05_r8 ! thermal conductivity of dry organic soil (Farouki, 1981)
    real(r8) :: om_b         = 2.7_r8  ! Clapp Hornberger paramater for oragnic soil (Letts, 2000)
    real(r8) :: organic_max  = 130._r8 ! organic matter (kg/m3) where soil is assumed to act like peat
    real(r8) :: csol_bedrock = 2.0e6_r8 ! vol. heat capacity of granite/sandstone  J/(m3 K)(Shabbir, 2000)
    real(r8) :: pc           = 0.5_r8   ! percolation threshold
    real(r8) :: pcbeta       = 0.139_r8 ! percolation exponent
    real(r8) :: perc_frac               ! "percolating" fraction of organic soil
    real(r8) :: perc_norm               ! normalize to 1 when 100% organic soil
    real(r8) :: uncon_hksat             ! series conductivity of mineral/organic soil
    real(r8) :: uncon_frac              ! fraction of "unconnected" soil
    integer  :: start(3),count(3)       ! netcdf start/count arrays

    real(r8) :: watsat_tmp, bsw_tmp, sucsat_tmp, press_tmp
    real(r8) :: bd, tkm, bsw2_tmp,psisat_tmp
    real(r8) :: vwcsat_tmp, xksat, hksat_tmp

    integer  :: ier                                ! error status
    character(len=256) :: locfn                    ! local filEname
    character(len= 32) :: subname = 'clm_pflotran_interfaceMod' ! subroutine name
    integer :: mxsoil_color                        ! maximum number of soil color classes

    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon

    logical :: readvar
    logical , pointer :: lakpoi(:)      ! true => landunit is a lake point

    integer, pointer :: clm_cell_ids_nindex(:)
    integer, pointer :: clm_surf_cell_ids_nindex(:)
    integer :: clm_npts
    integer :: clm_surf_npts

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

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! Assign local pointers to derived subtypes components (landunit-level)
    ltype           => lun%itype

    ! Assign local pointer to derived subtypes components (column-level)
    clandunit       => col%landunit
    cgridcell       => col%gridcell
    wtgcell         => col%wtgcell
    ctype           => col%itype
    hksat           => cps%hksat
    sucsat          => cps%sucsat
    watsat          => cps%watsat
    h2osoi_vol      => cws%h2osoi_vol
    h2osoi_liq      => cws%h2osoi_liq
    h2osoi_ice      => cws%h2osoi_ice
    t_soisno        => ces%t_soisno
    topo            => ldomain%topo
    zwt             => cws%zwt
    latdeg          => grc%latdeg
    londeg          => grc%londeg
    lakpoi          => lun%lakpoi
    dz              => cps%dz

    !------------------------------------------------------------------------
    allocate(pflotran_m)

    ! Create PFLOTRAN model
    pflotran_m => pflotranModelCreate(mpicom, pflotran_prefix)

    ! Initialize PETSc vector for data transfer between CLM and PFLOTRAN
    call CLMPFLOTRANIDataInit()

    ! Compute number of cells in CLM domain.
    ! Assumption-1: One column per CLM grid cell.
    ! Assumption-2: nz = nlevsoi and nz /= nlevgrnd. Need to add a flag in input
    !               file to differenticate if PFLOTRAN grid is 'nlevsoi' or
    !               'nlevgrnd' deep.
    clm_npts = (endg-begg+1)*nlevsoi
    clm_surf_npts = (endg-begg+1)
    allocate(clm_cell_ids_nindex( 1:clm_npts))
    allocate(clm_surf_cell_ids_nindex(1:clm_surf_npts))

    ! Save cell IDs of CLM grid
    clm_npts = 0
    clm_surf_npts = 0
    do g = begg, endg
       do j = 1,nlevsoi
          clm_npts = clm_npts + 1
          clm_cell_ids_nindex(clm_npts) = (ldecomp%gdc2glo(g)-1)*nlevsoi + j - 1
       enddo
       clm_surf_npts=clm_surf_npts + 1
       clm_surf_cell_ids_nindex(clm_surf_npts)=(ldecomp%gdc2glo(g)-1)*nlevsoi
    enddo

    ! CLM: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_3d = clm_npts
    clm_pf_idata%ngclm_3d = clm_npts

    ! CLM: Surface of subsurface domain (local and ghosted cells)
    clm_pf_idata%nlclm_surf_3d = (endg-begg+1)
    clm_pf_idata%ngclm_surf_3d = (endg-begg+1)
    ! For CLM: Same as surface of subsurface domain
    clm_pf_idata%nlclm_2d = clm_surf_npts
    clm_pf_idata%ngclm_2d = clm_surf_npts

    ! PFLOTRAN: Subsurface domain (local and ghosted cells)
    clm_pf_idata%nlpf_3d = pflotran_m%realization%patch%grid%nlmax
    clm_pf_idata%ngpf_3d = pflotran_m%realization%patch%grid%ngmax

    ! PFLOTRAN: Surface of subsurface domain (local and ghosted cells)
    if(pflotran_m%option%iflowmode == TH_MODE) then
      clm_pf_idata%nlpf_surf_3d = pflotranModelNSurfCells3DDomain(pflotran_m)
      clm_pf_idata%ngpf_surf_3d = pflotranModelNSurfCells3DDomain(pflotran_m)
    else
      clm_pf_idata%nlpf_surf_3d = 0
      clm_pf_idata%ngpf_surf_3d = 0
    endif
    
    ! PFLOTRAN: Surface domain (local and ghosted cells)
#ifdef SURFACE_FLOW
    if(pflotran_m%option%nsurfflowdof > 0) then
      clm_pf_idata%nlpf_2d = pflotran_m%simulation%surf_realization%patch%grid%nlmax
      clm_pf_idata%ngpf_2d = pflotran_m%simulation%surf_realization%patch%grid%ngmax
    else
      clm_pf_idata%nlpf_2d = 0
      clm_pf_idata%ngpf_2d = 0
    endif
#else
    clm_pf_idata%nlpf_2d = 0
    clm_pf_idata%ngpf_2d = 0
#endif

    ! Allocate vectors for data transfer between CLM and PFLOTRAN.
    call CLMPFLOTRANIDataCreateVec(MPI_COMM_WORLD)

    ! Initialize maps for transferring data between CLM and PFLOTRAN.
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, CLM2PF_FLUX_MAP_ID)
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, CLM2PF_SOIL_MAP_ID)
    call pflotranModelInitMapping(pflotran_m, clm_cell_ids_nindex, &
                                  clm_npts, PF2CLM_FLUX_MAP_ID)

    if(pflotran_m%option%iflowmode==TH_MODE) then
      call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, CLM2PF_GFLUX_MAP_ID)
    endif

#ifdef SURFACE_FLOW
    if(pflotran_m%option%nsurfflowdof > 0) then
      call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, PF2CLM_SURF_MAP_ID)
      call pflotranModelInitMapping(pflotran_m, clm_surf_cell_ids_nindex, &
                                    clm_surf_npts, CLM2PF_RFLUX_MAP_ID)
    endif
#endif

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

    allocate(sand3d(begg:endg,nlevsoi),clay3d(begg:endg,nlevsoi))
    allocate(organic3d(begg:endg,nlevsoi))

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

    ! --------------------------------------------------------------------
    ! If a organic matter dataset has been specified, read it
    ! --------------------------------------------------------------------

    call organicrd(organic3d)

    gcount = 0_r8 ! assumption that only 1 soil-column per grid cell
    do c = begc, endc

      ! Set gridcell and landunit indices
      g = cgridcell(c)
      l = clandunit(c)
      gcount = g - begg

      if (ltype(l)==istdlak .or. ltype(l)==istwet .or. ltype(l)==istice .or. ltype(l)==istice_mec) then
        write (iulog,*), 'WARNING: Land Unit type Lake/Wet/Ice/Ice_mec ... within the domain'
        write (iulog,*), 'CLM-CN -- PFLOTRAN does not support this land unit presently'
      else if (ltype(l)==isturb .and. (ctype(c) /= icol_road_perv) .and. (ctype(c) /= icol_road_imperv) )then
        ! Urban Roof, sunwall, shadewall properties set to special value
        write (iulog,*), 'WARNING: Land Unit type is Urban '
        write (iulog,*), 'CLM-CN -- PFLOTRAN does not support this land unit presently'
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
          if (ltype(l)==isturb) then
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

          press_tmp = 101325.0_r8 - 998.2_r8*9.81_r8*(zwt(g) - zsoi(lev))
          press_tmp = 101325.0_r8 - 998.2_r8*9.81_r8*(2.0_r8 - zsoi(lev))

          if (lev <= nlevsoi) then
            hksat_x_clm_loc(gcount*nlevsoi + lev ) = hksat_x_clm_loc(gcount*nlevsoi + lev ) + hksat_tmp*wtgcell(g)
            hksat_y_clm_loc(gcount*nlevsoi + lev ) = hksat_y_clm_loc(gcount*nlevsoi + lev ) + hksat_tmp*wtgcell(g)
            hksat_z_clm_loc(gcount*nlevsoi + lev ) = hksat_z_clm_loc(gcount*nlevsoi + lev ) + hksat(g,lev)*wtgcell(g)
            sucsat_clm_loc( gcount*nlevsoi + lev ) = sucsat_clm_loc( gcount*nlevsoi + lev ) + sucsat(g,lev)*wtgcell(g)
            watsat_clm_loc( gcount*nlevsoi + lev ) = watsat_clm_loc( gcount*nlevsoi + lev ) + watsat(g,lev)*wtgcell(g)
            bsw_clm_loc(    gcount*nlevsoi + lev ) = bsw_clm_loc(    gcount*nlevsoi + lev ) + bsw_tmp*wtgcell(g)
            press_clm_loc(  gcount*nlevsoi + lev ) = press_clm_loc(  gcount*nlevsoi + lev ) + press_tmp*wtgcell(g)

            if(pflotran_m%option%myrank.eq.-1) then
              write(*,'(I4,9F15.10)'), gcount*nlevsoi + lev, &
                                      sucsat_clm_loc( gcount*nlevsoi + lev ), &
                                      bsw_clm_loc(    gcount*nlevsoi + lev ), &
                                      watsat_clm_loc( gcount*nlevsoi + lev ), &
                                      hksat_x_clm_loc(gcount*nlevsoi + lev ), &
                                      hksat_y_clm_loc(gcount*nlevsoi + lev ), &
                                      hksat_z_clm_loc(gcount*nlevsoi + lev ), &
                                      998.2_r8*9.81_r8*(zwt(g) - zsoi(lev)), zwt(g), zsoi(lev)
              write(*,*),gcount*nlevsoi+lev, press_clm_loc(gcount*nlevsoi+lev)
            endif
          endif
        enddo 
      endif
    enddo ! do c = begc, endc

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
    call pflotranModelGetUpdatedStates(pflotran_m)

    ! Initialize soil temperature
    if(pflotran_m%option%iflowmode==TH_MODE) then
      call VecGetArrayF90(clm_pf_idata%temp_clm, temp_clm_loc, ierr)
      do c = begc,endc
        l = clandunit(c)
        if (.not. lakpoi(l)) then  !not lake
         g = cgridcell(c)
         gcount = g - begg
         do j = 1, nlevsoi
            t_soisno(c,j) = temp_clm_loc(gcount*nlevsoi+j)+273.15_r8
         enddo
         t_soisno(c,nlevsoi+1:nlevgrnd) = t_soisno(c,nlevsoi)
        endif
      enddo
      call VecRestoreArrayF90(clm_pf_idata%temp_clm, temp_clm_loc, ierr)
    endif

    ! Initialize soil moisture
    call VecGetArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr)
    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    do c = begc,endc
      l = clandunit(c)
      if (ltype(l) == istsoil) then
        g = cgridcell(c)
        gcount = g - begg
        do j = 1, nlevsoi
          h2osoi_liq(c,j) = sat_clm_loc(gcount*nlevsoi+j)*dz(c,j)*1.e3_r8
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/dz(c,j)/denh2o + &
                            h2osoi_ice(c,j)/dz(c,j)/denice
          h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
        enddo
      endif
    enddo
    call VecGetArrayF90(clm_pf_idata%sat_clm, sat_clm_loc, ierr)
    call VecGetArrayF90(clm_pf_idata%watsat_clm, watsat_clm_loc, ierr)

    deallocate(sand3d,clay3d,organic3d)

  end subroutine clm_pf_interface_init

end module clm_pflotran_interfaceMod

#endif
