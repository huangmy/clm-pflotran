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
  use clm_varpar      , only : nlevsoi, lsmlon, lsmlat, nlevgrnd
  use shr_kind_mod    , only: r8 => shr_kind_r8
  use decompMod       , only : get_proc_bounds, get_proc_global
  use domainMod       , only : ldomain

  use ncdio_pio       
  use fileutils       , only : getfil
  use spmdMod         , only : mpicom, MPI_INTEGER, masterproc
  use organicFileMod  , only : organicrd
  use clm_varcon      , only : istice, istdlak, istwet, isturb, istice_mec,  &
       icol_roof, icol_sunwall, icol_shadewall, &
       icol_road_perv, icol_road_imperv, zisoi, zsoi
  use abortutils      , only : endrun

  use clm_pflotran_interface_type
  use pflotran_model_module


  ! !PUBLIC TYPES:
  implicit none
  save

  type(pflotran_model_type),pointer,public :: pflotran_m

  !
  private    ! By default everything is private

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_pf_interface_init  ! Phase one initialization


contains

  subroutine clm_pf_interface_init()

    !
    ! !ARGUMENTS:

    implicit none
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

    real(r8), pointer :: qflx_sink_clmpf(:,:)       !
    real(r8), pointer :: sat_clmpf(:,:)             !
    real(r8), pointer :: hksat_x_clmpf(:,:)         ! hydraulic conductivity in x-dir at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: hksat_y_clmpf(:,:)         ! hydraulic conductivity in y-dir at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: hksat_z_clmpf(:,:)         ! hydraulic conductivity in z-dir at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: sucsat_clmpf(:,:)          ! minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: watsat_clmpf(:,:)          ! volumetric soil water at saturation (porosity) (nlevgrnd)
    real(r8), pointer :: bsw_clmpf(:,:)             ! Clapp and Hornberger "b" (nlevgrnd)
    real(r8), pointer :: topo_clmpf(:)              ! Topogrpahy
    real(r8), pointer :: zisoi_clmpf(:)             ! Depth at soil interfaces
    real(r8), pointer :: zwt_clmpf(:)               ! water table depth (m)

    integer , pointer :: ctype(:)                   ! column type index
    real(r8), pointer :: hksat(:,:)                 ! hydraulic conductivity at saturation (mm H2O /s) (nlevgrnd)
    real(r8), pointer :: sucsat(:,:)                ! minimum soil suction (mm) (nlevgrnd)
    real(r8), pointer :: watsat(:,:)                ! volumetric soil water at saturation (porosity) (nlevgrnd)
    integer , pointer :: cgridcell(:)               ! gridcell index of column
    real(r8), pointer :: wtgcell(:)                 ! weight (relative to gridcell)
    integer , pointer :: ltype(:)                   ! landunit type index
    real(r8), pointer :: h2osoi_vol(:,:)            ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
    real(r8), pointer :: topo(:)                    ! topography

    real(r8), pointer :: zwt(:)                     ! water table depth (m)



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
    integer  :: start(3),count(3)      ! netcdf start/count arrays

    real(r8) :: watsat_tmp, bsw_tmp, sucsat_tmp
    real(r8) :: bd, tkm, bsw2_tmp,psisat_tmp
    real(r8) :: vwcsat_tmp, xksat, hksat_tmp

    integer  :: ier                                ! error status
    character(len=256) :: locfn                    ! local filEname
    character(len= 32) :: subname = 'clm_pflotran_interfaceMod' ! subroutine name
    integer :: mxsoil_color                        ! maximum number of soil color classes

    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon

    logical :: readvar 

    !------------------------------------------------------------------------
    allocate(pflotran_m)

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! Assign local pointers to derived subtypes components (landunit-level)
    ltype           => clm3%g%l%itype

    ! Assign local pointer to derived subtypes components (column-level)
    cgridcell       => clm3%g%l%c%gridcell
    wtgcell         => clm3%g%l%c%wtgcell
    ctype           => clm3%g%l%c%itype
    hksat           => clm3%g%l%c%cps%hksat
    sucsat          => clm3%g%l%c%cps%sucsat
    watsat          => clm3%g%l%c%cps%watsat
    h2osoi_vol      => clm3%g%l%c%cws%h2osoi_vol
    topo            => ldomain%topo
    zwt             => clm3%g%l%c%cws%zwt

    write(iulog,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%          Within clm_pf_interface_init               %%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%                                                     %%'
    write(iulog,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    write(iulog,*) ' '
    write(iulog,*) 'begg=',begg, 'endg=',endg

    allocate(clm_pf_data%qflx_sink(begg:endg,1:nlevsoi))
    allocate(clm_pf_data%sat      (begg:endg,1:nlevsoi))

    allocate(clm_pf_data%hksat_x( begg:endg,1:nlevgrnd))
    allocate(clm_pf_data%hksat_y( begg:endg,1:nlevgrnd))
    allocate(clm_pf_data%hksat_z( begg:endg,1:nlevgrnd))
    allocate(clm_pf_data%sucsat ( begg:endg,1:nlevgrnd))
    allocate(clm_pf_data%watsat ( begg:endg,1:nlevgrnd))
    allocate(clm_pf_data%bsw    ( begg:endg,1:nlevgrnd))
    allocate(clm_pf_data%topo   ( nbeg:nend))
    allocate(clm_pf_data%zisoi  ( 0:nlevgrnd))
    allocate(clm_pf_data%zwt    ( begg:endg))


    qflx_sink_clmpf           => clm_pf_data%qflx_sink
    sat_clmpf                 => clm_pf_data%sat

    hksat_x_clmpf             => clm_pf_data%hksat_x
    hksat_y_clmpf             => clm_pf_data%hksat_y
    hksat_z_clmpf             => clm_pf_data%hksat_z
    sucsat_clmpf              => clm_pf_data%sucsat
    watsat_clmpf              => clm_pf_data%watsat
    bsw_clmpf                 => clm_pf_data%bsw
    topo_clmpf                => clm_pf_data%topo
    zisoi_clmpf               => clm_pf_data%zisoi
    zwt_clmpf                 => clm_pf_data%zwt

    clm_pf_data%begg       = begg
    clm_pf_data%endg       = endg
    clm_pf_data%nbeg       = nbeg
    clm_pf_data%nend       = nend
    clm_pf_data%nlevgrnd   = nlevgrnd
    clm_pf_data%nlevsoi    = nlevsoi
    clm_pf_data%ngrids     = endg - begg + 1


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
    
    ier = pio_inq_varid(ncid, 'mxsoil_color', varid)
    if (ier == PIO_NOERR) then
       ier = pio_inq_varid(ncid, 'mxsoil_color', varid)
       ier = pio_get_var(ncid, varid, mxsoil_color)
    else
       mxsoil_color = 8
    end if
 
    call ncd_io(ncid=ncid, varname='PCT_SAND', flag='read', data=sand3d, dim1name=grlnd, readvar=readvar)
    if(.not. readvar) call endrun( trim(subname)//' ERROR: PCT_SAND NOT on surfadata file' )
    
    call ncd_io(ncid=ncid, varname='PCT_CLAY', flag='read', data=clay3d, dim1name=grlnd,readvar=readvar)
    if(.not. readvar) call endrun( trim(subname)//' ERROR: PCT_CLAY NOT on surfadata file' )
    
    ! --------------------------------------------------------------------
    ! If a organic matter dataset has been specified, read it
    ! --------------------------------------------------------------------

    call organicrd(organic3d)

    ! Initialize the variables to ZERO
    do g = begg, endg
       zwt_clmpf(g)                = 0._r8
       do lev = 1,nlevgrnd
          qflx_sink_clmpf(g,lev)   = 0._r8
          hksat_x_clmpf(g,lev)     = 0._r8
          hksat_y_clmpf(g,lev)     = 0._r8
          hksat_z_clmpf(g,lev)     = 0._r8
          sucsat_clmpf(g,lev)      = 0._r8
          watsat_clmpf(g,lev)      = 0._r8
          bsw_clmpf(g,lev)         = 0._r8
       enddo
    end do


    do c = begc, endc

       ! Set gridcell and landunit indices
       g = cgridcell(c)

       if (ltype(l)==istdlak .or. ltype(l)==istwet .or. ltype(l)==istice .or. ltype(l)==istice_mec) then

          write (iulog,*), 'WARNING: Land Unit type Lake/Wet/Ice/Ice_mec ... within the domain'
          write (iulog,*), 'CLM-CN -- PFLOTRAN does not support this land unit presently'


       else if (ltype(l)==isturb .and. (ctype(c) /= icol_road_perv) .and. (ctype(c) /= icol_road_imperv) )then
          ! Urban Roof, sunwall, shadewall properties set to special value
          write (iulog,*), 'WARNING: Land Unit type is Urban '
          write (iulog,*), 'CLM-CN -- PFLOTRAN does not support this land unit presently'

      else  ! soil columns of both urban and non-urban types

          !write (iulog,*) 'g = ', g
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
                uncon_hksat=uncon_frac/((1._r8-om_frac)/xksat &
                     +((1._r8-perc_frac)*om_frac)/om_hksat)
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


             hksat_x_clmpf(g,lev)  = hksat_x_clmpf(g,lev)  + hksat_tmp     * wtgcell(g)
             hksat_y_clmpf(g,lev)  = hksat_x_clmpf(g,lev)
             hksat_z_clmpf(g,lev)  = hksat_z_clmpf(g,lev)  + hksat(g,lev)  * wtgcell(g)

             sucsat_clmpf(g,lev)   = sucsat_clmpf(g,lev)   + sucsat(g,lev) * wtgcell(g)
             watsat_clmpf(g,lev)   = watsat_clmpf(g,lev)   + watsat(g,lev) * wtgcell(g)
             bsw_clmpf(g,lev)      = bsw_clmpf(g,lev)      + bsw_tmp       * wtgcell(g)

             !write(iulog, *), ' clm-zsoi : ',zsoi(lev)
             !write(iulog, *), ' clm-hksat_x:', hksat_x_clmpf(g,lev)
             !write(iulog, *), ' clm-hksat_y:', hksat_y_clmpf(g,lev)
             !write(iulog, *), ' clm-hksat_z:', hksat_z_clmpf(g,lev)
             !write(iulog, *), ' clm-sucsat:', sucsat_clmpf(g,lev)
             !write(iulog, *), ' clm-watsat:', watsat_clmpf(g,lev)
             !write(iulog, *), ' clm-bsw:', bsw_clmpf(g,lev)
             

          enddo

          zwt_clmpf(g)      = zwt_clmpf(g)      + zwt(g)    * wtgcell(g)
          !write (iulog,*), 'zwt = ', zwt(c)

       endif

    enddo


    do j = 0, nlevgrnd
       zisoi_clmpf(j) = zisoi(j)                         !interface depths [m]
    enddo


    pflotran_m => pflotranModelCreate()
    
    call pflotranModelInitMapping(    pflotran_m )
    call pflotranModelSetSoilProp(    pflotran_m )
    call pflotranModelSetICs(         pflotran_m )
    call pflotranModelStepperRunInit( pflotran_m )

    deallocate(sand3d,clay3d,organic3d)

  end subroutine clm_pf_interface_init

end module clm_pflotran_interfaceMod
