module prep_ocn_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8 
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_atm, num_inst_rof, num_inst_ice 
  use seq_comm_mct,     only: num_inst_glc, num_inst_wav, num_inst_ocn
  use seq_comm_mct,     only: num_inst_xao, num_inst_frc
  use seq_comm_mct,     only: num_inst_max
  use seq_comm_mct,     only: CPLID, OCNID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs                               
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata  
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod, only: component_get_x2c_cx, component_get_c2x_cx
  use component_type_mod, only: ocn, atm, ice, rof, wav, glc

  implicit none
  save
  private

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_ocn_init
  public :: prep_ocn_mrg

  public :: prep_ocn_accum
  public :: prep_ocn_accum_avg

  public :: prep_ocn_calc_a2x_ox
  public :: prep_ocn_calc_i2x_ox
  public :: prep_ocn_calc_r2x_ox
  public :: prep_ocn_calc_g2x_ox
  public :: prep_ocn_calc_w2x_ox

  public :: prep_ocn_get_a2x_ox
  public :: prep_ocn_get_r2x_ox
  public :: prep_ocn_get_i2x_ox
  public :: prep_ocn_get_g2x_ox
  public :: prep_ocn_get_w2x_ox

  public :: prep_ocn_get_x2oacc_ox
  public :: prep_ocn_get_x2oacc_ox_cnt

  public :: prep_ocn_get_mapper_Sa2o  
  public :: prep_ocn_get_mapper_Va2o  
  public :: prep_ocn_get_mapper_Fa2o  
  public :: prep_ocn_get_mapper_Fr2o  
  public :: prep_ocn_get_mapper_Rr2o  
  public :: prep_ocn_get_mapper_SFi2o 
  public :: prep_ocn_get_mapper_Rg2o  
  public :: prep_ocn_get_mapper_Sw2o  

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_ocn_merge
  private :: getfld

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_Sa2o
  type(seq_map), pointer :: mapper_Va2o
  type(seq_map), pointer :: mapper_Fa2o
  type(seq_map), pointer :: mapper_Fr2o
  type(seq_map), pointer :: mapper_Rr2o
  type(seq_map), pointer :: mapper_SFi2o
  type(seq_map), pointer :: mapper_Rg2o
  type(seq_map), pointer :: mapper_Sw2o

  ! attribute vectors 
  type(mct_aVect), pointer :: a2x_ox(:) ! Atm export, ocn grid, cpl pes 
  type(mct_aVect), pointer :: r2x_ox(:) ! Rof export, ocn grid, cpl pes 
  type(mct_aVect), pointer :: i2x_ox(:) ! Ice export, ocn grid, cpl pes 
  type(mct_aVect), pointer :: g2x_ox(:) ! Glc export, ocn grid, cpl pes 
  type(mct_aVect), pointer :: w2x_ox(:) ! Wav export, ocn grid, cpl pes 

  type(mct_aVect), target  :: x2o_ox_inst  ! multi instance for averaging

  ! accumulation variables
  type(mct_aVect), pointer :: x2oacc_ox(:)  ! Ocn import, ocn grid, cpl pes 
  integer        , target  :: x2oacc_ox_cnt ! x2oacc_ox: number of time samples accumulated

  ! other module variables
  integer       :: mpicom_CPLID   ! MPI cpl communicator
  logical       :: flood_present  ! .true.  => rof is computing flood
  character(CS) :: vect_map       ! vector mapping type
  logical       :: x2o_average    ! logical for x2o averaging to 1 ocean instance from multi instances
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_ocn_init(infodata, atm_c2_ocn, atm_c2_ice, ice_c2_ocn, rof_c2_ocn, &
       wav_c2_ocn, glc_c2_ocn)
       
    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and all other non-mapping
    ! module variables except for accumulators
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    logical                 , intent(in)    :: atm_c2_ocn ! .true.=>atm to ocn coupling on
    logical                 , intent(in)    :: atm_c2_ice ! .true.=>atm to ice coupling on
    logical                 , intent(in)    :: ice_c2_ocn ! .true.=>ice to ocn coupling on
    logical                 , intent(in)    :: rof_c2_ocn ! .true.=>rof to ocn coupling on
    logical                 , intent(in)    :: wav_c2_ocn ! .true.=>wav to ocn coupling on
    logical                 , intent(in)    :: glc_c2_ocn ! .true.=>glc to ocn coupling on
    !
    ! Local Variables
    logical                  :: esmf_map_flag  ! .true. => use esmf for mapping
    logical                  :: ocn_present    ! .true.  => ocn is present
    logical                  :: atm_present    ! .true.  => atm is present
    logical                  :: ice_present    ! .true.  => ice is present
    logical                  :: iamroot_CPLID  ! .true. => CPLID masterproc
    logical                  :: samegrid_ao    ! samegrid atm and ocean
    logical                  :: samegrid_og    ! samegrid glc and ocean
    logical                  :: samegrid_ow    ! samegrid ocean and wave
    logical                  :: samegrid_ro    ! samegrid runoff and ocean
    integer                  :: atm_nx, atm_ny
    integer                  :: lsize_o
    integer                  :: eli, egi, eri 
    integer                  :: ewi, eai, eii, eoi 
    integer                  :: ka,km,k1,k2,k3 ! aVect field indices
    character(CL)            :: ocn_gnam       ! ocn grid
    character(CL)            :: atm_gnam       ! atm grid
    character(CL)            :: rof_gnam       ! rof grid
    character(CL)            :: wav_gnam       ! wav grid
    character(CL)            :: glc_gnam       ! glc grid
    type(mct_avect), pointer :: o2x_ox
    type(mct_avect), pointer :: x2o_ox
    character(*), parameter  :: subname = '(prep_ocn_init)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    character(*), parameter  :: F01 = "('"//subname//" : ', A, I8 )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata , &
         ocn_present=ocn_present       , &
         atm_present=atm_present       , &
         ice_present=ice_present       , &
         flood_present=flood_present   , &
         vect_map=vect_map             , &
         atm_gnam=atm_gnam             , &
         ocn_gnam=ocn_gnam             , &
         rof_gnam=rof_gnam             , &
         wav_gnam=wav_gnam             , &
         atm_nx=atm_nx                 , &
         atm_ny=atm_ny                 , &
         esmf_map_flag=esmf_map_flag   )

    allocate(mapper_Sa2o)
    allocate(mapper_Va2o)
    allocate(mapper_Fa2o)
    allocate(mapper_Fr2o)
    allocate(mapper_Rr2o)
    allocate(mapper_SFi2o)
    allocate(mapper_Rg2o)
    allocate(mapper_Sw2o)

    if (ocn_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       o2x_ox => component_get_c2x_cx(ocn(1)) 
       x2o_ox => component_get_x2c_cx(ocn(1)) 
       lsize_o = mct_aVect_lsize(o2x_ox)

       ! x2o_average setup logic
       if (num_inst_max == num_inst_ocn) then
          ! standard multi-instance merge
          x2o_average = .false.
       elseif (num_inst_max > 1 .and. num_inst_ocn == 1) then
          ! averaging ocean merge
          x2o_average = .true.
          if (iamroot_CPLID) then
             write(logunit,F01) 'x2o averaging on over instances =',num_inst_max
          end if
          call mct_aVect_init(x2o_ox_inst, x2o_ox, lsize_o)
          call mct_aVect_zero(x2o_ox_inst)
       else
          ! not allowed
          write(logunit,F00) ' ERROR in x2o_average setup logic '
          call shr_sys_abort(subname//' ERROR in x2o_average setup logic')
       endif

       allocate(a2x_ox(num_inst_atm))
       do eai = 1,num_inst_atm
          call mct_aVect_init(a2x_ox(eai), rList=seq_flds_a2x_fields, lsize=lsize_o)
          call mct_aVect_zero(a2x_ox(eai))
       enddo
       allocate(r2x_ox(num_inst_rof))
       do eri = 1,num_inst_rof
          call mct_aVect_init(r2x_ox(eri), rList=seq_flds_r2x_fields, lsize=lsize_o)
          call mct_aVect_zero(r2x_ox(eri))
       enddo
       allocate(g2x_ox(num_inst_glc))
       do egi = 1,num_inst_glc
          call mct_aVect_init(g2x_ox(egi), rList=seq_flds_g2x_fields, lsize=lsize_o)
          call mct_aVect_zero(g2x_ox(egi))
       end do
       allocate(w2x_ox(num_inst_wav))
       do ewi = 1,num_inst_wav
          call mct_aVect_init(w2x_ox(ewi), rList=seq_flds_w2x_fields, lsize=lsize_o)
          call mct_aVect_zero(w2x_ox(ewi))
       enddo
       allocate(i2x_ox(num_inst_ice))
       do eii = 1,num_inst_ice
          call mct_aVect_init(i2x_ox(eii), rList=seq_flds_i2x_fields, lsize=lsize_o)
          call mct_aVect_zero(i2x_ox(eii))
       enddo

       allocate(x2oacc_ox(num_inst_ocn))
       do eoi = 1,num_inst_ocn
          call mct_avect_init(x2oacc_ox(eoi), x2o_ox, lsize_o)
          call mct_aVect_zero(x2oacc_ox(eoi))
       end do
       x2oacc_ox_cnt = 0

       samegrid_ao = .true. 
       samegrid_ro = .true.
       samegrid_ow = .true.
       samegrid_og = .true.
       if (trim(atm_gnam) /= trim(ocn_gnam)) samegrid_ao = .false.
       if (trim(rof_gnam) /= trim(ocn_gnam)) samegrid_ro = .false.
       if (trim(ocn_gnam) /= trim(wav_gnam)) samegrid_ow = .false.
       if (trim(ocn_gnam) /= trim(glc_gnam)) samegrid_og = .false.

       if (atm_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fa2o'
          end if
          call seq_map_init_rcfile(mapper_Fa2o, atm(1), ocn(1), &
               'seq_maps.rc','atm2ocn_fmapname:','atm2ocn_fmaptype:',samegrid_ao, &
               'mapper_Fa2o initialization',esmf_map_flag)
          call shr_sys_flush(logunit)
       end if

       ! atm_c2_ice flag is here because ice and ocn are constrained to be on the same
       ! grid so the atm->ice mapping is set to the atm->ocn mapping to improve performance
       if (atm_c2_ocn .or. atm_c2_ice) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sa2o'
          end if
          call seq_map_init_rcfile(mapper_Sa2o, atm(1), ocn(1), &
               'seq_maps.rc','atm2ocn_smapname:','atm2ocn_smaptype:',samegrid_ao, &
               'mapper_Sa2o initialization',esmf_map_flag)

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Va2o'
          end if
          call seq_map_init_rcfile(mapper_Va2o, atm(1), ocn(1), &
               'seq_maps.rc','atm2ocn_vmapname:','atm2ocn_vmaptype:',samegrid_ao, &
               'mapper_Va2o initialization',esmf_map_flag)

          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Va2o vect'
          end if
          call seq_map_initvect(mapper_Va2o, vect_map, atm(1), ocn(1), string='mapper_Va2o initvect')
       endif
       call shr_sys_flush(logunit)

       ! needed for domain checking
       if (ice_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_SFi2o'
          end if
          call seq_map_init_rearrolap(mapper_SFi2o, ice(1), ocn(1), 'mapper_SFi2o')
       endif
       call shr_sys_flush(logunit)

       if (rof_c2_ocn) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rr2o'
          end if
          call seq_map_init_rcfile(mapper_Rr2o, rof(1), ocn(1), &
               'seq_maps.rc', 'rof2ocn_rmapname:', 'rof2ocn_rmaptype:',samegrid_ro, &
               'mapper_Rr2o initialization',esmf_map_flag)

          if (flood_present) then
             if (iamroot_CPLID) then
                write(logunit,*) ' '
                write(logunit,F00) 'Initializing mapper_Fr2o'
             end if
             call seq_map_init_rcfile( mapper_Fr2o, rof(1), ocn(1), &
                  'seq_maps.rc', 'rof2ocn_fmapname:', 'rof2ocn_fmaptype:',samegrid_ro, &
                  string='mapper_Fr2o initialization', esmf_map=esmf_map_flag)
          endif
       endif
       call shr_sys_flush(logunit)

       if (glc_c2_ocn) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Rg2o'
          end if
          call seq_map_init_rcfile(mapper_Rg2o, glc(1), ocn(1), &
               'seq_maps.rc', 'glc2ocn_rmapname:', 'glc2ocn_rmaptype:',samegrid_og, &
               'mapper_Rg2o initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)

       if (wav_c2_ocn) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sw2o'
          end if
          call seq_map_init_rcfile(mapper_Sw2o, wav(1), ocn(1), &
               'seq_maps.rc', 'wav2ocn_smapname:', 'wav2ocn_smaptype:',samegrid_ow, &
               'mapper_Sw2o initialization')
       endif
       call shr_sys_flush(logunit)

    end if

  end subroutine prep_ocn_init

  !================================================================================================

  subroutine prep_ocn_accum(timer)
    !---------------------------------------------------------------
    ! Description
    ! Accumulate ocn inputs
    ! Form partial sum of tavg ocn inputs (virtual "send" to ocn) 
    ! NOTE: this is done AFTER the call to the merge in prep_ocn_mrg
    !
    ! Arguments
    character(len=*)        , intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi
    type(mct_avect) , pointer   :: x2o_ox
    character(*)    , parameter :: subname = '(prep_ocn_accum)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer), barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       x2o_ox => component_get_x2c_cx(ocn(eoi))

       if (x2oacc_ox_cnt == 0) then
          call mct_avect_copy(x2o_ox, x2oacc_ox(eoi))
       else
          call mct_avect_accum(x2o_ox, x2oacc_ox(eoi))
       endif
    enddo
    x2oacc_ox_cnt = x2oacc_ox_cnt + 1
    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_accum

  !================================================================================================

  subroutine prep_ocn_accum_avg(timer_accum)
    !---------------------------------------------------------------
    ! Description
    ! Finish accumulation ocn inputs
    !
    ! Arguments
    character(len=*), intent(in)    :: timer_accum
    !
    ! Local Variables
    integer :: eoi
    type(mct_avect), pointer :: x2o_ox
    character(*), parameter  :: subname = '(prep_ocn_accum_avg)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_accum), barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       ! temporary formation of average
       if (x2oacc_ox_cnt > 1) then
          call mct_avect_avg(x2oacc_ox(eoi), x2oacc_ox_cnt)
       end if

       ! ***NOTE***THE FOLLOWING ACTUALLY MODIFIES x2o_ox
       x2o_ox   => component_get_x2c_cx(ocn(eoi)) 
       call mct_avect_copy(x2oacc_ox(eoi), x2o_ox)
    enddo
    x2oacc_ox_cnt = 0
    call t_drvstopf (trim(timer_accum))

  end subroutine prep_ocn_accum_avg

  !================================================================================================

  subroutine prep_ocn_mrg(infodata, fractions_ox, xao_ox, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Merge all ocn inputs
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(mct_aVect)         , intent(in)    :: fractions_ox(:)
    type(mct_aVect)         , intent(in)    :: xao_ox(:) ! Atm-ocn fluxes, ocn grid, cpl pes 
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eii, ewi, egi, eoi, eai, eri, exi, efi, emi
    real(r8)                 :: flux_epbalfact ! adjusted precip factor
    type(mct_avect), pointer :: x2o_ox
    integer                  :: cnt
    character(*), parameter  :: subname = '(prep_ocn_mrg)'
    !---------------------------------------------------------------

    call seq_infodata_GetData(infodata, &
         flux_epbalfact=flux_epbalfact)

    call t_drvstartf (trim(timer_mrg), barrier=mpicom_CPLID)

    ! Use emi here for instance averaging capability, num_inst_max = num_inst_ocn normally
    ! if NOT x2o_average, just fill each instance of component_get_x2c_cx(ocn(eoi))
    ! if     x2o_average, then computer merge into x2o_ox_inst and accumulate that to 
    !                     component_get_x2c_cx(ocn(1)) and then average it at the end

    if (x2o_average) then
       x2o_ox   => component_get_x2c_cx(ocn(1)) 
       call mct_aVect_zero(x2o_ox)
    endif

    cnt = 0
    do emi = 1,num_inst_max
       ! Use fortran mod to address ensembles in merge
       eoi = mod((emi-1),num_inst_ocn) + 1
       eai = mod((emi-1),num_inst_atm) + 1
       eii = mod((emi-1),num_inst_ice) + 1
       eri = mod((emi-1),num_inst_rof) + 1
       ewi = mod((emi-1),num_inst_wav) + 1
       egi = mod((emi-1),num_inst_glc) + 1
       exi = mod((emi-1),num_inst_xao) + 1
       efi = mod((emi-1),num_inst_frc) + 1

       if (x2o_average) then
          x2o_ox   => x2o_ox_inst
       else
          x2o_ox   => component_get_x2c_cx(ocn(eoi)) 
       endif

       call prep_ocn_merge( flux_epbalfact, a2x_ox(eai), i2x_ox(eii), r2x_ox(eri),  &
            w2x_ox(ewi), g2x_ox(egi), xao_ox(exi), fractions_ox(efi), x2o_ox )

       if (x2o_average) then
          x2o_ox   => component_get_x2c_cx(ocn(1)) 
          call mct_aVect_accum(x2o_ox_inst, x2o_ox)
          cnt = cnt + 1
       endif
    enddo

    if (x2o_average) then
       x2o_ox   => component_get_x2c_cx(ocn(1)) 
       call mct_avect_avg(x2o_ox,cnt)
    endif

    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_ocn_mrg

  !================================================================================================

  subroutine prep_ocn_merge( flux_epbalfact, a2x_o, i2x_o, r2x_o, w2x_o, g2x_o, xao_o, &
       fractions_o, x2o_o )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    real(r8)       , intent(in)    :: flux_epbalfact
    type(mct_aVect), intent(in)    :: a2x_o
    type(mct_aVect), intent(in)    :: i2x_o
    type(mct_aVect), intent(in)    :: r2x_o
    type(mct_aVect), intent(in)    :: w2x_o
    type(mct_aVect), intent(in)    :: g2x_o
    type(mct_aVect), intent(in)    :: xao_o
    type(mct_aVect), intent(in)    :: fractions_o
    type(mct_aVect), intent(inout) :: x2o_o
    !
    ! Local variables
    integer  :: n,ka,ki,ko,kir,kor
    integer  :: lsize
    integer  :: noflds,naflds,niflds,nxflds
    integer  :: kof,kaf,kif,kxf
    real(r8) :: ifrac,ifracr
    real(r8) :: afrac,afracr
    real(r8) :: frac_sum
    real(r8) :: avsdr, anidr, avsdf, anidf   ! albedos
    real(r8) :: fswabsv, fswabsi             ! sw
    character(CL) :: field_ocn   ! string converted to char
    character(CL) :: field_atm   ! string converted to char
    character(CL) :: field_ice   ! string converted to char
    character(CL) :: field_xao   ! string converted to char
    character(CL) :: itemc_ocn   ! string converted to char
    character(CL) :: itemc_atm   ! string converted to char
    character(CL) :: itemc_ice   ! string converted to char
    character(CL) :: itemc_xao   ! string converted to char
    integer, save :: index_a2x_Faxa_swvdr
    integer, save :: index_a2x_Faxa_swvdf
    integer, save :: index_a2x_Faxa_swndr
    integer, save :: index_a2x_Faxa_swndf
    integer, save :: index_i2x_Fioi_swpen
    integer, save :: index_xao_So_avsdr
    integer, save :: index_xao_So_anidr
    integer, save :: index_xao_So_avsdf
    integer, save :: index_xao_So_anidf
    integer, save :: index_a2x_Faxa_snowc
    integer, save :: index_a2x_Faxa_snowl
    integer, save :: index_a2x_Faxa_rainc
    integer, save :: index_a2x_Faxa_rainl
    integer, save :: index_r2x_Forr_rofl
    integer, save :: index_r2x_Forr_rofi
    integer, save :: index_r2x_Flrr_flood
    integer, save :: index_g2x_Fogg_rofl
    integer, save :: index_g2x_Fogg_rofi
    integer, save :: index_x2o_Foxx_swnet
    integer, save :: index_x2o_Faxa_snow 
    integer, save :: index_x2o_Faxa_rain 
    integer, save :: index_x2o_Faxa_prec  
    integer, save :: index_x2o_Foxx_rofl
    integer, save :: index_x2o_Foxx_rofi
    logical :: iamroot  
    logical, save, pointer :: amerge(:),imerge(:),xmerge(:)
    integer, save, pointer :: aindx(:), iindx(:), oindx(:), xindx(:)
    logical, save :: first_time = .true.
    character(*),parameter :: subName = '(prep_x2o_merge) '
    !----------------------------------------------------------------------- 

    call seq_comm_setptrs(CPLID, iamroot=iamroot)

    noflds = mct_aVect_nRattr(x2o_o)
    naflds = mct_aVect_nRattr(a2x_o)
    niflds = mct_aVect_nRattr(i2x_o)
    nxflds = mct_aVect_nRattr(xao_o)

    if (first_time) then
       index_a2x_Faxa_swvdr     = mct_aVect_indexRA(a2x_o,'Faxa_swvdr')
       index_a2x_Faxa_swvdf     = mct_aVect_indexRA(a2x_o,'Faxa_swvdf')
       index_a2x_Faxa_swndr     = mct_aVect_indexRA(a2x_o,'Faxa_swndr')
       index_a2x_Faxa_swndf     = mct_aVect_indexRA(a2x_o,'Faxa_swndf')
       index_i2x_Fioi_swpen     = mct_aVect_indexRA(i2x_o,'Fioi_swpen') 
       index_xao_So_avsdr       = mct_aVect_indexRA(xao_o,'So_avsdr')
       index_xao_So_anidr       = mct_aVect_indexRA(xao_o,'So_anidr')
       index_xao_So_avsdf       = mct_aVect_indexRA(xao_o,'So_avsdf')
       index_xao_So_anidf       = mct_aVect_indexRA(xao_o,'So_anidf')
       index_x2o_Foxx_swnet     = mct_aVect_indexRA(x2o_o,'Foxx_swnet')

       index_a2x_Faxa_snowc     = mct_aVect_indexRA(a2x_o,'Faxa_snowc')
       index_a2x_Faxa_snowl     = mct_aVect_indexRA(a2x_o,'Faxa_snowl')
       index_a2x_Faxa_rainc     = mct_aVect_indexRA(a2x_o,'Faxa_rainc')
       index_a2x_Faxa_rainl     = mct_aVect_indexRA(a2x_o,'Faxa_rainl')
       index_r2x_Forr_rofl      = mct_aVect_indexRA(r2x_o,'Forr_rofl') 
       index_r2x_Forr_rofi      = mct_aVect_indexRA(r2x_o,'Forr_rofi') 
       index_r2x_Flrr_flood     = mct_aVect_indexRA(r2x_o,'Flrr_flood') 
       index_g2x_Fogg_rofl      = mct_aVect_indexRA(g2x_o,'Fogg_rofl') 
       index_g2x_Fogg_rofi      = mct_aVect_indexRA(g2x_o,'Fogg_rofi') 
       index_x2o_Faxa_snow      = mct_aVect_indexRA(x2o_o,'Faxa_snow')
       index_x2o_Faxa_rain      = mct_aVect_indexRA(x2o_o,'Faxa_rain')
       index_x2o_Faxa_prec      = mct_aVect_indexRA(x2o_o,'Faxa_prec') 
       index_x2o_Foxx_rofl      = mct_aVect_indexRA(x2o_o,'Foxx_rofl') 
       index_x2o_Foxx_rofi      = mct_aVect_indexRA(x2o_o,'Foxx_rofi') 

       ! Compute all other quantities based on standardized naming convention (see below)
       ! Only ocn field states that have the name-prefix Sx_ will be merged
       ! Only field names have the same name-suffix (after the "_") will be merged
       !    (e.g. Si_fldname, Sa_fldname => merged to => Sx_fldname)
       ! All fluxes will be scaled by the corresponding afrac or ifrac 
       !   EXCEPT for 
       !    -- Faxa_snnet, Faxa_snow, Faxa_rain, Faxa_prec (derived)
       ! All i2x_o fluxes that have the name-suffix "Faii" (atm/ice fluxes) will be ignored 
       ! - only ice fluxes that are Fioi_... will be used in the ocean merges

       allocate(aindx(noflds), amerge(noflds))
       allocate(iindx(noflds), imerge(noflds))
       allocate(xindx(noflds), xmerge(noflds))
       aindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       amerge(:) = .true.
       imerge(:) = .true.
       xmerge(:) = .true.

       do kof = 1,noflds
          call getfld(kof, x2o_o, field_ocn, itemc_ocn)
          if (field_ocn(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through 
          end if
          if (field_ocn(1:1) == 'S' .and. field_ocn(2:2) /= 'x') then
             cycle ! ignore all ocn states that do not have a Sx_ prefix 
          end if
          if (trim(field_ocn) == 'Foxx_swnet'.or. &
              trim(field_ocn) == 'Faxa_snow' .or. &
              trim(field_ocn) == 'Faxa_rain' .or. &
              trim(field_ocn) == 'Faxa_prec') then
             cycle ! ignore swnet, snow, rain, prec - treated explicitly above
          end if
!          if (trim(field_ocn(1:5)) == 'Foxx_') then
!             cycle ! ignore runoff fields from land - treated in coupler
!          end if

          do kaf = 1,naflds
             call getfld(kaf, a2x_o, field_atm, itemc_atm)
             if (trim(itemc_ocn) == trim(itemc_atm)) then
                if ((trim(field_ocn) == trim(field_atm))) then
                   if (field_atm(1:1) == 'F') amerge(kof) = .false.
                end if
                aindx(kof) = kaf
                exit
             end if
          end do
          do kif = 1,niflds
             call getfld(kif, i2x_o, field_ice, itemc_ice)
             if (field_ice(1:1) == 'F' .and. field_ice(2:4) == 'aii') then
                cycle ! ignore all i2x_o fluxes that are ice/atm fluxes
             end if
             if (trim(itemc_ocn) == trim(itemc_ice)) then
                if ((trim(field_ocn) == trim(field_ice))) then
                   if (field_ice(1:1) == 'F') imerge(kof) = .false.
                end if
                iindx(kof) = kif
                exit
             end if
          end do
          do kxf = 1,nxflds
             call getfld(kxf, xao_o, field_xao, itemc_xao)
             if (trim(itemc_ocn) == trim(itemc_xao)) then
                if ((trim(field_ocn) == trim(field_xao))) then
                   if (field_xao(1:1) == 'F') xmerge(kof) = .false.
                end if
                xindx(kof) = kxf
                exit
             end if
          end do
          if (aindx(kof) == 0) itemc_atm = 'unset'
          if (iindx(kof) == 0) itemc_ice = 'unset'
          if (xindx(kof) == 0) itemc_xao = 'unset'

          if (iamroot) then
             write(logunit,10)trim(itemc_ocn),&
                  trim(itemc_xao),trim(itemc_ice),trim(itemc_atm)
10           format(' ',' ocn field: ',a15,', xao merge: ',a15, &
                  ', ice merge: ',a15,', atm merge: ',a15)
             write(logunit, *)'field_ocn,kof,imerge,amerge,xmerge= ',&
                  trim(field_ocn),kof,imerge(kof),xmerge(kof),amerge(kof) 
         end if
       end do

       first_time = .false.
    end if
    
    call mct_aVect_zero(x2o_o)

    call mct_aVect_copy(aVin=a2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=i2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=r2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=w2x_o, aVout=x2o_o, vector=mct_usevector)
    call mct_aVect_copy(aVin=xao_o, aVout=x2o_o, vector=mct_usevector)

    ! Compute input ocn state (note that this only applies to non-land portion of gridcell)

    ki  = mct_aVect_indexRa(fractions_o,"ifrac",perrWith=subName)
    ko  = mct_aVect_indexRa(fractions_o,"ofrac",perrWith=subName)
    kir = mct_aVect_indexRa(fractions_o,"ifrad",perrWith=subName)
    kor = mct_aVect_indexRa(fractions_o,"ofrad",perrWith=subName)
    lsize = mct_aVect_lsize(x2o_o)
    do n = 1,lsize

       ifrac = fractions_o%rAttr(ki,n)
       afrac = fractions_o%rAttr(ko,n)
       frac_sum = ifrac + afrac
       if ((frac_sum) /= 0._r8) then
          ifrac = ifrac / (frac_sum)
          afrac = afrac / (frac_sum)
       endif

       ifracr = fractions_o%rAttr(kir,n)
       afracr = fractions_o%rAttr(kor,n)
       frac_sum = ifracr + afracr
       if ((frac_sum) /= 0._r8) then
          ifracr = ifracr / (frac_sum)
          afracr = afracr / (frac_sum)
       endif

       ! Derived: compute net short-wave
       avsdr = xao_o%rAttr(index_xao_So_avsdr,n)  
       anidr = xao_o%rAttr(index_xao_So_anidr,n)  
       avsdf = xao_o%rAttr(index_xao_So_avsdf,n)  
       anidf = xao_o%rAttr(index_xao_So_anidf,n)  
       fswabsv  =  a2x_o%rAttr(index_a2x_Faxa_swvdr,n) * (1.0_R8 - avsdr) &
                 + a2x_o%rAttr(index_a2x_Faxa_swvdf,n) * (1.0_R8 - avsdf)
       fswabsi  =  a2x_o%rAttr(index_a2x_Faxa_swndr,n) * (1.0_R8 - anidr) &
                 + a2x_o%rAttr(index_a2x_Faxa_swndf,n) * (1.0_R8 - anidf)
       x2o_o%rAttr(index_x2o_Foxx_swnet,n) = (fswabsv + fswabsi)                 * afracr + &
                                             i2x_o%rAttr(index_i2x_Fioi_swpen,n) * ifrac

       ! Derived: compute total precipitation - scale total precip and runoff

       x2o_o%rAttr(index_x2o_Faxa_snow ,n) = a2x_o%rAttr(index_a2x_Faxa_snowc,n) * afrac + &
                                             a2x_o%rAttr(index_a2x_Faxa_snowl,n) * afrac 
       x2o_o%rAttr(index_x2o_Faxa_rain ,n) = a2x_o%rAttr(index_a2x_Faxa_rainc,n) * afrac + &
                                             a2x_o%rAttr(index_a2x_Faxa_rainl,n) * afrac

       x2o_o%rAttr(index_x2o_Faxa_snow ,n) = x2o_o%rAttr(index_x2o_Faxa_snow ,n) * flux_epbalfact
       x2o_o%rAttr(index_x2o_Faxa_rain ,n) = x2o_o%rAttr(index_x2o_Faxa_rain ,n) * flux_epbalfact

       x2o_o%rAttr(index_x2o_Faxa_prec ,n) = x2o_o%rAttr(index_x2o_Faxa_rain ,n) + &
                                             x2o_o%rAttr(index_x2o_Faxa_snow ,n)

       x2o_o%rAttr(index_x2o_Foxx_rofl, n) = (r2x_o%rAttr(index_r2x_Forr_rofl , n) + &
                                              r2x_o%rAttr(index_r2x_Flrr_flood, n) + &
                                              g2x_o%rAttr(index_g2x_Fogg_rofl , n)) * flux_epbalfact
       x2o_o%rAttr(index_x2o_Foxx_rofi, n) = (r2x_o%rAttr(index_r2x_Forr_rofi , n) + &
                                              g2x_o%rAttr(index_g2x_Fogg_rofi , n)) * flux_epbalfact
    end do

    do kof = 1,noflds
       do n = 1,lsize
          ifrac = fractions_o%rAttr(ki,n)
          afrac = fractions_o%rAttr(ko,n)
          frac_sum = ifrac + afrac
          if ((frac_sum) /= 0._r8) then
             ifrac = ifrac / (frac_sum)
             afrac = afrac / (frac_sum)
          endif
          if (iindx(kof) > 0) then
             if (imerge(kof)) then 
                x2o_o%rAttr(kof,n) = x2o_o%rAttr(kof,n) + i2x_o%rAttr(iindx(kof),n) * ifrac
             else
                x2o_o%rAttr(kof,n) = i2x_o%rAttr(iindx(kof),n) * ifrac
             end if
          end if
          if (aindx(kof) > 0) then
             if (amerge(kof)) then
                x2o_o%rAttr(kof,n) = x2o_o%rAttr(kof,n) + a2x_o%rAttr(aindx(kof),n) * afrac
             else
                x2o_o%rAttr(kof,n) = a2x_o%rAttr(aindx(kof),n) * afrac
             end if
          end if
          if (xindx(kof) > 0) then
             if (xmerge(kof)) then
                x2o_o%rAttr(kof,n) = x2o_o%rAttr(kof,n) + xao_o%rAttr(xindx(kof),n) * afrac
             else
                x2o_o%rAttr(kof,n) = xao_o%rAttr(xindx(kof),n) * afrac
             end if
          end if
       end do
    end do

  end subroutine prep_ocn_merge

  !================================================================================================

  subroutine getfld(n, av, field, suffix)
    integer         , intent(in)    :: n
    type(mct_aVect) , intent(in)    :: av 
    character(len=*), intent(out)   :: field
    character(len=*), intent(out)   :: suffix

    type(mct_string) :: mstring     ! mct char type

    call mct_aVect_getRList(mstring,n,av)
    field  = mct_string_toChar(mstring)
    suffix = trim(field(scan(field,'_'):))
    call mct_string_clean(mstring)

    if (field(1:1) /= 'S' .and. field(1:1) /= 'F' .and. field(1:2) /= 'PF') then
       write(6,*)'field attribute',trim(field),' must start with S or F or PF' 
       call shr_sys_abort()
    end if

  end subroutine getfld

  !================================================================================================

  subroutine prep_ocn_calc_a2x_ox(timer)
    !---------------------------------------------------------------
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: eai
    type(mct_avect), pointer :: a2x_ax
    character(*), parameter  :: subname = '(prep_ocn_calc_a2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       a2x_ax => component_get_c2x_cx(atm(eai))

       call seq_map_map(mapper_Sa2o, a2x_ax, a2x_ox(eai), fldlist=seq_flds_a2x_states, norm=.true.)

       call seq_map_map(mapper_Fa2o, a2x_ax, a2x_ox(eai), fldlist=seq_flds_a2x_fluxes, norm=.true.)

       !--- tcx the norm should be true below, it's false for bfb backwards compatability
       call seq_map_mapvect(mapper_Va2o, vect_map, a2x_ax, a2x_ox(eai), 'Sa_u', 'Sa_v', norm=.false.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_calc_a2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_i2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_ox (note that i2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*)     , intent(in) :: timer
    !
    ! Local Variables
    integer :: eii
    type(mct_avect), pointer :: i2x_ix
    character(*), parameter  :: subname = '(prep_ocn_calc_i2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eii = 1,num_inst_ice
       i2x_ix => component_get_c2x_cx(ice(eii))
       call seq_map_map(mapper_SFi2o, i2x_ix, i2x_ox(eii), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_ocn_calc_i2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_r2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create r2x_ox (note that r2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eri
    type(mct_avect), pointer :: r2x_rx
    character(*), parameter  :: subname = '(prep_ocn_calc_r2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eri = 1,num_inst_rof
       r2x_rx => component_get_c2x_cx(rof(eri))

       call seq_map_map(mapper_Rr2o, r2x_rx, r2x_ox(eri), &
            fldlist='Forr_rofl:Forr_rofi', norm=.false.)
       if (flood_present) then
          call seq_map_map(mapper_Fr2o, r2x_rx, r2x_ox(eri), &
               fldlist='Flrr_flood', norm=.true.)
       endif
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_ocn_calc_r2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_g2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create g2x_ox (note that g2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: egi
    type(mct_avect), pointer :: g2x_gx
    character(*),  parameter :: subname = '(prep_ocn_calc_g2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do egi = 1,num_inst_rof
       g2x_gx => component_get_c2x_cx(glc(egi))
       call seq_map_map(mapper_Rg2o, g2x_gx, g2x_ox(egi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_ocn_calc_g2x_ox

  !================================================================================================

  subroutine prep_ocn_calc_w2x_ox(timer)
    !---------------------------------------------------------------
    ! Description
    ! Create w2x_ox (note that w2x_ox is a local module variable)
    !
    ! Arguments
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: ewi
    type(mct_avect), pointer :: w2x_wx
    character(*), parameter  :: subname = '(prep_ocn_calc_w2x_ox)'
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do ewi = 1,num_inst_rof
       w2x_wx => component_get_c2x_cx(wav(ewi))
       call seq_map_map(mapper_Sw2o, w2x_wx, w2x_ox(ewi), norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))
  end subroutine prep_ocn_calc_w2x_ox

  !================================================================================================

  function prep_ocn_get_a2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_a2x_ox(:)
    prep_ocn_get_a2x_ox => a2x_ox(:)   
  end function prep_ocn_get_a2x_ox

  function prep_ocn_get_r2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_r2x_ox(:)
    prep_ocn_get_r2x_ox => r2x_ox(:)   
  end function prep_ocn_get_r2x_ox

  function prep_ocn_get_i2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_i2x_ox(:)
    prep_ocn_get_i2x_ox => i2x_ox(:)   
  end function prep_ocn_get_i2x_ox

  function prep_ocn_get_g2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_g2x_ox(:)
    prep_ocn_get_g2x_ox => g2x_ox(:)   
  end function prep_ocn_get_g2x_ox

  function prep_ocn_get_w2x_ox()
    type(mct_aVect), pointer :: prep_ocn_get_w2x_ox(:)
    prep_ocn_get_w2x_ox => w2x_ox(:)   
  end function prep_ocn_get_w2x_ox

  function prep_ocn_get_x2oacc_ox()
    type(mct_aVect), pointer :: prep_ocn_get_x2oacc_ox(:)
    prep_ocn_get_x2oacc_ox => x2oacc_ox(:)   
  end function prep_ocn_get_x2oacc_ox

  function prep_ocn_get_x2oacc_ox_cnt()
    integer, pointer :: prep_ocn_get_x2oacc_ox_cnt
    prep_ocn_get_x2oacc_ox_cnt => x2oacc_ox_cnt
  end function prep_ocn_get_x2oacc_ox_cnt

  function prep_ocn_get_mapper_Sa2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Sa2o
    prep_ocn_get_mapper_Sa2o => mapper_Sa2o  
  end function prep_ocn_get_mapper_Sa2o

  function prep_ocn_get_mapper_Va2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Va2o
    prep_ocn_get_mapper_Va2o => mapper_Va2o  
  end function prep_ocn_get_mapper_Va2o

  function prep_ocn_get_mapper_Fa2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Fa2o
    prep_ocn_get_mapper_Fa2o => mapper_Fa2o  
  end function prep_ocn_get_mapper_Fa2o

  function prep_ocn_get_mapper_Fr2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Fr2o
    prep_ocn_get_mapper_Fr2o => mapper_Fr2o  
  end function prep_ocn_get_mapper_Fr2o

  function prep_ocn_get_mapper_Rr2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Rr2o
    prep_ocn_get_mapper_Rr2o => mapper_Rr2o  
  end function prep_ocn_get_mapper_Rr2o

  function prep_ocn_get_mapper_SFi2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_SFi2o
    prep_ocn_get_mapper_SFi2o => mapper_SFi2o  
  end function prep_ocn_get_mapper_SFi2o

  function prep_ocn_get_mapper_Rg2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Rg2o
    prep_ocn_get_mapper_Rg2o => mapper_Rg2o  
  end function prep_ocn_get_mapper_Rg2o

  function prep_ocn_get_mapper_Sw2o()
    type(seq_map), pointer :: prep_ocn_get_mapper_Sw2o
    prep_ocn_get_mapper_Sw2o => mapper_Sw2o  
  end function prep_ocn_get_mapper_Sw2o

end module prep_ocn_mod

