module prep_atm_mod

  use shr_kind_mod,     only: r8 => SHR_KIND_R8 
  use shr_kind_mod,     only: cs => SHR_KIND_CS
  use shr_kind_mod,     only: cl => SHR_KIND_CL
  use shr_sys_mod,      only: shr_sys_abort, shr_sys_flush
  use seq_comm_mct,     only: num_inst_atm, num_inst_ocn, num_inst_ice, num_inst_lnd, num_inst_xao, &
                              num_inst_frc, CPLID, ATMID, logunit
  use seq_comm_mct,     only: seq_comm_getData=>seq_comm_setptrs                               
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_getdata  
  use seq_map_type_mod 
  use seq_map_mod
  use seq_flds_mod
  use t_drv_timers_mod
  use mct_mod
  use perf_mod
  use component_type_mod

  implicit none
  save
  PRIVATE

  !--------------------------------------------------------------------------
  ! Public interfaces
  !--------------------------------------------------------------------------

  public :: prep_atm_init
  public :: prep_atm_mrg

  public :: prep_atm_get_l2x_ax
  public :: prep_atm_get_i2x_ax
  public :: prep_atm_get_o2x_ax

  public :: prep_atm_calc_l2x_ax
  public :: prep_atm_calc_i2x_ax
  public :: prep_atm_calc_o2x_ax

  public :: prep_atm_get_mapper_So2a 
  public :: prep_atm_get_mapper_Fo2a 
  public :: prep_atm_get_mapper_Sl2a 
  public :: prep_atm_get_mapper_Fl2a 
  public :: prep_atm_get_mapper_Si2a 
  public :: prep_atm_get_mapper_Fi2a 

  !--------------------------------------------------------------------------
  ! Private interfaces
  !--------------------------------------------------------------------------

  private :: prep_atm_merge
  private :: getfld

  !--------------------------------------------------------------------------
  ! Private data
  !--------------------------------------------------------------------------

  ! mappers
  type(seq_map), pointer :: mapper_So2a 
  type(seq_map), pointer :: mapper_Sl2a 
  type(seq_map), pointer :: mapper_Si2a 
  type(seq_map), pointer :: mapper_Fo2a           ! needed for seq_frac_init
  type(seq_map), pointer :: mapper_Fl2a           ! needed for seq_frac_init
  type(seq_map), pointer :: mapper_Fi2a           ! needed for seq_frac_init
  
  ! attribute vectors 
  type(mct_aVect), target :: l2x_ax(num_inst_lnd) ! Lnd export, atm grid, cpl pes - allocated in driver
  type(mct_aVect), target :: i2x_ax(num_inst_ice) ! Ice export, atm grid, cpl pes - allocated in driver
  type(mct_aVect), target :: o2x_ax(num_inst_ocn) ! Ocn export, atm grid, cpl pes - allocated in driver
  
  ! other module variables
  integer :: mpicom_CPLID  ! MPI cpl communicator
  logical :: iamroot_CPLID ! .true. => CPLID masterproc
  !================================================================================================

contains

  !================================================================================================

  subroutine prep_atm_init(infodata, &
       atm, ocn, ocn_c2_atm, ice, ice_c2_atm, lnd, lnd_c2_atm)

    !---------------------------------------------------------------
    ! Description
    ! Initialize module attribute vectors and  mappers
    !
    ! Arguments
    type (seq_infodata_type) , intent(inout) :: infodata
    type(component_type)     , intent(in)    :: atm(:)
    type(component_type)     , intent(in)    :: ocn(:)
    logical                  , intent(in)    :: ocn_c2_atm ! .true.  => ocn to atm coupling on
    type(component_type)     , intent(in)    :: ice(:)
    logical                  , intent(in)    :: ice_c2_atm ! .true.  => ice to atm coupling on
    type(component_type)     , intent(in)    :: lnd(:)
    logical                  , intent(in)    :: lnd_c2_atm ! .true.  => lnd to atm coupling on
    !
    ! Local Variables
    integer                          :: lsize_a
    integer                          :: eli, eoi,  eii, eai
    integer                          :: ka,km,k1,k2,k3 ! aVect field indices
    logical                          :: samegrid_ao    ! samegrid atm and ocean
    logical                          :: samegrid_al    ! samegrid atm and land
    logical                          :: esmf_map_flag  ! .true. => use esmf for mapping
    logical                          :: atm_present    ! .true.  => atm is present
    logical                          :: ocn_present    ! .true.  => ocn is present
    logical                          :: ice_present    ! .true.  => ice is present
    logical                          :: lnd_present    ! .true.  => lnd is prsent
    character(CL)                    :: ocn_gnam       ! ocn grid
    character(CL)                    :: atm_gnam       ! atm grid
    character(CL)                    :: lnd_gnam       ! lnd grid
    type(mct_avect), pointer         :: a2x_ax
    type(mct_gsMap), pointer         :: gsMap_ax
    type(mct_gsMap), pointer         :: gsMap_ox
    type(mct_gsMap), pointer         :: gsMap_ix
    type(mct_gsMap), pointer         :: gsMap_lx
    character(*), parameter          :: subname = '(prep_atm_init)'
    character(*), parameter          :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call seq_infodata_getData(infodata, &
         atm_present=atm_present,       &
         ocn_present=ocn_present,       &
         ice_present=ice_present,       &
         lnd_present=lnd_present,       &
         atm_gnam=atm_gnam,             &
         ocn_gnam=ocn_gnam,             &
         lnd_gnam=lnd_gnam,             &
         esmf_map_flag=esmf_map_flag)

    allocate(mapper_So2a)
    allocate(mapper_Sl2a) 
    allocate(mapper_Si2a) 
    allocate(mapper_Fo2a) 
    allocate(mapper_Fl2a) 
    allocate(mapper_Fi2a) 

    if (atm_present) then

       call seq_comm_getData(CPLID, &
            mpicom=mpicom_CPLID, iamroot=iamroot_CPLID)

       a2x_ax => component_get_c2x_cx(atm(1)) 
       lsize_a = mct_aVect_lsize(a2x_ax)

       do eli = 1,num_inst_lnd
          call mct_aVect_init(l2x_ax(eli), rList=seq_flds_l2x_fields, lsize=lsize_a)
          call mct_aVect_zero(l2x_ax(eli))
       end do
       do eoi = 1,num_inst_ocn
          call mct_aVect_init(o2x_ax(eoi), rList=seq_flds_o2x_fields, lsize=lsize_a)
          call mct_aVect_zero(o2x_ax(eoi))
       enddo
       do eii = 1,num_inst_ice
          call mct_aVect_init(i2x_ax(eii), rList=seq_flds_i2x_fields, lsize=lsize_a)
          call mct_aVect_zero(i2x_ax(eii))
       enddo
   
       samegrid_al = .true. 
       samegrid_ao = .true.
       if (trim(atm_gnam) /= trim(lnd_gnam)) samegrid_al = .false.
       if (trim(atm_gnam) /= trim(ocn_gnam)) samegrid_ao = .false.
       
       if (ocn_c2_atm) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_So2a'
          end if
          gsmap_ox => component_get_gsmap_cx(ocn(1)) 
          gsmap_ax => component_get_gsmap_cx(atm(1)) 
          call seq_map_init_rcfile(mapper_So2a, gsmap_ox, gsmap_ax, mpicom_CPLID, &
               'seq_maps.rc','ocn2atm_smapname:','ocn2atm_smaptype:',samegrid_ao, &
               'mapper_So2a initialization',esmf_map_flag)
       end if
       
       ! needed for domain checking
       if (ocn_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fo2a'
          end if
          gsmap_ox => component_get_gsmap_cx(ocn(1)) 
          gsmap_ax => component_get_gsmap_cx(atm(1)) 
          call seq_map_init_rcfile(mapper_Fo2a, gsmap_ox, gsmap_ax, mpicom_CPLID, &
               'seq_maps.rc','ocn2atm_fmapname:','ocn2atm_fmaptype:',samegrid_ao, &
               'mapper_Fo2a initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)
       
       if (ice_c2_atm) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Si2a'
          end if
          gsmap_ix => component_get_gsmap_cx(ice(1)) 
          gsmap_ax => component_get_gsmap_cx(atm(1)) 
          call seq_map_init_rcfile(mapper_Si2a, gsmap_ix, gsmap_ax, mpicom_CPLID, &
               'seq_maps.rc','ice2atm_smapname:','ice2atm_smaptype:',samegrid_ao, &
               'mapper_Si2a initialization',esmf_map_flag)
       end if
       
       ! needed for domain checking
       if (ice_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fi2a'
          end if
          gsmap_ix => component_get_gsmap_cx(ice(1)) 
          gsmap_ax => component_get_gsmap_cx(atm(1)) 
          call seq_map_init_rcfile(mapper_Fi2a, gsmap_ix, gsmap_ax, mpicom_CPLID, &
               'seq_maps.rc','ice2atm_fmapname:','ice2atm_fmaptype:',samegrid_ao, &
               'mapper_Fi2a initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)
       
       ! needed for domain checking
       if (lnd_present) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Fl2a'
          end if
          gsmap_lx => component_get_gsmap_cx(lnd(1)) 
          gsmap_ax => component_get_gsmap_cx(atm(1)) 
          call seq_map_init_rcfile(mapper_Fl2a, gsmap_lx, gsmap_ax,mpicom_CPLID, &
               'seq_maps.rc','lnd2atm_fmapname:','lnd2atm_fmaptype:',samegrid_al, &
               'mapper_Fl2a initialization',esmf_map_flag)
       endif
       call shr_sys_flush(logunit)

       if (lnd_c2_atm) then
          if (iamroot_CPLID) then
             write(logunit,*) ' '
             write(logunit,F00) 'Initializing mapper_Sl2a'
          end if
          gsmap_lx => component_get_gsmap_cx(lnd(1)) 
          gsmap_ax => component_get_gsmap_cx(atm(1)) 
          call seq_map_init_rcfile(mapper_Sl2a, gsmap_lx, gsmap_ax,mpicom_CPLID, &
               'seq_maps.rc','lnd2atm_smapname:','lnd2atm_smaptype:',samegrid_al, &
               'mapper_Sl2a initialization',esmf_map_flag)
       end if

       
    end if

  end subroutine prep_atm_init

  !================================================================================================

  subroutine prep_atm_mrg(infodata, atm, fractions_ax, xao_ax, timer_mrg)

    !---------------------------------------------------------------
    ! Description
    ! Prepare run phase, including running the merge
    !
    ! Arguments
    type(seq_infodata_type) , intent(in)    :: infodata
    type(component_type)    , intent(inout) :: atm(:)
    type(mct_aVect)         , intent(in)    :: fractions_ax(:)
    type(mct_aVect)         , intent(in)    :: xao_ax(:) 
    character(len=*)        , intent(in)    :: timer_mrg
    !
    ! Local Variables
    integer                  :: eli, eoi, eii, exi, efi, eai
    type(mct_avect), pointer :: x2a_ax
    character(*), parameter  :: subname = '(prep_atm_mrg)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer_mrg),barrier=mpicom_CPLID)
    do eai = 1,num_inst_atm
       ! Use fortran mod to address ensembles in merge
       eli = mod((eai-1),num_inst_lnd) + 1
       eoi = mod((eai-1),num_inst_ocn) + 1
       eii = mod((eai-1),num_inst_ice) + 1
       exi = mod((eai-1),num_inst_xao) + 1
       efi = mod((eai-1),num_inst_frc) + 1

       x2a_ax => component_get_x2c_cx(atm(eai)) ! This is actually modifying x2a_ax
       call prep_atm_merge(l2x_ax(eli), o2x_ax(eoi), xao_ax(exi), i2x_ax(eii), &
            fractions_ax(efi), x2a_ax)
    enddo
    call t_drvstopf  (trim(timer_mrg))

  end subroutine prep_atm_mrg

  !================================================================================================

  subroutine prep_atm_merge( l2x_a, o2x_a, xao_a, i2x_a, fractions_a, x2a_a )

    !----------------------------------------------------------------------- 
    !
    ! Arguments
    type(mct_aVect), intent(in)    :: l2x_a
    type(mct_aVect), intent(in)    :: o2x_a
    type(mct_aVect), intent(in)    :: xao_a
    type(mct_aVect), intent(in)    :: i2x_a
    type(mct_aVect), intent(in)    :: fractions_a
    type(mct_aVect), intent(inout) :: x2a_a
    !
    ! Local workspace
    real(r8) :: fracl, fraci, fraco
    integer  :: n,ka,ki,kl,ko,kx,kof,kif,klf
    integer  :: lsize       
    integer  :: index_x2a_Sf_lfrac
    integer  :: index_x2a_Sf_ifrac
    integer  :: index_x2a_Sf_ofrac
    character(CL) :: field_atm   ! string converted to char
    character(CL) :: field_lnd   ! string converted to char
    character(CL) :: field_ice   ! string converted to char
    character(CL) :: field_xao   ! string converted to char
    character(CL) :: field_ocn   ! string converted to char
    character(CL) :: itemc_atm   ! string converted to char
    character(CL) :: itemc_lnd   ! string converted to char
    character(CL) :: itemc_ice   ! string converted to char
    character(CL) :: itemc_xao   ! string converted to char
    character(CL) :: itemc_ocn   ! string converted to char
    logical :: iamroot  
    logical :: first_time = .true.
    logical, pointer, save :: lmerge(:),imerge(:),xmerge(:),omerge(:)
    integer, pointer, save :: lindx(:), iindx(:), oindx(:),xindx(:)
    integer, save          :: naflds, klflds,niflds,noflds,nxflds
    !-----------------------------------------------------------------------
    !
    call seq_comm_setptrs(CPLID, iamroot=iamroot)

    if (first_time) then
          
       naflds = mct_aVect_nRattr(x2a_a)
       klflds = mct_aVect_nRattr(l2x_a)
       niflds = mct_aVect_nRattr(i2x_a)
       noflds = mct_aVect_nRattr(o2x_a)
       nxflds = mct_aVect_nRattr(xao_a)

       allocate(lindx(naflds), lmerge(naflds))
       allocate(iindx(naflds), imerge(naflds))
       allocate(xindx(naflds), xmerge(naflds))
       allocate(oindx(naflds), omerge(naflds))

       lindx(:) = 0
       iindx(:) = 0
       xindx(:) = 0
       oindx(:) = 0
       lmerge(:)  = .true.
       imerge(:)  = .true.
       xmerge(:)  = .true.
       omerge(:)  = .true.

       ! Field naming rules
       ! Only atm states that are Sx_... will be merged
       ! Only fluxes that are F??x_... will be merged 
       ! All fluxes will be multiplied by corresponding component fraction

       do ka = 1,naflds
          call getfld(ka, x2a_a, field_atm, itemc_atm)
          if (field_atm(1:2) == 'PF') then
             cycle ! if flux has first character as P, pass straight through 
          end if
          if (field_atm(1:1) == 'S' .and. field_atm(2:2) /= 'x') then
             cycle ! any state fields that are not Sx_ will just be copied
          end if

          do kl = 1,klflds
             call getfld(kl, l2x_a, field_lnd, itemc_lnd)
             if (trim(itemc_atm) == trim(itemc_lnd)) then
                if ((trim(field_atm) == trim(field_lnd))) then
                   if (field_lnd(1:1) == 'F') lmerge(ka) = .false.
                end if
                lindx(ka) = kl
                exit 
             end if
          end do
          do ki = 1,niflds
             call getfld(ki, i2x_a, field_ice, itemc_ice)
             if (field_ice(1:1) == 'F' .and. field_ice(2:4) == 'ioi') then
                cycle ! ignore all fluxes that are ice/ocn fluxes
             end if
             if (trim(itemc_atm) == trim(itemc_ice)) then
                if ((trim(field_atm) == trim(field_ice))) then
                   if (field_ice(1:1) == 'F') imerge(ka) = .false.
                end if
                iindx(ka) = ki
                exit 
             end if
          end do
          do kx = 1,nxflds
             call getfld(kx, xao_a, field_xao, itemc_xao)
             if (trim(itemc_atm) == trim(itemc_xao)) then
                if ((trim(field_atm) == trim(field_xao))) then
                   if (field_xao(1:1) == 'F') xmerge(ka) = .false.
                end if
                xindx(ka) = kx
                exit 
             end if
          end do
          do ko = 1,noflds
             call getfld(ko, o2x_a, field_ocn, itemc_ocn)
             if (trim(itemc_atm) == trim(itemc_ocn)) then
                if ((trim(field_atm) == trim(field_ocn))) then
                   if (field_ocn(1:1) == 'F') omerge(ka) = .false.
                end if
                oindx(ka) = ko
                exit 
             end if
          end do
          if (lindx(ka) == 0) itemc_lnd = 'unset'
          if (iindx(ka) == 0) itemc_ice = 'unset'
          if (xindx(ka) == 0) itemc_xao = 'unset'
          if (oindx(ka) == 0) itemc_ocn = 'unset'

          if (iamroot) then
             write(logunit,10)trim(itemc_atm),trim(itemc_lnd),&
                  trim(itemc_ice),trim(itemc_xao),trim(itemc_ocn)
10           format(' ',' atm field: ',a15,', lnd merge: ',a15, &
                  ', ice merge: ',a15,', xao merge: ',a15,', ocn merge: ',a15)
             write(logunit, *)'field_atm,lmerge, imerge, xmerge, omerge= ',&
                  trim(field_atm),lmerge(ka),imerge(ka),xmerge(ka),omerge(ka)
         end if
       end do
       first_time = .false.
    end if

    ! Zero attribute vector

    call mct_avect_zero(x2a_a)

    ! Update surface fractions

    kif=mct_aVect_indexRA(fractions_a,"ifrac")
    klf=mct_aVect_indexRA(fractions_a,"lfrac")
    kof=mct_aVect_indexRA(fractions_a,"ofrac")
    lsize = mct_avect_lsize(x2a_a)

    index_x2a_Sf_lfrac = mct_aVect_indexRA(x2a_a,'Sf_lfrac')
    index_x2a_Sf_ifrac = mct_aVect_indexRA(x2a_a,'Sf_ifrac')
    index_x2a_Sf_ofrac = mct_aVect_indexRA(x2a_a,'Sf_ofrac')
    do n = 1,lsize
       x2a_a%rAttr(index_x2a_Sf_lfrac,n) = fractions_a%Rattr(klf,n)
       x2a_a%rAttr(index_x2a_Sf_ifrac,n) = fractions_a%Rattr(kif,n)
       x2a_a%rAttr(index_x2a_Sf_ofrac,n) = fractions_a%Rattr(kof,n)
    end do

    ! Copy attributes that do not need to be merged
    ! These are assumed to have the same name in 
    ! (o2x_a and x2a_a) and in (l2x_a and x2a_a), etc.

    call mct_aVect_copy(aVin=l2x_a, aVout=x2a_a, vector=mct_usevector)
    call mct_aVect_copy(aVin=o2x_a, aVout=x2a_a, vector=mct_usevector)
    call mct_aVect_copy(aVin=i2x_a, aVout=x2a_a, vector=mct_usevector) 
    call mct_aVect_copy(aVin=xao_a, aVout=x2a_a, vector=mct_usevector)

    ! If flux to atm is coming only from the ocean (based on field being in o2x_a) - 
    ! -- then scale by both ocean and ice fraction
    ! If flux to atm is coming only from the land or ice or coupler
    ! -- then do scale by fraction above
    
    do ka = 1,naflds
       do n = 1,lsize
          fracl = fractions_a%Rattr(klf,n)
          fraci = fractions_a%Rattr(kif,n)
          fraco = fractions_a%Rattr(kof,n)
          if (lindx(ka) > 0 .and. fracl > 0._r8) then
             if (lmerge(ka)) then 
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + l2x_a%rAttr(lindx(ka),n) * fracl
             else
                x2a_a%rAttr(ka,n) = l2x_a%rAttr(lindx(ka),n) * fracl
             end if
          end if
          if (iindx(ka) > 0 .and. fraci > 0._r8) then
             if (imerge(ka)) then 
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + i2x_a%rAttr(iindx(ka),n) * fraci
             else
                x2a_a%rAttr(ka,n) = i2x_a%rAttr(iindx(ka),n) * fraci
             end if
          end if
          if (xindx(ka) > 0 .and. fraco > 0._r8) then
             if (xmerge(ka)) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + xao_a%rAttr(xindx(ka),n) * fraco
             else
                x2a_a%rAttr(ka,n) = xao_a%rAttr(xindx(ka),n) * fraco
             end if
          end if
          if (oindx(ka) > 0) then
             if (omerge(ka) .and. fraco > 0._r8) then
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco
             end if
             if (.not. omerge(ka)) then
                !--- NOTE: This IS using the ocean fields and ice fraction !! ---
                x2a_a%rAttr(ka,n) = o2x_a%rAttr(oindx(ka),n) * fraci
                x2a_a%rAttr(ka,n) = x2a_a%rAttr(ka,n) + o2x_a%rAttr(oindx(ka),n) * fraco 
             end if
          end if
       end do
    end do

  end subroutine prep_atm_merge

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

  subroutine prep_atm_calc_o2x_ax(ocn, fractions_ox, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create o2x_ax (note that o2x_ax is a local module variable)
    !
    ! Arguments
    type(component_type) , intent(in) :: ocn(:)
    type(mct_aVect) , intent(in) :: fractions_ox(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eoi, efi 
    type(mct_aVect) , pointer :: o2x_ox
    character(*), parameter   :: subname = '(prep_atm_calc_o2x_ax)'
    character(*), parameter   :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eoi = 1,num_inst_ocn
       efi = mod((eoi-1),num_inst_frc) + 1

       o2x_ox => component_get_c2x_cx(ocn(eoi))
       call seq_map_map(mapper_So2a, o2x_ox, o2x_ax(eoi),&
            fldlist=seq_flds_o2x_states,norm=.true., &
            avwts_s=fractions_ox(efi),avwtsfld_s='ofrac')
       call seq_map_map(mapper_Fo2a, o2x_ox, o2x_ax(eoi),&
            fldlist=seq_flds_o2x_fluxes,norm=.true.)
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_atm_calc_o2x_ax

  !================================================================================================

  subroutine prep_atm_calc_i2x_ax(ice, fractions_ix, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create i2x_ax (note that i2x_ax is a local module variable)
    !
    ! Arguments
    type(component_type) , intent(in) :: ice(:)
    type(mct_aVect) , intent(in) :: fractions_ix(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eii, efi
    type(mct_aVect) , pointer :: i2x_ix
    character(*), parameter   :: subname = '(prep_atm_calc_i2x_ax)'
    character(*), parameter   :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eii = 1,num_inst_ice
       efi = mod((eii-1),num_inst_frc) + 1

       i2x_ix => component_get_c2x_cx(ice(eii))
       call seq_map_map(mapper_Si2a, i2x_ix, i2x_ax(eii), &
            fldlist=seq_flds_i2x_states, &
            avwts_s=fractions_ix(eii), avwtsfld_s='ifrac')
       call seq_map_map(mapper_Fi2a, i2x_ix, i2x_ax(eii), &
            fldlist=seq_flds_i2x_fluxes, &
            avwts_s=fractions_ix(eii), avwtsfld_s='ifrac')
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_atm_calc_i2x_ax

  !================================================================================================

  subroutine prep_atm_calc_l2x_ax(lnd, fractions_lx, timer)
    !---------------------------------------------------------------
    ! Description
    ! Create l2x_ax (note that l2x_ax is a local module variable)
    !
    ! Arguments
    type(component_type) , intent(in) :: lnd(:)
    type(mct_aVect) , intent(in) :: fractions_lx(:)
    character(len=*), intent(in) :: timer
    !
    ! Local Variables
    integer :: eli, efi 
    type(mct_avect), pointer :: l2x_lx
    character(*), parameter  :: subname = '(prep_atm_calc_l2x_ax)'
    character(*), parameter  :: F00 = "('"//subname//" : ', 4A )"
    !---------------------------------------------------------------

    call t_drvstartf (trim(timer),barrier=mpicom_CPLID)
    do eli = 1,num_inst_lnd
       efi = mod((eli-1),num_inst_frc) + 1

       l2x_lx => component_get_c2x_cx(lnd(eli))
       call seq_map_map(mapper_Sl2a, l2x_lx, l2x_ax(eli), &
            fldlist=seq_flds_l2x_states, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
       call seq_map_map(mapper_Fl2a, l2x_lx, l2x_ax(eli), &
            fldlist=seq_flds_l2x_fluxes, norm=.true., &
            avwts_s=fractions_lx(efi), avwtsfld_s='lfrin')
    enddo
    call t_drvstopf  (trim(timer))

  end subroutine prep_atm_calc_l2x_ax

  !================================================================================================

  function prep_atm_get_l2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_l2x_ax(:)
    prep_atm_get_l2x_ax => l2x_ax(:)   
  end function prep_atm_get_l2x_ax

  function prep_atm_get_i2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_i2x_ax(:)
    prep_atm_get_i2x_ax => i2x_ax(:)   
  end function prep_atm_get_i2x_ax

  function prep_atm_get_o2x_ax()
    type(mct_aVect), pointer :: prep_atm_get_o2x_ax(:)
    prep_atm_get_o2x_ax => o2x_ax(:)   
  end function prep_atm_get_o2x_ax

  function prep_atm_get_mapper_So2a()
    type(seq_map), pointer :: prep_atm_get_mapper_So2a
    prep_atm_get_mapper_So2a => mapper_So2a  
  end function prep_atm_get_mapper_So2a

  function prep_atm_get_mapper_Fo2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fo2a
    prep_atm_get_mapper_Fo2a => mapper_Fo2a  
  end function prep_atm_get_mapper_Fo2a

  function prep_atm_get_mapper_Sl2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Sl2a
    prep_atm_get_mapper_Sl2a => mapper_Sl2a  
  end function prep_atm_get_mapper_Sl2a

  function prep_atm_get_mapper_Fl2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fl2a
    prep_atm_get_mapper_Fl2a => mapper_Fl2a  
  end function prep_atm_get_mapper_Fl2a

  function prep_atm_get_mapper_Si2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Si2a
    prep_atm_get_mapper_Si2a => mapper_Si2a  
  end function prep_atm_get_mapper_Si2a

  function prep_atm_get_mapper_Fi2a()
    type(seq_map), pointer :: prep_atm_get_mapper_Fi2a
    prep_atm_get_mapper_Fi2a => mapper_Fi2a  
  end function prep_atm_get_mapper_Fi2a

  !================================================================================================

end module prep_atm_mod

