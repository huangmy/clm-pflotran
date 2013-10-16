module seq_map_type_mod

  use shr_kind_mod , only: R8 => SHR_KIND_R8, IN=>SHR_KIND_IN
  use shr_kind_mod , only: CL => SHR_KIND_CL, CX => SHR_KIND_CX
  use shr_mct_mod  , only: shr_mct_sMatPInitnc, shr_mct_queryConfigFile
  use shr_sys_mod
  use shr_const_mod
  use mct_mod
#ifdef USE_ESMF_LIB
  use esmf
  use esmfshr_mod
  use seq_map_esmf
#endif

  type seq_map
    logical                 :: copy_only
    logical                 :: rearrange_only
    logical                 :: esmf_map
    type(mct_rearr)         :: rearr
    type(mct_sMatp)         :: sMatp
    !
    !---- for comparing
    integer(IN)             :: counter   ! indicates which seq_maps this mapper points to
    character(CL)           :: strategy  ! indicates the strategy for this mapper, (copy, rearrange, X, Y)
    character(CX)           :: mapfile   ! indicates the mapping file used
    type(mct_gsMap),pointer :: gsmap_s
    type(mct_gsMap),pointer :: gsmap_d
    !
    !---- for cart3d
    character(CL)           :: cart3d_init
    real(R8), pointer       :: slon_s(:)
    real(R8), pointer       :: clon_s(:)
    real(R8), pointer       :: slat_s(:)
    real(R8), pointer       :: clat_s(:)
    real(R8), pointer       :: slon_d(:)
    real(R8), pointer       :: clon_d(:)
    real(R8), pointer       :: slat_d(:)
    real(R8), pointer       :: clat_d(:)
    integer(IN)             :: mpicom    ! mpicom
    !
#ifdef USE_ESMF_LIB
    !---- import and export States for this mapper object, 
    !---- routehandle is stored in the exp_state for repeated remapping use
    type(ESMF_State)        :: imp_state
    type(ESMF_State)        :: exp_state
#endif
  end type seq_map
  public seq_map

end module seq_map_type_mod
