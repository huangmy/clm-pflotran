!#include <misc.h>
!#include <preproc.h>

module clm_pflotran_interface_type

  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  !
  ! !PUBLIC TYPES:
  implicit none

  private

  type, public :: clm_pflotran_type

     ! Parameters
     real(r8), pointer :: hksat_x(:,:)         ! hydraulic conductivity in x-dir at saturation (mm H2O /s) (nlevgrnd)
     real(r8), pointer :: hksat_y(:,:)         ! hydraulic conductivity in y-dir at saturation (mm H2O /s) (nlevgrnd)
     real(r8), pointer :: hksat_z(:,:)         ! hydraulic conductivity in z-dir at saturation (mm H2O /s) (nlevgrnd)
     real(r8), pointer :: sucsat(:,:)          ! minimum soil suction (mm) (nlevgrnd)
     real(r8), pointer :: watsat(:,:)          ! volumetric soil water at saturation (porosity) (nlevgrnd)
     real(r8), pointer :: bsw(:,:)             ! Clapp and Hornberger "b" (nlevgrnd)  
     real(r8), pointer :: topo(:)              ! Topography    
     real(r8), pointer :: zisoi(:)             ! Depth at interface

     ! Initial Conditions
     real(r8), pointer :: zwt(:)               ! Water table depth [m]

     ! From CLM to PFLOTRAN
     real(r8), pointer :: qflx_sink(:,:)       ! [kg/sec], assuming top surf-area = 1 m^2 

     ! From PFLOTRAN to CLM
     real(r8), pointer :: sat(:,:)             ! (Theta/Theta_sat)

     integer :: begg, endg
     integer :: nbeg, nend
     integer :: nlevgrnd, nlevsoi, ngrids
     
  end type clm_pflotran_type

  type(clm_pflotran_type)    , public, target     , save :: clm_pf_data



end module clm_pflotran_interface_type
