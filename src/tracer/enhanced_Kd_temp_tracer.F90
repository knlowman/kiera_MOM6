!> A tracer package that mimics salinity
module enhanced_Kd_temp_tracer

! This file is part of MOM6. See LICENSE.md for the license.

use MOM_debugging,     only : hchksum
use MOM_diag_mediator, only : post_data, register_diag_field, safe_alloc_ptr
use MOM_diag_mediator, only : diag_ctrl
use MOM_error_handler, only : MOM_error, FATAL, WARNING, MOM_mesg
use MOM_file_parser, only : get_param, log_param, log_version, param_file_type
use MOM_forcing_type, only : forcing
use MOM_grid, only : ocean_grid_type
use MOM_hor_index, only : hor_index_type
use MOM_io, only : file_exists, read_data, slasher, vardesc, var_desc, query_vardesc
use MOM_open_boundary, only : ocean_OBC_type
use MOM_restart, only : query_initialized, MOM_restart_CS
use MOM_sponge, only : set_up_sponge_field, sponge_CS
use MOM_time_manager, only : time_type
use MOM_tracer_registry, only : register_tracer, tracer_registry_type
use MOM_tracer_diabatic, only : tracer_vertdiff, applyTracerBoundaryFluxesInOut
use MOM_tracer_Z_init, only : tracer_Z_init
use MOM_unit_scaling, only : unit_scale_type
use MOM_variables, only : surface
use MOM_variables, only : thermo_var_ptrs
use MOM_verticalGrid, only : verticalGrid_type

use coupler_types_mod, only : coupler_type_set_data, ind_csurf
use atmos_ocean_fluxes_mod, only : aof_set_coupler_flux

implicit none ; private

#include <MOM_memory.h>

public register_enhanced_Kd_temp_tracer, initialize_enhanced_Kd_temp_tracer
public enhanced_Kd_temp_tracer_column_physics
public enhanced_Kd_temp_tracer_end

!> The control structure for the diffusive temperature change tracer
type, public :: enhanced_Kd_temp_tracer_CS ; private
  type(time_type), pointer :: Time => NULL() !< A pointer to the ocean model's clock.
  type(tracer_registry_type), pointer :: tr_Reg => NULL() !< A pointer to the MOM tracer registry
  real, pointer :: extra_dT(:,:,:) => NULL()   !< The array of diffusive heating tracer used in this
                                               !! subroutine [degC]
  real, pointer :: Kd_int_tuned(:,:,:) => NULL() !< The added diffusivity [m^2/s]
  logical :: enhanced_Kd_temp_may_reinit = .true. !< Hard coding since this should not matter

  integer :: id_encd_Kd_tracer = -1   !< A diagnostic ID

  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                                   !! regulate the timing of diagnostic output.
  type(MOM_restart_CS), pointer :: restart_CSp => NULL() !< A pointer to the restart control structure

  type(vardesc) :: tr_desc !< A description and metadata for the enhanced Kd temperature change tracer
end type enhanced_Kd_temp_tracer_CS

contains

!> Register the enhanced Kd temperature change tracer with MOM6
function register_enhanced_Kd_temp_tracer(HI, GV, param_file, CS, tr_Reg, restart_CS)
  type(hor_index_type),       intent(in) :: HI   !< A horizontal index type structure
  type(verticalGrid_type),    intent(in) :: GV   !< The ocean's vertical grid structure
  type(param_file_type),      intent(in) :: param_file !< A structure to parse for run-time parameters
  type(enhanced_Kd_temp_tracer_CS),  pointer  :: CS !< The control structure returned by a previous
                                               !! call to register_enhanced_Kd_temp_tracer.
  type(tracer_registry_type), pointer    :: tr_Reg !< A pointer that is set to point to the control
                                                  !! structure for the tracer advection and
                                                  !! diffusion module
  type(MOM_restart_CS),       pointer    :: restart_CS !< A pointer to the restart control structure
! This subroutine is used to register tracer fields and subroutines
! to be used with MOM.

  ! Local variables
  character(len=40)  :: mdl = "enhanced_Kd_temp_tracer" ! This module's name.
  character(len=200) :: inputdir ! The directory where the input files are.
  character(len=48)  :: var_name ! The variable's name.
  character(len=3)   :: name_tag ! String for creating identifying enhanced_Kd_temp
! This include declares and sets the variable "version".
#include "version_variable.h"
  real, pointer :: tr_ptr(:,:,:) => NULL()
  logical :: register_enhanced_Kd_temp_tracer
  integer :: isd, ied, jsd, jed, nz, i, j
  isd = HI%isd ; ied = HI%ied ; jsd = HI%jsd ; jed = HI%jed ; nz = GV%ke

  if (associated(CS)) then
    call MOM_error(WARNING, "register_enhanced_Kd_temp_tracer called with an "// &
                             "associated control structure.")
    return
  endif
  allocate(CS)

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "")

  allocate(CS%extra_dT(isd:ied,jsd:jed,nz)) ; CS%extra_dT(:,:,:) = 0.0

  CS%tr_desc = var_desc(trim("enhanced_Kd_temp"), "degC", &
                     "Enhanced Kd temperature change passive tracer", caller=mdl)

  tr_ptr => CS%extra_dT(:,:,:)
  call query_vardesc(CS%tr_desc, name=var_name, caller="register_enhanced_Kd_temp_tracer")
  ! Register the tracer for horizontal advection, diffusion, and restarts.
  call register_tracer(tr_ptr, tr_Reg, param_file, HI, GV, name="enhanced_Kd_temp", &
                       longname="Enhanced Kd temperature change passive tracer", units="degC", &
                       registry_diags=.true., restart_CS=restart_CS, &
                       mandatory=.not.CS%enhanced_Kd_temp_may_reinit)

  CS%tr_Reg => tr_Reg
  CS%restart_CSp => restart_CS
  register_enhanced_Kd_temp_tracer = .true.

end function register_enhanced_Kd_temp_tracer

!> Initialize the enhanced Kd temperature change tracer
subroutine initialize_enhanced_Kd_temp_tracer(restart, day, G, GV, h, diag, OBC, CS, &
                                  sponge_CSp, tv)
  logical,                            intent(in) :: restart !< .true. if the fields have already
                                                         !! been read from a restart file.
  type(time_type),            target, intent(in) :: day  !< Time of the start of the run.
  type(ocean_grid_type),              intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type),            intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                                      intent(in) :: h    !< Layer thicknesses [H ~> m or kg m-2]
  type(diag_ctrl),            target, intent(in) :: diag !< A structure that is used to regulate
                                                         !! diagnostic output.
  type(ocean_OBC_type),               pointer    :: OBC  !< This open boundary condition type specifies
                                                         !! whether, where, and what open boundary
                                                         !! conditions are used.
  type(enhanced_Kd_temp_tracer_CS),        pointer    :: CS !< The control structure returned by a previous
                                                       !! call to register_enhanced_Kd_temp_tracer.
  type(sponge_CS),                    pointer    :: sponge_CSp !< Pointer to the control structure for the sponges.
  type(thermo_var_ptrs),              intent(in) :: tv   !< A structure pointing to various thermodynamic variables
!   This subroutine initializes the tracer fields in CS%extra_dT(:,:,:).

  ! Local variables
  character(len=16) :: name     ! A variable's name in a NetCDF file.
  character(len=72) :: longname ! The long name of that variable.
  character(len=48) :: units    ! The dimensions of the variable.
  character(len=48) :: flux_units ! The units for age tracer fluxes, either
                                ! years m3 s-1 or years kg s-1.
  logical :: OK
  integer :: i, j, k, is, ie, js, je, isd, ied, jsd, jed, nz
  integer :: IsdB, IedB, JsdB, JedB

  if (.not.associated(CS)) return
  if (.not.associated(CS%extra_dT)) return

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed
  IsdB = G%IsdB ; IedB = G%IedB ; JsdB = G%JsdB ; JedB = G%JedB

  CS%Time => day
  CS%diag => diag
  name = "enhanced_Kd_temp"

  call query_vardesc(CS%tr_desc, name=name, caller="initialize_enhanced_Kd_temp_tracer")
  if ((.not.restart) .or. (.not.query_initialized(CS%extra_dT, name, CS%restart_CSp))) then
    do k=1,nz ; do j=jsd,jed ; do i=isd,ied
      CS%extra_dT(i,j,k) = 0.0
    enddo ; enddo ; enddo
  endif

  if (associated(OBC)) then
  ! Steal from updated DOME in the fullness of time.
  endif

  CS%id_encd_Kd_tracer = register_diag_field("ocean_model", "enhanced_Kd_temp_diff", CS%diag%axesTL, &
        day, "Temperature change due to enhanced diapycnal diffusivity.", "degC")

end subroutine initialize_enhanced_Kd_temp_tracer

!> Apply sources, sinks and diapycnal diffusion to the tracers in this package.
subroutine enhanced_Kd_temp_tracer_column_physics(h_old, h_new, ea, eb, fluxes, dt, G, GV, US, CS, tv, debug, &
              evap_CFL_limit, minimum_forcing_depth)
  type(ocean_grid_type),   intent(in) :: G    !< The ocean's grid structure
  type(verticalGrid_type), intent(in) :: GV   !< The ocean's vertical grid structure
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_old !< Layer thickness before entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: h_new !< Layer thickness after entrainment [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: ea   !< an array to which the amount of fluid entrained
                                              !! from the layer above during this call will be
                                              !! added [H ~> m or kg m-2].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), &
                           intent(in) :: eb   !< an array to which the amount of fluid entrained
                                              !! from the layer below during this call will be
                                              !! added [H ~> m or kg m-2].
  type(forcing),           intent(in) :: fluxes !< A structure containing pointers to thermodynamic
                                              !! and tracer forcing fields.  Unused fields have NULL ptrs.
  real,                    intent(in) :: dt   !< The amount of time covered by this call [T ~> s]
  type(unit_scale_type),   intent(in) :: US   !< A dimensional unit scaling type
  type(enhanced_Kd_temp_tracer_CS), pointer :: CS  !< The control structure returned by a previous
                                              !! call to register_enhanced_Kd_temp_tracer.
  type(thermo_var_ptrs),   intent(in) :: tv   !< A structure pointing to various thermodynamic variables
  logical,                 intent(in) :: debug !< If true calculate checksums
  real,          optional, intent(in) :: evap_CFL_limit !< Limit on the fraction of the water that can
                                              !! be fluxed out of the top layer in a timestep [nondim]
  real,          optional, intent(in) :: minimum_forcing_depth !< The smallest depth over which
                                              !! fluxes can be applied [H ~> m or kg m-2]

!   This subroutine applies diapycnal diffusion and any other column
! tracer physics or chemistry to the tracers from this file.

! The arguments to this subroutine are redundant in that
!     h_new(k) = h_old(k) + ea(k) - eb(k-1) + eb(k) - ea(k+1)

  ! Local variables
  real :: year, h_total, scale, htot, Ih_limit
  integer :: secs, days
  integer :: i, j, k, is, ie, js, je, nz, k_max
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)) :: h_work ! Used so that h can be modified
  real :: dTdz, dTdz_1, dTdz_2, num, denom
  logical :: error_flag
  character(len=300)  ::  output_str
  
  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = GV%ke

  if (.not.associated(CS)) return
  if (.not.associated(CS%extra_dT)) return

  call tracer_vertdiff(h_old, ea, eb, dt, CS%extra_dT, G, GV)

  error_flag = .false.

  do k=1,nz ; do j=js,je ; do i=is,ie
    if (k==1) then
      dTdz = (tv%T(i,j,k+1)-tv%T(i,j,k))/(0.5*(max(0.1,h_new(i,j,k+1)+h_new(i,j,k))))
      num = tv%Kd_int_tuned(i,j,K+1)*dTdz - 0.0
    else if (k==nz) then
      dTdz = (tv%T(i,j,k)-tv%T(i,j,k-1))/(0.5*(max(0.1,h_new(i,j,k)+h_new(i,j,k-1))))
      num = 0.0 - tv%Kd_int_tuned(i,j,K)*dTdz
    else
      dTdz_1 = (tv%T(i,j,k+1)-tv%T(i,j,k))/(0.5*(max(0.1,h_new(i,j,k+1)+h_new(i,j,k))))
      dTdz_2 = (tv%T(i,j,k)-tv%T(i,j,k-1))/(0.5*(max(0.1,h_new(i,j,k)+h_new(i,j,k-1))))
      num = tv%Kd_int_tuned(i,j,K+1)*dTdz_1 - tv%Kd_int_tuned(i,j,K)*dTdz_2
    endif
    denom = max(0.1,h_new(i,j,k))
    CS%extra_dT(i,j,k) = CS%extra_dT(i,j,k) + dt*num/denom

    if (denom < 0.01) then
      write (output_str, '(a, es10.3, a, i0)') 'denom: ', denom, 'k = ', k
      call MOM_mesg(trim(output_str))
    endif
    if ((dt*num/denom) > 1E-4) then
      write (output_str, '(a, es10.3, a, es10.3, a, es10.3, a, i0)') 'dt*num/denom: ', dt*num/denom, &
        'num: ', num, 'denom: ', denom, 'k = ', k
      call MOM_mesg(trim(output_str))
!      error_flag = .true.
    endif

  enddo ; enddo ; enddo
!  if (denom < 0.01) then
!    write (output_str, '(a, es10.3)'), 'denom: ', denom
!    call MOM_mesg(trim(output_str))
!  endif
!  if (error_flag) call MOM_error(FATAL, "enhanced_Kd_temp_tracer: "// &
!          "dt*num/denom > 0.01.")

  if (CS%id_encd_Kd_tracer>0) call post_data(CS%id_encd_Kd_tracer, CS%extra_dT, CS%diag)

end subroutine enhanced_Kd_temp_tracer_column_physics

!> Deallocate memory associated with this tracer package
subroutine enhanced_Kd_temp_tracer_end(CS)
  type(enhanced_Kd_temp_tracer_CS), pointer :: CS !< The control structure returned by a previous
                                              !! call to register_enhanced_Kd_temp_tracer.
  integer :: m

  if (associated(CS)) then
    if (associated(CS%extra_dT)) deallocate(CS%extra_dT)
    deallocate(CS)
  endif
end subroutine enhanced_Kd_temp_tracer_end

!> \namespace enhanced_Kd_temp_tracer
!!
!!  By Kiera Lowman, 2025
!!
!!  This file contains the routines necessary to model a passive

end module enhanced_Kd_temp_tracer
