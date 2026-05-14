!> Calculates the energy requirements of mixing for control case and case with added diffusivity.
module MOM_diapyc_energy_tuning

! This file is part of MOM6. See LICENSE.md for the license.

!! \author By Kiera Lowman, July 2024

use MOM_diag_mediator,      only : diag_ctrl, Time_type, post_data, register_diag_field, register_scalar_field
use MOM_error_handler,      only : MOM_error, FATAL, WARNING, MOM_mesg, is_root_pe
use MOM_error_handler,      only : callTree_enter, callTree_leave, callTree_waypoint, callTree_showQuery
use MOM_file_parser,        only : get_param, log_version, param_file_type, read_param
use MOM_grid,               only : ocean_grid_type
use MOM_unit_scaling,       only : unit_scale_type
use MOM_variables,          only : thermo_var_ptrs
use MOM_verticalGrid,       only : verticalGrid_type
use MOM_EOS,                only : calculate_specific_vol_derivs, calculate_density, EOS_domain
use MOM_spatial_means,      only : global_area_integral
use MOM_string_functions,   only : uppercase, lowercase

use MOM_diapyc_energy_req,  only: diapyc_energy_req_calc, diapyc_energy_req_CS
use MOM_diapyc_energy_req,  only: diapyc_energy_req_init, diapyc_energy_req_end

use MOM_time_manager,       only : time_type, time_type_to_real, operator(//)
use time_manager_mod,       only : length_of_year

use, intrinsic :: ieee_arithmetic

implicit none ; private

#include <MOM_memory.h>

public diapyc_energy_tuning_init, diapyc_energy_tuning_calc, tuning_get_added_diff, diapyc_energy_tuning_end 

! A note on unit descriptions in comments: MOM6 uses units that can be rescaled for dimensional
! consistency testing. These are noted in comments with units like Z, H, L, and T, along with
! their mks counterparts with notation like "a velocity [Z T-1 ~> m s-1]".  If the units
! vary with the Boussinesq approximation, the Boussinesq variant is given first.

!> This control structure holds parameters for the MOM_diapyc_energy_tuning module
type, public :: diapyc_energy_tuning_CS ; private
  logical :: initialized = .false. !< A variable that is here because empty
                               !! structures are not permitted by some compilers.
  type(diag_ctrl), pointer :: diag => NULL() !< A structure that is used to
                               !! regulate the timing of diagnostic output.

  type(time_type), pointer :: Time => NULL() !< Pointer to model time (needed for ramping up added mixing energy)
!  type(time_type), pointer :: Time_init => NULL() !< Pointer to model initialization time (needed for ramping up added mixing energy)
  
  real :: Kd_add        !< The scale of diffusivity that was added on the previous timestep [Z2 T-1 ~> m2 s-1].
  real :: energy_target !< Target change in global power requirements for diapycnal mixing [W].

  real :: ramp_time !< Time period over which to linearly ramp up energy input [years].
  real :: year_offset !< Start year of ensemble member, to subtract when calculating elapsed time.
  
  character(len=20) :: Kd_prof  !< The name of the density-dependent diffusivity profile.

  !< 4 values that define the coordinate potential density range over which a diffusivity scaled by
  !! Kd_add is added [R ~> kg m-3].
  real, dimension(4)  :: rho_range
  real, dimension(4)  :: lat_range !! 4 values that define latitude range where Kd is enhanced [degLat].
  real, dimension(4)  :: lon_range !! 4 values that define longitude range where Kd is enhanced [degLon].
  logical :: use_abs_lat  !< If true, use the absolute value of latitude when setting lat_fn.

  !>@{ Diagnostic IDs
  integer :: id_EnKdTuned=-1,id_Kd_scaling=-1,id_Kd_int_added=-1,id_EnChangeTuned=-1,id_Kd_int_base=-1, &
             id_Kd_lay_added=-1, id_Kd_lay_base=-1
  !>@}

  type(diapyc_energy_req_CS),   pointer :: diapyc_en_rec_CSp     => NULL() !< Control structure for a child module

end type diapyc_energy_tuning_CS

!>@{ Coded parmameters for specifying mixing schemes
!character(len=20), parameter :: SURF_STRING        = "SURF"
!character(len=20), parameter :: THERM_STRING       = "THERM"
!character(len=20), parameter :: MID_STRING         = "MID"
!character(len=20), parameter :: BOT_STRING         = "BOT"
!>@}

contains

!> This subroutine tunes the amplitude of the user-specified profile of added diffusivity such that the globally integrated change
!! in energy required for mixing is equal to a set value. Note that it is assumed that the user seeks an increase in the energy
!! required for diapycnal mixing.
subroutine diapyc_energy_tuning_calc(h_3d, dt, tv, G, GV, US, CS, T_f, S_f, Kd_int, Kd_int_base, Kd_int_added, Kd_lay, &
                                      Kd_lay_base, Kd_lay_added)
  type(ocean_grid_type),          intent(in)    :: G    !< The ocean's grid structure.
  type(verticalGrid_type),        intent(in)    :: GV   !< The ocean's vertical grid structure.
  type(unit_scale_type),          intent(in)    :: US   !< A dimensional unit scaling type
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
                                  intent(in)    :: h_3d !< Layer thickness before entrainment [H ~> m or kg m-2].
  type(thermo_var_ptrs),          intent(inout) :: tv   !< A structure containing pointers to any
                                                        !! available thermodynamic fields.
                                                        !! Absent fields have NULL ptrs.
  real,                           intent(in)    :: dt   !< The amount of time covered by the req_calc call [T ~> s].
  type(diapyc_energy_tuning_CS),     pointer       :: CS   !< This module's control structure.

  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
                                  intent(in)    :: T_f  !< Temperature with massless layers filled in vertically [degC].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
                                  intent(in)    :: S_f  !< Salinity with massless layers filled in vertically [ppt].

  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke+1), &
    intent(inout)   :: Kd_int       !< TOTAL INTERFACE DIFFUSIVITY PROFILE [Z2 T-1 ~> m2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke+1), &
    intent(in)    :: Kd_int_base !<  Base interface diffusivity profile [Z2 T-1 ~> m2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke+1), &
    intent(out)   :: Kd_int_added !< Added interface diffusivity returned after tuning [Z2 T-1 ~> m2 s-1].

  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
          intent(inout)   :: Kd_lay      !< TOTAL LAYER DIFFUSIVITY PROFILE [Z2 T-1 ~> m2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
          intent(in)    :: Kd_lay_base !<  Base layer diffusivity profile [Z2 T-1 ~> m2 s-1].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke), &
          intent(out)   :: Kd_lay_added !< Added layer diffusivity returned after tuning [Z2 T-1 ~> m2 s-1].

  ! Local variables
  real :: energy_change, global_target  ! Actual and target global energy change for tuned added diffusity. [W]
  real :: base_energy, tot_energy  ! Global energy required for base diffusivity and with added diffusivity [W].
  real :: energy_error, error_tol  ! Difference between desired and actual energy change and error tolerance [W].
  real :: Kd_lower, Kd_upper ! bounds to use for bisection method
  real :: PE_err_lower, PE_err_upper   ! Energy change error associated with Kd_lower and Kd_upper [W].
  integer :: num_iter, max_iter       ! Iteration counter and max number of allowed iterations.

  ! --- SAFETY GUARDS FOR CS%Kd_add ----------------------------
  real, parameter :: Kd_add_min = 0.0
  real, parameter :: Kd_add_max = 1.0e2
  real, parameter :: max_step_factor = 10.0

  real :: Kd_prev

  real :: elapsed_years  ! Current year of the simulation.
  
  real, dimension(G%isd:G%ied,G%jsd:G%jed) :: &
      energy_Kd ! 2D array of the column-integrated energy used by diapycnal mixing after tuning
                                          ! [R Z L2 T-3 ~> W m-2].

  real, dimension(GV%ke) :: &
    T0, S0, &   ! T0 & S0 are columns of initial temperatures and salinities [degC] and g/kg.
    h_col       ! h_col is a column of thicknesses h at tracer points [H ~> m or kg m-2].
  real, dimension(GV%ke+1) :: &
    Kd_col       ! A column of diapycnal diffusivities at interfaces [Z2 T-1 ~> m2 s-1].
  
  integer :: i, j, k, is, ie, js, je, nz, itt
!  logical :: may_print

  logical :: showCallTree
  character(len=300)  ::  output_str

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke

  if (.not. associated(CS)) call MOM_error(FATAL, "diapyc_energy_tuning: "// &
         "Module must be initialized before it is used.")
  if (CS%Kd_add == 0) call MOM_error(FATAL, "diapyc_energy_tuning: "// &
          "TUNE_KD_ADD = 0.")
  if (CS%energy_target == 0) call MOM_error(FATAL, "diapyc_energy_tuning: "// &
              "KD_POWER_CHANGE = 0; must specify non-zero value.")

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("Beginning diapyc_energy_tuning_calc()")

  num_iter = 0
  max_iter = 200

!  elapsed_years = time_type_to_real(CS%Time - CS%Time_init) / time_type_to_real(length_of_year())
  elapsed_years = time_type_to_real(CS%Time)/time_type_to_real(length_of_year()) - CS%year_offset

  if (CS%ramp_time /= 0.0) then
    if (elapsed_years < CS%ramp_time) then
      global_target = CS%energy_target * (elapsed_years/CS%ramp_time)
    else
      global_target = CS%energy_target
    endif
  else ! if ramp-up time not supplied, turn on at full-force immediately
    global_target = CS%energy_target
  endif

  error_tol = global_target * 0.001 ! using 1% error tolerance
  Kd_lower = 0.0
  Kd_upper = 0.0

  energy_Kd(:,:) = 0.0 ! previously set to 1.88 for debugging
 
!  write (output_str, '(A, ES14.7)') "From tuning subroutine, max val of initial Kd profile", maxval(Kd_int_added)
!  call MOM_mesg(''//output_str)

  ! Calculating global energy requirement for diapycnal mixing for base diffusivity case.
  do j=js,je ; do i=is,ie
    if (G%mask2dT(i,j) > 0.5) then
      do k=1,nz 
        T0(K) = tv%T(i,j,K) ; S0(K) = tv%S(i,j,K)
        h_col(K) = max(h_3d(i,j,K),GV%H_subroundoff)
      enddo

      do k=1,nz+1
        Kd_col(K) = Kd_int_base(i,j,K)
      enddo

      ! may_print = is_root_PE() .and. (i==ie) .and. (j==je)
      call diapyc_energy_req_calc(h_col, T0, S0, Kd_col, energy_Kd(i,j), dt, tv, G, GV, US, CS=CS%diapyc_en_rec_CSp)
    endif
  enddo ; enddo

  base_energy = global_area_integral(energy_Kd, G, scale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  PE_err_lower = 0
  PE_err_upper = 0
  energy_error = error_tol*2.0

  if (showCallTree) call callTree_waypoint("Entering tuning loop of diapyc_energy_tuning_calc()")

  ! Tuning added diffusivity profile.
  do while ((abs(energy_error) > error_tol) .and. (num_iter < max_iter))

     num_iter = num_iter + 1

     ! Save previous value for safety checks
     Kd_prev = CS%Kd_add

     ! I shouldn't need to add another condition for if Kd_add_scaling greater/less than Kd bounds
     if (num_iter /= 1) then

       if (energy_error > 0.0) then
         Kd_upper = CS%Kd_add
         PE_err_upper = energy_error

         if (Kd_lower == 0.0) then
           CS%Kd_add = 0.1*CS%Kd_add
         else
           CS%Kd_add = Kd_upper - PE_err_upper * (Kd_upper - Kd_lower)/(PE_err_upper-PE_err_lower)
         endif

       else if (energy_error < 0.0) then
         Kd_lower = CS%Kd_add
         PE_err_lower = energy_error

         if (Kd_upper == 0.0) then
           CS%Kd_add = 10.0*CS%Kd_add
         else
           CS%Kd_add = Kd_upper - PE_err_upper * (Kd_upper - Kd_lower)/(PE_err_upper-PE_err_lower)
         endif

       endif
     endif

     ! NaN / Inf -> revert to previous
     if (.not. ieee_is_finite(CS%Kd_add)) then
       CS%Kd_add = Kd_prev
     endif

     ! Enforce physical bounds
     CS%Kd_add = min(max(CS%Kd_add, Kd_add_min), Kd_add_max)

     ! If both bounds exist, force iterate to remain bracketed
     if (Kd_lower > 0.0 .and. Kd_upper > 0.0) then
       CS%Kd_add = min(max(CS%Kd_add, Kd_lower), Kd_upper)
     endif

     ! Limit per-iteration jump size
     if (Kd_prev > 0.0) then
       if (CS%Kd_add > max_step_factor*Kd_prev) CS%Kd_add = max_step_factor*Kd_prev
       if (CS%Kd_add < Kd_prev/max_step_factor) CS%Kd_add = Kd_prev/max_step_factor
     endif
     ! ----------------------------------------------------------------------

!     write (output_str, '(a, es12.5, a, es12.5, a, es12.5)') 'Kd_lower: ', Kd_lower, '  Kd_upper: ', &
!         Kd_upper, '  Kd_add_scaling', CS%Kd_add
!     call MOM_mesg(''//output_str)
    
!     if (Kd_add_scaling == 0) then ! for debugging
!         Kd_add_scaling = 1.0E-05
!     endif

!     if (present(Kd_lay)) then
!       if (.not. present(Kd_lay_base)) then
!         call MOM_error(FATAL, "diapyc_energy_tuning: "// &
!         "Kd_lay_base not provided but Kd_lay given as argument.")
!       else
!         call tuning_get_added_diff(tv, G, GV, US, CS, Kd_int_added, Kd_lay_added=Kd_lay_added, T_f=T_f, S_f=S_f)
!         Kd_lay(:,:,:) = Kd_lay_base(:,:,:)+Kd_lay_added(:,:,:)
!         write (output_str, '(A)') "Called tuning_get_added_diff with Kd_lay_added"
!         call MOM_mesg(''//output_str)
!       endif
!     else 
!       call tuning_get_added_diff(tv, G, GV, US, CS, Kd_int_added, T_f=T_f, S_f=S_f)
!     endif


     call tuning_get_added_diff(tv, G, GV, US, CS, Kd_int_added, Kd_lay_added, T_f=T_f, S_f=S_f)
     Kd_lay(:,:,:) = Kd_lay_base(:,:,:)+Kd_lay_added(:,:,:)
     Kd_int(:,:,:) = Kd_int_base(:,:,:)+Kd_int_added(:,:,:)

     if(.not.associated(tv%Kd_int_tuned)) then
       allocate(tv%Kd_int_tuned(G%isd:G%ied,G%jsd:G%jed,GV%ke+1), source=0.0)
     endif
     tv%Kd_int_tuned(:,:,:) = Kd_int_added(:,:,:)

     energy_Kd(:,:) = 0.0
     
     do j=js,je ; do i=is,ie
       if (G%mask2dT(i,j) > 0.5) then
         do k=1,nz
           T0(K) = tv%T(i,j,K) ; S0(K) = tv%S(i,j,K)
           h_col(K) = max(h_3d(i,j,K),GV%H_subroundoff)
         enddo
         do k=1,nz+1
           Kd_col(K) = Kd_int(i,j,K)
         enddo
         ! may_print = is_root_PE() .and. (i==ie) .and. (j==je)
         call diapyc_energy_req_calc(h_col, T0, S0, Kd_col, energy_Kd(i,j), dt, tv, G, GV, US, CS=CS%diapyc_en_rec_CSp)
       endif
     enddo ; enddo

     tot_energy = global_area_integral(energy_Kd, G, scale=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
     energy_change = tot_energy - base_energy
     energy_error = energy_change - global_target

  enddo

  write (output_str, '(a, i3, a, f0.3, a, es12.4, a, es12.4, a, es12.4)') 'Tot. iter.: ', num_iter, ' Elapsed yrs: ', elapsed_years, &
             '  Final Kd_add: ', CS%Kd_add, '  Energy change: ', energy_change, '  Energy error: ', energy_error
  call MOM_mesg(''//output_str) 
!  if (showCallTree) call callTree_waypoint("Finished ALL iterations of second do loop of diapyc_tuning_calc()")

  if (num_iter == max_iter .and. (abs(energy_error) > error_tol)) then
        call MOM_error(FATAL, "diapyc_energy_tuning_calc: "// &
         "Maximum number of iterations reached.")
  endif

  if (CS%id_EnKdTuned>0) call post_data(CS%id_EnKdTuned, energy_Kd, CS%diag)
  if (CS%id_EnChangeTuned>0) call post_data(CS%id_EnChangeTuned, energy_change, CS%diag)
  if (CS%id_Kd_scaling>0) call post_data(CS%id_Kd_scaling, CS%Kd_add, CS%diag)
  if (CS%id_Kd_int_added>0) call post_data(CS%id_Kd_int_added, Kd_int_added, CS%diag)
  if (CS%id_Kd_int_base>0) call post_data(CS%id_Kd_int_base, Kd_int_base, CS%diag)
  if (CS%id_Kd_lay_added>0) call post_data(CS%id_Kd_lay_added, Kd_lay_added, CS%diag)
  if (CS%id_Kd_lay_base>0) call post_data(CS%id_Kd_lay_base, Kd_lay_base, CS%diag)
 if (showCallTree) call callTree_leave("Leaving diapyc_energy_tuning_calc()")

end subroutine diapyc_energy_tuning_calc

subroutine tuning_get_added_diff(tv, G, GV, US, CS, Kd_int_added, Kd_lay_added, T_f, S_f)
  type(ocean_grid_type),                    intent(in)    :: G   !< The ocean's grid structure.
  type(verticalGrid_type),                  intent(in)    :: GV  !< The ocean's vertical grid structure
!  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(in)    :: h   !< Layer thickness [H ~> m or kg m-2].
  type(thermo_var_ptrs),                    intent(in)    :: tv  !< A structure containing pointers
                                                                 !! to any available thermodynamic
                                                                 !! fields. Absent fields have NULL ptrs.
  type(unit_scale_type),                    intent(in)    :: US  !< A dimensional unit scaling type
  type(diapyc_energy_tuning_CS),            pointer       :: CS  !< This module's control structure.
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)+1), intent(out)   :: Kd_int_added !< The diapycnal
                                                                 !! diffusivity that is being added at
                                                                 !! each interface [Z2 T-1 ~> m2 s-1].
  real, dimension(SZI_(G),SZJ_(G),SZK_(G)), intent(out)   :: Kd_lay_added !< The diapycnal
                                                                 !! diffusivity that is being added at
                                                                 !! each interface [Z2 T-1 ~> m2 s-1].
 real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke),   optional, intent(in)    :: T_f !< Temperature with massless
                                                                  !! layers filled in vertically [degC].
  real, dimension(G%isd:G%ied,G%jsd:G%jed,GV%ke),   optional, intent(in)    :: S_f !< Salinity with massless
                                                                  !! layers filled in vertically [ppt].
  ! Local variables
  real :: Rcv(SZI_(G),SZK_(G)) ! The coordinate density in layers [R ~> kg m-3].
  real :: p_ref(SZI_(G))       ! An array of tv%P_Ref pressures [R L2 T-2 ~> Pa].
  real :: rho_fn      ! The density dependence of the input function, 0-1 [nondim].
  real :: lat_fn      ! The latitude dependence of the input function, 0-1 [nondim].
  real :: lon_fn      ! The longitude dependence of the input function, 0-1 [nondim].
  logical :: use_EOS  ! If true, density is calculated from T & S using an
                      ! equation of state.
  integer, dimension(2) :: EOSdom ! The i-computational domain for the equation of state
  integer :: i, j, k, is, ie, js, je, nz
  integer :: isd, ied, jsd, jed

  logical :: showCallTree

!  real :: kappa_fill  ! diffusivity used to fill massless layers
!  real :: dt_fill     ! timestep used to fill massless layers
!  character(len=200) :: mesg

  is = G%isc ; ie = G%iec ; js = G%jsc ; je = G%jec ; nz = G%ke
  isd = G%isd ; ied = G%ied ; jsd = G%jsd ; jed = G%jed

  showCallTree = callTree_showQuery()
  if (showCallTree) call callTree_enter("Calling tuning_get_added_diff()")

  if (.not.associated(CS)) call MOM_error(FATAL,"diapyc_energy_tuning: "//&
         "Module must be initialized before it is used.")

  use_EOS = associated(tv%eqn_of_state)
  if (.not. use_EOS) return

  Kd_lay_added(:,:,:) = 0.0
  Kd_int_added(:,:,:) = 0.0

  do i=is,ie ; p_ref(i) = tv%P_Ref ; enddo
  EOSdom(:) = EOS_domain(G%HI)
  
  do j=js,je
    if (present(T_f) .and. present(S_f)) then
      do k=1,nz
        call calculate_density(T_f(:,j,k), S_f(:,j,k), p_ref, Rcv(:,k), tv%eqn_of_state, EOSdom)
      enddo
    else
      do k=1,nz
        call calculate_density(tv%T(:,j,k), tv%S(:,j,k), p_ref, Rcv(:,k), tv%eqn_of_state, EOSdom)
      enddo
    endif
    
    do k=1,nz ; do i=is,ie
      rho_fn = val_weights(Rcv(i,k), CS%rho_range)
      if (any(CS%lat_range /= -1.0e9)) then
        if (CS%use_abs_lat) then
          lat_fn = val_weights(abs(G%geoLatT(i,j)), CS%lat_range)
        else
          lat_fn = val_weights(G%geoLatT(i,j), CS%lat_range)
        endif
      else
        lat_fn = 1.0
      endif
      if (any(CS%lon_range /= -1.0e9)) then
        lon_fn = val_weights(G%geoLonT(i,j), CS%lon_range)
      else
        lon_fn = 1.0
      endif
      if (rho_fn * lat_fn * lon_fn > 0.0) then
        Kd_lay_added(i,j,k) = CS%Kd_add * rho_fn * lat_fn * lon_fn
      endif
    enddo ; enddo

    do K=2,nz ; do i=is,ie
      rho_fn = val_weights( 0.5*(Rcv(i,k-1) + Rcv(i,k)), CS%rho_range)
      if (any(CS%lat_range /= -1.0e9)) then
        if (CS%use_abs_lat) then
          lat_fn = val_weights(abs(G%geoLatT(i,j)), CS%lat_range)
        else
          lat_fn = val_weights(G%geoLatT(i,j), CS%lat_range)
        endif
      else
        lat_fn = 1.0
      endif
      if (any(CS%lon_range /= -1.0e9)) then
        lon_fn = val_weights(G%geoLonT(i,j), CS%lon_range)
      else
        lon_fn = 1.0
      endif
      if (rho_fn * lat_fn * lon_fn > 0.0) then
        Kd_int_added(i,j,K) = CS%Kd_add * rho_fn * lat_fn * lon_fn
      endif
    enddo ; enddo
  enddo

  if (showCallTree) call callTree_leave("Leaving tuning_get_added_diff()")

end subroutine tuning_get_added_diff

!> This subroutine checks whether the 4 values of range are in ascending order.
function range_OK(range) result(OK)
  real, dimension(4), intent(in) :: range  !< Four values to check.
  logical                        :: OK     !< Return value.

  OK = ((range(1) <= range(2)) .and. (range(2) <= range(3)) .and. &
        (range(3) <= range(4)))

end function range_OK

!> This subroutine returns a value that goes smoothly from 0 to 1, stays
!! at 1, and then goes smoothly back to 0 at the four values of range.  The
!! transitions are cubic, and have zero first derivatives where the curves
!! hit 0 and 1.  The values in range must be in ascending order, as can be
!! checked by calling range_OK.
function val_weights(val, range) result(ans)
  real,               intent(in) :: val    !< Value for which we need an answer [arbitrary units].
  real, dimension(4), intent(in) :: range  !< Range over which the answer is non-zero [arbitrary units].
  real                           :: ans    !< Return value [nondim].
  ! Local variables
  real :: x   ! A nondimensional number between 0 and 1.

  ans = 0.0
  if ((val > range(1)) .and. (val < range(4))) then
    if (val < range(2)) then
      ! x goes from 0 to 1; ans goes from 0 to 1, with 0 derivatives at the ends.
      x = (val - range(1)) / (range(2) - range(1))
      ans = x**2 * (3.0 - 2.0 * x)
    elseif (val > range(3)) then
      ! x goes from 0 to 1; ans goes from 0 to 1, with 0 derivatives at the ends.
      x = (range(4) - val) / (range(4) - range(3))
      ans = x**2 * (3.0 - 2.0 * x)
    else
      ans = 1.0
    endif
  endif

end function val_weights

!> Initialize parameters and allocate memory associated with the diapycnal energy requirement module.
subroutine diapyc_energy_tuning_init(Time, G, GV, US, param_file, diag, CS)
  type(time_type), target,    intent(in)    :: Time        !< model time
!  type(time_type), target,    intent(in)    :: Time_init   !< model start time, added Feb 5, 2026
  type(ocean_grid_type),      intent(in)    :: G           !< model grid structure
  type(verticalGrid_type),    intent(in)    :: GV          !< ocean vertical grid structure
  type(unit_scale_type),      intent(in)    :: US          !< A dimensional unit scaling type
  type(param_file_type),      intent(in)    :: param_file  !< file to parse for parameter values
  type(diag_ctrl),    target, intent(inout) :: diag        !< structure to regulate diagnostic output
  type(diapyc_energy_tuning_CS), pointer       :: CS          !< module control structure

  integer, save :: init_calls = 0
  
  ! This include declares and sets the variable "version".
# include "version_variable.h"

  character(len=40)  :: mdl = "MOM_diapyc_energy_tuning" ! This module's name.
  character(len=256) :: mesg    ! Message for error messages.

  if (associated(CS)) then   
    call MOM_error(WARNING, "diapyc_energy_tuning_init called with an "// &
                             "associated control structure.")
    return
  else ; allocate(CS) ; endif

  CS%initialized = .true.
  CS%diag => diag
  CS%Time => Time !< added May 28
!  CS%Time_init => Time_init !< added Feb 5, 2026

  ! Read all relevant parameters and write them to the model log.
  call log_version(param_file, mdl, version, "The following parameters are used for energy-tuning the added diffusivity.")

  !< Four successive values that define a range of potential densities over which the extra diffusivity is applied.  The four values
  !! specify the density at which the extra diffusivity starts to increase from 0, "hits its full value, starts to decrease again, and
  !! is back to 0.
  call get_param(param_file, mdl, "TUNING_RHO_V1", CS%rho_range(1), &
                 units="kg m-3", default=-1.0e9, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "TUNING_RHO_V2", CS%rho_range(2), &
                 units="kg m-3", default=-1.0e9, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "TUNING_RHO_V3", CS%rho_range(3), &
                 units="kg m-3", default=-1.0e9, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "TUNING_RHO_V4", CS%rho_range(4), &
                 units="kg m-3", default=-1.0e9, scale=US%kg_m3_to_R)
  call get_param(param_file, mdl, "TUNING_LON_RANGE", CS%lon_range(:), &
                 units="degree", default=-1.0e9)
  call get_param(param_file, mdl, "TUNING_LAT_RANGE", CS%lat_range(:), &
                 units="degree", default=-1.0e9)
  call get_param(param_file, mdl, "TUNING_USE_ABS_LAT", CS%use_abs_lat, &
                 "If true, use the absolute value of latitude when "//&
                 "checking whether a point fits into range of latitudes.", &
                 default=.false.)

  call get_param(param_file, mdl, "TUNE_KD_ADD", CS%Kd_add, &
                 "A user-specified guess for the amplitude of additional diffusivity over a range of "//&
                 "density.", default=0.0, units="m2 s-1", &
                 scale=US%m2_s_to_Z2_T)
  call get_param(param_file, mdl, "KD_POWER_CHANGE", CS%energy_target, &
                 "Target change in power associated with diapycnal mixing.", units="W")
  call get_param(param_file, mdl, "TUNING_RAMP_YRS", CS%ramp_time, &
                 "Time period (in years) over which to linearly ramp up additional diapycnal &
                 mixing energy input.", units="years", default=0.0)
  call get_param(param_file, mdl, "TUNING_OFFSET_YRS", CS%year_offset, &
                 "Ensemble member start time (in years) to use for calculating elapased time.", &
                 units="years", default=0.0)

  if (any( CS%lat_range /= -1.0e9 ) .and. (.not.range_OK(CS%lat_range)) ) then
    write(mesg, '(4(1pe15.6))') CS%lat_range(1:4)
    call MOM_error(FATAL, "diapyc_energy_tuning: bad latitude range: \n  "//&
                    trim(mesg))
  endif
  if (any( CS%lon_range /= -1.0e9 ) .and. (.not.range_OK(CS%lon_range)) ) then
    write(mesg, '(4(1pe15.6))') CS%lon_range(1:4)
    call MOM_error(FATAL, "diapyc_energy_tuning: bad longitude range: \n  "//&
                    trim(mesg))
  endif
  if (.not.range_OK(CS%rho_range)) then
    write(mesg, '(4(1pe15.6))') CS%rho_range(1:4)
    call MOM_error(FATAL, "diapyc_energy_tuning: bad density range: \n  "//&
                    trim(mesg))
  endif

  CS%id_EnKdTuned = register_diag_field('ocean_model', 'energy_Kd_tuned', diag%axesT1, Time, &
                 "Column-integrated rate of energy consumption by diapycnal mixing after tuning (base + added Kd).", units="W m-2", &
                 conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  CS%id_EnChangeTuned = register_scalar_field('ocean_model', 'PE_chg_Kd_tuned', Time, diag, &
          'Globally integrated change in energy due to diapycnal mixing.', units='W') !, conversion=US%RZ3_T3_to_W_m2*US%L_to_Z**2)
  ! unsure of conversion factor for PE_chg_Kd_tuned
  CS%id_Kd_scaling = register_scalar_field('ocean_model', 'Kd_scaling_tuned', Time, diag, &
                 "Energy-tuned scaling of added diapycnal diffusivity profile.", units="m2 s-1")
  CS%id_Kd_int_added = register_diag_field('ocean_model', 'Kd_int_tuned', diag%axesTi, Time, &
          'Added diapycnal diffusivity at interfaces after tuning.', units='m2 s-1', conversion=US%Z2_T_to_m2_s)
  CS%id_Kd_int_base = register_diag_field('ocean_model', 'Kd_int_base', diag%axesTi, Time, &
          'Base diapycnal diffusivity at interfaces without added profile.', units='m2 s-1', conversion=US%Z2_T_to_m2_s)
  CS%id_Kd_lay_added = register_diag_field('ocean_model', 'Kd_lay_tuned', diag%axesTL, Time, &
          'Added layer diapycnal diffusivity after tuning.', units='m2 s-1', conversion=US%Z2_T_to_m2_s)
  CS%id_Kd_lay_base = register_diag_field('ocean_model', 'Kd_lay_base', diag%axesTL, Time, &
          'Base layer diapycnal diffusivity without added profile.', units='m2 s-1', conversion=US%Z2_T_to_m2_s)

  call diapyc_energy_req_init(Time, G, GV, US, param_file, diag, CS%diapyc_en_rec_CSp)

end subroutine diapyc_energy_tuning_init

!> Clean up and deallocate memory associated with the diapycnal energy requirement module.
subroutine diapyc_energy_tuning_end(CS)
  type(diapyc_energy_tuning_CS), pointer :: CS !< Diapycnal energy requirement control structure that
                                            !! will be deallocated in this subroutine.
  if (associated(CS%diapyc_en_rec_CSp)) then
    call diapyc_energy_req_end(CS%diapyc_en_rec_CSp)
  endif

  if (associated(CS)) deallocate(CS)
end subroutine diapyc_energy_tuning_end

end module MOM_diapyc_energy_tuning
