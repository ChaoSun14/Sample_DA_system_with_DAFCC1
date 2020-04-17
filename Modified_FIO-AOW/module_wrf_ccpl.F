MODULE module_wrf_ccpl
 
USE module_domain_type
USE module_domain, ONLY : domain
USE module_timing
USE CCPL_interface_mod
USE module_wrf_error
type, private :: wrf_state_var
   integer               :: comp_id
   integer               :: parent_comp_id
   integer               :: time_step
   integer               :: grid_H2D_id
   integer               :: decomp_id
   integer               :: num_local_cells
   integer               :: num_global_cells
   integer,  allocatable :: local_cell_global_indexes(:)
   real,     allocatable :: center_lats(:), center_lons(:)
   integer,  allocatable :: mask(:)
   logical               :: initialization
   real, allocatable :: U10_buf(:,:)
   real, allocatable :: V10_buf(:,:)
   real, allocatable :: RHOA_buf(:,:)
   real, allocatable :: UST_buf(:,:)
   real, allocatable :: T2_buf(:,:)
   real, allocatable :: Q2_buf(:,:)
   real, allocatable :: QFX_buf(:,:)
   real, allocatable :: LH_buf(:,:)
   real, allocatable :: SH_buf(:,:)
   real, allocatable :: RAIN_buf(:,:)
   real, allocatable :: SWDOWN_buf(:,:)
   real, allocatable :: LWDOWN_buf(:,:)
   real, allocatable :: MSLP_buf(:,:)  

   real, allocatable :: SST_buf(:,:)
   real, allocatable :: SSS_buf(:,:)
   real, allocatable :: SSU_buf(:,:)
   real, allocatable :: SSV_buf(:,:)
   
   real, allocatable :: whitecap_fraction_buf(:,:)
   real, allocatable :: u_orbital_srf_buf(:,:)
   real, allocatable :: v_orbital_srf_buf(:,:)
   real, allocatable :: u_stokes_srf_buf(:,:)
   real, allocatable :: v_stokes_srf_buf(:,:)
   real, allocatable :: break_volume_buf(:,:)
   real, allocatable :: tau_in_buf(:,:)
   real, allocatable :: significant_wave_height_buf(:,:)
   real, allocatable :: peak_wave_length_buf(:,:)
   real, allocatable :: peak_wave_period_buf(:,:)

end type wrf_state_var

!! CCPL
#ifdef CCPL_DA
type, private :: wrf_da_var
   integer               :: ensemble_comp_id,ensemble_comm
   integer               :: wrf_frame_id,parent_comp_id,parent_comm
   integer               :: wrf_ensemble_comp_id,wrf_ensemble_comm
   integer               :: comp_id,comm
   integer               :: grid_H2D_id,grid_vert_id1,grid_vert_id2,grid_regional_time_id,grid_3d_id1,grid_3d_id2
   integer               :: decomp_id,time_step,da_timer_id
   integer               :: ensemble_nums,ensemble_mem_id,nprocs,nlat,nlon,nsig,nsig_soil
   integer               :: gsi_ens_member_inst_id,gsi_ens_ave_inst_id,enkf_inst_id
   integer               :: da_num_local_cells
   integer,  allocatable :: local_cells_global_index_stagger(:)
   integer,  allocatable :: mask(:)
   integer,  allocatable :: da_fields_id(:)
   integer,  allocatable :: control_vars(:)
   integer,  allocatable :: comp_or_grid_ids(:)
   integer,  allocatable :: decomp_ids(:)
   integer,  allocatable :: field_inst_ids(:)
   integer,  allocatable :: timers_id(:)
   integer, dimension(:),  pointer :: regional_time=>null()
   real,                   pointer :: ccpl_pt=>null()
   real, dimension(:),     pointer :: aeta1=>null(),aeta2=>null(),eta1=>null(),eta2=>null()
   real, dimension(:,:),   pointer :: ccpl_dx_mc=>null(),ccpl_dy_mc=>null(),ccpl_xlat=>null(),ccpl_xlong=>null(),ccpl_psfc0=>null(),ccpl_phb=>null(),ccpl_landmask=>null(),ccpl_seaice=>null(),ccpl_sst=>null(),ccpl_ivgtyp=>null(),ccpl_isltyp=>null(),ccpl_vegfra=>null(),ccpl_snow=>null(),ccpl_u10=>null(),ccpl_v10=>null(),ccpl_smois=>null(),ccpl_tslb=>null(),ccpl_tsk=>null(),ccpl_mu=>null(),ccpl_mub=>null()
    real, dimension(:,:,:), pointer :: ccpl_u=>null(),ccpl_v=>null(),ccpl_t=>null(),ccpl_q=>null(),ccpl_w=>null(),ccpl_ph=>null()

end type wrf_da_var
#endif
!! CCPL
   type(wrf_state_var), public :: wrf_state_variables(100)
   integer                     :: wrf_frame_id,local_communicator
   character(len=800), private :: message
   logical                     :: wave_coupled,ocean_coupled

   integer                     :: ensemble_comp_id,ensemble_comm,ensemble_nums,ensemble_mem_id 
   character(256)              :: ensemble_mem_id_str
!! CCPL
#ifdef CCPL_DA
   type(wrf_da_var), public    :: wrf_da_variables(100)
   integer                     :: set_parent_id,wrf_ensemble_comp_id,wrf_ensemble_comm
#endif
!! CCPL
contains 

!! CCPL
SUBROUTINE register_component(comm)
implicit none
    integer, intent(inout) :: comm

    call CCPL_get_ensemble_info(ensemble_nums,ensemble_mem_id)
    write(ensemble_mem_id_str,*) ensemble_mem_id
    ensemble_comm = CCPL_NULL_COMM
    ensemble_comp_id = CCPL_register_component(-1, "WMP_member"//trim(adjustl(ensemble_mem_id_str)),"active_coupled_system", ensemble_comm, change_dir=.false., annotation = "register ensemble member "//trim(adjustl(ensemble_mem_id_str))//" of wrf-masnum-pom to C-Coupler")
    comm = CCPL_NULL_COMM
    wrf_frame_id = CCPL_register_component(ensemble_comp_id, "wrf","atm", comm, change_dir=.false., annotation = "register WRF FRAME to C-Coupler")
    call CCPL_change_to_ensemble_dir(wrf_frame_id,ensemble_mem_id)
    local_communicator = comm
END SUBROUTINE register_component

!SUBROUTINE register_wrf_frame_coupling_configuration
!implicit none
!     integer             :: parent_comp_id
!     parent_comp_id     = -1
!     local_communicator = CCPL_NULL_COMM
!     wrf_frame_id       = CCPL_register_component(parent_comp_id, "wrf", "atm", local_communicator,.false.,change_dir=.true., annotation= "register WRF FRAME to c-coupler")
!     write(message,*) "zhaobiao, local_communicator:", local_communicator
!     CALL wrf_message(message)

!END SUBROUTINE register_wrf_frame_coupling_configuration
!! CCPL

SUBROUTINE register_component_coupling_configuration (wrf_grid)
implicit none
INCLUDE "mpif.h"
  type(domain), pointer :: wrf_grid
  integer               :: parent_comp_id
  integer               :: wrf_grid_id
  integer               :: comp_id
  character(len=1024)   :: annotation
  character(len=20)     :: nestid
  character(len=80)     :: comp_name
  integer               :: grid_H2D_id
  integer               :: decomp_id
  integer               :: field_mark_ocn, field_mark_wave, field_mark_atm 
  integer               :: field_id_sss, field_id_sst, field_id_ssu, field_id_ssv
  integer               :: field_id_ustokes, field_id_vstokes, field_id_uorb, field_id_vorb, field_id_whitecap, field_id_breakrate, field_id_hwave, field_id_pwave, field_id_lwavep, field_id_tau_in
  integer               :: field_id_u10, field_id_v10, field_id_ust, field_id_rhoa, field_id_qfx, field_id_sh, field_id_lh, field_id_t2, field_id_q2, field_id_rain, field_id_mslp, field_id_swdown, field_id_lwdown
  integer               :: timer_id, import_interface_id, export_interface_id
  integer, allocatable  :: fields_id(:)
  integer               :: num_comps, individual_or_family(100)
  character(len=1024)   :: comps_full_names(12)
!! CCPL
#ifdef CCPL_DA
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  integer               :: grid_vert_id1,grid_vert_id2,grid_regional_time_id,grid_3d_id1,grid_3d_id2
  integer               :: da_field_id_regional_time,da_field_id_pt,da_field_id_aeta1,da_field_id_aeta2,da_field_id_eta1,da_field_id_eta2
  integer               :: da_field_id_dx_mc,da_field_id_dy_mc,da_field_id_xlat,da_field_id_xlong,da_field_id_psfc0,da_field_id_phb,da_field_id_mu,da_field_id_mub,da_field_id_landmask
  integer               :: da_field_id_seaice,da_field_id_sst,da_field_id_ivgtyp,da_field_id_isltyp,da_field_id_vegfra
  integer               :: da_field_id_snow,da_field_id_u10,da_field_id_v10,da_field_id_smois,da_field_id_tslb,da_field_id_tsk
  integer               :: da_field_id_u,da_field_id_v,da_field_id_t,da_field_id_q,da_field_id_w,da_field_id_ph
  integer               :: da_field_mark,myrank,nprocs,MyError,da_timer_id
  integer, allocatable  :: da_fields_id(:)
  real                  :: regional_time_coord_values(6)
  character(len=1024)   :: ens_comp_full_name,ens_comps_full_names(1)
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#endif
!! CCPL
  field_mark_ocn  = 0
  field_mark_wave = 1
  field_mark_atm  = 2

!==============================================Register the component model using CCPL_register_component ================== ===================================================================
  IF(wrf_grid%id .eq. 1)THEN
     write(nestid,"(I2)") wrf_grid%id
     annotation        =   "component "//"WRF_TOP_d0"//trim(adjustl(nestid))//" start registration"
     parent_comp_id    =   wrf_frame_id
     comp_name         =   "WRF_TOP_d0"//trim(adjustl(nestid))
  ELSE
    write(nestid,"(I2)") wrf_grid%id
    annotation        =   "component "//"WRF_NEST_d0"//trim(adjustl(nestid))//" start registration"
    parent_comp_id    =   wrf_state_variables(wrf_grid%parent_id)%comp_id      
    comp_name         =   "WRF_NEST_d0"//trim(adjustl(nestid))   
  END IF
  !if (.false.) then
!! CCPL
#ifdef CCPL_DA
  comp_id     =   CCPL_register_component(parent_comp_id, comp_name, "atm", local_communicator,change_dir=.false., annotation=annotation)
#else  
  comp_id     =   CCPL_register_component(parent_comp_id, comp_name, "atm", local_communicator,change_dir=.true., annotation=annotation)
#endif
!! CCPL

  write(message,*) trim(comp_name)," comp_id: ", comp_id, " domain:", wrf_grid%id, " parent comp_id: ",parent_comp_id
  write(*,*) trim(comp_name)," comp_id: ", comp_id, " domain:", wrf_grid%id, " parent comp_id: ",parent_comp_id
  CALL wrf_debug(100,message)
  wrf_grid_id = wrf_grid%id
  wrf_state_variables(wrf_grid_id)%comp_id        = comp_id
  wrf_state_variables(wrf_grid_id)%time_step      = wrf_grid%dt
  wrf_state_variables(wrf_grid_id)%parent_comp_id      = parent_comp_id
  wrf_state_variables(wrf_grid_id)%initialization = .true.
  CALL CCPL_set_normal_time_step(comp_id, wrf_state_variables(wrf_grid_id)%time_step)

!! CCPL
#ifdef CCPL_DA
!-------------------------------------------------------------------------------------------------------------------------
  wrf_da_variables(wrf_grid_id)%time_step            = wrf_grid%dt
  set_parent_id=CCPL_get_parent_id_of_component(ensemble_comp_id,"get parent id of ensemble_comp")
  ens_comp_full_name="WMP_member*@wrf@"//trim(adjustl(comp_name))
  ens_comps_full_names(1)=ens_comp_full_name  
  wrf_ensemble_comp_id=CCPL_register_ensemble_set_of_components("wrf_ensemble_set",comp_id,1,ens_comps_full_names,ensemble_nums,"register wrf ensemble set to C-Coupler")
  wrf_ensemble_comm=CCPL_get_component_comm(wrf_ensemble_comp_id,"get wrf_ensemble_comm")
  CALL mpi_comm_rank (local_communicator, myrank, MyError)
  CALL mpi_comm_size (local_communicator, nprocs, MyError)
  wrf_da_variables(wrf_grid_id)%nlon                 = wrf_grid%ed31-1         !nlon
  wrf_da_variables(wrf_grid_id)%nlat                 = wrf_grid%ed33-1         !nlat
  wrf_da_variables(wrf_grid_id)%nsig                 = wrf_grid%ed32-1         !nsig
  wrf_da_variables(wrf_grid_id)%nsig_soil            = size(wrf_grid%smois,2)  !nsig_soil 
  wrf_da_variables(wrf_grid_id)%nprocs               = nprocs
  wrf_da_variables(wrf_grid_id)%comp_id              = comp_id
  wrf_da_variables(wrf_grid_id)%wrf_frame_id         = wrf_frame_id
  wrf_da_variables(wrf_grid_id)%parent_comp_id       = parent_comp_id
  wrf_da_variables(wrf_grid_id)%ensemble_comp_id     = ensemble_comp_id
  wrf_da_variables(wrf_grid_id)%wrf_ensemble_comp_id = wrf_ensemble_comp_id
  wrf_da_variables(wrf_grid_id)%comm                 = local_communicator
  wrf_da_variables(wrf_grid_id)%parent_comm          = local_communicator
  wrf_da_variables(wrf_grid_id)%ensemble_comm        = ensemble_comm
  wrf_da_variables(wrf_grid_id)%wrf_ensemble_comm    = wrf_ensemble_comm
  wrf_da_variables(wrf_grid_id)%ensemble_mem_id      = ensemble_mem_id
  wrf_da_variables(wrf_grid_id)%ensemble_nums        = ensemble_nums
  call CCPL_set_normal_time_step(wrf_ensemble_comp_id, wrf_da_variables(wrf_grid_id)%time_step, annotation="setting the time step for wrf ensemble")
!------------------------------------------------------------------------------------------------------------------------------
#endif
!! CCPL
!===============================================================================================================================================================================================

 
!=================================================Register wrf H2D grid using CCPL_register_H2D_grid_via_local_data =============================================================================
  CALL wrf_grid_parallel_decomposition(wrf_grid)
  grid_H2D_id    =   CCPL_register_H2D_grid_via_local_data(comp_id, "wrf_grid_via_local", "LON_LAT", "degrees", "acyclic", wrf_state_variables(wrf_grid_id)%num_global_cells, &
                                                           wrf_state_variables(wrf_grid_id)%num_local_cells, wrf_state_variables(wrf_grid_id)%local_cell_global_indexes,     &
                                                           -999999., -999999., -999999., -999999., wrf_state_variables(wrf_grid_id)%center_lons, &
                                                           wrf_state_variables(wrf_grid_id)%center_lats, wrf_state_variables(wrf_grid_id)%mask, annotation="register wrf  H2D  grid")
  wrf_state_variables(wrf_grid_id)%grid_H2D_id    = grid_H2D_id
!================================================================================================================================================================================================



!================================================== Register a new parallel decomposition of an H2D grid ========================================================================================
  decomp_id      =   CCPL_register_normal_parallel_decomp("decomp_wrf_grid", wrf_state_variables(wrf_grid_id)%grid_H2D_id, wrf_state_variables(wrf_grid_id)%num_local_cells, &
                                                   wrf_state_variables(wrf_grid_id)%local_cell_global_indexes, annotation="decompose wrf grid")
  wrf_state_variables(wrf_grid_id)%decomp_id      = decomp_id
!================================================================================================================================================================================================

!! CCPL
#ifdef CCPL_DA
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  wrf_state_variables(wrf_grid_id)%mask=1
  wrf_da_variables(wrf_grid_id)%grid_H2D_id    =   CCPL_register_H2D_grid_via_local_data(comp_id, "da_wrf_grid_via_local", "LON_LAT", "degrees", "acyclic", wrf_state_variables(wrf_grid_id)%num_global_cells, &
                                                           wrf_state_variables(wrf_grid_id)%num_local_cells, wrf_state_variables(wrf_grid_id)%local_cell_global_indexes,     &
                                                           -999999., -999999., -999999., -999999., wrf_state_variables(wrf_grid_id)%center_lons, &
                                                           wrf_state_variables(wrf_grid_id)%center_lats, wrf_state_variables(wrf_grid_id)%mask, annotation="register wrf  H2D  grid for ccpl_da")
  wrf_da_variables(wrf_grid_id)%decomp_id      =   CCPL_register_normal_parallel_decomp("da_decomp_wrf_grid", wrf_da_variables(wrf_grid_id)%grid_H2D_id, wrf_da_variables(wrf_grid_id)%da_num_local_cells, &
                                                   wrf_da_variables(wrf_grid_id)%local_cells_global_index_stagger, annotation="decompose wrf grid for ccpl_da")
  
  grid_vert_id1=CCPL_register_V1D_Z_grid_via_model_data(comp_id, "wrf_vertical_grid_nsig", "Pa", wrf_grid%znu(1:wrf_da_variables(wrf_grid_id)%nsig), annotation="register vertical grid znu for wrf")
  grid_vert_id2=CCPL_register_V1D_Z_grid_via_model_data(comp_id, "wrf_vertical_grid_nsig_1", "Pa", wrf_grid%znw(1:wrf_da_variables(wrf_grid_id)%nsig+1), annotation="register vertical grid znw for wrf")
  regional_time_coord_values=(/1.,2.,3.,4.,5.,6./)
  grid_regional_time_id=CCPL_register_V1D_Z_grid_via_model_data(comp_id, "grid_for_regional_time", "virtual_gird", regional_time_coord_values, annotation="register simplified grid for regional time")
  grid_3d_id1=CCPL_register_MD_grid_via_multi_grids(comp_id, "wrf_3d_grid_1", wrf_da_variables(wrf_grid_id)%grid_H2D_id, grid_vert_id1, annotation="register wrf 3D-grid with nsig vertical levels")
  grid_3d_id2=CCPL_register_MD_grid_via_multi_grids(comp_id, "wrf_3d_grid_2", wrf_da_variables(wrf_grid_id)%grid_H2D_id, grid_vert_id2, annotation="register wrf 3D-grid with nsig+1 vertical levels") 
  wrf_da_variables(wrf_grid_id)%grid_vert_id1         = grid_vert_id1
  wrf_da_variables(wrf_grid_id)%grid_vert_id2         = grid_vert_id2
  wrf_da_variables(wrf_grid_id)%grid_regional_time_id = grid_regional_time_id
  wrf_da_variables(wrf_grid_id)%grid_3d_id1           = grid_3d_id1
  wrf_da_variables(wrf_grid_id)%grid_3d_id2           = grid_3d_id2
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#endif
!! CCPL

!================================================== Register fields using  CCPL_register_field_instance =========================================================================================
  CALL allocate_coupling_buf_fields(wrf_grid)
  field_id_sss      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%SSS_buf, "sss", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_ocn, &
                                                      usage_tag=CCPL_TAG_CPL_REST, field_unit="Psu", annotation="register field instance of SSS") 
  field_id_sst      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%SST_buf, "sst", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_ocn, &
                                                      usage_tag=CCPL_TAG_CPL_REST, field_unit="Celsius", annotation="register field instance of SST")
  field_id_ssu      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%SSU_buf, "ssu", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_ocn, &
                                                      usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of SSU")
  field_id_ssv      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%SSV_buf, "ssv", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_ocn, &
                                                      usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of SSV")

  field_id_ustokes  =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%u_stokes_srf_buf, "u_stokes_srf", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_wave,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of u component of stokes drift")
  field_id_vstokes  =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%v_stokes_srf_buf, "v_stokes_srf", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_wave,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of v component of stokes drift")
  field_id_uorb     =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%u_orbital_srf_buf, "u_orbital_srf", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_wave,     &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of u component of orbital speed")
  field_id_vorb     =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%v_orbital_srf_buf, "v_orbital_srf", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_wave,     &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of v component of orbital speed")
  field_id_whitecap =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%whitecap_fraction_buf,"whitecap_fraction",wrf_state_variables(wrf_grid_id)%decomp_id,wrf_state_variables(wrf_grid_id)%grid_H2D_id,field_mark_wave, &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="unitless", annotation="register field instance of whitecape fraction")
  field_id_breakrate=   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%break_volume_buf, "break_volume", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_wave,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m3 m-2 s-1 ", annotation="register field instance of breaking rate")
  field_id_hwave    =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%significant_wave_height_buf, "significant_wave_height", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id,  &
                                                     field_mark_wave, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of significant wave height")
  field_id_lwavep   =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%peak_wave_length_buf, "peak_wave_length", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id,                &
                                                     field_mark_wave, usage_tag=CCPL_TAG_CPL_REST, field_unit="m", annotation="register field instance of peak wave length")
  field_id_pwave    =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%peak_wave_period_buf, "peak_wave_period", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id,                &
                                                     field_mark_wave, usage_tag=CCPL_TAG_CPL_REST, field_unit="second", annotation="register field instance of peak wave period")
  field_id_tau_in   =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%tau_in_buf, "tau_in", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_wave,                   &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="N", annotation="register field instance of wave induced stress")
  
  
  field_id_u10      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%U10_buf, "u10", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of u component wind speed at 10 meters")
  field_id_v10      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%V10_buf, "v10", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of v component wind speed at 10 meters")
  field_id_ust      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%UST_buf, "ust", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="m s-1", annotation="register field instance of friction velocity")
  field_id_rhoa     =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%RHOA_buf, "rhoa", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,     &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="kg m-3", annotation="register field instance of surface air density")
  field_id_t2       =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%T2_buf, "t2", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,         &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="K", annotation="register field instance of air temperature at 2 meters")
  field_id_q2       =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%Q2_buf, "q2", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,         &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="kg kg-1", annotation="register field instance of air specific humidity at 2 meters")
  field_id_qfx      =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%QFX_buf, "qfx", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,       &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="kg m-2 s-1", annotation="register field instance of surface moisture flux")
  field_id_sh       =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%SH_buf, "sh", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,         &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="W m-2 ", annotation="register field instance of sensible heat flux")
  field_id_lh       =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%LH_buf, "lh", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,         &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="W m-2", annotation="register field instance of latent heat flux")
  field_id_rain     =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%RAIN_buf, "rain", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,     &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="mm s-1", annotation="register field instance of  precipitation")
  field_id_mslp     =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%MSLP_buf, "mslp", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm,     &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="Pasca", annotation="register field instance of mean sea level pressure")
  field_id_swdown   =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%SWDOWN_buf, "swdown", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm, &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="W m-2", annotation="register field instance of net short wave flux ")
  field_id_lwdown   =   CCPL_register_field_instance(wrf_state_variables(wrf_grid_id)%LWDOWN_buf, "lwdown", wrf_state_variables(wrf_grid_id)%decomp_id, wrf_state_variables(wrf_grid_id)%grid_H2D_id, field_mark_atm, &
                                                     usage_tag=CCPL_TAG_CPL_REST, field_unit="W m-2", annotation="register field instance of downward  long wave flux ")
!================================================================================================================================================================================================================

!! CCPL
#ifdef CCPL_DA
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  CALL allocate_da_buf_fields(wrf_grid)
  call send_da_fields(wrf_grid)
  da_field_mark  = 3
  da_field_id_regional_time=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%regional_time,"regional_time",-1,&
    wrf_da_variables(wrf_grid_id)%grid_regional_time_id,da_field_mark,CCPL_TAG_CPL_REST,annotation="register time for GSI")
  da_field_id_pt=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_pt,"P_TOP",-1,&
    wrf_da_variables(wrf_grid_id)%comp_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register pressure top of the model")
  da_field_id_aeta1=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%aeta1,"ZNU_1",-1,&
    wrf_da_variables(wrf_grid_id)%grid_vert_id1,da_field_mark,CCPL_TAG_CPL_REST, annotation="register eta values on half (mass) levels")
  da_field_id_aeta2=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%aeta2,"ZNU_2",-1,&
    wrf_da_variables(wrf_grid_id)%grid_vert_id1,da_field_mark,CCPL_TAG_CPL_REST, annotation="register eta values on half (mass) levels")
  da_field_id_eta1=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%eta1,"ZNW_1",-1,&
    wrf_da_variables(wrf_grid_id)%grid_vert_id2,da_field_mark,CCPL_TAG_CPL_REST, annotation="register eta values on full (w) levels")
  da_field_id_eta2=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%eta2,"ZNW_2",-1,&
    wrf_da_variables(wrf_grid_id)%grid_vert_id2,da_field_mark,CCPL_TAG_CPL_REST, annotation="register eta values on full (w) levels")
  da_field_id_dx_mc=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_dx_mc,"MAPFAC_MX",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register x direction map scale factor")
  da_field_id_dy_mc=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_dy_mc,"MAPFAC_MY",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register y direction map scale factor")
  da_field_id_xlat=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_xlat,"XLAT",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register latitudes")
  da_field_id_xlong=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_xlong,"XLONG",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register longitudes")
  da_field_id_mu=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_mu,"mu",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register perturbation dry air mass in column")
  da_field_id_mub=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_mub,"mub",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register base state dry air mass in column")
  da_field_id_psfc0=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_psfc0,"psfc0",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register surface pressure")
  da_field_id_phb=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_phb,"PHB",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="m2 s-2",annotation="register base-state geopotential")
  da_field_id_landmask=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_landmask,"landmask",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register land mask (1 for land, 0 for water)")
  da_field_id_seaice=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_seaice,"seaice",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register sea ice flag")
  da_field_id_sst=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_sst,"ccpl_sst",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="K",annotation="register skin sea surface temperature")
  da_field_id_ivgtyp=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_ivgtyp,"ivgtyp",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register domain vegetation category")
  da_field_id_isltyp=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_isltyp,"isltyp",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register domain soil category")
  da_field_id_vegfra=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_vegfra,"vegfra",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, annotation="register vegetation fraction")
  da_field_id_snow=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_snow,"snow",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST,field_unit="kg m-2",annotation="register accumulated melted snow") 
  da_field_id_u10=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_u10,"ccpl_u10",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST,field_unit="m s-1",annotation="register field u10")
  da_field_id_v10=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_v10,"ccpl_v10",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST,field_unit="m s-1",annotation="register field v10")
  da_field_id_smois=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_smois,"smois",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="m3 m-3",annotation="register soil moistures (k=1)")
  da_field_id_tslb=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_tslb,"tslb",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="K",annotation="register soil temperature (k=1)")
  da_field_id_tsk=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_tsk,"tsk",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_H2D_id,da_field_mark,CCPL_TAG_CPL_REST, field_unit="K",annotation="register surface skin temperature")
  da_field_id_t=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_t,"ccpl_t",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_3d_id1,da_field_mark,CCPL_TAG_CPL_REST, field_unit="K",annotation="register perturbation potential temperature theta-t0")
  da_field_id_u=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_u,"ccpl_u",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_3d_id1,da_field_mark,CCPL_TAG_CPL_REST, field_unit="m s-1",annotation="register x-wind component at mass point")
  da_field_id_v=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_v,"ccpl_v",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_3d_id1,da_field_mark,CCPL_TAG_CPL_REST, field_unit="m s-1",annotation="register y-wind component at mass point")
  da_field_id_q=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_q,"ccpl_q",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_3d_id1,da_field_mark,CCPL_TAG_CPL_REST, field_unit="kg kg-1",annotation="register water vapor mixing ratio")
  da_field_id_w=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_w,"ccpl_w",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_3d_id1,da_field_mark,CCPL_TAG_CPL_REST, field_unit="m s-1",annotation="register z-wind component")
  da_field_id_ph=CCPL_register_field_instance(wrf_da_variables(wrf_grid_id)%ccpl_ph,"ccpl_ph",wrf_da_variables(wrf_grid_id)%decomp_id,&
    wrf_da_variables(wrf_grid_id)%grid_3d_id1,da_field_mark,CCPL_TAG_CPL_REST, field_unit="m2 s-2",annotation="register perturbation geopotential")
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  da_timer_id = CCPL_define_single_timer(comp_id, "seconds", 21600, 86400, 0, annotation="define a single timer for comp_id")
#endif
!! CCPL

!                                                    Define a single timer that is a periodic timer for the given component model
  timer_id = CCPL_define_single_timer(comp_id, "seconds", 300, 0, 0, annotation="define a single timer for comp_id") 

!=================================================== Register  import and export interface for a number of field instances ======================================================================================
  allocate(fields_id(13))
  fields_id(1)  = field_id_sss
  fields_id(2)  = field_id_sst
  fields_id(3)  = field_id_ssu
  fields_id(4)  = field_id_ssv 
  !import_interface_id = CCPL_register_import_interface("ATM_receive_from_OCN_"//trim(adjustl(ensemble_mem_id_str)), 4, fields_id, timer_id, 1, annotation="register interface for receiving data from OCN")
  import_interface_id = CCPL_register_import_interface("ATM_receive_from_OCN", 4, fields_id, timer_id, 1, annotation="register interface for receiving data from OCN")
  fields_id(1)  = field_id_ustokes
  fields_id(2)  = field_id_vstokes
  fields_id(3)  = field_id_uorb
  fields_id(4)  = field_id_vorb
  fields_id(5)  = field_id_whitecap
  fields_id(6)  = field_id_breakrate
  fields_id(7)  = field_id_hwave
  fields_id(8)  = field_id_lwavep
  fields_id(9)  = field_id_pwave
  fields_id(10) = field_id_tau_in
  import_interface_id = CCPL_register_import_interface("ATM_receive_from_WAVE", 10, fields_id, timer_id, 1, annotation="register interface for receiving data from WAVE")

  fields_id(1)  = field_id_u10
  fields_id(2)  = field_id_v10
  fields_id(3)  = field_id_ust
  fields_id(4)  = field_id_rain
  fields_id(5)  = field_id_t2
  fields_id(6)  = field_id_q2
  fields_id(7)  = field_id_qfx
  fields_id(8)  = field_id_sh
  fields_id(9)  = field_id_lh
  fields_id(10) = field_id_rhoa
  fields_id(11) = field_id_mslp
  fields_id(12) = field_id_swdown
  fields_id(13) = field_id_lwdown
  export_interface_id = CCPL_register_export_interface("ATM_send_to_OCN", 13, fields_id, timer_id, annotation="register interface for sending ATM data to OCN")
  export_interface_id = CCPL_register_export_interface("ATM_send_to_WAVE", 4, fields_id, timer_id, annotation="register interface for sending ATM data to WAVE")
  deallocate(fields_id)
!=========================================== coupling generation ============================================================================================================================================
!  Disable external coupling genration used for coupled nesting
!  CALL CCPL_get_configurable_comps_full_names(comp_id, "external_comps_for_coupling_generation", num_comps, comps_full_names, individual_or_family, annotation="test CCPL_get_configurable_comps_full_names")
!  CALL CCPL_do_external_coupling_generation(num_comps, comps_full_names, individual_or_family)
  CALL CCPL_do_family_coupling_generation(ensemble_comp_id, ensemble_comp_id, "Do family coupling generation of the coupled model")
!=========================================== check wich model will be coupled ===============================================================================================================================
  wave_coupled  = CCPL_is_comp_type_coupled(comp_id,"wave", annotation="atm is coupling with wave")
  ocean_coupled = CCPL_is_comp_type_coupled(comp_id,"ocn", annotation="atm is coupling with ocean") 
!================================================================================================================================================================================================================

!! CCPL
#ifdef CCPL_DA
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  allocate(wrf_da_variables(wrf_grid_id)%da_fields_id(1:32))
  wrf_da_variables(wrf_grid_id)%da_fields_id(1)=da_field_id_pt
  wrf_da_variables(wrf_grid_id)%da_fields_id(2)=da_field_id_aeta1
  wrf_da_variables(wrf_grid_id)%da_fields_id(3)=da_field_id_aeta2
  wrf_da_variables(wrf_grid_id)%da_fields_id(4)=da_field_id_eta1
  wrf_da_variables(wrf_grid_id)%da_fields_id(5)=da_field_id_eta2
  wrf_da_variables(wrf_grid_id)%da_fields_id(6)=da_field_id_dx_mc
  wrf_da_variables(wrf_grid_id)%da_fields_id(7)=da_field_id_dy_mc
  wrf_da_variables(wrf_grid_id)%da_fields_id(8)=da_field_id_xlat
  wrf_da_variables(wrf_grid_id)%da_fields_id(9)=da_field_id_xlong
  wrf_da_variables(wrf_grid_id)%da_fields_id(10)=da_field_id_psfc0
  wrf_da_variables(wrf_grid_id)%da_fields_id(11)=da_field_id_phb
  wrf_da_variables(wrf_grid_id)%da_fields_id(12)=da_field_id_landmask
  wrf_da_variables(wrf_grid_id)%da_fields_id(13)=da_field_id_seaice
  wrf_da_variables(wrf_grid_id)%da_fields_id(14)=da_field_id_sst
  wrf_da_variables(wrf_grid_id)%da_fields_id(15)=da_field_id_ivgtyp
  wrf_da_variables(wrf_grid_id)%da_fields_id(16)=da_field_id_isltyp
  wrf_da_variables(wrf_grid_id)%da_fields_id(17)=da_field_id_vegfra
  wrf_da_variables(wrf_grid_id)%da_fields_id(18)=da_field_id_snow
  wrf_da_variables(wrf_grid_id)%da_fields_id(19)=da_field_id_u10
  wrf_da_variables(wrf_grid_id)%da_fields_id(20)=da_field_id_v10
  wrf_da_variables(wrf_grid_id)%da_fields_id(21)=da_field_id_smois
  wrf_da_variables(wrf_grid_id)%da_fields_id(22)=da_field_id_tslb
  wrf_da_variables(wrf_grid_id)%da_fields_id(23)=da_field_id_tsk
  wrf_da_variables(wrf_grid_id)%da_fields_id(24)=da_field_id_regional_time
  wrf_da_variables(wrf_grid_id)%da_fields_id(25)=da_field_id_mu
  wrf_da_variables(wrf_grid_id)%da_fields_id(26)=da_field_id_mub
  wrf_da_variables(wrf_grid_id)%da_fields_id(27)=da_field_id_t
  wrf_da_variables(wrf_grid_id)%da_fields_id(28)=da_field_id_q
  wrf_da_variables(wrf_grid_id)%da_fields_id(29)=da_field_id_u
  wrf_da_variables(wrf_grid_id)%da_fields_id(30)=da_field_id_v
  wrf_da_variables(wrf_grid_id)%da_fields_id(31)=da_field_id_w
  wrf_da_variables(wrf_grid_id)%da_fields_id(32)=da_field_id_ph 
  allocate(wrf_da_variables(wrf_grid_id)%control_vars(1:5))
  wrf_da_variables(wrf_grid_id)%control_vars(1) = wrf_da_variables(wrf_grid_id)%nlat         !nlat
  wrf_da_variables(wrf_grid_id)%control_vars(2) = wrf_da_variables(wrf_grid_id)%nlon         !nlon
  wrf_da_variables(wrf_grid_id)%control_vars(3) = wrf_da_variables(wrf_grid_id)%nsig         !nsig
  wrf_da_variables(wrf_grid_id)%control_vars(4) = wrf_da_variables(wrf_grid_id)%nsig_soil    !nsig_soil
  wrf_da_variables(wrf_grid_id)%control_vars(5) = wrf_da_variables(wrf_grid_id)%nprocs
  allocate(wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(1:6))
  wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(1) = wrf_da_variables(wrf_grid_id)%grid_H2D_id
  wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(2) = wrf_da_variables(wrf_grid_id)%grid_vert_id1
  wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(3) = wrf_da_variables(wrf_grid_id)%grid_vert_id2
  wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(4) = wrf_da_variables(wrf_grid_id)%grid_3d_id1
  wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(5) = wrf_da_variables(wrf_grid_id)%grid_3d_id2
  wrf_da_variables(wrf_grid_id)%comp_or_grid_ids(6) = wrf_da_variables(wrf_grid_id)%grid_regional_time_id
  allocate(wrf_da_variables(wrf_grid_id)%decomp_ids(1:6))
  wrf_da_variables(wrf_grid_id)%decomp_ids     = -1
  wrf_da_variables(wrf_grid_id)%decomp_ids(1)  = wrf_da_variables(wrf_grid_id)%decomp_id
  allocate(wrf_da_variables(wrf_grid_id)%field_inst_ids(1:32))
  wrf_da_variables(wrf_grid_id)%field_inst_ids = wrf_da_variables(wrf_grid_id)%da_fields_id(1:32)
  allocate(wrf_da_variables(wrf_grid_id)%timers_id(1))
  wrf_da_variables(wrf_grid_id)%timers_id   = timer_id
  wrf_da_variables(wrf_grid_id)%da_timer_id   = da_timer_id
  write(*,*) "==============================================================================="
  write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_MEMBER initialize"
  call start_timing
  wrf_da_variables(wrf_grid_id)%gsi_ens_member_inst_id=CCPL_ensemble_procedures_inst_init(wrf_da_variables(wrf_grid_id)%wrf_ensemble_comp_id, wrf_da_variables(wrf_grid_id)%comp_id, "gsi_ens_member", wrf_da_variables(wrf_grid_id)%field_inst_ids, wrf_da_variables(wrf_grid_id)%comp_or_grid_ids, wrf_da_variables(wrf_grid_id)%decomp_ids, wrf_da_variables(wrf_grid_id)%timers_id, wrf_da_variables(wrf_grid_id)%control_vars, annotation="do ensemble procedure GSI_ENS_MEMBER initialize")
  call end_timing(" CCPL_OUT: GSI_ENS_MEMBER initialize ")
  write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_MEMBER initialize"
  write(*,*) "==============================================================================="

  write(*,*) "==============================================================================="
  write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_AVE initialize"
  call start_timing
  wrf_da_variables(wrf_grid_id)%gsi_ens_ave_inst_id=CCPL_ensemble_procedures_inst_init(wrf_da_variables(wrf_grid_id)%wrf_ensemble_comp_id, wrf_da_variables(wrf_grid_id)%comp_id, "gsi_ens_aver", wrf_da_variables(wrf_grid_id)%field_inst_ids, wrf_da_variables(wrf_grid_id)%comp_or_grid_ids, wrf_da_variables(wrf_grid_id)%decomp_ids, wrf_da_variables(wrf_grid_id)%timers_id, wrf_da_variables(wrf_grid_id)%control_vars, annotation="do ensemble procedure GSI_ENS_AVE initialize")
  call end_timing(" CCPL_OUT: GSI_ENS_AVE initialize ")
  write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_AVE initialize"
  write(*,*) "==============================================================================="

  write(*,*) "==============================================================================="
  write(*,*) "[CCPL <WRF>] Start ensemble procedure ENKF initialize"
  call start_timing
  wrf_da_variables(wrf_grid_id)%enkf_inst_id=CCPL_ensemble_procedures_inst_init(wrf_da_variables(wrf_grid_id)%wrf_ensemble_comp_id, wrf_da_variables(wrf_grid_id)%comp_id, "enkf", wrf_da_variables(wrf_grid_id)%field_inst_ids, wrf_da_variables(wrf_grid_id)%comp_or_grid_ids, wrf_da_variables(wrf_grid_id)%decomp_ids, wrf_da_variables(wrf_grid_id)%timers_id, wrf_da_variables(wrf_grid_id)%control_vars, annotation="do ensemble procedure ENKF initialize")
  call end_timing(" CCPL_OUT: ENKF initialize ")
  write(*,*) "[CCPL <WRF>] Finish ensemble procedure ENKF initialize"
  write(*,*) "==============================================================================="
  call CCPL_change_to_ensemble_dir(wrf_da_variables(wrf_grid_id)%wrf_frame_id, wrf_da_variables(wrf_grid_id)%ensemble_mem_id)
!------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#endif
!! CCPL

END SUBROUTINE register_component_coupling_configuration 





SUBROUTINE run_coupling(wrf_grid)
  implicit none
  type(domain), pointer :: wrf_grid
  integer               :: wrf_grid_id
  integer               :: i
  integer               :: field_update_status(13)
  logical               :: interface_status
  wrf_grid_id = wrf_grid%id
  IF(wrf_state_variables(wrf_grid_id)%initialization)THEN
       CALL send_coupling_fields(wrf_grid)
       IF(ocean_coupled) interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_send_to_OCN", .true., annotation="initialization,execute ATM_send_to_OCN")
       IF(wave_coupled)  interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_send_to_WAVE", .true., annotation="initialization,execute ATM_send_to_WAVE")
       
       IF(ocean_coupled) interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_receive_from_OCN", .true., field_update_status, annotation="initialization,execute ATM_receive_from_OCN")
       IF(wave_coupled)  interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_receive_from_WAVE", .true., field_update_status, annotation="initialization,execute ATM_receive_from_WAVE")
       CALL receive_coupling_fields(wrf_grid) 
       
       write(message,*)"zhaobiao,run_coupling,initialization, domain ", wrf_grid%id
       CALL wrf_debug(100,message)
  ELSE
       CALL send_coupling_fields(wrf_grid)
       IF(ocean_coupled) interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_send_to_OCN", .false., annotation="integrate,execute ATM_send_to_OCN")
       IF(wave_coupled)  interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_send_to_WAVE", .false., annotation="integrate,execute ATM_send_to_WAVE")
       
       IF(ocean_coupled) interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_receive_from_OCN", .false., field_update_status,annotation="integrate,execute ATM_receive_from_OCN")
       IF(wave_coupled)  interface_status = CCPL_execute_interface_using_name(wrf_state_variables(wrf_grid_id)%comp_id, "ATM_receive_from_WAVE", .false., field_update_status,annotation="integrate,execute ATM_receive_from_WAVE")
       CALL receive_coupling_fields(wrf_grid)
!! CCPL
#ifndef CCPL_DA       
       CALL CCPL_advance_time(wrf_state_variables(wrf_grid_id)%comp_id)
#endif       
!! CCPL
       !CALL CCPL_do_restart_write_IO(wrf_state_variables(wrf_grid_id)%comp_id,.false.)     
       write(message,*)"zhaobiao,run_coupling,integrate,domain ", wrf_grid%id
       CALL wrf_debug(100,message)
  END IF
  wrf_state_variables(wrf_grid_id)%initialization = .false.

END SUBROUTINE run_coupling

!! CCPL
#ifdef CCPL_DA
SUBROUTINE run_ccpl_da(wrf_grid)   

    USE ESMF_mod, ONLY : ESMF_TimeGet
    USE module_domain, ONLY : domain_clock_get
    USE module_utility, ONLY : WRFU_Time

    IMPLICIT NONE  
INCLUDE "mpif.h"      
    TYPE(domain) , POINTER      :: wrf_grid
    INTEGER                     :: err,wrf_grid_id        
    TYPE(WRFU_Time)             :: currentTime
    INTEGER                     :: iyear,imonth,iday,ihour,iminute,isecond
    LOGICAL                     :: prepare_da

    wrf_grid_id = wrf_grid%id    
    call domain_clock_get(wrf_grid, current_time=currentTime)
    call ESMF_TimeGet(currentTime, YY=iyear,MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond)        
    call CCPL_check_current_time(wrf_da_variables(wrf_grid_id)%comp_id, iyear*10000+imonth*100+iday, ihour*3600+iminute*60+isecond, "check time") 
    if (CCPL_is_timer_on(wrf_da_variables(wrf_grid_id)%da_timer_id)) then
        prepare_da=.true.
    else
        prepare_da=.false.
    end if
    if (prepare_da) then    
    call send_da_fields(wrf_grid)        
    write(*,*) "==============================================================================="
    write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_AVE run"
    call start_timing
    call CCPL_ensemble_procedures_inst_run(wrf_da_variables(wrf_grid_id)%gsi_ens_ave_inst_id,0,-1,annotation="do GSI_ENS_AVE run")  
    call end_timing(" CCPL_OUT: GSI_ENS_AVE run ")
    write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_AVE run"
    write(*,*) "==============================================================================="   
    write(*,*) "==============================================================================="
    write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_MEMBER run"   
    call start_timing
    call CCPL_ensemble_procedures_inst_run(wrf_da_variables(wrf_grid_id)%gsi_ens_member_inst_id,0,-1,annotation="do GSI_ENS_MEMBER run")
    call end_timing(" CCPL_OUT: GSI_ENS_MEMBER run ")
    write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_MEMBER run"
    write(*,*) "==============================================================================="
    write(*,*) "==============================================================================="
    write(*,*) "[CCPL <WRF>] Start ensemble procedure ENKF run"
    call start_timing
    call CCPL_ensemble_procedures_inst_run(wrf_da_variables(wrf_grid_id)%enkf_inst_id,1,-1,annotation="do ENKF run")
    call end_timing(" CCPL_OUT: ENKF run ")
    write(*,*) "[CCPL <WRF>] Finish ensemble procedure ENKF run"
    write(*,*) "==============================================================================="
    call CCPL_change_to_ensemble_dir(wrf_da_variables(wrf_grid_id)%wrf_frame_id, wrf_da_variables(wrf_grid_id)%ensemble_mem_id)
    call receive_da_fields(wrf_grid)
    CALL CCPL_do_restart_write_IO(wrf_da_variables(wrf_grid_id)%comp_id, .false.,annotation="WRF do restart write")  

END SUBROUTINE run_ccpl_da
#endif
!! CCPL





SUBROUTINE wrf_grid_parallel_decomposition(wrf_grid_temp)
   implicit none
   TYPE(domain) , POINTER :: wrf_grid_temp
   integer                :: wrf_grid_temp_id
   integer                :: num_local_rows, num_local_cols, num_local_cells
   integer                :: num_global_rows, num_global_cols, num_global_cells
   integer                :: i, j, n
   integer                :: iep, isp, jep, jsp, ied, isd, jed, jsd
!! CCPL
#ifdef CCPL_DA
  !-------------------------------------------------------------------------------------------------------------------------------------------------
   integer                :: dm_is,dm_ie,dm_js,dm_je,Isize,Jsize,da_num_local_cells
  !-------------------------------------------------------------------------------------------------------------------------------------------------
#endif
!! CCPL
   isd = wrf_grid_temp%sd31
   ied = wrf_grid_temp%ed31-1
   isp = wrf_grid_temp%sp31
   iep = MIN(wrf_grid_temp%ep31,ied)
   jsd = wrf_grid_temp%sd33
   jed = wrf_grid_temp%ed33-1
   jsp = wrf_grid_temp%sp33
   jep = MIN(wrf_grid_temp%ep33,jed)
   wrf_grid_temp_id = wrf_grid_temp%id
   
   num_local_rows   = jep - jsp + 1
   num_local_cols   = iep - isp + 1
   num_global_rows  = jed - jsd + 1
   num_global_cols  = ied - isd + 1
   num_local_cells  = num_local_rows  * num_local_cols
   num_global_cells = num_global_rows * num_global_cols
   wrf_state_variables(wrf_grid_temp_id)%num_local_cells   = num_local_cells
   wrf_state_variables(wrf_grid_temp_id)%num_global_cells  = num_global_cells
   write(message,*)"zhaobiao, number of cells for row and column: ",num_global_rows,num_global_cols,num_local_rows,num_local_cols
   CALL wrf_debug(100,message)


   IF (wrf_grid_temp%sd32 .ne. wrf_grid_temp%sp32 .or. wrf_grid_temp%ed32 .ne. wrf_grid_temp%ep32)THEN
      CALL CCPL_report_error(wrf_state_variables(wrf_grid_temp_id)%comp_id, .false., "the parallel decomposition on vertical direction is used in WRF")  
   END IF

   IF (num_local_cells .gt. 0) THEN
      allocate(wrf_state_variables(wrf_grid_temp_id)%local_cell_global_indexes(num_local_cells))
      allocate(wrf_state_variables(wrf_grid_temp_id)%center_lats(num_local_cells))
      allocate(wrf_state_variables(wrf_grid_temp_id)%center_lons(num_local_cells))
      allocate(wrf_state_variables(wrf_grid_temp_id)%mask(num_local_cells))
      wrf_state_variables(wrf_grid_temp_id)%mask = 0
      n = 0
      DO j = jsp, jep
         DO i = isp, iep
            n = n + 1
            wrf_state_variables(wrf_grid_temp_id)%local_cell_global_indexes(n) = (j-1)*num_global_cols+i
            wrf_state_variables(wrf_grid_temp_id)%center_lats(n)               = wrf_grid_temp%xlat(i,j)
            wrf_state_variables(wrf_grid_temp_id)%center_lons(n)               = wrf_grid_temp%xlong(i,j)
            IF(wrf_grid_temp%landmask(i,j).eq.0)THEN
              wrf_state_variables(wrf_grid_temp_id)%mask(n) = 1
            END IF
         END DO
      END DO
   END IF
!! CCPL
#ifdef CCPL_DA
   dm_is = wrf_grid_temp%sm31
   dm_ie = wrf_grid_temp%em31
   dm_js = wrf_grid_temp%sm33
   dm_je = wrf_grid_temp%em33
   Isize=dm_ie-dm_is+1
   Jsize=dm_je-dm_js+1
   da_num_local_cells=Isize*Jsize
   wrf_da_variables(wrf_grid_temp_id)%da_num_local_cells=da_num_local_cells
   if (da_num_local_cells .gt. 0) then
     allocate(wrf_da_variables(wrf_grid_temp_id)%local_cells_global_index_stagger(da_num_local_cells))
     allocate(wrf_da_variables(wrf_grid_temp_id)%mask(da_num_local_cells))  
     wrf_da_variables(wrf_grid_temp_id)%local_cells_global_index_stagger = CCPL_NULL_INT
     n=1
     do j=dm_js,dm_je
         do i=dm_is,dm_ie
             if (j.ge.wrf_grid_temp%sp33.and.j.le.wrf_grid_temp%ep33.and.&
                 i.ge.wrf_grid_temp%sp31.and.i.le.wrf_grid_temp%ep31) then
               if (j.ne.wrf_grid_temp%ed33.and.i.ne.wrf_grid_temp%ed31) then
                   wrf_da_variables(wrf_grid_temp_id)%local_cells_global_index_stagger(n) = (j-1)*(wrf_grid_temp%ed31-1)+i
               end if
             end if
             wrf_da_variables(wrf_grid_temp_id)%mask(n) = 1
             n=n+1
         end do
     end do
   end if
#endif
!! CCPL
   write(*,*) "Finish wrf_grid_parallel_decomposition"
END SUBROUTINE wrf_grid_parallel_decomposition




SUBROUTINE allocate_coupling_buf_fields(wrf_grid_temp)
   implicit none
   TYPE(domain) , POINTER :: wrf_grid_temp
   integer                :: wrf_grid_temp_id
   integer                :: i, j
   integer                :: iep, isp, jep, jsp, ied, isd, jed, jsd
   isd = wrf_grid_temp%sd31
   ied = wrf_grid_temp%ed31-1
   isp = wrf_grid_temp%sp31
   iep = MIN(wrf_grid_temp%ep31,ied)
   jsd = wrf_grid_temp%sd33
   jed = wrf_grid_temp%ed33-1
   jsp = wrf_grid_temp%sp33
   jep = MIN(wrf_grid_temp%ep33,jed)
   wrf_grid_temp_id = wrf_grid_temp%id  
 
      IF (wrf_state_variables(wrf_grid_temp_id)%num_local_cells .gt. 0) THEN
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%U10_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%V10_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%RHOA_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%Q2_buf    (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%RAIN_buf    (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%MSLP_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%T2_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%UST_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%SH_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%LH_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%QFX_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%LWDOWN_buf  (isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%SWDOWN_buf  (isp:iep,jsp:jep))

        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%SST_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%SSS_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%SSU_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%SSV_buf(isp:iep,jsp:jep))

        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%u_stokes_srf_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%v_stokes_srf_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%u_orbital_srf_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%v_orbital_srf_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%whitecap_fraction_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%break_volume_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%significant_wave_height_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%peak_wave_length_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%peak_wave_period_buf(isp:iep,jsp:jep))
        ALLOCATE(wrf_state_variables(wrf_grid_temp_id)%tau_in_buf(isp:iep,jsp:jep))
        
        wrf_state_variables(wrf_grid_temp_id)%U10_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%V10_buf = 0.0 
        wrf_state_variables(wrf_grid_temp_id)%RHOA_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%Q2_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%RAIN_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%MSLP_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%T2_buf = 0.0 
        wrf_state_variables(wrf_grid_temp_id)%UST_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%SH_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%LH_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%QFX_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%LWDOWN_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%SWDOWN_buf = 0.0 
        
        wrf_state_variables(wrf_grid_temp_id)%SSS_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%SST_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%SSU_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%SSV_buf = 0.0
        
        wrf_state_variables(wrf_grid_temp_id)%u_stokes_srf_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%v_stokes_srf_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%u_orbital_srf_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%v_orbital_srf_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%whitecap_fraction_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%break_volume_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%significant_wave_height_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%peak_wave_length_buf = 0.0
        wrf_state_variables(wrf_grid_temp_id)%peak_wave_period_buf = 0.0 
        wrf_state_variables(wrf_grid_temp_id)%tau_in_buf = 0.0
       END IF
      
END SUBROUTINE allocate_coupling_buf_fields



SUBROUTINE send_coupling_fields(wrf_grid_temp)
   implicit none
   TYPE(domain) , POINTER :: wrf_grid_temp
   integer                :: wrf_grid_temp_id
   integer                :: i, j
   integer                :: iep, isp, jep, jsp, ied, isd, jed, jsd
   real,parameter         :: eps=1.0e-10
   isd = wrf_grid_temp%sd31
   ied = wrf_grid_temp%ed31-1
   isp = wrf_grid_temp%sp31
   iep = MIN(wrf_grid_temp%ep31,ied)
   jsd = wrf_grid_temp%sd33
   jed = wrf_grid_temp%ed33-1
   jsp = wrf_grid_temp%sp33
   jep = MIN(wrf_grid_temp%ep33,jed)
   wrf_grid_temp_id = wrf_grid_temp%id
        DO j = jsp, jep
         DO i = isp, iep
            !tansfer u10 and v10 from wrf grid to earth grid when wrf use lambert or polar projection
            wrf_state_variables(wrf_grid_temp_id)%U10_buf(i,j)    = wrf_grid_temp%u10(i,j)*wrf_grid_temp%cosa(i,j)-wrf_grid_temp%v10(i,j)*wrf_grid_temp%sina(i,j)
            wrf_state_variables(wrf_grid_temp_id)%V10_buf(i,j)    = wrf_grid_temp%v10(i,j)*wrf_grid_temp%cosa(i,j)+wrf_grid_temp%u10(i,j)*wrf_grid_temp%sina(i,j)
            wrf_state_variables(wrf_grid_temp_id)%UST_buf(i,j)    = wrf_grid_temp%ust(i,j)
            wrf_state_variables(wrf_grid_temp_id)%RHOA_buf(i,j)   = 1.0/(wrf_grid_temp%alt(i,1,j)+eps)
            wrf_state_variables(wrf_grid_temp_id)%T2_buf(i,j)     = wrf_grid_temp%t2(i,j)
            wrf_state_variables(wrf_grid_temp_id)%Q2_buf(i,j)     = wrf_grid_temp%q2(i,j)/(wrf_grid_temp%q2(i,j)+1) !convert water mixing ratio to specific humidity
            wrf_state_variables(wrf_grid_temp_id)%QFX_buf(i,j)    = wrf_grid_temp%qfx(i,j)
            wrf_state_variables(wrf_grid_temp_id)%LH_buf(i,j)     = -1.0*wrf_grid_temp%lh(i,j)
            wrf_state_variables(wrf_grid_temp_id)%SH_buf(i,j)     = -1.0*wrf_grid_temp%hfx(i,j)
            wrf_state_variables(wrf_grid_temp_id)%SWDOWN_buf(i,j) = wrf_grid_temp%gsw(i,j)
            wrf_state_variables(wrf_grid_temp_id)%LWDOWN_buf(i,j) = wrf_grid_temp%glw(i,j)
            wrf_state_variables(wrf_grid_temp_id)%MSLP_buf(i,j)   = wrf_grid_temp%psfc(i,j)*exp((9.81*wrf_grid_temp%ht(i,j))/(287.0*wrf_grid_temp%t2(i,j)*(1.0+0.61*wrf_grid_temp%q2(i,j))))
            wrf_state_variables(wrf_grid_temp_id)%RAIN_buf(i,j)   =(wrf_grid_temp%raincv(i,j)+wrf_grid_temp%rainncv(i,j))/wrf_grid_temp%dt       !time step averaged precipitation
         END DO
        END DO
       write(message,*)"zhaobiao,send_coupling_fields,domain",wrf_grid_temp%id
       CALL wrf_debug(100,message)
END SUBROUTINE send_coupling_fields



SUBROUTINE receive_coupling_fields(wrf_grid_temp)
   implicit none
   TYPE(domain) , POINTER :: wrf_grid_temp
   integer                :: wrf_grid_temp_id
   integer                :: i, j
   integer                :: iep, isp, jep, jsp, ied, isd, jed, jsd
   real,parameter         :: eps=1.0e-10 
   isd = wrf_grid_temp%sd31
   ied = wrf_grid_temp%ed31-1
   isp = wrf_grid_temp%sp31
   iep = MIN(wrf_grid_temp%ep31,ied)
   jsd = wrf_grid_temp%sd33
   jed = wrf_grid_temp%ed33-1
   jsp = wrf_grid_temp%sp33
   jep = MIN(wrf_grid_temp%ep33,jed)
   wrf_grid_temp_id = wrf_grid_temp%id  

      IF (wrf_state_variables(wrf_grid_temp_id)%num_local_cells .gt. 0) THEN

        DO j = jsp, jep
         DO i = isp, iep
            wrf_grid_temp%sst(i,j)          =    (wrf_state_variables(wrf_grid_temp_id)%SST_buf(i,j)+273.15)*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%sss(i,j)          =    wrf_state_variables(wrf_grid_temp_id)%SSS_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))       
            wrf_grid_temp%uoce(i,j)         =    (wrf_state_variables(wrf_grid_temp_id)%SSU_buf(i,j)*wrf_grid_temp%cosa(i,j)+wrf_state_variables(wrf_grid_temp_id)%SSV_buf(i,j)*wrf_grid_temp%sina(i,j))*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%voce(i,j)         =    (wrf_state_variables(wrf_grid_temp_id)%SSV_buf(i,j)*wrf_grid_temp%cosa(i,j)-wrf_state_variables(wrf_grid_temp_id)%SSU_buf(i,j)*wrf_grid_temp%sina(i,j))*(1-wrf_grid_temp%landmask(i,j))
            
            wrf_grid_temp%ustokes(i,j)      =    (wrf_state_variables(wrf_grid_temp_id)%u_stokes_srf_buf(i,j)*wrf_grid_temp%cosa(i,j)+wrf_state_variables(wrf_grid_temp_id)%v_stokes_srf_buf(i,j)*wrf_grid_temp%sina(i,j))*(1-wrf_grid_temp%landmask(i,j)) 
            wrf_grid_temp%vstokes(i,j)      =    (wrf_state_variables(wrf_grid_temp_id)%v_stokes_srf_buf(i,j)*wrf_grid_temp%cosa(i,j)-wrf_state_variables(wrf_grid_temp_id)%u_stokes_srf_buf(i,j)*wrf_grid_temp%sina(i,j))*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%uorb(i,j)         =    (wrf_state_variables(wrf_grid_temp_id)%u_orbital_srf_buf(i,j)*wrf_grid_temp%cosa(i,j)+wrf_state_variables(wrf_grid_temp_id)%v_orbital_srf_buf(i,j)*wrf_grid_temp%sina(i,j))*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%vorb(i,j)         =    (wrf_state_variables(wrf_grid_temp_id)%v_orbital_srf_buf(i,j)*wrf_grid_temp%cosa(i,j)-wrf_state_variables(wrf_grid_temp_id)%u_orbital_srf_buf(i,j)*wrf_grid_temp%sina(i,j))*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%whitecap(i,j)     =    wrf_state_variables(wrf_grid_temp_id)%whitecap_fraction_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%breakrate(i,j)    =    wrf_state_variables(wrf_grid_temp_id)%break_volume_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%hwave(i,j)        =    wrf_state_variables(wrf_grid_temp_id)%significant_wave_height_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%lwavep(i,j)       =    wrf_state_variables(wrf_grid_temp_id)%peak_wave_length_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%pwave(i,j)        =    wrf_state_variables(wrf_grid_temp_id)%peak_wave_period_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))
            wrf_grid_temp%tau_in(i,j)       =    wrf_state_variables(wrf_grid_temp_id)%tau_in_buf(i,j)*(1-wrf_grid_temp%landmask(i,j))

            wrf_grid_temp%mslp(i,j)         =    wrf_grid_temp%psfc(i,j)*exp((9.81*wrf_grid_temp%ht(i,j))/(287.0*wrf_grid_temp%t2(i,j)*(1.0+0.61*wrf_grid_temp%q2(i,j)))) 
            wrf_grid_temp%rhoa(i,j)         =    1.0/(wrf_grid_temp%alt(i,1,j)+eps)                  
         END DO
        END DO
       END IF
       write(message,*)"zhaobiao,receive_coupling_fields,domain ",wrf_grid_temp%id
       CALL wrf_debug(100,message)

END SUBROUTINE receive_coupling_fields



SUBROUTINE release_wrf_buf
USE module_configure, ONLY : model_config_rec 
implicit none
integer  :: i
 DO i=1,model_config_rec%max_dom
   IF(ALLOCATED(wrf_state_variables(i)%U10_buf)) DEALLOCATE(wrf_state_variables(i)%U10_buf)
   IF(ALLOCATED(wrf_state_variables(i)%V10_buf)) DEALLOCATE(wrf_state_variables(i)%V10_buf)
   IF(ALLOCATED(wrf_state_variables(i)%RHOA_buf)) DEALLOCATE(wrf_state_variables(i)%RHOA_buf)
   IF(ALLOCATED(wrf_state_variables(i)%UST_buf)) DEALLOCATE(wrf_state_variables(i)%UST_buf)
   IF(ALLOCATED(wrf_state_variables(i)%T2_buf)) DEALLOCATE(wrf_state_variables(i)%T2_buf)
   IF(ALLOCATED(wrf_state_variables(i)%Q2_buf)) DEALLOCATE(wrf_state_variables(i)%Q2_buf)
   IF(ALLOCATED(wrf_state_variables(i)%QFX_buf)) DEALLOCATE(wrf_state_variables(i)%QFX_buf)
   IF(ALLOCATED(wrf_state_variables(i)%LH_buf)) DEALLOCATE(wrf_state_variables(i)%LH_buf)
   IF(ALLOCATED(wrf_state_variables(i)%SH_buf)) DEALLOCATE(wrf_state_variables(i)%SH_buf)
   IF(ALLOCATED(wrf_state_variables(i)%RAIN_buf)) DEALLOCATE(wrf_state_variables(i)%RAIN_buf)
   IF(ALLOCATED(wrf_state_variables(i)%SWDOWN_buf)) DEALLOCATE(wrf_state_variables(i)%SWDOWN_buf)
   IF(ALLOCATED(wrf_state_variables(i)%LWDOWN_buf)) DEALLOCATE(wrf_state_variables(i)%LWDOWN_buf)
   IF(ALLOCATED(wrf_state_variables(i)%MSLP_buf)) DEALLOCATE(wrf_state_variables(i)%MSLP_buf)
   
   IF(ALLOCATED(wrf_state_variables(i)%SST_buf)) DEALLOCATE(wrf_state_variables(i)%SST_buf)
   IF(ALLOCATED(wrf_state_variables(i)%SSS_buf)) DEALLOCATE(wrf_state_variables(i)%SSS_buf)
   IF(ALLOCATED(wrf_state_variables(i)%SSU_buf)) DEALLOCATE(wrf_state_variables(i)%SSU_buf)
   IF(ALLOCATED(wrf_state_variables(i)%SSV_buf)) DEALLOCATE(wrf_state_variables(i)%SSV_buf)
   
   IF(ALLOCATED(wrf_state_variables(i)%whitecap_fraction_buf)) DEALLOCATE(wrf_state_variables(i)%whitecap_fraction_buf)
   IF(ALLOCATED(wrf_state_variables(i)%u_orbital_srf_buf)) DEALLOCATE(wrf_state_variables(i)%u_orbital_srf_buf)
   IF(ALLOCATED(wrf_state_variables(i)%v_orbital_srf_buf)) DEALLOCATE(wrf_state_variables(i)%v_orbital_srf_buf)
   IF(ALLOCATED(wrf_state_variables(i)%u_stokes_srf_buf)) DEALLOCATE(wrf_state_variables(i)%u_stokes_srf_buf)
   IF(ALLOCATED(wrf_state_variables(i)%v_stokes_srf_buf)) DEALLOCATE(wrf_state_variables(i)%v_stokes_srf_buf)
   IF(ALLOCATED(wrf_state_variables(i)%break_volume_buf)) DEALLOCATE(wrf_state_variables(i)%break_volume_buf)
   IF(ALLOCATED(wrf_state_variables(i)%tau_in_buf)) DEALLOCATE(wrf_state_variables(i)%tau_in_buf)
   IF(ALLOCATED(wrf_state_variables(i)%significant_wave_height_buf)) DEALLOCATE(wrf_state_variables(i)%significant_wave_height_buf)
   IF(ALLOCATED(wrf_state_variables(i)%peak_wave_length_buf)) DEALLOCATE(wrf_state_variables(i)%peak_wave_length_buf)
   IF(ALLOCATED(wrf_state_variables(i)%peak_wave_period_buf)) DEALLOCATE(wrf_state_variables(i)%peak_wave_period_buf)
  
   IF(ALLOCATED(wrf_state_variables(i)%local_cell_global_indexes)) DEALLOCATE(wrf_state_variables(i)%local_cell_global_indexes)
   IF(ALLOCATED(wrf_state_variables(i)%center_lats)) DEALLOCATE(wrf_state_variables(i)%center_lats)
   IF(ALLOCATED(wrf_state_variables(i)%center_lons)) DEALLOCATE(wrf_state_variables(i)%center_lons)
   IF(ALLOCATED(wrf_state_variables(i)%mask)) DEALLOCATE(wrf_state_variables(i)%mask)
 END DO
END SUBROUTINE release_wrf_buf


!! CCPL
#ifdef CCPL_DA
SUBROUTINE allocate_da_buf_fields(wrf_grid_temp)

   implicit none
   TYPE(domain) , POINTER :: wrf_grid_temp
   integer                :: wrf_grid_temp_id
   integer                :: ims,ime,jms,jme,nsig,nsig_soil

   ims = wrf_grid_temp%sm31
   ime = wrf_grid_temp%em31
   jms = wrf_grid_temp%sm33
   jme = wrf_grid_temp%em33
   wrf_grid_temp_id = wrf_grid_temp%id
   nsig = wrf_grid_temp%ed32-1  
   nsig_soil = size(wrf_grid_temp%smois,2)
   
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_pt)
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%regional_time(6))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%aeta1(1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%aeta2(1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%eta1(1:nsig+1))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%eta2(1:nsig+1))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_dx_mc(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_dy_mc(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_xlat(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_xlong(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_mu(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_mub(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_psfc0(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_phb(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_landmask(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_seaice(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_sst(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_ivgtyp(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_isltyp(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_vegfra(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_snow(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_u10(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_v10(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_smois(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_tslb(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_tsk(ims:ime,jms:jme))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_t(ims:ime,jms:jme,1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_u(ims:ime,jms:jme,1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_v(ims:ime,jms:jme,1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_q(ims:ime,jms:jme,1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_w(ims:ime,jms:jme,1:nsig))
   ALLOCATE(wrf_da_variables(wrf_grid_temp_id)%ccpl_ph(ims:ime,jms:jme,1:nsig))
   wrf_da_variables(wrf_grid_temp_id)%ccpl_pt = wrf_grid_temp%p_top
   wrf_da_variables(wrf_grid_temp_id)%regional_time = 0.0
   wrf_da_variables(wrf_grid_temp_id)%aeta1 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%aeta2 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%eta1 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%eta2 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_dx_mc = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_dy_mc = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_xlat = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_xlong = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_mu = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_mub = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_psfc0 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_phb = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_landmask = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_seaice = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_sst = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_ivgtyp = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_isltyp = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_vegfra = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_snow = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_u10 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_v10 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_smois = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_tslb = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_tsk = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_t = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_u = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_v = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_q = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_w = 0.0
   wrf_da_variables(wrf_grid_temp_id)%ccpl_ph = 0.0
      
END SUBROUTINE allocate_da_buf_fields

SUBROUTINE send_da_fields(wrf_grid_temp)

   USE ESMF_mod, ONLY : ESMF_TimeGet
   USE module_domain, ONLY : domain_clock_get
   USE module_utility, ONLY : WRFU_Time
   use module_state_description, ONLY: P_QV
   implicit none

   TYPE(domain) , POINTER :: wrf_grid_temp
   TYPE(WRFU_Time)        :: currentTime
   integer                :: wrf_grid_temp_id
   integer                :: i,j,k
   integer                :: iyear,imonth,iday,ihour,iminute,isecond
   integer                :: ims,ime,jms,jme,nsig,nsig_soil

   ims = wrf_grid_temp%sm31
   ime = wrf_grid_temp%em31
   jms = wrf_grid_temp%sm33
   jme = wrf_grid_temp%em33
   wrf_grid_temp_id = wrf_grid_temp%id
   nsig = wrf_grid_temp%ed32-1  
   nsig_soil = size(wrf_grid_temp%smois,2)

   call domain_clock_get(wrf_grid_temp, current_time=currentTime)
   call ESMF_TimeGet(currentTime, YY=iyear,MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond) 
        
   wrf_da_variables(wrf_grid_temp_id)%ccpl_pt = wrf_grid_temp%p_top
   wrf_da_variables(wrf_grid_temp_id)%regional_time(1) = iyear
   wrf_da_variables(wrf_grid_temp_id)%regional_time(2) = imonth
   wrf_da_variables(wrf_grid_temp_id)%regional_time(3) = iday
   wrf_da_variables(wrf_grid_temp_id)%regional_time(4) = ihour
   wrf_da_variables(wrf_grid_temp_id)%regional_time(5) = iminute
   wrf_da_variables(wrf_grid_temp_id)%regional_time(6) = isecond
   wrf_da_variables(wrf_grid_temp_id)%aeta1 = wrf_grid_temp%znu(1:nsig)
   wrf_da_variables(wrf_grid_temp_id)%aeta2 = 0.0
   wrf_da_variables(wrf_grid_temp_id)%eta1 = wrf_grid_temp%znw(1:nsig+1)
   wrf_da_variables(wrf_grid_temp_id)%eta2 = 0.0

   wrf_da_variables(wrf_grid_temp_id)%ccpl_dx_mc(ims:ime,jms:jme)    = 1.0/(wrf_grid_temp%msft(ims:ime,jms:jme)*wrf_grid_temp%rdx)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_dy_mc(ims:ime,jms:jme)    = 1.0/(wrf_grid_temp%msft(ims:ime,jms:jme)*wrf_grid_temp%rdy)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_xlat(ims:ime,jms:jme)     = wrf_grid_temp%xlat(ims:ime,jms:jme)  !XLAT: LATITUDE, SOUTH IS NEGATIVE
   wrf_da_variables(wrf_grid_temp_id)%ccpl_xlong(ims:ime,jms:jme)    = wrf_grid_temp%xlong(ims:ime,jms:jme)  !XLONG: LONGITUDE, WEST IS NEGATIVE
   wrf_da_variables(wrf_grid_temp_id)%ccpl_mu(ims:ime,jms:jme)       = wrf_grid_temp%mu_2(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_mub(ims:ime,jms:jme)      = wrf_grid_temp%mub(ims:ime,jms:jme) 
   wrf_da_variables(wrf_grid_temp_id)%ccpl_psfc0(ims:ime,jms:jme)    = wrf_grid_temp%mu_2(ims:ime,jms:jme)+wrf_grid_temp%mub(ims:ime,jms:jme)+wrf_grid_temp%p_top
   wrf_da_variables(wrf_grid_temp_id)%ccpl_phb(ims:ime,jms:jme)      = wrf_grid_temp%phb(ims:ime,1,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_seaice(ims:ime,jms:jme)   = wrf_grid_temp%xice(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_sst(ims:ime,jms:jme)      = wrf_grid_temp%sst(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_ivgtyp(ims:ime,jms:jme)   = wrf_grid_temp%ivgtyp(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_isltyp(ims:ime,jms:jme)   = wrf_grid_temp%isltyp(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_vegfra(ims:ime,jms:jme)   = wrf_grid_temp%vegfra(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_snow(ims:ime,jms:jme)     = wrf_grid_temp%snow(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_u10(ims:ime,jms:jme)      = wrf_grid_temp%u10(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_v10(ims:ime,jms:jme)      = wrf_grid_temp%v10(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_smois(ims:ime,jms:jme)    = wrf_grid_temp%smois(ims:ime,1,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_tslb(ims:ime,jms:jme)     = wrf_grid_temp%tslb(ims:ime,1,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_tsk(ims:ime,jms:jme)      = wrf_grid_temp%tsk(ims:ime,jms:jme)
   wrf_da_variables(wrf_grid_temp_id)%ccpl_landmask(ims:ime,jms:jme) = wrf_grid_temp%xland(ims:ime,jms:jme)
   where(abs(wrf_da_variables(wrf_grid_temp_id)%ccpl_landmask-2.0)<0.01) wrf_da_variables(wrf_grid_temp_id)%ccpl_landmask = 0.0
   do k=1, nsig
     !wrf_da_variables(wrf_grid_temp_id)%ccpl_t(ims:ime,jms:jme,k) = wrf_grid_temp%th_phy_m_t0(ims:ime,k,jms:jme)+300.0
     wrf_da_variables(wrf_grid_temp_id)%ccpl_t(ims:ime,jms:jme,k)   = wrf_grid_temp%t_2(ims:ime,k,jms:jme)+300.0
     wrf_da_variables(wrf_grid_temp_id)%ccpl_u(ims:ime-1,jms:jme,k) = 0.5*(wrf_grid_temp%u_2(ims:ime-1,k,jms:jme)+wrf_grid_temp%u_2(ims+1:ime,k,jms:jme))
     !wrf_da_variables(wrf_grid_temp_id)%ccpl_u(ime,jms:jme,k)       = 1.5*wrf_grid_temp%u_2(ime,k,jms:jme)-0.5*wrf_grid_temp%u_2(ime-1,k,jms:jme)
     wrf_da_variables(wrf_grid_temp_id)%ccpl_v(ims:ime,jms:jme-1,k) = 0.5*(wrf_grid_temp%v_2(ims:ime,k,jms:jme-1)+wrf_grid_temp%v_2(ims:ime,k,jms+1:jme))
     !wrf_da_variables(wrf_grid_temp_id)%ccpl_v(ims:ime,jme,k)       = 1.5*wrf_grid_temp%v_2(ims:ime,k,jme)-0.5*wrf_grid_temp%v_2(ims:ime,k,jme-1)
     wrf_da_variables(wrf_grid_temp_id)%ccpl_q(ims:ime,jms:jme,k)   = wrf_grid_temp%moist(ims:ime,k,jms:jme,P_QV)
     wrf_da_variables(wrf_grid_temp_id)%ccpl_w(ims:ime,jms:jme,k)   = 0.5*(wrf_grid_temp%w_2(ims:ime,k,jms:jme)+wrf_grid_temp%w_2(ims:ime,k+1,jms:jme))
     wrf_da_variables(wrf_grid_temp_id)%ccpl_ph(ims:ime,jms:jme,k)  = 0.5*(wrf_grid_temp%ph_2(ims:ime,k,jms:jme)+wrf_grid_temp%ph_2(ims:ime,k+1,jms:jme))
   end do
  
END SUBROUTINE send_da_fields

SUBROUTINE receive_da_fields(wrf_grid_temp)

   USE ESMF_mod, ONLY : ESMF_TimeGet
   USE module_domain, ONLY : domain_clock_get
   USE module_utility, ONLY : WRFU_Time
   use module_state_description, ONLY: P_QV
   implicit none

   TYPE(domain) , POINTER :: wrf_grid_temp
   TYPE(WRFU_Time)        :: currentTime
   integer                :: wrf_grid_temp_id
   integer                :: i,j,k
   integer                :: ims,ime,jms,jme,nsig,nsig_soil
   integer                :: ipe,ips,jpe,jps,ide,ids,jde,jds

   ims = wrf_grid_temp%sm31
   ime = wrf_grid_temp%em31
   jms = wrf_grid_temp%sm33
   jme = wrf_grid_temp%em33
     
   ips = wrf_grid_temp%sp31
   ipe = wrf_grid_temp%ep31  
   jps = wrf_grid_temp%sp33
   jpe = wrf_grid_temp%ep33

   ide = wrf_grid_temp%ed31
   jde = wrf_grid_temp%ed33

   wrf_grid_temp_id = wrf_grid_temp%id
   nsig = wrf_grid_temp%ed32-1  
   nsig_soil = size(wrf_grid_temp%smois,2)
   
   !U, V, T
   do k=1, nsig
     wrf_grid_temp%t_2(ims:ime,k,jms:jme) = wrf_da_variables(wrf_grid_temp_id)%ccpl_t(ims:ime,jms:jme,k)-300.0
     !wrf_grid_temp%moist(ims:ime,k,jms:jme,P_QV) = wrf_da_variables(wrf_grid_temp_id)%ccpl_q(ims:ime,jms:jme,k)
     wrf_grid_temp%u_2(ips:ipe-1,k,jps:jpe) = 0.5*(wrf_da_variables(wrf_grid_temp_id)%ccpl_u(ips-1:ipe-2,jps:jpe,k) + wrf_da_variables(wrf_grid_temp_id)%ccpl_u(ips:ipe,jps:jpe,k))
     wrf_grid_temp%v_2(ips:ipe,k,jps:jpe-1) = 0.5*(wrf_da_variables(wrf_grid_temp_id)%ccpl_v(ips:ipe,jps-1:jpe-2,k) + wrf_da_variables(wrf_grid_temp_id)%ccpl_v(ips:ipe,jps:jpe-1,k))
   end do


END SUBROUTINE receive_da_fields

SUBROUTINE release_wrf_da_buf

   USE module_configure, ONLY : model_config_rec 
   implicit none

   integer  :: i
   DO i=1,model_config_rec%max_dom
     IF(ALLOCATED(wrf_da_variables(i)%regional_time)) DEALLOCATE(wrf_da_variables(i)%regional_time)
     IF(ALLOCATED(wrf_da_variables(i)%aeta1)) DEALLOCATE(wrf_da_variables(i)%aeta1)
     IF(ALLOCATED(wrf_da_variables(i)%aeta2)) DEALLOCATE(wrf_da_variables(i)%aeta2)
     IF(ALLOCATED(wrf_da_variables(i)%eta1)) DEALLOCATE(wrf_da_variables(i)%eta1)
     IF(ALLOCATED(wrf_da_variables(i)%eta2)) DEALLOCATE(wrf_da_variables(i)%eta2)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_dx_mc)) DEALLOCATE(wrf_da_variables(i)%ccpl_dx_mc)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_dy_mc)) DEALLOCATE(wrf_da_variables(i)%ccpl_dy_mc)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_xlat)) DEALLOCATE(wrf_da_variables(i)%ccpl_xlat)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_xlong)) DEALLOCATE(wrf_da_variables(i)%ccpl_xlong)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_mu)) DEALLOCATE(wrf_da_variables(i)%ccpl_mu)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_mub)) DEALLOCATE(wrf_da_variables(i)%ccpl_mub)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_psfc0)) DEALLOCATE(wrf_da_variables(i)%ccpl_psfc0)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_phb)) DEALLOCATE(wrf_da_variables(i)%ccpl_phb)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_landmask)) DEALLOCATE(wrf_da_variables(i)%ccpl_landmask)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_seaice)) DEALLOCATE(wrf_da_variables(i)%ccpl_seaice)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_sst)) DEALLOCATE(wrf_da_variables(i)%ccpl_sst)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_ivgtyp)) DEALLOCATE(wrf_da_variables(i)%ccpl_ivgtyp)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_isltyp)) DEALLOCATE(wrf_da_variables(i)%ccpl_isltyp)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_vegfra)) DEALLOCATE(wrf_da_variables(i)%ccpl_vegfra)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_snow)) DEALLOCATE(wrf_da_variables(i)%ccpl_snow)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_u10)) DEALLOCATE(wrf_da_variables(i)%ccpl_u10)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_v10)) DEALLOCATE(wrf_da_variables(i)%ccpl_v10)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_smois)) DEALLOCATE(wrf_da_variables(i)%ccpl_smois)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_tslb)) DEALLOCATE(wrf_da_variables(i)%ccpl_tslb)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_tsk)) DEALLOCATE(wrf_da_variables(i)%ccpl_tsk)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_t)) DEALLOCATE(wrf_da_variables(i)%ccpl_t)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_u)) DEALLOCATE(wrf_da_variables(i)%ccpl_u)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_v)) DEALLOCATE(wrf_da_variables(i)%ccpl_v)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_q)) DEALLOCATE(wrf_da_variables(i)%ccpl_q)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_w)) DEALLOCATE(wrf_da_variables(i)%ccpl_w)
     IF(ALLOCATED(wrf_da_variables(i)%ccpl_ph)) DEALLOCATE(wrf_da_variables(i)%ccpl_ph)
   END DO

END SUBROUTINE release_wrf_da_buf
#endif
!! CCPL

END MODULE module_wrf_ccpl 
