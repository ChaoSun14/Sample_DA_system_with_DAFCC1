! sunchao 2019-03-21 
! 
! abstract: This module hold specification kinds for variable declaration 
!           when gsi is used as external link library for C-Coupler. 
 


module gsi_ccpl_kinds

    use kinds, only: i_kind,r_kind, r_single
    
    implicit none
    
    logical, public  :: ccpl_run = .False. ! if use online data transmission
    logical, public  :: io_read  = .True. ! if read files
    logical, public  :: io_check = .False. ! if check io_read and online transmission
    !=====================================================================
    !exe_run : ccpl_run = F, io_read = T, io_check = F 
    !lib.so  : 
    !   1) ccpl_run = T <-- use online data transmission 
    !      io_read  = F <-- do not read files
    !      io_check = F <-- do not check io_read and online transmission
    !
    !   2) ccpl_run = F <-- do not use online data transmission 
    !      io_read  = T <-- read files
    !      io_check = T <-- check io_read and online transmission
    !=====================================================================

    integer, public  :: gsi_ccpl_comm
    integer, public  :: gsi_ccpl_comm_all, gsi_grid_ccpl_comm, color, gsi_grid_instance_id, gsi_mem_instance_id, gsi_ensmean_instance_id, gsi_comp_id, gsi_grid_id, gsi_grid_id_input, gsi_vert_grid1, gsi_vert_grid2, gsi_3d_grid1, gsi_3d_grid2, gsi_regional_time_grid_id, gsi_global_decomp_id, gsi_instance_id, gsi_decomp_id, ccpl_npe_all, ccpl_mype_all, ccpl_npe_wrf, ccpl_npe_grid, ccpl_mype_grid, ccpl_mem_id
    integer, public  :: ccpl_nlon_regional, ccpl_nlat_regional, ccpl_nsig, ccpl_nsig_soil 
    integer, pointer, public        :: ccpl_regional_time(:)=>null(), ccpl_regional_time_grid(:)=>null()
    real(r_kind), pointer, public   :: ccpl_pt=>null()
    real(r_single), pointer, public :: ccpl_aeta1(:)=>null(), ccpl_eta1(:)=>null(), ccpl_aeta2(:)=>null(), ccpl_eta2(:)=>null()
    real(r_single), pointer, public :: ccpl_glat(:,:)=>null(), ccpl_glon(:,:)=>null(), ccpl_dx_mc(:,:)=>null(), ccpl_dy_mc(:,:)=>null()
    real(r_kind), pointer, public   :: ccpl_psfc(:,:)=>null(), ccpl_snow(:,:)=>null(), ccpl_fis(:,:)=>null(), ccpl_landmask(:,:)=>null(), ccpl_xice(:,:)=>null(), ccpl_sst(:,:)=>null(), ccpl_tsk(:,:)=>null(), ccpl_tslb(:,:)=>null(), ccpl_smois(:,:)=>null(), ccpl_ivgtyp(:,:)=>null(), ccpl_vegfrac(:,:)=>null(), ccpl_isltyp(:,:)=>null()
    real(r_kind), pointer, public   :: ccpl_u(:,:,:)=>null(), ccpl_v(:,:,:)=>null(), ccpl_q(:,:,:)=>null(), ccpl_pot(:,:,:)=>null()
    integer, allocatable, public    :: local_cells_global_index(:), local_cells_global_index_with_buffer(:)

    contains

    subroutine gsi_ccpl_kinds_cleanup

        implicit none
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        write(*,*) "[CCPL <GSI>] Start deallocate fields pointer"
        if (associated(ccpl_psfc))      deallocate(ccpl_psfc)
        if (associated(ccpl_snow))      deallocate(ccpl_snow)
        if (associated(ccpl_fis))       deallocate(ccpl_fis)
        if (associated(ccpl_landmask))  deallocate(ccpl_landmask)
        if (associated(ccpl_xice))      deallocate(ccpl_xice)
        if (associated(ccpl_sst))       deallocate(ccpl_sst)
        if (associated(ccpl_tsk))       deallocate(ccpl_tsk)
        if (associated(ccpl_tslb))      deallocate(ccpl_tslb)
        if (associated(ccpl_smois))     deallocate(ccpl_smois)
        if (associated(ccpl_ivgtyp))    deallocate(ccpl_ivgtyp)
        if (associated(ccpl_vegfrac))   deallocate(ccpl_vegfrac)
        if (associated(ccpl_isltyp))    deallocate(ccpl_isltyp)
        if (associated(ccpl_u))         deallocate(ccpl_u)
        if (associated(ccpl_v))         deallocate(ccpl_v)
        if (associated(ccpl_q))         deallocate(ccpl_q)
        if (associated(ccpl_pot))       deallocate(ccpl_pot)
        if (allocated(local_cells_global_index))  deallocate(local_cells_global_index)
        if (allocated(local_cells_global_index_with_buffer))  deallocate(local_cells_global_index_with_buffer)
        write(*,*) "[CCPL <GSI>] Finish deallocate fields pointer"
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"

    end subroutine gsi_ccpl_kinds_cleanup

end module
