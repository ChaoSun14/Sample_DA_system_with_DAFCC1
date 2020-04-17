module enkf_ccpl_kinds

    use kinds,     only: i_kind, r_kind, r_single, i_long, r_double
    
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
    
    integer, public  :: enkf_ccpl_comm, enkf_ccpl_comm_grid, enkf_npe_grid, enkf_mype_grid
    integer, public  :: ccpl_ens_num, ccpl_nlon, ccpl_nlat, ccpl_nsig
    integer, public  :: enkf_grid_instance_id, enkf_grid_comp_id, enkf_h2d_grid_id_input, enkf_h2d_md_grid_id_input, enkf_h2d_md_grid_id, enkf_vert_md_grid_id1, enkf_global_decomp_id
    integer, public  :: enkf_instance_id, enkf_comp_id, enkf_h2d_id, enkf_vert_grid1, enkf_v3d_grid1, enkf_ens_num_grid, enkf_h2d_md_id, enkf_v3d_md_id, enkf_decomp_id
    integer, allocatable, public    :: local_cells_global_index(:)
    integer, pointer, public        :: ccpl_nlon_ens(:)=>null(), ccpl_nlat_ens(:)=>null(), ccpl_nsig_ens(:)=>null()
    real, pointer, public			:: global_ccpl_xlat(:,:,:)=>null(), global_ccpl_xlon(:,:,:)=>null(), global_ccpl_mu(:,:,:)=>null(), global_ccpl_mub(:,:,:)=>null(), ccpl_znu(:,:)=>null()
    real(r_single), pointer, public :: ccpl_pt_ens(:)=>null(), ccpl_pt=>null() 
    real(r_single), pointer, public :: ccpl_u(:,:,:)=>null(), ccpl_v(:,:,:)=>null(), ccpl_t(:,:,:)=>null(), ccpl_q(:,:,:)=>null(), ccpl_w(:,:,:)=>null(), ccpl_ph(:,:,:)=>null(), ccpl_mu(:,:)=>null(), ccpl_mub(:,:)=>null(), ccpl_psfc(:,:)=>null()
    real(r_single), allocatable, public    :: ccpl_anal_chunk(:,:,:)
    real(r_double), allocatable, public    :: ccpl_qsat_chunk(:,:,:)
    !real(r_single), public :: ccpl_pt

    contains
    
    subroutine enkf_ccpl_kinds_cleanup

        implicit none
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        write(*,*) "[CCPL <ENKF>] Start deallocate fields pointer"
        if (associated(ccpl_pt_ens))      deallocate(ccpl_pt_ens)
        if (associated(ccpl_pt))          deallocate(ccpl_pt)
        if (associated(ccpl_u))           deallocate(ccpl_u)
        if (associated(ccpl_v))           deallocate(ccpl_v)
        if (associated(ccpl_t))           deallocate(ccpl_t)
        if (associated(ccpl_t))           deallocate(ccpl_t)
        if (associated(ccpl_q))           deallocate(ccpl_q)
        if (associated(ccpl_w))           deallocate(ccpl_w)
        if (associated(ccpl_ph))          deallocate(ccpl_ph)
        if (associated(ccpl_mu))          deallocate(ccpl_mu)
        if (associated(ccpl_mub))         deallocate(ccpl_mub)
        if (associated(ccpl_psfc))        deallocate(ccpl_psfc)

        if (allocated(local_cells_global_index))  deallocate(local_cells_global_index)
        if (allocated(ccpl_anal_chunk))           deallocate(ccpl_anal_chunk)
        if (allocated(ccpl_qsat_chunk))           deallocate(ccpl_qsat_chunk)
        write(*,*) "[CCPL <ENKF>] Finish deallocate fields pointer"
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"

    end subroutine enkf_ccpl_kinds_cleanup

end module
