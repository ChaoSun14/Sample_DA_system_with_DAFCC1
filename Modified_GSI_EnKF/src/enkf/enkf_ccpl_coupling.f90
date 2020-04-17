! sunchao 2019-03-21 
! 
! abstract: This module is 

module enkf_ccpl_coupling_mod
    
    use CCPL_interface_mod
    use enkf_ccpl_kinds, only: enkf_ccpl_comm, ccpl_ens_num
    use enkf_ccpl_kinds, only: enkf_instance_id, enkf_comp_id, enkf_h2d_id, enkf_vert_grid1, enkf_v3d_grid1, enkf_ens_num_grid, enkf_h2d_md_id, enkf_v3d_md_id, enkf_decomp_id
    use enkf_ccpl_kinds, only: ccpl_nlon, ccpl_nlat, ccpl_nsig, global_ccpl_xlat, global_ccpl_xlon, global_ccpl_mu, global_ccpl_mub, ccpl_znu, ccpl_pt
    use enkf_ccpl_kinds, only: ccpl_u, ccpl_v, ccpl_t, ccpl_q, ccpl_w, ccpl_ph, ccpl_mu, ccpl_mub, ccpl_psfc
    use enkf_ccpl_kinds, only: local_cells_global_index, ccpl_anal_chunk, ccpl_qsat_chunk
    use enkf_ccpl_kinds, only: enkf_ccpl_comm_grid, enkf_npe_grid, enkf_mype_grid
    use kinds,    only: r_double, r_kind, r_single, i_kind 
    !use mpisetup, only: numproc, nproc
    !use loadbal, only: npts_max, indxproc, numptsperproc  
    
    implicit none
    integer, private :: numproc, nproc, ierror
    integer(i_kind), private, allocatable, dimension(:)   :: numptsperproc
    integer(i_kind), private, allocatable, dimension(:,:) :: indxproc
contains
    subroutine register_grids_decomps

        use mpi
        implicit none
         
        integer              :: n,np,npts,npts_max,i,j,kk
        integer, allocatable :: mask(:)
        real(r_kind)         :: max_lon, min_lon, max_lat, min_lat
        real(r_kind), allocatable :: xlat(:), xlon(:), mype_lat(:), mype_lon(:), x_global_ccpl_mu(:,:), x_global_ccpl_mub(:,:)
    
        call mpi_bcast(global_ccpl_xlat, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        call mpi_bcast(global_ccpl_xlon, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        call mpi_bcast(global_ccpl_mu, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        call mpi_bcast(global_ccpl_mub, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        !write(*,*) "global_ccpl_mu(1,1,1), (ccpl_nlon,ccpl_nlat,1): ", global_ccpl_mu(1,1,1), global_ccpl_mu(ccpl_nlon,ccpl_nlat,1)
        
        enkf_comp_id = CCPL_external_procedures_get_comp_ID(enkf_instance_id, annotation="get enkf comp id")
        !enkf_h2d_md_id_input = CCPL_external_procedures_para_get_grid_ID(enkf_instance_id, 1)
        !enkf_v3d_md_id1 = CCPL_external_procedures_para_get_grid_ID(enkf_instance_id, 4)
        !enkf_h2d_id_input = CCPL_external_procedures_para_get_grid_ID(enkf_instance_id, 7)
        npts=ccpl_nlon*ccpl_nlat
        call mpi_comm_size(enkf_ccpl_comm,numproc,ierror)
        call mpi_comm_rank(enkf_ccpl_comm,nproc,ierror)
        allocate(numptsperproc(numproc))
        numptsperproc = 0
        np = 0
        do n=1,npts
            np = np + 1
            if (np > numproc) np = 1
            numptsperproc(np) = numptsperproc(np)+1
        end do
        npts_max = maxval(numptsperproc)
        allocate(indxproc(numproc,npts_max))
        numptsperproc = 0
        np = 0
        do n=1,npts
            np = np + 1
            if (np > numproc) np = 1
            numptsperproc(np) = numptsperproc(np)+1
            indxproc(np,numptsperproc(np)) = n
        end do

        allocate(mask(numptsperproc(nproc+1)))
        allocate(mype_lat(numptsperproc(nproc+1)),mype_lon(numptsperproc(nproc+1)))
        allocate(xlat(ccpl_nlon*ccpl_nlat), xlon(ccpl_nlon*ccpl_nlat))
        allocate(x_global_ccpl_mu(ccpl_nlon*ccpl_nlat,ccpl_ens_num), x_global_ccpl_mub(ccpl_nlon*ccpl_nlat,ccpl_ens_num))
        allocate(ccpl_mu(numptsperproc(nproc+1),ccpl_ens_num), ccpl_mub(numptsperproc(nproc+1),ccpl_ens_num))
        kk = 1
        do j = 1, ccpl_nlat
            do i = 1, ccpl_nlon
               xlat(kk) = global_ccpl_xlat(i,j,1)
               xlon(kk) = global_ccpl_xlon(i,j,1)
               x_global_ccpl_mu(kk,:) = global_ccpl_mu(i,j,:)
               x_global_ccpl_mub(kk,:) = global_ccpl_mub(i,j,:)
               kk = kk+1
            end do
        end do
        mask = 1
        do i = 1, numptsperproc(nproc+1)
            mype_lat(i)=xlat(indxproc(nproc+1,i))
            mype_lon(i)=xlon(indxproc(nproc+1,i))
            ccpl_mu(i,:)=x_global_ccpl_mu(indxproc(nproc+1,i),:)
            ccpl_mub(i,:)=x_global_ccpl_mub(indxproc(nproc+1,i),:)
        end do
        deallocate(xlat, xlon)
        deallocate(x_global_ccpl_mu, x_global_ccpl_mub)
        max_lon = -999999.!2*PI
        min_lon = -999999.!0
        max_lat = -999999.!PI/2.0
        min_lat = -999999.!-PI/2.0
        enkf_h2d_id = CCPL_register_H2D_grid_via_local_data(enkf_comp_id,"enkf_H2D_grid", "LON_LAT",  "degrees", "acyclic", ccpl_nlon*ccpl_nlat, numptsperproc(nproc+1), indxproc(nproc+1,:), min_lon, max_lon, min_lat, max_lat, mype_lon, mype_lat, mask, annotation="register H2D grid for enkf" )
        enkf_decomp_id = CCPL_register_normal_parallel_decomp("enkf_decomp", enkf_h2d_id, numptsperproc(nproc+1), indxproc(nproc+1,:), annotation = "register parallel decomp for enkf")
        enkf_vert_grid1 = CCPL_register_V1D_Z_grid_via_model_data(enkf_comp_id, "enkf_vertical_grid_nsig", "Pa", ccpl_znu(1:ccpl_nsig,1), annotation="register simplified vertical grid 1 for enkf")
        enkf_v3d_grid1 = CCPL_register_MD_grid_via_multi_grids(enkf_comp_id, "enkf_3d_grid_1", enkf_h2d_id, enkf_vert_grid1, annotation="register enkf 3D-grid with nsig vertical levels")
        enkf_ens_num_grid = CCPL_register_1D_normal_grid_without_data(enkf_comp_id, "enkf_ensemble_number_grid", ccpl_ens_num, annotation="register virtual ensemble number grid for enkf")
        enkf_h2d_md_id = CCPL_register_MD_grid_via_multi_grids(enkf_comp_id, "enkf_h2d_md_grid", enkf_h2d_id, enkf_ens_num_grid, annotation="register enkf h2d grid with ensemble number grid")
        enkf_v3d_md_id = CCPL_register_MD_grid_via_multi_grids(enkf_comp_id, "enkf_v3d_md_grid", enkf_v3d_grid1, enkf_ens_num_grid, annotation="register enkf v2d grid with ensemble number grid")
        !enkf_h2d_md_id = CCPL_register_MD_grid_via_multi_grids(enkf_comp_id, "enkf_h2d_md_grid", enkf_ens_num_grid, enkf_h2d_id, annotation="register enkf h2d grid with ensemble number grid")
        !enkf_v3d_md_id = CCPL_register_MD_grid_via_multi_grids(enkf_comp_id, "enkf_v3d_md_grid", enkf_ens_num_grid, enkf_v3d_grid1, annotation="register enkf v2d grid with ensemble number grid")
    end subroutine register_grids_decomps

    subroutine enkf_fields_declare

        use constants, only: zero, one, three, deg2rad, half, two, r0_01
        implicit none
        integer :: field_id_u, field_id_v, field_id_t, field_id_q, field_id_w, field_id_ph, field_id_mu, field_id_mub, field_id_psfc
        integer :: dims_size_2d(2), dims_size_3d(3)

        write(*,*) "+++++++ start declare fields +++++++++++"
        write(*,*) "numptsperproc(nproc+1): ", numptsperproc(nproc+1)
        write(*,*) "ccpl_ens_num: ", ccpl_ens_num

        
        dims_size_2d(1)= numptsperproc(nproc+1)
        dims_size_2d(2)= ccpl_ens_num
        
        dims_size_3d(1)= numptsperproc(nproc+1)
        dims_size_3d(2)= ccpl_nsig
        dims_size_3d(3)= ccpl_ens_num
        
        
        
        field_id_psfc = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_psfc, "psfc0", CCPL_PARA_TYPE_IN, enkf_decomp_id, enkf_h2d_md_id, dims_size_2d, annotation="declare surface pressure")
        !field_id_mu = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_mu, "mu", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_h2d_md_id, dims_size_2d, annotation="declare perturbation dry air mass in column") 
        !field_id_mub = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_mub, "mub", CCPL_PARA_TYPE_IN, enkf_decomp_id, enkf_h2d_md_id, dims_size_2d, annotation="declare base state dry air mass in column") 
        field_id_u = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_u, "ccpl_u", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_v3d_md_id, dims_size_3d, annotation="register x-wind component at mass point")
        field_id_v = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_v, "ccpl_v", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_v3d_md_id, dims_size_3d, annotation="register y-wind component at mass point")
        field_id_q = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_q, "ccpl_q", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_v3d_md_id, dims_size_3d, annotation="register water vapor mixing ratio")
        field_id_t = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_t, "ccpl_t", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_v3d_md_id, dims_size_3d, annotation="declare perturbation potential temperature theta-t0")
        field_id_w = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_w, "ccpl_w", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_v3d_md_id, dims_size_3d, annotation="declare z-wind component")
        field_id_ph = CCPL_external_procedures_para_declare_field(enkf_instance_id, ccpl_ph, "ccpl_ph", CCPL_PARA_TYPE_INOUT, enkf_decomp_id, enkf_v3d_md_id, dims_size_3d, annotation="declare perturbation geopotential")
         
        write(*,*) "++++++++ finish declare fields +++++++ "
!        write(*,*) "glat= ", ccpl_nlat_regional

    end subroutine enkf_fields_declare

    subroutine enkf_fields_transfer_in

        use params, only: nlevs,nvars,ndim,nbackgrounds,nanals,pseudo_rh
        use constants, only: zero,one,cp,fv,rd,grav,rad2deg
        use mpi
        implicit none
    
        integer :: i,j,k,kk,n
        logical :: ice
        real(r_single), dimension(:,:), allocatable              :: enkf_virttemp
        real(r_single), dimension(:,:), allocatable              :: enkf_pressure
        real(r_single), dimension(:,:), allocatable              :: enkf_spechumd
        real(r_kind), allocatable ::  x_global_ccpl_mu(:,:)

        call mpi_bcast(global_ccpl_xlat, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        call mpi_bcast(global_ccpl_xlon, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        call mpi_bcast(global_ccpl_mu, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        call mpi_bcast(global_ccpl_mub, ccpl_nlon*ccpl_nlat*ccpl_ens_num, MPI_REAL4, 0, enkf_ccpl_comm_grid, ierror) 
        allocate(x_global_ccpl_mu(ccpl_nlon*ccpl_nlat,ccpl_ens_num))
        kk = 1
        do j = 1, ccpl_nlat
            do i = 1, ccpl_nlon
               x_global_ccpl_mu(kk,:) = global_ccpl_mu(i,j,:)
               kk = kk+1
            end do
        end do
        do i = 1, numptsperproc(nproc+1)
            ccpl_mu(i,:)=x_global_ccpl_mu(indxproc(nproc+1,i),:)
        end do
        deallocate(x_global_ccpl_mu)
        !global_ccpl_xlat = global_ccpl_xlat*rad2deg
        !global_ccpl_xlon = global_ccpl_xlon*rad2deg
        ccpl_t = ccpl_t-300. 
        if (allocated(ccpl_anal_chunk)) deallocate(ccpl_anal_chunk)
        allocate(ccpl_anal_chunk(1:nanals,1:numptsperproc(nproc+1),1:ndim))
        
        ! nvars=3 Updating U, V, T, and MU
        ! nvars=4 Updating U, V, T, QVAPOR, and MU
        ! nvars=5 Updating U, V, T, QVAPOR, PH, and MU
        ! nvars=6 Updating U, V, T, QVAPOR, W, PH, and MU
        if (nvars .eq. 3) then 
            do i = 1, nanals              
                ccpl_anal_chunk(i,:,1:ccpl_nsig) = ccpl_u(:,:,i)
                ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) = ccpl_v(:,:,i)
                ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig) = ccpl_t(:,:,i)
                ccpl_anal_chunk(i,:,3*ccpl_nsig+1) = ccpl_mu(:,i)
            end do
        else if (nvars .eq. 4) then 
            do i = 1, nanals              
                ccpl_anal_chunk(i,:,1:ccpl_nsig) = ccpl_u(:,:,i)
                ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) = ccpl_v(:,:,i)
                ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig) = ccpl_t(:,:,i)
                ccpl_anal_chunk(i,:,3*ccpl_nsig+1:4*ccpl_nsig) = ccpl_q(:,:,i)
                ccpl_anal_chunk(i,:,4*ccpl_nsig+1) = ccpl_mu(:,i)
            end do
        else if (nvars .eq. 5) then 
            do i = 1, nanals   
                ccpl_anal_chunk(i,:,1:ccpl_nsig) = ccpl_u(:,:,i)
                ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) = ccpl_v(:,:,i)
                ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig) = ccpl_t(:,:,i)
                ccpl_anal_chunk(i,:,3*ccpl_nsig+1:4*ccpl_nsig) = ccpl_q(:,:,i)
                ccpl_anal_chunk(i,:,4*ccpl_nsig+1:5*ccpl_nsig) = ccpl_ph(:,:,i)
                ccpl_anal_chunk(i,:,5*ccpl_nsig+1) = ccpl_mu(:,i)
            end do   
        else if (nvars .eq. 6) then 
            do i = 1, nanals              
                ccpl_anal_chunk(i,:,1:ccpl_nsig) = ccpl_u(:,:,i)
                ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) = ccpl_v(:,:,i)
                ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig) = ccpl_t(:,:,i)
                ccpl_anal_chunk(i,:,3*ccpl_nsig+1:4*ccpl_nsig) = ccpl_q(:,:,i)
                ccpl_anal_chunk(i,:,4*ccpl_nsig+1:5*ccpl_nsig) = ccpl_w(:,:,i)
                ccpl_anal_chunk(i,:,5*ccpl_nsig+1:6*ccpl_nsig) = ccpl_ph(:,:,i)
                ccpl_anal_chunk(i,:,6*ccpl_nsig+1) = ccpl_mu(:,i)
            end do 
        end if 
        

        if (allocated(ccpl_qsat_chunk)) deallocate(ccpl_qsat_chunk)
        allocate(ccpl_qsat_chunk(nanals,numptsperproc(nproc+1),ccpl_nsig))
        if(.not. allocated(enkf_virttemp)) allocate(enkf_virttemp(numptsperproc(nproc+1),ccpl_nsig))
        if(.not. allocated(enkf_pressure)) allocate(enkf_pressure(numptsperproc(nproc+1),ccpl_nsig))
        if(.not. allocated(enkf_spechumd)) allocate(enkf_spechumd(numptsperproc(nproc+1),ccpl_nsig))

        do n = 1, nanals
        do k = 1, ccpl_nsig   
            do i = 1, numptsperproc(nproc+1)
                ! Compute the dry hydrostatic pressure at the respective
                ! grid coordinate; This is dry pressure not full
                ! pressure, ignore this difference, since we are only
                ! using this to compute qsat, which in turn is only used
                ! to compute normalized humidity analysis variable
                enkf_pressure(i,k) = ccpl_znu(k,n)*(ccpl_mu(i,n)   &
                  & + ccpl_mub(i,n)) + ccpl_pt

                ! Compute mixing ratio from specific humidity.  
                enkf_spechumd(i,k) = (ccpl_q(i,k,n))/(1.0 +       &
                  & ccpl_q(i,k,n))

                ! Compute virtual temp (this is only used to compute
                ! saturation specific humidity (call genqsat1)
                enkf_virttemp(i,k) = ((ccpl_t(i,k,n) +        &
                  & 300.0)/((1000.0/(enkf_pressure(i,k)/100.0))    &
                  & **(rd/cp))) * (1. + fv*enkf_spechumd(i,k))
            end do ! 
        end do ! do k = 1, ccpl_nsig
        ! Compute the saturation specific humidity
        if (pseudo_rh) then
            ice = .false.
            call genqsat1(enkf_spechumd, ccpl_qsat_chunk(n,:,:), enkf_pressure/100.0, enkf_virttemp, ice, numptsperproc(nproc+1), ccpl_nsig)
        else
            ccpl_qsat_chunk(n,:,:) = 1._r_double
        endif
       
        end do

        if(allocated(enkf_virttemp))       deallocate(enkf_virttemp)
        if(allocated(enkf_pressure))       deallocate(enkf_pressure)
        if(allocated(enkf_spechumd))       deallocate(enkf_spechumd)
        !write(*,*) "============================= U ============================="
        ! ccpl_anal_chunk(i,:,1:ccpl_nsig) = ccpl_u(:,:,i)
        write(*,*) "ccpl_anal_chunk(1,1:10,1): ", ccpl_anal_chunk(1,1:10,1)
        write(*,*) "ccpl_u(1:10,1,1): ", ccpl_u(1:10,1,1)
        !write(*,*) "============================= T ============================="
        write(*,*) "ccpl_anal_chunk(3,1:10,141): ", ccpl_anal_chunk(3,1:10,141)
        write(*,*) "ccpl_t(1:10,141,3): ", ccpl_t(1:10,141,3)
    end subroutine enkf_fields_transfer_in

    subroutine enkf_fields_transfer_out

        use params, only: nlevs,nvars,ndim,nbackgrounds,nanals,pseudo_rh,cliptracers
        use constants, only: zero,one,cp,fv,rd,grav,zero
        implicit none
    
        integer :: i,j,k,n
        real    :: clip
        real(r_kind), allocatable ::  x_global_ccpl_mu(:,:)
        

        ! nvars=3 Updating U, V, T, and MU
        ! nvars=4 Updating U, V, T, QVAPOR, and MU
        ! nvars=5 Updating U, V, T, QVAPOR, PH, and MU
        ! nvars=6 Updating U, V, T, QVAPOR, W, PH, and MU
         if (nvars .eq. 3) then 
            do i = 1, nanals              
                ccpl_u(:,:,i) = ccpl_u(:,:,i) + ccpl_anal_chunk(i,:,1:ccpl_nsig) 
                ccpl_v(:,:,i) = ccpl_v(:,:,i) + ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) 
                ccpl_t(:,:,i) = ccpl_t(:,:,i) + ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig)
                ccpl_mu(:,i) = ccpl_mu(:,i) + ccpl_anal_chunk(i,:,3*ccpl_nsig+1)
            end do
        else if (nvars .eq. 4) then 
            do i = 1, nanals              
                ccpl_u(:,:,i) = ccpl_u(:,:,i) + ccpl_anal_chunk(i,:,1:ccpl_nsig) 
                ccpl_v(:,:,i) = ccpl_v(:,:,i) + ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) 
                ccpl_t(:,:,i) = ccpl_t(:,:,i) + ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig)
                ccpl_q(:,:,i) = ccpl_q(:,:,i) + ccpl_anal_chunk(i,:,3*ccpl_nsig+1:4*ccpl_nsig) 
                ccpl_mu(:,i) = ccpl_mu(:,i) + ccpl_anal_chunk(i,:,4*ccpl_nsig+1)
            end do
        else if (nvars .eq. 5) then 
            do i = 1, nanals              
                ccpl_u(:,:,i) = ccpl_u(:,:,i) + ccpl_anal_chunk(i,:,1:ccpl_nsig) 
                ccpl_v(:,:,i) = ccpl_v(:,:,i) + ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) 
                ccpl_t(:,:,i) = ccpl_t(:,:,i) + ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig)
                ccpl_q(:,:,i) = ccpl_q(:,:,i) + ccpl_anal_chunk(i,:,3*ccpl_nsig+1:4*ccpl_nsig) 
                ccpl_ph(:,:,i) = ccpl_ph(:,:,i) + ccpl_anal_chunk(i,:,4*ccpl_nsig+1:5*ccpl_nsig)
                ccpl_mu(:,i) = ccpl_mu(:,i) + ccpl_anal_chunk(i,:,5*ccpl_nsig+1)
            end do   
        else if (nvars .eq. 6) then 
            do i = 1, nanals              
                ccpl_u(:,:,i) = ccpl_u(:,:,i) + ccpl_anal_chunk(i,:,1:ccpl_nsig) 
                ccpl_v(:,:,i) = ccpl_v(:,:,i) + ccpl_anal_chunk(i,:,ccpl_nsig+1:2*ccpl_nsig) 
                ccpl_t(:,:,i) = ccpl_t(:,:,i) + ccpl_anal_chunk(i,:,2*ccpl_nsig+1:3*ccpl_nsig)
                ccpl_q(:,:,i) = ccpl_q(:,:,i) + ccpl_anal_chunk(i,:,3*ccpl_nsig+1:4*ccpl_nsig)
                ccpl_w(:,:,i) = ccpl_w(:,:,i) + ccpl_anal_chunk(i,:,4*ccpl_nsig+1:5*ccpl_nsig) 
                ccpl_ph(:,:,i) = ccpl_ph(:,:,i) + ccpl_anal_chunk(i,:,5*ccpl_nsig+1:6*ccpl_nsig)
                ccpl_mu(:,i) = ccpl_mu(:,i) + ccpl_anal_chunk(i,:,6*ccpl_nsig+1)
            end do 
        end if 
       ccpl_t=ccpl_t+300.
        
       allocate(x_global_ccpl_mu(ccpl_nlon*ccpl_nlat,ccpl_ens_num))
        do i = 1, numptsperproc(nproc+1)
            x_global_ccpl_mu(indxproc(nproc+1,i),:) = ccpl_mu(i,:)
        end do
        k = 1
        do j = 1, ccpl_nlat
            do i = 1, ccpl_nlon
               global_ccpl_mu(i,j,:) = x_global_ccpl_mu(k,:)
               k = k+1
            end do
        end do
        deallocate(x_global_ccpl_mu)
        
        ! Clip all tracers (assume names start with 'Q')

       if (cliptracers .and. nvars .ge. 4) then

          clip = tiny(ccpl_q(1,1,1))
          where (ccpl_q < clip) ccpl_q = clip

       end if 
       
        if(allocated(ccpl_anal_chunk))       deallocate(ccpl_anal_chunk)

    end subroutine enkf_fields_transfer_out

end module enkf_ccpl_coupling_mod


