module enkf_ccpl_grid_vars_mod

    use CCPL_interface_mod
    use enkf_ccpl_kinds, only: enkf_grid_instance_id, enkf_grid_comp_id, enkf_h2d_grid_id_input, enkf_h2d_md_grid_id_input, enkf_h2d_md_grid_id, enkf_ens_num_grid, enkf_vert_md_grid_id1, enkf_global_decomp_id
    use enkf_ccpl_kinds, only: ccpl_ens_num, ccpl_nlon, ccpl_nlat, ccpl_nsig, global_ccpl_xlat, global_ccpl_xlon, global_ccpl_mu, global_ccpl_mub, ccpl_znu, ccpl_pt, ccpl_nlon_ens, ccpl_nlat_ens, ccpl_nsig_ens, ccpl_pt_ens
    use enkf_ccpl_kinds, only: enkf_ccpl_comm_grid, enkf_npe_grid, enkf_mype_grid 
    implicit none


contains
   subroutine register_grids_decomps_from_ccpl

   		implicit none

   		integer :: i
   		real, allocatable :: ens_num_index(:)


        enkf_grid_comp_id = CCPL_external_procedures_get_comp_ID(enkf_grid_instance_id, annotation="get enkf comp id")
        enkf_h2d_md_grid_id_input = CCPL_external_procedures_para_get_grid_ID(enkf_grid_instance_id, 1)
        enkf_h2d_md_grid_id=enkf_h2d_md_grid_id_input
        enkf_vert_md_grid_id1=CCPL_external_procedures_para_get_grid_ID(enkf_grid_instance_id, 2)
        !enkf_vert_grid2=CCPL_external_procedures_para_get_grid_ID(enkf_grid_instance_id, 3)
        enkf_h2d_grid_id_input = CCPL_external_procedures_para_get_grid_ID(enkf_grid_instance_id, 7)
        ccpl_ens_num = CCPL_external_procedures_para_get_control_var(enkf_grid_instance_id, 1, annotation="get the number of ensemble members")
        enkf_ens_num_grid = CCPL_register_1D_normal_grid_without_data(enkf_grid_comp_id, "enkf_grid_ensemble_number_grid", ccpl_ens_num, annotation="register virtual ensemble number grid for enkf")
        !allocate(ens_num_index(ccpl_ens_num))
        !do i = 1, ccpl_ens_num
        	!ens_num_index(i) = i
        !end do
        !enkf_ens_num_grid = CCPL_register_V1D_Z_grid_via_model_data(enkf_grid_comp_id, "enkf_ensemble_number_grid", "num", ens_num_index, annotation="register virtual ensemble number grid for enkf")
        ccpl_nlat = CCPL_external_procedures_para_get_control_var(enkf_grid_instance_id, 2, annotation="get the number of latitudes")
        ccpl_nlon = CCPL_external_procedures_para_get_control_var(enkf_grid_instance_id, 3, annotation="get the number of longitudes")
        ccpl_nsig = CCPL_external_procedures_para_get_control_var(enkf_grid_instance_id, 4, annotation="get the number of levels")

        write(*,*) "+++++++++++++++++++++++++++"
        write(*,*) "ccpl_ens_num : ", ccpl_ens_num
        write(*,*) "nlat : ", ccpl_nlat
        write(*,*) "nlon : ", ccpl_nlon
        write(*,*) "nsig : ", ccpl_nsig


    end subroutine register_grids_decomps_from_ccpl

    subroutine enkf_ccpl_grid_vars_transfer

        implicit none
        integer :: field_id_nlat, field_id_nlon, field_id_nsig, field_id_xlat, field_id_xlon, field_id_pt, field_id_znu, field_id_mu, field_id_mub 
        integer :: dims_size(3), i, j, k, ierror
        integer, allocatable :: local_cells_global_index_all(:)
        real :: pt

        call mpi_comm_size(enkf_ccpl_comm_grid,enkf_npe_grid,ierror)
        call mpi_comm_rank(enkf_ccpl_comm_grid,enkf_mype_grid,ierror)
      
        !field_id_nlat=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_nlat, "nlat", CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of latitudes")
        !field_id_nlon=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_nlon, "nlon", CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of longitudes")
        !field_id_nsig=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_nsig, "nsig", CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of levels")
        !field_id_pt=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_pt, "P_TOP", CCPL_PARA_TYPE_IN, -1, -1, annotation="declare pressure top of the model")

        !field_id_nlat=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_nlat_ens, "nlat", CCPL_PARA_TYPE_IN, -1, enkf_ens_num_grid, (/ccpl_ens_num/), annotation="declare regional no. of latitudes")
        !field_id_nlon=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_nlon_ens, "nlon", CCPL_PARA_TYPE_IN, -1, enkf_ens_num_grid, (/ccpl_ens_num/), annotation="declare regional no. of longitudes")
        !field_id_nsig=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_nsig_ens, "nsig", CCPL_PARA_TYPE_IN, -1, enkf_ens_num_grid, (/ccpl_ens_num/), annotation="declare regional no. of levels")
        field_id_pt=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_pt_ens, "P_TOP", CCPL_PARA_TYPE_IN, -1, enkf_ens_num_grid, (/ccpl_ens_num/), annotation="declare pressure top of the model")

        !ccpl_nlat => ccpl_nlat_ens(1)
        !ccpl_nlon => ccpl_nlon_ens(1)
        !ccpl_nsig => ccpl_nsig_ens(1)
        ccpl_pt => ccpl_pt_ens(1)
        !ccpl_pt = 100.

        field_id_znu=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, ccpl_znu, "ZNU_1", CCPL_PARA_TYPE_IN, -1, enkf_vert_md_grid_id1, (/ccpl_nsig, ccpl_ens_num/), annotation="declare eta values on half (mass) levels")
        allocate(local_cells_global_index_all(1:ccpl_nlat*ccpl_nlon))
        local_cells_global_index_all=CCPL_NULL_INT
        if (enkf_mype_grid==0) then
        k=0
        do j=1, ccpl_nlat
            do i=1, ccpl_nlon
               k=k+1
               local_cells_global_index_all(k)=(j-1)*ccpl_nlon+i
            end do 
        end do
        end if
        enkf_global_decomp_id = CCPL_register_normal_parallel_decomp("enkf_global_decomp", enkf_h2d_grid_id_input, ccpl_nlat*ccpl_nlon, local_cells_global_index_all, annotation = "register grid decomp for enkf global grid vars")
        dims_size(1)=ccpl_nlon
        dims_size(2)=ccpl_nlat
        dims_size(3)=ccpl_ens_num
        !allocate(global_ccpl_xlat(ccpl_nlon,ccpl_nlat,ccpl_ens_num),global_ccpl_xlon(ccpl_nlon,ccpl_nlat,ccpl_ens_num),global_ccpl_mu(ccpl_nlon,ccpl_nlat,ccpl_ens_num),global_ccpl_mub(ccpl_nlon,ccpl_nlat,ccpl_ens_num))
        field_id_xlat=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, global_ccpl_xlat, "XLAT", CCPL_PARA_TYPE_IN, enkf_global_decomp_id, enkf_h2d_md_grid_id, dims_size, annotation="declare latitudes")
        field_id_xlon=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, global_ccpl_xlon, "XLONG", CCPL_PARA_TYPE_IN, enkf_global_decomp_id, enkf_h2d_md_grid_id, dims_size, annotation="declare longitudes")
        field_id_mu=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, global_ccpl_mu, "mu", CCPL_PARA_TYPE_IN, enkf_global_decomp_id, enkf_h2d_md_grid_id, dims_size, annotation="declare mu")
        field_id_mub=CCPL_external_procedures_para_declare_field(enkf_grid_instance_id, global_ccpl_mub, "mub", CCPL_PARA_TYPE_IN, enkf_global_decomp_id, enkf_h2d_md_grid_id, dims_size, annotation="declare mub")

     end subroutine enkf_ccpl_grid_vars_transfer

end module enkf_ccpl_grid_vars_mod
