module gsi_ccpl_grid_vars_mod

    use CCPL_interface_mod
    use gsi_ccpl_kinds, only: gsi_grid_instance_id, gsi_grid_id_input, gsi_comp_id, gsi_grid_id, gsi_vert_grid1, gsi_vert_grid2, gsi_3d_grid1, gsi_3d_grid2, gsi_regional_time_grid_id, gsi_global_decomp_id, gsi_grid_ccpl_comm
    use gsi_ccpl_kinds, only: ccpl_nlon_regional, ccpl_nlat_regional, ccpl_nsig, ccpl_nsig_soil, ccpl_pt, ccpl_aeta1, ccpl_eta1, ccpl_aeta2, ccpl_eta2, ccpl_glat, ccpl_glon, ccpl_dx_mc, ccpl_dy_mc, ccpl_npe_wrf  
    use kinds, only: r_kind
    implicit none

contains
   subroutine register_grids_decomps_from_wrf
        use kinds, only: r_kind
        gsi_comp_id = CCPL_external_procedures_get_comp_ID(gsi_grid_instance_id, annotation="get gsi comp id")
        gsi_grid_id_input = CCPL_external_procedures_para_get_grid_ID(gsi_grid_instance_id,1)
        gsi_grid_id=gsi_grid_id_input
        gsi_vert_grid1=CCPL_external_procedures_para_get_grid_ID(gsi_grid_instance_id,2)
        gsi_vert_grid2=CCPL_external_procedures_para_get_grid_ID(gsi_grid_instance_id,3)
        write(*,*) "finish register grids from wrf"
    end subroutine register_grids_decomps_from_wrf

    subroutine gsi_ccpl_grid_vars_transfer

        use kinds, only: i_kind,r_kind, r_single
        use constants, only: zero, one, three, deg2rad,half, two, r0_01
        implicit none
        integer :: field_id_nlat,field_id_nlon,field_id_nsig,field_id_nsig_soil,field_id_pt,field_id_aeta1,field_id_aeta2,field_id_eta1,field_id_eta2
        integer :: field_id_dx_mc,field_id_dy_mc,field_id_glat,field_id_glon
        integer :: dims_size(2),i,j,k,mype_grid,ierror
        integer, allocatable :: local_cells_global_index_all(:)
        !field_id_nlat=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_nlat_regional,"nlat",CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of latitudes")
        !field_id_nlon=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_nlon_regional,"nlon",CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of longitudes")
        !field_id_nsig=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_nsig,"nsig",CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of levels")
        !field_id_nsig_soil=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_nsig_soil,"nsig_soil",CCPL_PARA_TYPE_IN, -1, -1, annotation="declare regional no. of soil levels")
        ccpl_nlat_regional = CCPL_external_procedures_para_get_control_var(gsi_grid_instance_id, 2, annotation="get the number of latitudes")
        ccpl_nlon_regional = CCPL_external_procedures_para_get_control_var(gsi_grid_instance_id, 3, annotation="get the number of longitudes")
        ccpl_nsig          = CCPL_external_procedures_para_get_control_var(gsi_grid_instance_id, 4, annotation="get the number of levels")
        ccpl_nsig_soil     = CCPL_external_procedures_para_get_control_var(gsi_grid_instance_id, 5, annotation="get the number of soil levels")
        field_id_pt=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_pt,"P_TOP",CCPL_PARA_TYPE_IN, -1, -1, annotation="declare pressure top of the model")
        write(*,*) "+++++++++++++++++++++++++++"
        write(*,*) "nlat,nlon,nsig,nsig_soil ",ccpl_nlat_regional,ccpl_nlon_regional,ccpl_nsig,ccpl_nsig_soil
        write(*,*) "ccpl_nsig, pt", ccpl_nsig, ccpl_pt
        field_id_aeta1=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_aeta1,"ZNU_1",CCPL_PARA_TYPE_IN, -1, gsi_vert_grid1, (/ccpl_nsig/), annotation="declare eta values on half (mass) levels")
        field_id_aeta2=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_aeta2,"ZNU_2",CCPL_PARA_TYPE_IN, -1, gsi_vert_grid1, (/ccpl_nsig/), annotation="declare eta values on half (mass) levels")
        field_id_eta1=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_eta1,"ZNW_1",CCPL_PARA_TYPE_IN, -1, gsi_vert_grid2, (/ccpl_nsig+1/), annotation="declare eta values on full (w) levels")
        field_id_eta2=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_eta2,"ZNW_2",CCPL_PARA_TYPE_IN, -1, gsi_vert_grid2, (/ccpl_nsig+1/), annotation="declare eta values on full (w) levels")
        ccpl_npe_wrf = CCPL_external_procedures_para_get_control_var(gsi_grid_instance_id, 6, annotation="get the number process of model")
        call mpi_comm_rank(gsi_grid_ccpl_comm,mype_grid,ierror)
        allocate(local_cells_global_index_all(1:ccpl_nlat_regional*ccpl_nlon_regional))
        local_cells_global_index_all=CCPL_NULL_INT
        if(mype_grid<1) then
        k=0
        do j=1, ccpl_nlat_regional
            do i=1, ccpl_nlon_regional
               k=k+1
               local_cells_global_index_all(k)=(j-1)*ccpl_nlon_regional+i
            end do 
        end do
        end if
        gsi_global_decomp_id = CCPL_register_normal_parallel_decomp("gsi_global_decomp",gsi_grid_id, ccpl_nlat_regional*ccpl_nlon_regional, local_cells_global_index_all, annotation = "register grid decomp for gsi global grid vars")
        dims_size(1)=ccpl_nlon_regional
        dims_size(2)=ccpl_nlat_regional
        field_id_glat=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_glat,"XLAT",CCPL_PARA_TYPE_IN, gsi_global_decomp_id, gsi_grid_id, dims_size, annotation="declare latitudes")
        field_id_glon=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_glon,"XLONG",CCPL_PARA_TYPE_IN, gsi_global_decomp_id, gsi_grid_id, dims_size, annotation="declare longitudes")
        field_id_dx_mc=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_dx_mc,"MAPFAC_MX",CCPL_PARA_TYPE_IN, gsi_global_decomp_id, gsi_grid_id, dims_size, annotation="declare x direction map scale factor")
        field_id_dy_mc=CCPL_external_procedures_para_declare_field(gsi_grid_instance_id,ccpl_dy_mc,"MAPFAC_MY",CCPL_PARA_TYPE_IN, gsi_global_decomp_id, gsi_grid_id, dims_size, annotation="declare y direction map scale factor")

        deallocate(local_cells_global_index_all)

     end subroutine gsi_ccpl_grid_vars_transfer

end module gsi_ccpl_grid_vars_mod














