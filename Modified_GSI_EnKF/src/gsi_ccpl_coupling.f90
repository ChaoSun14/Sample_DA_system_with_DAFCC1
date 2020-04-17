! sunchao 2019-03-21 
! 
! abstract: This module is 

module gsi_ccpl_coupling_mod
    
    use CCPL_interface_mod
    use gsi_ccpl_kinds, only: gsi_instance_id, gsi_grid_id_input, gsi_comp_id, gsi_decomp_id, gsi_grid_id, gsi_vert_grid1, gsi_vert_grid2, gsi_3d_grid1, gsi_3d_grid2, gsi_regional_time_grid_id, gsi_ccpl_comm, gsi_ccpl_comm_all, color, gsi_grid_ccpl_comm
    use gsi_ccpl_kinds, only: ccpl_nlat_regional, ccpl_nlon_regional, ccpl_nsig, ccpl_glat, ccpl_glon, ccpl_aeta1, ccpl_eta1, ccpl_dx_mc, ccpl_dy_mc
    use gsi_ccpl_kinds, only: ccpl_psfc, ccpl_snow, ccpl_fis, ccpl_landmask, ccpl_xice, ccpl_sst, ccpl_tsk, ccpl_tslb, ccpl_smois, ccpl_ivgtyp, ccpl_vegfrac, ccpl_isltyp, ccpl_u, ccpl_v, ccpl_q, ccpl_pot, ccpl_regional_time
    use gsi_ccpl_kinds, only: local_cells_global_index, local_cells_global_index_with_buffer, gsi_ccpl_kinds_cleanup
    use gsi_ccpl_kinds, only: ccpl_npe_all, ccpl_mype_all, ccpl_npe_wrf, ccpl_mem_id
    use kinds, only: r_kind, r_single    
    !use gridmod, only: lat1, lon1, lat2, lon2, ijn, istart, jstart, nlat, nlon, region_lat, region_lon, nsig
    !use mpimod, only: npe, mype
    implicit none
    integer, private :: lon1,lon2,lat1,lat2,nlat,nlon,npe,mype,nsig,ierror,split_comm
    integer, allocatable, private :: ilat1(:),istart(:),jlon1(:),jstart(:)
contains

    subroutine register_grids_decomps

        use mpi
        use constants, only: pi
        use general_sub2grid_mod, only: general_deter_subdomain
        use general_commvars_mod, only: ltosi, ltosj, ltosi_s, ltosj_s
        implicit none
         
        integer              :: i,j,k,ii,jj,kk
        integer              :: start_index_mype, lat_start_index, lon_start_index
        integer, allocatable :: mask(:)
        real(r_single)         :: max_lon, min_lon, max_lat, min_lat, regional_time_coord_values(6)
        real(r_single), allocatable :: mype_lat(:), mype_lon(:)
        logical              :: buff_flag1,buff_flag2
        
        logical :: periodic
        logical, allocatable :: periodic_s(:)
        
        call mpi_bcast(ccpl_glat, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        call mpi_bcast(ccpl_glon, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        !call mpi_bcast(ccpl_dx_mc, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        !call mpi_bcast(ccpl_dy_mc, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        !do j=1, ccpl_nlat_regional
        !    do i=1, ccpl_nlon_regional
        !        write(*,*) "[CCPL <GSI>] ", ccpl_glat(i,j), ccpl_glon(i,j), ccpl_dx_mc(i,j), ccpl_dy_mc(i,j), i, j
        !    end do
        !end do
        !call mpi_barrier(gsi_grid_ccpl_comm, ierror)
        write(*,*) "[CCPL <GSI>] ", ccpl_glat(1,1:10), ccpl_glon(1,1:10)
        ccpl_npe_wrf = CCPL_external_procedures_para_get_control_var(gsi_instance_id, 6, annotation="get the number process of model")
        call mpi_comm_size(gsi_ccpl_comm_all,ccpl_npe_all,ierror)
        call mpi_comm_rank(gsi_ccpl_comm_all,ccpl_mype_all,ierror)
        if (ccpl_npe_wrf.ne.ccpl_npe_all) then
            color = ccpl_mype_all / ccpl_npe_wrf
            call mpi_comm_split(gsi_ccpl_comm_all, color, ccpl_mype_all, split_comm, ierror)
            gsi_ccpl_comm = split_comm
            ccpl_mem_id = 0
        else
            color = -1
            gsi_ccpl_comm = gsi_ccpl_comm_all
            ccpl_mem_id = CCPL_external_procedures_para_get_control_var(gsi_instance_id, 1, annotation="get the ensemble member id of current process")    
        end if
        call mpi_comm_size(gsi_ccpl_comm,npe,ierror)
        call mpi_comm_rank(gsi_ccpl_comm,mype,ierror)
    
        write(*,*) "[CCPL <GSI>]  gsi_ccpl_comm_all = ", gsi_ccpl_comm_all
        write(*,*) "[CCPL <GSI>]  ccpl_npe_wrf, ccpl_mem_id = ", ccpl_npe_wrf, ccpl_mem_id
        nlat=ccpl_nlat_regional
        nlon=ccpl_nlon_regional
        nsig=ccpl_nsig
        if (allocated(ilat1))  deallocate(ilat1)
        if (allocated(istart))  deallocate(istart)
        if (allocated(jlon1))  deallocate(jlon1)
        if (allocated(jstart))  deallocate(jstart)
        if (allocated(periodic_s))  deallocate(periodic_s)
        allocate(ilat1(npe),istart(npe),jlon1(npe),jstart(npe),periodic_s(npe))
        call general_deter_subdomain(npe,mype,ccpl_nlat_regional,ccpl_nlon_regional,.true., &
                                periodic,periodic_s,lon1,lon2,lat1,lat2,ilat1,istart,jlon1,jstart)

        gsi_comp_id = CCPL_external_procedures_get_comp_ID(gsi_instance_id, annotation="get gsi comp id") 
        !gsi_grid_id_input   = CCPL_external_procedures_para_get_grid_ID(gsi_instance_id,1)
        
        allocate(local_cells_global_index(lat1*lon1), mask(lat1*lon1))
        allocate(mype_lat(lat1*lon1),mype_lon(lat1*lon1))
        mask = 1
        k = 0
        do j=1, lon1
            do i=1, lat1
                k = k+1
                local_cells_global_index(k)=(jstart(mype+1)+j-2)*nlat+istart(mype+1)+i-1
                mype_lon(k)=ccpl_glon(jstart(mype+1)+j-1, istart(mype+1)+i-1)
                mype_lat(k)=ccpl_glat(jstart(mype+1)+j-1, istart(mype+1)+i-1)
                !mype_lon(k)=region_lon(istart(mype+1)+i-1, jstart(mype+1)+j-1)
                !mype_lat(k)=region_lat(istart(mype+1)+i-1, jstart(mype+1)+j-1)
            end do
        end do
        allocate(local_cells_global_index_with_buffer(lat2*lon2))
        k=0
        do j=1, lon2
            buff_flag1=.false.
            if(jstart(mype+1).eq.1) then
                lon_start_index=jstart(mype+1)
                if(j.eq.1) then
                    buff_flag1=.true.
                else
                    jj=j-1
                end if
            else
                lon_start_index=jstart(mype+1)-1
                jj=j
            end if
            if((lon_start_index+jj-2).eq.nlon) buff_flag1=.true.
            do i=1, lat2
                buff_flag2=.false.
                if(istart(mype+1).eq.1) then
                    lat_start_index=istart(mype+1)
                    if(i.eq.1) then
                        buff_flag2=.true.
                    else
                        ii=i-1
                    end if
                else
                    lat_start_index=istart(mype+1)-1
                    ii=i
                end if
                if((lat_start_index+ii-2).eq.nlat) buff_flag2=.true.
                k=k+1
                if(buff_flag1.or.buff_flag2) then
                    local_cells_global_index_with_buffer(k)=CCPL_NULL_INT
                else
                    local_cells_global_index_with_buffer(k)=(lon_start_index+jj-2)*nlat+lat_start_index+ii-1
                end if
            end do
        end do
        if (color.ne.0.and.color.ne.-1) then
            local_cells_global_index_with_buffer=CCPL_NULL_INT
        end if

        !write(*,*) "local_cells_global_index_with_buffer: ", local_cells_global_index_with_buffer
        
        max_lon = -999999.!2*PI
        min_lon = -999999.!0
        max_lat = -999999.!PI/2.0
        min_lat = -999999.!-PI/2.0
        gsi_grid_id = CCPL_register_H2D_grid_via_local_data(gsi_comp_id,"gsi_H2D_grid", "LON_LAT",  "degrees","acyclic", nlat*nlon, lat1*lon1, local_cells_global_index, min_lon, max_lon, min_lat, max_lat, mype_lon, mype_lat, mask, annotation="register H2D grid for gsi" )
        !gsi_grid_id = CCPL_external_procedures_para_get_grid_ID(gsi_instance_id,1)
        gsi_decomp_id = CCPL_register_normal_parallel_decomp("gsi_decomp",gsi_grid_id, lat2*lon2, local_cells_global_index_with_buffer, annotation = "register parallel decomp (with buffer) for gsi")
        regional_time_coord_values=(/1.,2.,3.,4.,5.,6./)
        gsi_regional_time_grid_id=CCPL_register_V1D_Z_grid_via_model_data(gsi_comp_id, "gsi_regional_time_grid", "virtual_gird", regional_time_coord_values, annotation="register simplified regional_time_grid for gsi")
        gsi_vert_grid1 = CCPL_register_V1D_Z_grid_via_model_data(gsi_comp_id, "gsi_vertical_grid_nsig", "Pa", ccpl_aeta1(1:nsig), annotation="register simplified vertical grid 1 for gsi")
        gsi_vert_grid2 = CCPL_register_V1D_Z_grid_via_model_data(gsi_comp_id, "gsi_vertical_grid_nsig_1", "Pa", ccpl_eta1(1:nsig), annotation="register simplified vertical grid 2 for gsi")
        gsi_3d_grid1=CCPL_register_MD_grid_via_multi_grids(gsi_comp_id,"gsi_3d_grid_1",gsi_grid_id,gsi_vert_grid1,annotation="register gsi 3D-grid with nsig vertical levels")
        gsi_3d_grid2=CCPL_register_MD_grid_via_multi_grids(gsi_comp_id,"gsi_3d_grid_2",gsi_grid_id,gsi_vert_grid2,annotation="register gsi 3D-grid with nsig+1 vertical levels")
        deallocate(local_cells_global_index, mask, local_cells_global_index_with_buffer)
        deallocate(mype_lat,mype_lon)
   end subroutine register_grids_decomps

   subroutine gsi_fields_declare
        use constants, only: zero, one, three, deg2rad,half, two, r0_01
        implicit none
        integer :: field_id_regional_time,field_id_nlat,field_id_nlon,field_id_nsig,field_id_nsig_soil,field_id_pt,field_id_aeta1,field_id_aeta2,field_id_eta1,field_id_eta2
        integer :: field_id_dx_mc,field_id_dy_mc,field_id_glat,field_id_glon,field_id_psfc,field_id_fis
        integer :: field_id_landmask,field_id_xice,field_id_sst,field_id_ivgtyp,field_id_isltyp,field_id_vegfrac,field_id_snow,field_id_smois,field_id_tslb,field_id_tsk
        integer :: field_id_u,field_id_v,field_id_pot,field_id_q
        integer :: dims_size_2d(2),dims_size_3d(3)

        dims_size_2d(1)=lat2
        dims_size_2d(2)=lon2
        dims_size_3d(1)=lat2
        dims_size_3d(2)=lon2
        dims_size_3d(3)=nsig
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        write(*,*) "[CCPL <GSI>] Start declare fields"
        write(*,*) "[CCPL <GSI>] dims_size_2d: ", dims_size_2d
        write(*,*) "[CCPL <GSI>] dims_size_3d: ", dims_size_3d

        field_id_regional_time=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_regional_time,"regional_time",CCPL_PARA_TYPE_IN, -1, gsi_regional_time_grid_id, (/6/),annotation="declare regional time")
        field_id_smois=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_smois,"smois",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare soil moistures (k=1)")
        field_id_tslb=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_tslb,"tslb",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare soil temperature (k=1)")
        field_id_psfc=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_psfc,"psfc0",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare surface pressure")
         field_id_fis=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_fis,"PHB",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare base-state geopotential")
         field_id_snow=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_snow,"snow",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare accumulated melted snow")

         field_id_pot=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_pot,"ccpl_t",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_3d_grid1, dims_size_3d, annotation="declare perturbation potential temperature theta-t0") 
         field_id_u=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_u,"ccpl_u",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_3d_grid1, dims_size_3d, annotation="register x-wind component at mass point")
         field_id_v=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_v,"ccpl_v",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_3d_grid1, dims_size_3d, annotation="register y-wind component at mass point")
         field_id_q=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_q,"ccpl_q",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_3d_grid1, dims_size_3d, annotation="register water vapor mixing ratio")
         field_id_landmask=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_landmask,"landmask",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare land mask (1 for land, 0 for water)")
         field_id_xice=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_xice,"seaice",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare sea ice flag")
         field_id_sst=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_sst,"ccpl_sst",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare skin sea surface temperature")
         field_id_tsk=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_tsk,"tsk",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare surface skin temperature")
         field_id_ivgtyp=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_ivgtyp,"ivgtyp",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare domain vegetation category")
         field_id_isltyp=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_isltyp,"isltyp",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare domain soil category")
         field_id_vegfrac=CCPL_external_procedures_para_declare_field(gsi_instance_id,ccpl_vegfrac,"vegfra",CCPL_PARA_TYPE_IN, gsi_decomp_id, gsi_grid_id, dims_size_2d, annotation="declare vegetation fraction")

        write(*,*) "[CCPL <GSI>] Finish declare fields"
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"
   
   end subroutine gsi_fields_declare
   
   subroutine boundary_fill_all_fields
        
        use mpi
        implicit none

        real(r_single),parameter:: one_single = 1.0_r_single
        real(r_single),parameter:: r45 = 45.0_r_single
        real(r_single) rad2deg_single
        call mpi_bcast(ccpl_glat, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        call mpi_bcast(ccpl_glon, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        call mpi_bcast(ccpl_dx_mc, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        call mpi_bcast(ccpl_dy_mc, ccpl_nlat_regional*ccpl_nlon_regional, MPI_REAL4, 0, gsi_grid_ccpl_comm, ierror)
        rad2deg_single=r45/atan(one_single) 
        ccpl_glat = ccpl_glat/rad2deg_single
        ccpl_glon = ccpl_glon/rad2deg_single
        write(*,*) "[CCPL <GSI>] ", ccpl_glat(1,1:10), ccpl_glon(1,1:10)
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        write(*,*) "[CCPL <GSI>] Start boundary_fill_fields"

        call boundary_fill_3D(ccpl_u)
        call boundary_fill_3D(ccpl_v)
        call boundary_fill_3D(ccpl_q)
        call boundary_fill_3D(ccpl_pot)
        call boundary_fill_2D(ccpl_psfc) 
        call boundary_fill_2D(ccpl_snow)
        call boundary_fill_2D(ccpl_fis)
        call boundary_fill_2D(ccpl_landmask)
        call boundary_fill_2D(ccpl_xice)
        call boundary_fill_2D(ccpl_sst)
        call boundary_fill_2D(ccpl_tsk)
        call boundary_fill_2D(ccpl_tslb)
        call boundary_fill_2D(ccpl_smois)
        call boundary_fill_2D(ccpl_ivgtyp)
        call boundary_fill_2D(ccpl_vegfrac)
        call boundary_fill_2D(ccpl_isltyp)

        write(*,*) "[CCPL <GSI>] Finish boundary_fill_fields"
        write(*,*) "++++++++++++++++++++++++++++++++++++++++"
        
   end subroutine boundary_fill_all_fields

   subroutine boundary_fill_2D(field_inout)
        
        implicit none

        real(r_kind), pointer, INTENT(INOUT) :: field_inout(:,:)
        real(r_kind), allocatable            :: field_out(:,:)
        integer              :: i,j,k,ii,jj,kk
        integer              :: lat_start_index, lon_start_index
        logical              :: buff_flag1,buff_flag2
        
        allocate(field_out(lat2,lon2))
        do j=1, lon2
            if(jstart(mype+1).eq.1) then
                lon_start_index=jstart(mype+1)
                if(j.eq.1) then
                    jj=j+1
                else
                    jj=j
                end if
            else
                lon_start_index=jstart(mype+1)-1
                jj=j
            end if
            if((lon_start_index+jj-2).eq.nlon) then
                jj=j-1
            end if
            do i=1, lat2
                if(istart(mype+1).eq.1) then
                    lat_start_index=istart(mype+1)
                    if(i.eq.1) then
                        ii=i+1
                    else
                        ii=i
                    end if
                else
                    lat_start_index=istart(mype+1)-1
                    ii=i
                end if
                if((lat_start_index+ii-2).eq.nlat) then
                    ii=ii-1
                end if
                k=k+1
                field_out(i,j)=field_inout(ii,jj)
                !write(*,*) "lon_start_index ", lon_start_index+jj-1
            end do
        end do
        field_inout=field_out
        deallocate(field_out)
   end subroutine boundary_fill_2D

   subroutine boundary_fill_3D(field_inout)
        
        implicit none

        real(r_kind), pointer, INTENT(INOUT) :: field_inout(:,:,:)
        real(r_kind), allocatable            :: field_out(:,:,:)
        integer              :: i,j,k,ii,jj,kk
        integer              :: lat_start_index, lon_start_index
        logical              :: buff_flag1,buff_flag2

        allocate(field_out(lat2,lon2,nsig))
        write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write(*,*) "ni,nj,nk= ",size(field_out,1),size(field_out,2),size(field_out,3)
        do k=1, nsig
        do j=1, lon2
            if(jstart(mype+1).eq.1) then
                lon_start_index=jstart(mype+1)
                if(j.eq.1) then
                    jj=j+1
                else
                    jj=j
                end if
            else
                lon_start_index=jstart(mype+1)-1
                jj=j
            end if
            if((jstart(mype+1)+jj-2).gt.nlon) then
                jj=j-1
            end if
            do i=1, lat2
                if(istart(mype+1).eq.1) then
                    lat_start_index=istart(mype+1)
                    if(i.eq.1) then
                        ii=i+1
                    else
                        ii=i
                    end if
                else
                    lat_start_index=istart(mype+1)-1
                    ii=i
                end if
                if((istart(mype+1)+ii-2).gt.nlat) then
                    ii=ii-1
                end if
                field_out(i,j,k)=field_inout(ii,jj,k)
                !write(*,*) "lon_start_index ", lon_start_index+jj-1
            end do
        end do
        end do
        field_inout=field_out
        deallocate(field_out)
   end subroutine boundary_fill_3D
        
end module gsi_ccpl_coupling_mod


