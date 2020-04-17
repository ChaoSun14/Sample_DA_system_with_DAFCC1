    module module_wrf_DA
        USE ccpl_interface_mod
        USE module_timing
        integer :: root_comm,root_comp_id,my_comm,my_comp_id,proc_id
        integer :: ensemble_nums, ensemble_mem_id, gsi_grid_instance_id, gsi_ens_member_inst_id, gsi_ens_ave_inst_id, enkf_inst_id
        integer :: export_interface_id,import_interface_id
        integer, pointer :: regional_time(:)=>null(),nlat=>null(),nlon=>null(),nsig=>null(),nsig_soil=>null()
        integer :: control_vars(5), comp_or_grid_ids(6), decomp_ids(6), field_inst_ids(32), da_timers_id(1)
        real, pointer :: ccpl_pt=>null()
        real, dimension(:),  pointer    :: aeta1=>null(),aeta2=>null(),eta1=>null(),eta2=>null()
        real, dimension(:,:), pointer   :: ccpl_dx_mc=>null(),ccpl_dy_mc=>null(),ccpl_xlat=>null(),ccpl_xlong=>null(),ccpl_psfc0=>null(),ccpl_phb=>null(),ccpl_landmask=>null(),ccpl_seaice=>null(),ccpl_sst=>null(),ccpl_ivgtyp=>null(),ccpl_isltyp=>null(),ccpl_vegfra=>null(),ccpl_snow=>null(),ccpl_u10=>null(),ccpl_v10=>null(),ccpl_smois=>null(),ccpl_tslb=>null(),ccpl_tsk=>null(),ccpl_mu=>null(),ccpl_mub=>null()
        real, dimension(:,:,:), pointer :: ccpl_u=>null(),ccpl_v=>null(),ccpl_t=>null(),ccpl_q=>null(),ccpl_w=>null(),ccpl_ph=>null()
        
        integer :: grid_H2D_id,grid_vert_id1,grid_vert_id2,grid_regional_time_id,grid_3d_id1,grid_3d_id2,decomp_id
        integer :: field_id_regional_time,field_id_nlat,field_id_nlon,field_id_nsig,field_id_nsig_soil,field_id_pt,field_id_aeta1,field_id_aeta2,field_id_eta1,field_id_eta2
        integer :: field_id_dx_mc,field_id_dy_mc,field_id_xlat,field_id_xlong,field_id_psfc0,field_id_phb,field_id_mu,field_id_mub
        integer :: field_id_landmask,field_id_seaice,field_id_sst,field_id_ivgtyp,field_id_isltyp,field_id_vegfra,field_id_snow,field_id_u10,field_id_v10,field_id_smois,field_id_tslb,field_id_tsk
        integer :: field_id_u,field_id_v,field_id_t,field_id_q,field_id_w, field_id_ph
        integer, dimension(:), allocatable  :: fields_id       
        integer :: myrank,nprocs
        integer, private :: time_step,timer_id,ims,ime,jms,jme
        character*14  :: string_time

    contains 
            
     SUBROUTINE register_component(comm)
        IMPLICIT NONE
        integer, intent(inout) :: comm
        character(256)         :: ensemble_mem_id_str

        call CCPL_get_ensemble_info(ensemble_nums,ensemble_mem_id)
        write(ensemble_mem_id_str,*) ensemble_mem_id
        comm = CCPL_NULL_COMM
        root_comp_id = CCPL_register_component(-1, "wrf","atm", comm, change_dir=.false., annotation = "register root wrf to C-Coupler")
        root_comm  = comm
        comm = CCPL_NULL_COMM
        my_comp_id = CCPL_register_component(root_comp_id, "wrf_member"//trim(adjustl(ensemble_mem_id_str)),"atm", comm, change_dir=.false., annotation = "register wrf ensemble member "//trim(adjustl(ensemble_mem_id_str))//" to C-Coupler")
        call CCPL_change_to_ensemble_dir(root_comp_id,ensemble_mem_id)
        my_comm = comm
     END SUBROUTINE register_component
     
    SUBROUTINE initialize_DA_instances(grid)
        USE ESMF_mod, ONLY : ESMF_TimeIntervalGet, ESMF_TimeGet
        USE module_domain, ONLY : domain, domain_get_time_step, domain_clock_get
        USE module_utility, ONLY : WRFU_Time
        use module_state_description, ONLY: P_QV
        USE ccpl_interface_mod
        IMPLICIT NONE
INCLUDE "mpif.h"
        ! Input data 
        TYPE(domain) , POINTER             :: grid
        ! Local data
        TYPE(WRFU_Time) :: currentTime
        integer :: dm_is,dm_ie,dm_js,dm_je
        integer :: i,j,k,ip,jc,ierr,Isize,Jsize,num_local_cells,num_vars,num_comps
        integer :: iyear,imonth,iday,ihour,iminute,isecond
        integer :: MyError,num_sum
        real    :: min_lat,max_lat,min_lon,max_lon,rad2deg,regional_time_coord_values(6)
        character*1024  :: comps_full_names(12)
        integer, dimension(:), allocatable  :: start,length,num_array,mask,local_cells_global_index,local_cells_global_index_stagger
        real, dimension(:), allocatable     :: ylong,ylat
        
        proc_id = CCPL_get_current_process_id_in_component(my_comp_id)
        CALL mpi_comm_rank (my_comm, myrank, MyError)
        CALL mpi_comm_size (my_comm, nprocs, MyError)
!Find the real part of grid data                
        dm_is = grid%sp31
        dm_ie = grid%ep31
        dm_js = grid%sp33
        dm_je = grid%ep33

        IF (grid%ed31.eq.dm_ie) THEN
            dm_ie=dm_ie-1
        END IF
        IF (grid%ed33.eq.dm_je) THEN
            dm_je=dm_je-1
        END IF
        Isize=dm_ie-dm_is+1
        Jsize=dm_je-dm_js+1
        allocate(start(Jsize),length(Jsize))
        jc=0
        DO j=dm_js,dm_je
           jc=jc+1
           start(jc)=(j-1)*(grid%ed31-1)+dm_is
           length(jc)=Isize
        END DO
        num_sum = 0
        num_local_cells = Isize*Jsize
        allocate(local_cells_global_index(0:num_local_cells-1))
        allocate(num_array(nprocs))
        CALL MPI_ALLGATHER(num_local_cells,1,MPI_INTEGER,num_array,1,MPI_INTEGER,my_comm,ierr)
        k=0
        j=1
        do i=0,num_local_cells-1
            local_cells_global_index(i)= start(j)+k
            k=k+1
            if(k.eq.length(j)) then
                k=0
                j=j+1
            endif
        enddo
!Choose the exact part of xlong and xlat,change these to ylong and ylat
        allocate(ylong(Isize*Jsize),ylat(Isize*Jsize),mask(Isize*Jsize))
        k=1
        do j=dm_js,dm_je
            do i=dm_is,dm_ie
                ylong(k) = grid%xlong(i,j)
                if(ylong(k)<0.0) then
                    ylong(k)=360.0+ylong(k)
                endif
                ylat(k) = grid%xlat(i,j)
                mask(k) = 1
                k = k + 1
            enddo
        enddo

        min_lat=-999999.!-90.0
        max_lat=-999999.!90.0
        min_lon=-999999.!0.0
        max_lon=-999999.!360.0
        
        grid_H2D_id = CCPL_register_H2D_grid_via_local_data(my_comp_id,"wrf_H2D_grid", "LON_LAT", "degrees","acyclic",sum(num_array),num_local_cells,local_cells_global_index,min_lon,max_lon, min_lat, max_lat, ylong, ylat,mask,annotation="register H2D grid for wrf")
! register decomposition  for the staggered domains     
        dm_is = grid%sm31
        dm_ie = grid%em31
        dm_js = grid%sm33
        dm_je = grid%em33
        ims=dm_is
        ime=dm_ie
        jms=dm_js
        jme=dm_je
        Isize=dm_ie-dm_is+1
        Jsize=dm_je-dm_js+1
        num_local_cells=Isize*Jsize
        allocate(local_cells_global_index_stagger(num_local_cells))
        local_cells_global_index_stagger = CCPL_NULL_INT
        k=1
        jc=0
        do j=dm_js,dm_je
            do i=dm_is,dm_ie
               if (j.ge.grid%sp33.and.j.le.grid%ep33.and.&
                   i.ge.grid%sp31.and.i.le.grid%ep31) then
                  if (j.ne.grid%ed33.and.i.ne.grid%ed31) then
                     local_cells_global_index_stagger(k) = (j-1)*(grid%ed31-1)+i
                     jc=jc+1
                     if (local_cells_global_index_stagger(k)==122) write(*,*) "special location ", i, j, grid%ed31, grid%sp31, grid%ep31, grid%sp33, grid%ep33
                  end if
               end if
               k = k+1
             end do
         end do
        !write(*,*) "local_cells_global_index_stagger nums= ", jc
        decomp_id = CCPL_register_normal_parallel_decomp("wrf_decomp",grid_H2D_id,num_local_cells,local_cells_global_index_stagger,annotation = "register parallel decomp for wrf")
        allocate(nlat,nlon,nsig,nsig_soil,ccpl_pt)
        nsig=grid%ed32-1
        grid_vert_id1=CCPL_register_V1D_Z_grid_via_model_data(my_comp_id, "wrf_vertical_grid_nsig", "Pa", grid%znu(1:nsig), annotation="register vertical grid znu for wrf")
        grid_vert_id2=CCPL_register_V1D_Z_grid_via_model_data(my_comp_id, "wrf_vertical_grid_nsig_1", "Pa", grid%znw(1:nsig+1), annotation="register vertical grid znw for wrf")
        call ESMF_TimeIntervalGet(domain_get_time_step(grid),S=time_step)
        call CCPL_set_normal_time_step(my_comp_id, time_step,annotation="setting the time step for wrf")
        call CCPL_set_normal_time_step(root_comp_id, time_step,annotation="setting the time step for wrf_root")

        call domain_clock_get(grid, current_time=currentTime)
        call ESMF_TimeGet(currentTime, YY=iyear,MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond)
        allocate(regional_time(6))
        regional_time(1)=iyear
        regional_time(2)=imonth
        regional_time(3)=iday
        regional_time(4)=ihour
        regional_time(5)=iminute
        regional_time(6)=isecond
        !read(regional_time,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)') iyear,imonth,iday,ihour,iminute,isecond
        regional_time_coord_values=(/1.,2.,3.,4.,5.,6./)
        grid_regional_time_id=CCPL_register_V1D_Z_grid_via_model_data(my_comp_id, "grid_for_regional_time", "virtual_gird", regional_time_coord_values, annotation="register simplified grid for regional time")
        field_id_regional_time=CCPL_register_field_instance(regional_time,"regional_time",-1,grid_regional_time_id,0,CCPL_TAG_CPL_REST,annotation="register time for GSI")
        nlon = grid%ed31-1
        nlat = grid%ed33-1
        nsig_soil=size(grid%smois,2)
        ccpl_pt=grid%p_top
        field_id_pt=CCPL_register_field_instance(ccpl_pt,"P_TOP",-1,my_comp_id,0,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register pressure top of the model")
        allocate(aeta1(1:nsig),aeta2(1:nsig),eta1(1:nsig+1),eta2(1:nsig+1))
        aeta1=grid%znu(1:nsig)
        aeta2=0.0
        eta1=grid%znw(1:nsig+1)
        eta2=0.0
        field_id_aeta1=CCPL_register_field_instance(aeta1,"ZNU_1",-1,grid_vert_id1,0,CCPL_TAG_CPL_REST, annotation="register eta values on half (mass) levels")
        field_id_aeta2=CCPL_register_field_instance(aeta2,"ZNU_2",-1,grid_vert_id1,0,CCPL_TAG_CPL_REST, annotation="register eta values on half (mass) levels")
        field_id_eta1=CCPL_register_field_instance(eta1,"ZNW_1",-1,grid_vert_id2,0,CCPL_TAG_CPL_REST, annotation="register eta values on full (w) levels")
        field_id_eta2=CCPL_register_field_instance(eta2,"ZNW_2",-1,grid_vert_id2,0,CCPL_TAG_CPL_REST, annotation="register eta values on full (w) levels")
        !DX_MC=one_single/(MAPFAC_M*rdx)    msft: Map scale factor on mass grid
        !DY_MC=one_single/(MAPFAC_M*rdy)    rdx/rdy: INVERSE X/Y GRID LENGTH
        allocate(ccpl_dx_mc(ims:ime,jms:jme),ccpl_dy_mc(ims:ime,jms:jme))
        ccpl_dx_mc=1.0/(grid%msft*grid%rdx) !Map scale factor on mass grid, x direction
        ccpl_dy_mc=1.0/(grid%msft*grid%rdy) !Map scale factor on mass grid, y direction
        field_id_dx_mc=CCPL_register_field_instance(ccpl_dx_mc,"MAPFAC_MX",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register x direction map scale factor")
        field_id_dy_mc=CCPL_register_field_instance(ccpl_dy_mc,"MAPFAC_MY",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register y direction map scale factor")
        allocate(ccpl_xlat(ims:ime,jms:jme),ccpl_xlong(ims:ime,jms:jme))
        !rad2deg=r45/atan(one_single)
        !rad2deg=1.0/(acos(-1.0)/180.0)
        !ccpl_xlat=grid%xlat/rad2deg   !XLAT: LATITUDE, SOUTH IS NEGATIVE
        !ccpl_xlong=grid%xlong/rad2deg  !XLONG: LONGITUDE, WEST IS NEGATIVE
        ccpl_xlat=grid%xlat   !XLAT: LATITUDE, SOUTH IS NEGATIVE
        ccpl_xlong=grid%xlong  !XLONG: LONGITUDE, WEST IS NEGATIVE
        field_id_xlat=CCPL_register_field_instance(ccpl_xlat,"XLAT",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register latitudes")
        field_id_xlong=CCPL_register_field_instance(ccpl_xlong,"XLONG",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register longitudes")      
        !mub  : base state dry air mass in column (Pa)
        !mu_2 : perturbation dry air mass in column (Pa)
        allocate(ccpl_mu(ims:ime,jms:jme),ccpl_mub(ims:ime,jms:jme))
        ccpl_mu = grid%mu_2(ims:ime,jms:jme)
        ccpl_mub = grid%mub(ims:ime,jms:jme)
        field_id_mu=CCPL_register_field_instance(ccpl_mu,"mu",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register perturbation dry air mass in column")
        field_id_mub=CCPL_register_field_instance(ccpl_mub,"mub",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register base state dry air mass in column")
        !psfc0=MU+MUB+pt
        allocate(ccpl_psfc0(ims:ime,jms:jme),ccpl_phb(ims:ime,jms:jme))
        ccpl_psfc0=grid%mu_2(ims:ime,jms:jme)+grid%mub(ims:ime,jms:jme)+grid%p_top
        field_id_psfc0=CCPL_register_field_instance(ccpl_psfc0,"psfc0",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="Pa",annotation="register surface pressure")
        ! PHB (zsfc*g) : base-state geopotential (m2 s-2)
        ! PHB=PHB(i,j,1) (field3(i,j,k),i=1,nlon_regional),j=1,nlat_regional)
        ccpl_phb=grid%phb(ims:ime,1,jms:jme)
        field_id_phb=CCPL_register_field_instance(ccpl_phb,"PHB",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="m2 s-2",annotation="register base-state geopotential")
        grid_3d_id1=CCPL_register_MD_grid_via_multi_grids(my_comp_id,"wrf_3d_grid_1",grid_H2D_id,grid_vert_id1,annotation="register wrf 3D-grid with nsig vertical levels")
        grid_3d_id2=CCPL_register_MD_grid_via_multi_grids(my_comp_id,"wrf_3d_grid_2",grid_H2D_id,grid_vert_id2,annotation="register wrf 3D-grid with nsig+1 vertical levels")    
        !th_phy_m_t0:  perturbation potential temperature theta-t0
        !u_phy: x-wind component at mass point
        !v_phy: y-wind component at mass point
        !moist(:,:,:,P_QV): Water vapor mixing ratio
        allocate(ccpl_t(ims:ime,jms:jme,1:nsig),ccpl_u(ims:ime,jms:jme,1:nsig),ccpl_v(ims:ime,jms:jme,1:nsig),ccpl_q(ims:ime,jms:jme,1:nsig))
        allocate(ccpl_w(ims:ime,jms:jme,1:nsig),ccpl_ph(ims:ime,jms:jme,1:nsig))
        do k=1, nsig
            ccpl_t(ims:ime,jms:jme,k)=grid%th_phy_m_t0(ims:ime,k,jms:jme)+300.0
            ccpl_u(ims:ime,jms:jme,k)=0.5*(grid%u_2(ims:ime,k,jms:jme)+grid%u_2(ims+1:ime+1,k,jms:jme))
            ccpl_u(ims:ime,jms:jme,k)=0.5*(grid%u_2(ims:ime,k,jms:jme)+grid%u_2(ims+1:ime+1,k,jms:jme))
            ccpl_v(ims:ime,jms:jme,k)=0.5*(grid%v_2(ims:ime,k,jms:jme)+grid%v_2(ims:ime,k,jms+1:jme+1))
            ccpl_w(ims:ime,jms:jme,k)=0.5*(grid%w_2(ims:ime,k,jms:jme)+grid%w_2(ims:ime,k+1,jms:jme))             
            ccpl_ph(ims:ime,jms:jme,k)=0.5*(grid%ph_2(ims:ime,k,:)+grid%ph_2(ims:ime,k+1,jms:jme))
            ccpl_q(ims:ime,jms:jme,k)=grid%moist(ims:ime,k,jms:jme,P_QV)
        end do
        field_id_t=CCPL_register_field_instance(ccpl_t,"ccpl_t",decomp_id,grid_3d_id1,0,CCPL_TAG_CPL_REST, field_unit="K",annotation="register perturbation potential temperature theta-t0")
        field_id_u=CCPL_register_field_instance(ccpl_u,"ccpl_u",decomp_id,grid_3d_id1,0,CCPL_TAG_CPL_REST, field_unit="m s-1",annotation="register x-wind component at mass point")
        field_id_v=CCPL_register_field_instance(ccpl_v,"ccpl_v",decomp_id,grid_3d_id1,0,CCPL_TAG_CPL_REST, field_unit="m s-1",annotation="register y-wind component at mass point")
        field_id_q=CCPL_register_field_instance(ccpl_q,"ccpl_q",decomp_id,grid_3d_id1,0,CCPL_TAG_CPL_REST, field_unit="kg kg-1",annotation="register water vapor mixing ratio")
        field_id_w=CCPL_register_field_instance(ccpl_w,"ccpl_w",decomp_id,grid_3d_id1,0,CCPL_TAG_CPL_REST, field_unit="m s-1",annotation="register z-wind component")
        field_id_ph=CCPL_register_field_instance(ccpl_ph,"ccpl_ph",decomp_id,grid_3d_id1,0,CCPL_TAG_CPL_REST, field_unit="m2 s-2",annotation="register perturbation geopotential")
        allocate(ccpl_landmask(ims:ime,jms:jme),ccpl_seaice(ims:ime,jms:jme),ccpl_sst(ims:ime,jms:jme))
        allocate(ccpl_ivgtyp(ims:ime,jms:jme),ccpl_isltyp(ims:ime,jms:jme),ccpl_vegfra(ims:ime,jms:jme))
        !xland: LAND MASK (1 FOR LAND, 2 FOR WATER)
        !landmask: LAND MASK (1 FOR LAND, 0 FOR WATER)
        !xice: SEA ICE FLAG
        !sstsk: SKIN SEA SURFACE TEMPERATURE
        !ivgtyp: DOMINANT VEGETATION CATEGORY
        !isltyp: DOMINANT SOIL CATEGORY
        !vegfra: VEGETATION FRACTION
        !ccpl_landmask=grid%landmask(ims:ime,jms:jme)
        ccpl_landmask=grid%xland(ims:ime,jms:jme)
        where(abs(ccpl_landmask-2.0)<0.01) ccpl_landmask=0.0
        ccpl_seaice=grid%xice(ims:ime,jms:jme)
        ccpl_sst=grid%sst(ims:ime,jms:jme)
        ccpl_ivgtyp=grid%ivgtyp(ims:ime,jms:jme)
        ccpl_isltyp=grid%isltyp(ims:ime,jms:jme)
        ccpl_vegfra=grid%vegfra(ims:ime,jms:jme)
        field_id_landmask=CCPL_register_field_instance(ccpl_landmask,"landmask",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register land mask (1 for land, 0 for water)")
        field_id_seaice=CCPL_register_field_instance(ccpl_seaice,"seaice",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register sea ice flag")
        field_id_sst=CCPL_register_field_instance(ccpl_sst,"ccpl_sst",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="K",annotation="register skin sea surface temperature")
        field_id_ivgtyp=CCPL_register_field_instance(ccpl_ivgtyp,"ivgtyp",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register domain vegetation category")
        field_id_isltyp=CCPL_register_field_instance(ccpl_isltyp,"isltyp",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register domain soil category")
        field_id_vegfra=CCPL_register_field_instance(ccpl_vegfra,"vegfra",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, annotation="register vegetation fraction") 
        allocate(ccpl_snow(ims:ime,jms:jme),ccpl_u10(ims:ime,jms:jme),ccpl_v10(ims:ime,jms:jme))
        allocate(ccpl_smois(ims:ime,jms:jme),ccpl_tslb(ims:ime,jms:jme),ccpl_tsk(ims:ime,jms:jme))
        !acsnom: ACCUMULATED MELTED SNOW (kg m-2)
        !u10: U at 10 M (m s-1)
        !v20: V at 10 M (m s-1)
        !sm: SOIL MOISTURES (m3 m-3)
        !tslb: SOIL TEMPERATURE (K)
        !tsk: SURFACE SKIN TEMPERATURE (K)
        ccpl_snow=grid%snow(ims:ime,jms:jme)
        ccpl_u10=grid%u10(ims:ime,jms:jme)
        ccpl_v10=grid%v10(ims:ime,jms:jme)
        ccpl_smois=grid%smois(ims:ime,1,jms:jme)
        ccpl_tslb=grid%tslb(ims:ime,1,jms:jme)
        ccpl_tsk=grid%tsk(ims:ime,jms:jme)
        field_id_snow=CCPL_register_field_instance(ccpl_snow,"snow",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST,field_unit="kg m-2",annotation="register accumulated melted snow") 
        field_id_u10=CCPL_register_field_instance(ccpl_u10,"ccpl_u10",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST,field_unit="m s-1",annotation="register field u10")
        field_id_v10=CCPL_register_field_instance(ccpl_v10,"ccpl_v10",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST,field_unit="m s-1",annotation="register field v10")
        field_id_smois=CCPL_register_field_instance(ccpl_smois,"smois",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="m3 m-3",annotation="register soil moistures (k=1)")
        field_id_tslb=CCPL_register_field_instance(ccpl_tslb,"tslb",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="K",annotation="register soil temperature (k=1)")
        field_id_tsk=CCPL_register_field_instance(ccpl_tsk,"tsk",decomp_id,grid_H2D_id,0,CCPL_TAG_CPL_REST, field_unit="K",annotation="register surface skin temperature")
        
        timer_id = CCPL_define_single_timer(my_comp_id, "seconds", 21600, 0, 0, annotation="define a single timer for wrf")      

        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Finish fields declare"
        write(*,*) "[CCPL <WRF>] wrf_set_comm, wrf_comm = ", root_comm, my_comm        
        write(*,*) "==============================================================================="
                 
        allocate(fields_id(1:32))
        fields_id(1)=field_id_pt
        fields_id(2)=field_id_aeta1
        fields_id(3)=field_id_aeta2
        fields_id(4)=field_id_eta1
        fields_id(5)=field_id_eta2
        fields_id(6)=field_id_dx_mc
        fields_id(7)=field_id_dy_mc
        fields_id(8)=field_id_xlat
        fields_id(9)=field_id_xlong
        fields_id(10)=field_id_psfc0
        fields_id(11)=field_id_phb
        fields_id(12)=field_id_landmask
        fields_id(13)=field_id_seaice
        fields_id(14)=field_id_sst
        fields_id(15)=field_id_ivgtyp
        fields_id(16)=field_id_isltyp
        fields_id(17)=field_id_vegfra
        fields_id(18)=field_id_snow
        fields_id(19)=field_id_u10
        fields_id(20)=field_id_v10
        fields_id(21)=field_id_smois
        fields_id(22)=field_id_tslb
        fields_id(23)=field_id_tsk
        fields_id(24)=field_id_regional_time
        fields_id(25)=field_id_mu
        fields_id(26)=field_id_mub
        fields_id(27)=field_id_t
        fields_id(28)=field_id_q
        fields_id(29)=field_id_u
        fields_id(30)=field_id_v
        fields_id(31)=field_id_w
        fields_id(32)=field_id_ph 


        control_vars(1) = nlat
        control_vars(2) = nlon
        control_vars(3) = nsig
        control_vars(4) = nsig_soil
        control_vars(5) = nprocs
        comp_or_grid_ids(1) = grid_H2D_id
        comp_or_grid_ids(2) = grid_vert_id1
        comp_or_grid_ids(3) = grid_vert_id2
        comp_or_grid_ids(4) = grid_3d_id1
        comp_or_grid_ids(5) = grid_3d_id2
        comp_or_grid_ids(6) = grid_regional_time_id
        decomp_ids = -1
        decomp_ids(1) = decomp_id
        field_inst_ids = fields_id(1:32)
        da_timers_id = timer_id

        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_MEMBER initialize"
        call start_timing
        gsi_ens_member_inst_id=CCPL_ensemble_procedures_inst_init(root_comp_id, my_comp_id, "gsi_ens_member", field_inst_ids, comp_or_grid_ids, decomp_ids, da_timers_id, control_vars, annotation="do ensemble procedure GSI_ENS_MEMBER initialize")
        call end_timing(" CCPL_OUT: GSI_ENS_MEMBER initialize ")
        write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_MEMBER initialize"
        write(*,*) "==============================================================================="
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_AVE initialize"
        call start_timing
        gsi_ens_ave_inst_id=CCPL_ensemble_procedures_inst_init(root_comp_id, my_comp_id, "gsi_ens_aver", field_inst_ids, comp_or_grid_ids, decomp_ids, da_timers_id, control_vars, annotation="do ensemble procedure GSI_ENS_AVE initialize")
        call end_timing(" CCPL_OUT: GSI_ENS_AVE initialize ")
        write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_AVE initialize"
        write(*,*) "==============================================================================="
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Start ensemble procedure ENKF initialize"
        call start_timing
        enkf_inst_id=CCPL_ensemble_procedures_inst_init(root_comp_id, my_comp_id, "enkf", field_inst_ids, comp_or_grid_ids, decomp_ids, da_timers_id, control_vars, annotation="do ensemble procedure ENKF initialize")
        call end_timing(" CCPL_OUT: ENKF initialize ")
        write(*,*) "[CCPL <WRF>] Finish ensemble procedure ENKF initialize"
        write(*,*) "==============================================================================="
        call CCPL_change_to_ensemble_dir(root_comp_id,ensemble_mem_id)

        
    END SUBROUTINE initialize_DA_instances
    
    SUBROUTINE run_DA_instances(grid)   
        
        USE ESMF_mod, ONLY : ESMF_TimeGet
        USE module_domain, ONLY : domain_clock_get
        USE module_utility, ONLY : WRFU_Time
        USE module_domain, ONLY : domain
        IMPLICIT NONE
        
        TYPE(domain) , POINTER             :: grid
        INTEGER                            :: err
        
        TYPE(WRFU_Time)             :: currentTime
        INTEGER                     :: iyear,imonth,iday,ihour,iminute,isecond
        LOGICAL                     :: prepare_da

        call domain_clock_get(grid, current_time=currentTime)
        call ESMF_TimeGet(currentTime, YY=iyear,MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond)
        
        call CCPL_check_current_time(my_comp_id, iyear*10000+imonth*100+iday, ihour*3600+iminute*60+isecond, "check time")
        if (CCPL_is_timer_on(timer_id)) then
            prepare_da = .true.
        else
            prepare_da = .false.
        end if
        if (prepare_da) then
            write(*,*) "==============================================================================="
            write(*,*) "[CCPL <WRF>] do CCPL fields copy in"
            call ccpl_fields_copy_in(grid)
        end if
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_AVE run"
        call start_timing
        call CCPL_ensemble_procedures_inst_run(gsi_ens_ave_inst_id,0,-1,annotation="do GSI_ENS_AVE run")
        call end_timing(" CCPL_OUT: GSI_ENS_AVE run ")
        write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_AVE run"
        write(*,*) "==============================================================================="
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Start ensemble procedure GSI_ENS_MEMBER run"
        call start_timing
        call CCPL_ensemble_procedures_inst_run(gsi_ens_member_inst_id,0,-1,annotation="do GSI_ENS_MEMBER run")
        call end_timing(" CCPL_OUT: GSI_ENS_MEMBER run ")
        write(*,*) "[CCPL <WRF>] Finish ensemble procedure GSI_ENS_MEMBER run"
        write(*,*) "==============================================================================="
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] Start ensemble procedure ENKF run"
        call start_timing
        call CCPL_ensemble_procedures_inst_run(enkf_inst_id,1,-1,annotation="do ENKF run")
        call end_timing(" CCPL_OUT: ENKF run ")
        write(*,*) "[CCPL <WRF>] Finish ensemble procedure ENKF run"
        write(*,*) "==============================================================================="

        call CCPL_change_to_ensemble_dir(root_comp_id,ensemble_mem_id)
        write(*,*) "==============================================================================="
        write(*,*) "[CCPL <WRF>] do CCPL fields copy out"
        call ccpl_fields_copy_out(grid)
        write(*,*) "==============================================================================="
        CALL CCPL_do_restart_write_IO(my_comp_id, .false.,annotation="WRF do restart write")
    

    END SUBROUTINE run_DA_instances

    SUBROUTINE ccpl_fields_copy_in(grid)


        USE ESMF_mod, ONLY : ESMF_TimeGet
        USE module_domain, ONLY : domain, domain_clock_get
        USE module_utility, ONLY : WRFU_Time
        use module_state_description, ONLY: P_QV
        IMPLICIT NONE
        
        TYPE(domain) , POINTER      :: grid
        INTEGER                     :: i,j,k,ip,err
        TYPE(WRFU_Time)             :: currentTime
        INTEGER                     :: iyear,imonth,iday,ihour,iminute,isecond

        call domain_clock_get(grid, current_time=currentTime)
        call ESMF_TimeGet(currentTime, YY=iyear,MM=imonth, DD=iday, H=ihour, M=iminute, S=isecond)
        
        regional_time(1)=iyear
        regional_time(2)=imonth
        regional_time(3)=iday
        regional_time(4)=ihour
        regional_time(5)=iminute
        regional_time(6)=isecond
        nsig = grid%ed32-1
        nlon = grid%ed31-1
        nlat = grid%ed33-1
        nsig_soil=size(grid%smois,2)
        ccpl_pt=grid%p_top
        
        aeta1=grid%znu(1:nsig)
        aeta2=0.0
        eta1=grid%znw(1:nsig+1)
        eta2=0.0
        ccpl_dx_mc=1.0/(grid%msft*grid%rdx) !Map scale factor on mass grid, x direction
        ccpl_dy_mc=1.0/(grid%msft*grid%rdy) !Map scale factor on mass grid, y direction
        
        !rad2deg=45.0/atan(1.0)
        !ccpl_xlat=grid%xlat/(45.0/atan(1.0))   !XLAT: LATITUDE, SOUTH IS NEGATIVE
        !ccpl_xlong=grid%xlong/(45.0/atan(1.0))  !XLONG: LONGITUDE, WEST IS NEGATIVE
        ccpl_xlat=grid%xlat  !XLAT: LATITUDE, SOUTH IS NEGATIVE
        ccpl_xlong=grid%xlong  !XLONG: LONGITUDE, WEST IS NEGATIVE
        ccpl_mu = grid%mu_2(ims:ime,jms:jme)
        ccpl_mub = grid%mub(ims:ime,jms:jme)       
        ccpl_psfc0=grid%mu_2(ims:ime,jms:jme)+grid%mub(ims:ime,jms:jme)+grid%p_top
        ccpl_phb=grid%phb(ims:ime,1,jms:jme)
        
        do k=1, nsig
            ccpl_t(ims:ime,jms:jme,k)=grid%th_phy_m_t0(ims:ime,k,jms:jme)+300.0
            ccpl_u(ims:ime,jms:jme,k)=0.5*(grid%u_2(ims:ime,k,jms:jme)+grid%u_2(ims+1:ime+1,k,jms:jme))
            ccpl_v(ims:ime,jms:jme,k)=0.5*(grid%v_2(ims:ime,k,jms:jme)+grid%v_2(ims:ime,k,jms+1:jme+1))
            ccpl_w(ims:ime,jms:jme,k)=0.5*(grid%w_2(ims:ime,k,jms:jme)+grid%w_2(ims:ime,k+1,jms:jme))             
            ccpl_ph(ims:ime,jms:jme,k)=0.5*(grid%ph_2(ims:ime,k,:)+grid%ph_2(ims:ime,k+1,jms:jme))
            ccpl_q(ims:ime,jms:jme,k)=grid%moist(ims:ime,k,jms:jme,P_QV)
        end do

        !ccpl_landmask=grid%landmask(ims:ime,jms:jme)
        ccpl_landmask=grid%xland(ims:ime,jms:jme)
        where(abs(ccpl_landmask-2.0)<0.01) ccpl_landmask=0.0
        ccpl_seaice=grid%xice(ims:ime,jms:jme)
        ccpl_sst=grid%sst(ims:ime,jms:jme)
        ccpl_ivgtyp=grid%ivgtyp(ims:ime,jms:jme)
        ccpl_isltyp=grid%isltyp(ims:ime,jms:jme)
        ccpl_vegfra=grid%vegfra(ims:ime,jms:jme)
        ccpl_snow=grid%snow(ims:ime,jms:jme)
        ccpl_u10=grid%u10(ims:ime,jms:jme)
        ccpl_v10=grid%v10(ims:ime,jms:jme)
        ccpl_smois(ims:ime,jms:jme)=grid%smois(ims:ime,1,jms:jme)
        ccpl_tslb=grid%tslb(ims:ime,1,jms:jme)
        ccpl_tsk=grid%tsk(ims:ime,jms:jme)

    END SUBROUTINE ccpl_fields_copy_in

    SUBROUTINE ccpl_fields_copy_out(grid)


        USE ESMF_mod, ONLY : ESMF_TimeGet
        USE module_domain, ONLY : domain, domain_clock_get
        USE module_utility, ONLY : WRFU_Time
        use module_state_description, ONLY: P_QV
        IMPLICIT NONE
        
        TYPE(domain) , POINTER      :: grid
        INTEGER                     :: i,j,k,ip,err

        nsig = grid%ed32-1
        nlon = grid%ed31-1
        nlat = grid%ed33-1
        nsig_soil=size(grid%smois,2)

        grid%mu_2(ims:ime,jms:jme) = ccpl_mu
        grid%phb(ims:ime,1,jms:jme) = ccpl_phb
        
        do k=1, nsig
            if (k==1) then
                !grid%w_2(ims:ime,k,jms:jme) = 1.5*(ccpl_w(:,:,k) - 0.5*ccpl_w(:,:,k+1))
                !grid%ph_2(ims:ime,k,jms:jme) = 1.5*(ccpl_ph(:,:,k) - 0.5*ccpl_ph(:,:,k+1))
            else ! k=2,nsig
                grid%w_2(ims:ime,k,jms:jme) = 0.5*(ccpl_w(:,:,k-1) + ccpl_w(:,:,k))
                grid%ph_2(ims:ime,k,jms:jme) = 0.5*(ccpl_ph(:,:,k-1) + ccpl_ph(:,:,k))
            end if
            grid%moist(ims:ime,k,jms:jme,P_QV) = ccpl_q(:,:,k)
            grid%th_phy_m_t0(ims:ime,k,jms:jme) = ccpl_t(:,:,k)-300.0    
            grid%u_2(ims:ime-1,k,jms:jme) = 0.5*(ccpl_u(ims:ime-1,:,k) + ccpl_u(ims+1:ime,:,k))
            grid%u_2(ime,:,:) = 1.5*ccpl_u(ime,:,:) -0.5* ccpl_u(ime-1,:,:)
            grid%v_2(ims:ime,k,jme) = 1.5*ccpl_v(:,jme,k) - 0.5*ccpl_v(:,jme-1,k)
            grid%v_2(ims:ime,k,jms:jme-1) = 0.5*(ccpl_v(:,jms:jme-1,k) + ccpl_v(:,jms+1:jme,k))
        end do
        grid%w_2(ims:ime,nsig+1,jms:jme) = 1.5*(ccpl_w(:,:,nsig) - 0.5*ccpl_w(:,:,nsig-1))
        grid%ph_2(ims:ime,nsig+1,jms:jme) = 1.5*(ccpl_ph(:,:,nsig) - 0.5*ccpl_ph(:,:,nsig-1))

    END SUBROUTINE ccpl_fields_copy_out
    
    end module

