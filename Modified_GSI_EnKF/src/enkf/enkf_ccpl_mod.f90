! sunchao 2019-03-08 
! 
! abstract: This module is used as an external link library for C-Coupler, 
! which contains the main module of enkf_wrf and introduces the link library of 
! C-Coupler.
  
  module enkf_ccpl_mod
  
  use CCPL_interface_mod ! c-coupler module
  use kinds, only: r_kind,r_double,i_kind
  use params, only: read_namelist,letkf_flag,readin_localization,lupd_satbiasc,&
                     numiter, nanals, lupd_obspace_serial
  use mpisetup, only: mpi_initialize, mpi_initialize_io, mpi_cleanup, nproc, &
                      numproc, mpi_wtime
  use enkf_obsmod, only: readobs, obfit_prior, obsprd_prior, &
                          deltapredx, nobs_sat, obfit_post, obsprd_post, &
                          obsmod_cleanup, biasprednorminv
  use innovstats, only: print_innovstats
  use gridinfo, only: getgridinfo, gridinfo_cleanup, npts,lonsgrd,latsgrd
  use statevec, only: read_ensemble, write_ensemble, statevec_cleanup
  use loadbal, only: load_balance, loadbal_cleanup
  use enkf, only: enkf_update
  use letkf, only: letkf_update
  use radinfo, only: radinfo_write, predx, jpch_rad, npred
  use inflation, only: inflate_ens
  use radinfo, only: init_rad, init_rad_vars, final_rad_vars
  use omp_lib, only: omp_get_max_threads
  use enkf_ccpl_kinds 
  use enkf_ccpl_grid_vars_mod
  use enkf_ccpl_coupling_mod
  use convinfo, only: convinfo_destroy
  use state_vectors, only: final_anasv
  implicit none

  integer(i_kind) j,n,nth
  real(r_double) t1,t2,tt1,tt2

  public :: enkf_grid_ccpl_init
  public :: enkf_grid_ccpl_run
  public :: enkf_ccpl_init
  public :: enkf_ccpl_run
  public :: enkf_ccpl_finalize

  contains
!==============================================================================  
! DESCRIPTION: This module contains the enkf_wrf main program, which has been split
! in initialize/run/finalize segments, and subroutines created for these steps:
! enkf_ccpl_init(), enkf_ccpl_run() and enkf_ccpl_finalize(). For more details
! about the enkf main program, refer to enkf_main.f90 in the same directory.
!==============================================================================
    
    subroutine enkf_grid_ccpl_init(ccpl_procedure_inst_id) bind(c)

       implicit none
       integer, intent(in)    ::  ccpl_procedure_inst_id
       tt1 = mpi_wtime()
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"
       write(*,*) "[CCPL <ENKF>] Start enkf grid vars initialization"
       write(*, *) "[CCPL <ENKF>] ccpl_procedure_inst_id is ", ccpl_procedure_inst_id
       t1 = mpi_wtime()
       enkf_grid_instance_id = ccpl_procedure_inst_id
       enkf_ccpl_comm_grid = CCPL_external_procedures_get_local_comm(ccpl_procedure_inst_id, annotation="get local comm")
       t2 = mpi_wtime()
       if (nproc == 0) print *,'time in enkf_grid_ccpl_init::get_enkf_grid_instance&comm =',t2-t1,'on proc',nproc
       t1 = mpi_wtime()
       call register_grids_decomps_from_ccpl
       call enkf_ccpl_grid_vars_transfer
       t2 = mpi_wtime()
       if (nproc == 0) print *,'time in enkf_grid_ccpl_init::grids_decomps&grid_vars_declare =',t2-t1,'on proc',nproc
       tt2 = mpi_wtime()
       if (nproc == 0) print *,'time in enkf_grid_ccpl_init =',tt2-tt1,'on proc',nproc
       write(*,*) "[CCPL <ENKF>] Finish enfk grid vars initialization"
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"
    
    end subroutine enkf_grid_ccpl_init 

    subroutine enkf_ccpl_init(ccpl_procedure_inst_id) bind(c)

      implicit none
      integer, intent(in)    ::  ccpl_procedure_inst_id
       
      ccpl_run = .True.
      io_read  = .False.
      io_check = .False.
      if(.false.) then ! ture for io_check
        ccpl_run = .False.
        io_read  = .True. 
        io_check = .True.
      end if
      
      tt1 = mpi_wtime()
      write(*,*) "++++++++++++++++++++++++++++++++++++++++"
      write(*,*) "[CCPL <ENKF>] Start enkf initialization"
      !write(*, *) "[CCPL <ENKF>] ccpl_procedure_inst_id is ", ccpl_procedure_inst_id
      !write(*,*) "ccpl_znu", ccpl_znu 
      t1 = mpi_wtime()
      enkf_ccpl_comm=CCPL_external_procedures_get_local_comm(ccpl_procedure_inst_id, annotation="get local comm")
      enkf_instance_id = ccpl_procedure_inst_id
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_ccpl_init::get_enkf_instance&comm =',t2-t1,'on proc',nproc
      t1 = mpi_wtime()
      call register_grids_decomps
      call enkf_fields_declare
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_ccpl_init::grids_decomps&fields_declare =',t2-t1,'on proc',nproc
      tt2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_ccpl_init =',tt2-tt1,'on proc',nproc
      
      write(*,*) "[CCPL <ENKF>] Finish enkf initialization"
      write(*,*) "++++++++++++++++++++++++++++++++++++++++"
    
    end subroutine enkf_ccpl_init

    subroutine enkf_grid_ccpl_run() bind(c)

       implicit none
       !tt1 = mpi_wtime()
       !write(*,*) "-------enkf_grid_run------"
       !write(*,*) "ccpl_znu", ccpl_znu   
       !write(*,*) "-------finish enkf_grid_run------"
       !tt2 = mpi_wtime()
       !if (nproc == 0) print *,'time in enkf_grid_ccpl_run =',tt2-tt1,'on proc',nproc
   
    end subroutine enkf_grid_ccpl_run

    subroutine enkf_ccpl_run() bind(c)

      use enkf_ccpl_coupling_mod, only: enkf_fields_transfer_in, enkf_fields_transfer_out
      implicit none
      
      write(*,*) "++++++++++++++++++++++++++++++++++++++++"
      write(*,*) "[CCPL <ENKF>] Start enkf run"
      
      !write(*,*) "ccpl_psfc x*y,z zise: ", size(ccpl_psfc(:,1)), size(ccpl_psfc(1,:))
      !write(*,*) "ccpl_psfc(1:10,1): ", ccpl_psfc(1:10,1)
      !write(*,*) "ccpl_t(1:10,1,1): ", ccpl_t(1:10,1,1)
      
      ! initialize MPI.      
      call mpi_initialize()
      tt1 = mpi_wtime()
      t1 = mpi_wtime()
      if (nproc==0) call w3tagb('ENKF_ANL',2011,0319,0055,'NP25')
      ! Initial radinfo variables (some flags may be over-ridden in enkf namelist)
      call init_rad()

      ! read namelist.
      call read_namelist()

      ! initialize MPI communicator for IO tasks.
      call mpi_initialize_io(nanals)

      ! Initialize derived radinfo variables
      call init_rad_vars()

      nth= omp_get_max_threads()
      if(nproc== 0)write(6,*) 'enkf_main:  number of threads ',nth
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf initialize (MPI, radinfo variables, namelist, MPI communicator for IO tasks, derived radinfo variables) =',t2-t1,'on proc',nproc
      t1 = mpi_wtime()
      call enkf_fields_transfer_in
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_ccpl_run::enkf_fields_transfer_out =',t2-t1,'on proc',nproc
      ! read horizontal grid information and pressure fields from
      ! 6-h forecast ensemble mean file.
      t1 = mpi_wtime()
      call getgridinfo()
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in getgridinfo =',t2-t1,'on proc',nproc

      ! read obs, initial screening.
      t1 = mpi_wtime()
      call readobs()
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in read_obs =',t2-t1,'on proc',nproc

      ! print innovation statistics for prior on root task.
      if (nproc == 0) then
         print *,'innovation statistics for prior:'
         call print_innovstats(obfit_prior, obsprd_prior)
      end if
      ! read in vertical profile of horizontal and  vertical localization length
      ! scales, set values for each ob.
      if (readin_localization) call read_locinfo()
      ! do load balancing (partitioning of grid points, observations among processors)
      t1 = mpi_wtime()
      call load_balance()
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in load_balance =',t2-t1,'on proc',nproc


      ! read in ensemble members, distribute pieces to each task.
      t1 = mpi_wtime()
      call read_ensemble()
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in read_ensemble =',t2-t1,'on proc',nproc

      t1 = mpi_wtime()
      ! state and bias correction coefficient update iteration.
      if(letkf_flag) then
         ! do ob space update using serial filter if desired
         if (lupd_obspace_serial) call enkf_update()
         call letkf_update()
      else
         call enkf_update()
      end if
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_update =',t2-t1,'on proc',nproc

      ! posterior inflation.
      t1 = mpi_wtime()
      call inflate_ens()
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in inflate_ens =',t2-t1,'on proc',nproc 
      
      ! print innovation statistics for posterior on root task.
      if (nproc == 0 .and. numiter > 0) then
         print *,'innovation statistics for posterior:'
         call print_innovstats(obfit_post, obsprd_post)
      ! write out bias coeffs on root.
         if (nobs_sat > 0 .and. lupd_satbiasc) then
            ! re-scale bias coefficients.
            do j=1,jpch_rad
                do n=1,npred
                   predx(n,j) = predx(n,j)*biasprednorminv(n)
                enddo
            enddo
            call radinfo_write()
         end if
      end if
     ! free memory (radinfo memory freed in radinfo_write)
     ! and write out analysis ensemble.
      call obsmod_cleanup()

      t1 = mpi_wtime()
      call write_ensemble()
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in write_ensemble =',t2-t1,'on proc',nproc

      t1 = mpi_wtime()
      call enkf_fields_transfer_out
      t2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_ccpl_run::enkf_fields_transfer_out =',t2-t1,'on proc',nproc

      call gridinfo_cleanup()
      call statevec_cleanup()
      call loadbal_cleanup()

      call convinfo_destroy()
      call final_anasv()
      call final_rad_vars()

      ! write log file (which script can check to verify completion).
      if (nproc .eq. 0) then
         call write_logfile()
      endif
      ! finalize MPI.
      if (nproc==0) call w3tage('ENKF_ANL')
      call mpi_cleanup()
      
      tt2 = mpi_wtime()
      if (nproc == 0) print *,'time in enkf_ccpl_run =',tt2-tt1,'on proc',nproc
      
      write(*,*) "[CCPL <ENKF>] Finish enkf run"
      write(*,*) "++++++++++++++++++++++++++++++++++++++++"

    end subroutine enkf_ccpl_run

    subroutine enkf_ccpl_finalize() bind(c)

      implicit none

    end subroutine enkf_ccpl_finalize



  end module enkf_ccpl_mod

 
