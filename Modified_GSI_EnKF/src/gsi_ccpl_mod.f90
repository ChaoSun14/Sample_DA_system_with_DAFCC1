! sunchao 2019-03-08 
! 
! abstract: This module is used as an external link library for C-Coupler, 
! which contains the main module of GSI and introduces the link library of 
! C-Coupler.
  
  module gsi_ccpl_mod

  use gsimod, only: gsimain_initialize,gsimain_run,gsimain_finalize
  use gsi_4dvar, only: l4dvar
  use gsi_4dcouplermod, only: gsi_4dcoupler_init_traj
  use gsi_4dcouplermod, only: gsi_4dcoupler_final_traj
  use timermod, only: timer_pri
  use kinds, only: i_kind,r_double
  use mpeu_util, only: die
  use CCPL_interface_mod ! c-coupler module
  use gsi_ccpl_kinds
  use gsi_ccpl_grid_vars_mod
  use gsi_ccpl_coupling_mod
  use mpimod, only: npe, mype
  use mpi
  implicit none
  
  public :: gsi_grid_ccpl_init
  public :: gsi_grid_ccpl_run
  public :: gsi_ccpl_init
  public :: gsi_ccpl_run
  public :: gsi_ccpl_finalize

  real(r_double) :: t1,t2,tt1,tt2
  
  contains

!==============================================================================  
! DESCRIPTION: This module contains the GSI main program, which has been split
! in initialize/run/finalize segments, and subroutines created for these steps:
! gsimain_initialize(), gsimain_run() and gsimain_finalize(). For more details
! about the GSI main program, refer to gsimain.f90 in the same directory.
!==============================================================================
   
  subroutine gsi_grid_ccpl_init(ccpl_procedure_inst_id) bind(c)

       implicit none
       integer, intent(in)    ::  ccpl_procedure_inst_id
       tt1 = mpi_wtime()
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"
       write(*,*) "[CCPL <GSI>] Start gsi grid vars initialization"
       write(*,*) "[CCPL <GSI>] ccpl_procedure_inst_id is ", ccpl_procedure_inst_id
       t1 = mpi_wtime()
       gsi_grid_instance_id = ccpl_procedure_inst_id
       gsi_grid_ccpl_comm = CCPL_external_procedures_get_local_comm(ccpl_procedure_inst_id, annotation="get local comm")
       t2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_grid_ccpl_init::get_gsi_grid_instance&comm =',t2-t1,'on proc', mype
       t1 = mpi_wtime()
       call register_grids_decomps_from_wrf
       call gsi_ccpl_grid_vars_transfer
       t2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_grid_ccpl_init::register_grids_decomps&grid_vars_declare =',t2-t1,'on proc', mype
       tt2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_grid_ccpl_init =',tt2-tt1,'on proc', mype
       write(*,*) "[CCPL <GSI>] Finish gsi grid vars initialization"
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"

   end subroutine gsi_grid_ccpl_init   


   subroutine gsi_ccpl_init(ccpl_procedure_inst_id) bind(c)
       !use gridmod
       implicit none
       integer, intent(in)    ::  ccpl_procedure_inst_id

       ccpl_run = .True.
       io_read  = .False. 
       io_check = .False.
       
       if (.false.) then !true: Enable io_read and check io_read with online_transn
           ccpl_run = .False.
           io_read  = .True. 
           io_check = .True.
       end if
       

       tt1 = mpi_wtime()
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"
       write(*,*) "[CCPL <GSI>] Start gsi initialization"
       !write(*,*) "[CCPL <GSI>] ccpl_procedure_inst_id is ", ccpl_procedure_inst_id
       !write(*, *) "ccpl_nlon_regional", ccpl_nlon_regional
       !write(*, *) "ccpl_glat(1,1:10)", ccpl_glat(1,1:10)
       t1 = mpi_wtime()
       gsi_instance_id = ccpl_procedure_inst_id 
       gsi_ccpl_comm_all = CCPL_external_procedures_get_local_comm(ccpl_procedure_inst_id, annotation="get local comm")
       t2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_ccpl_init::get_gsi_instance&comm =',t2-t1,'on proc', mype
       t1 = mpi_wtime()
       call register_grids_decomps
       call gsi_fields_declare
       t2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_ccpl_init::register_grids_decomps&fields_declare =',t2-t1,'on proc', mype
      
       !write(*,*) "[CCPL <GSI>] ccpl_regional_time:", ccpl_regional_time
       write(*,*) "[CCPL <GSI>] Finish gsi initialization"
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"
       tt2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_ccpl_init =',tt2-tt1,'on proc', mype

   end subroutine gsi_ccpl_init

   
   subroutine gsi_grid_ccpl_run() bind(c)
       
       implicit none
       if (color.eq.0.or.color.eq.-1) then
       !t1 = mpi_wtime()
       !write(*,*) "[CCPL <GSI>] gsi grid run"
       !write(*,*) "[CCPL <GSI>] size(ccpl_glat)", size(ccpl_glat,1), size(ccpl_glat,2) 
       !t2 = mpi_wtime()
       !if (mype == 0) print *,'time in gsi_grid_ccpl_run =',t2-t1,'on proc', mype
       end if

   end subroutine gsi_grid_ccpl_run

   subroutine gsi_ccpl_run() bind(c)
       
       implicit none
       integer(i_kind):: ier, ierror
       character(len=*),parameter:: myname='gsimain'
       
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"
       write(*,*) "[CCPL <GSI>] Start gsi run"
       write(*,*) "[CCPL <GSI>] ccpl_regional_time: ", ccpl_regional_time

       tt1 = mpi_wtime()
       t1 = mpi_wtime()
       call boundary_fill_all_fields
       t2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_ccpl_run::boundary_fill_all_fields =',t2-t1,'on proc', mype
       !write(*,*) "[CCPL <GSI>] ccpl_regional_time:", ccpl_regional_time
       !write(*,*) "[CCPL <GSI>]  ccpl_glat(1,50:70), ccpl_glon(1,50:70) ", ccpl_glat(1,50:70), ccpl_glon(1,50:70)
       !write(*,*) "[CCPL <GSI>] ccpl_landmask:", ccpl_landmask
       if (color.eq.0.or.color.eq.-1) then
       
       call gsimain_initialize
       
       call gsimain_run(init_pass=.true.,last_pass=.true.)
       
       ! Finalize atmospheric AD and TL model trajectory
       if(l4dvar) then
           call gsi_4dcoupler_final_traj(rc=ier)
           if(ier/=0) call die(myname,'gsi_4dcoupler_final_traj(), rc =',ier)
       endif
       call timer_pri(6)
       
       call gsimain_finalize
       
       tt2 = mpi_wtime()
       if (mype == 0) print *,'time in gsi_ccpl_run =',tt2-tt1,'on proc', mype
       end if
       
       write(*,*) "[CCPL <GSI>] Finish gsi run"
       write(*,*) "++++++++++++++++++++++++++++++++++++++++"

    end subroutine gsi_ccpl_run
    
    subroutine gsi_ccpl_finalize() bind(c)
        
        implicit none

    end subroutine gsi_ccpl_finalize
     
  end module gsi_ccpl_mod


