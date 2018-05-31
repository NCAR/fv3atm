module IPD_CCPP_driver

  use IPD_typedefs,       only: IPD_init_type,                       &
                                IPD_control_type,  IPD_data_type,    &
                                IPD_diag_type,     IPD_restart_type, &
                                IPD_interstitial_type, IPD_fastphys_type

  use ccpp_api,           only: ccpp_t,                              &
                                ccpp_init,                           &
                                ccpp_finalize,                       &
                                ccpp_physics_init,                   &
                                ccpp_physics_run,                    &
                                ccpp_physics_finalize,               &
                                ccpp_field_add

  use CCPP_data,          only: cdata_domain,                        &
                                cdata_block
  
! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules.inc"
! End include auto-generated list of modules for ccpp

  use iso_c_binding,      only: c_loc

  implicit none

!------------------------------------------------------!
!  Pointer to CCPP containers defined in CCPP_data     !
!------------------------------------------------------!
  type(ccpp_t), pointer :: cdata => null()

!----------------
! Public Entities
!----------------
! functions
  public IPD_CCPP_step

  CONTAINS
!*******************************************************************************************

!-------------------------------
!  IPD step generalized for CCPP
!-------------------------------
  subroutine IPD_CCPP_step (step, IPD_Control, IPD_Data, IPD_Diag, IPD_Restart, &
                            IPD_Interstitial, IPD_Fastphys, nblks, ccpp_suite, ierr)

#ifdef OPENMP
    use omp_lib
#endif

    implicit none

    character(len=*),                    intent(in)              :: step
    type(IPD_control_type),      target, intent(inout), optional :: IPD_Control
    type(IPD_data_type),         target, intent(inout), optional :: IPD_Data(:)
    type(IPD_diag_type),         target, intent(inout), optional :: IPD_Diag(:)
    type(IPD_restart_type),      target, intent(inout), optional :: IPD_Restart
    type(IPD_interstitial_type), target, intent(inout), optional :: IPD_Interstitial(:)
    type(IPD_fastphys_type),     target, intent(inout), optional :: IPD_Fastphys
    integer,                     target, intent(in),    optional :: nblks
    character(len=*),                    intent(in),    optional :: ccpp_suite
    integer,                             intent(out)             :: ierr
    ! Local variables
    integer :: nb
    integer :: nthrds, nt

    ierr = 0

#ifdef OPENMP
    nthrds = omp_get_max_threads()
#else
    nthrds = 1
#endif

    if (trim(step)=="init") then

      if (.not.present(IPD_Control)) then
        write(0,*) 'Optional argument IPD_Control required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      if (.not.present(IPD_Data)) then
        write(0,*) 'Optional argument IPD_Data required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      if (.not.present(IPD_Diag)) then
        write(0,*) 'Optional argument IPD_Diag required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      if (.not.present(IPD_Restart)) then
        write(0,*) 'Optional argument IPD_Restart required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      if (.not.present(IPD_Interstitial)) then
          write(0,*) 'Optional argument IPD_Interstitial required for IPD-CCPP init step'
          ierr = 1
          return
      end if

      if (.not.present(IPD_Fastphys)) then
        write(0,*) 'Optional argument IPD_Fastphys required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      if (.not.present(ccpp_suite)) then
        write(0,*) 'Optional argument ccpp_suite required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      ! Associate cdata with the global/domain cdata structure;
      ! this is needed because ccpp_fields.inc uses 'cdata' in
      ! its ccpp_field_add statements.
      associate_domain: associate (cdata => cdata_domain)
        !--- Initialize CCPP framework
        call ccpp_init(trim(ccpp_suite), cdata, ierr)
        if (ierr/=0) then
          write(0,*) 'An error occurred in ccpp_init'
          return
        end if

        ! ccpp_fields.inc contains block- and/or thread-dependent data, which
        ! is not used when physics are run over the entire domain of an MPI task;
        ! need to set block number nb and thread number nt to safe values
        nb = 1
        nt = 1
! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields.inc"
! End include auto-generated list of calls to ccpp_field_add

        !--- Initialize CCPP physics
        call ccpp_physics_init(cdata, ierr)
        if (ierr/=0) then
          write(0,*) 'An error occurred in ccpp_physics_init'
          return
        end if
      end associate associate_domain

      ! Allocate cdata structures for blocks and threads
      allocate(cdata_block(1:nblks,1:nthrds))

! For some reason, this didn't work with PGI - DH* 20180517 - check if this is still true
#ifndef __PGI
      ! Loop over blocks for each of the threads
!$OMP parallel num_threads (nthrds) &
!$OMP          default (shared) &
!$OMP          private (nb,nt,cdata) &
!$OMP          reduction (+:ierr)
#ifdef OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
#else
      do nt=1,nthrds
#endif
      do nb = 1,nblks
        ! Associate cdata with the cdata structure for this block
        ! and thread; this is needed because ccpp_fields.inc uses
        ! 'cdata' in its ccpp_field_add statements.
        associate_block: associate (cdata => cdata_block(nb,nt))
! Again, for some reason, this didn't work with PGI - DH* 20180517 - check if this is still true
#ifndef __PGI
          !--- Initialize CCPP framework for blocks/threads, use suite from scalar cdata to avoid reading the SDF multiple times
          call ccpp_init(ccpp_suite, cdata, ierr, suite=cdata_domain%suite)
#else
          !--- Initialize CCPP framework for blocks/threads, cannot use suite from scalar cdata with PGI (crashes)
          call ccpp_init(ccpp_suite, cdata, ierr)
#endif
          if (ierr/=0) then
            write(0,'(2(a,i4))') "An error occurred in ccpp_init for block ", nb, " and thread ", nt
            exit
          end if
! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields.inc"
! End include auto-generated list of calls to ccpp_field_add
        end associate associate_block
      end do
#ifndef __PGI
!$OMP end parallel
#else
      end do
#endif
      if (ierr/=0) return

    else if (trim(step)=="fast_physics") then

      call ccpp_physics_run(cdata_domain, group_name='fast_physics', ierr=ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_physics_run for group fast_physics"
         return
      end if

    ! Finalize
    else if (trim(step)=="finalize") then

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP finalize step'
        ierr = 1
        return
      end if

      !--- Finalize CCPP physics
      call ccpp_physics_finalize(cdata_domain, ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_physics_finalize"
         return
      end if

!$OMP parallel num_threads (nthrds) &
!$OMP          default (shared) &
!$OMP          private (nb,nt) &
!$OMP          reduction (+:ierr)
#ifdef OPENMP
      nt = omp_get_thread_num()+1
#else
      nt = 1
#endif
      do nb = 1,nblks
        !--- Finalize CCPP framework for blocks/threads
        call ccpp_finalize(cdata_block(nb,nt), ierr)
        if (ierr/=0) then
           write(0,'(a,i4,a,i4)') "An error occurred in ccpp_finalize for block ", nb, " and thread ", nt
           exit
        end if
      end do
!$OMP end parallel
      if (ierr/=0) return

      ! Deallocate cdata structure for blocks and threads
      deallocate(cdata_block)

      !--- Finalize CCPP framework
      call ccpp_finalize(cdata_domain, ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_finalize"
         return
      end if

    else

      write(0,'(2a)') 'Error, undefined IPD step ', trim(step)
      ierr = 1
      return

    end if

  end subroutine IPD_CCPP_step

end module IPD_CCPP_driver
