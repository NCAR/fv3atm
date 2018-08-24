module IPD_CCPP_driver

  use IPD_typedefs,       only: IPD_init_type,                       &
                                IPD_control_type,  IPD_data_type,    &
                                IPD_diag_type,     IPD_restart_type, &
                                IPD_interstitial_type

  use ccpp_api,           only: ccpp_t,                              &
                                ccpp_init,                           &
                                ccpp_finalize,                       &
                                ccpp_physics_init,                   &
                                ccpp_physics_run,                    &
                                ccpp_physics_finalize,               &
                                ccpp_field_add,                      &
                                ccpp_initialized

  use CCPP_data,          only: cdata_tile,                          &
                                cdata_domain,                        &
                                cdata_block,                         &
                                ccpp_suite,                          &
                                CCPP_shared
  
! Begin include auto-generated list of modules for ccpp
#include "ccpp_modules_slow_physics.inc"
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
                            IPD_Interstitial, nblks, ierr)

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
    integer,                     target, intent(in),    optional :: nblks
    integer,                             intent(out)             :: ierr
    ! Local variables
    integer :: nb
    integer :: nthrds, nt
    integer :: ierr2

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

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP init step'
        ierr = 1
        return
      end if

      ! Associate cdata with the global/domain cdata structure;
      ! this is needed because ccpp_fields.inc uses 'cdata' in
      ! its ccpp_field_add statements.
      associate_domain: associate (cdata => cdata_domain)
        !--- Initialize CCPP framework, if cdata_tile for fast physics
        !    is already initialized, use its suite to avoid reading the
        !    SDF multiple times.
        if (ccpp_initialized(cdata_tile)) then
          call ccpp_init(trim(ccpp_suite), cdata, ierr=ierr, cdata_target=cdata_tile)
        else
          call ccpp_init(trim(ccpp_suite), cdata, ierr=ierr)
        end if
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
#include "ccpp_fields_slow_physics.inc"
! End include auto-generated list of calls to ccpp_field_add

      end associate associate_domain

      ! Allocate cdata structures for blocks and threads
      allocate(cdata_block(1:nblks,1:nthrds))

      ! Loop over all blocks and threads
      do nt=1,nthrds
        do nb=1,nblks
          ! Associate cdata with the cdata structure for this block
          ! and thread; this is needed because ccpp_fields.inc uses
          ! 'cdata' in its ccpp_field_add statements.
          associate_block: associate (cdata => cdata_block(nb,nt))
            !--- Initialize CCPP framework for blocks/threads, use suite from scalar cdata
            !    to avoid reading the SDF multiple times. If cdata_tile is initialized, use
            !    this version (since cdata_domain is just a copy), otherwise use cdata_domain
            if (ccpp_initialized(cdata_tile)) then
              call ccpp_init(trim(ccpp_suite), cdata, ierr=ierr, cdata_target=cdata_tile)
            else
              call ccpp_init(trim(ccpp_suite), cdata, ierr=ierr, cdata_target=cdata_domain)
            end if
            if (ierr/=0) then
              write(0,'(2(a,i4))') "An error occurred in ccpp_init for block ", nb, " and thread ", nt
              return
            end if
! Begin include auto-generated list of calls to ccpp_field_add
#include "ccpp_fields_slow_physics.inc"
! End include auto-generated list of calls to ccpp_field_add
          end associate associate_block
        end do
      end do

   else if (trim(step)=="physics_init") then

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP physics_init step'
        ierr = 1
        return
      end if

      ! Loop over blocks, don't use threading on the outside but allowing threading
      ! inside the initialization, because physics initialization code does a lot of
      ! reading/writing data from disk, allocating fields, etc.
      CCPP_shared(:)%nthreads = nthrds

      ! DH* TODO - I believe that they way this works in FV3, physics init can only
      ! affect block- and thread-independent data - hence, it would be sufficient to
      ! run physics_init once over cdata_domain?!? Update notes accordingly if true. *DH
      ! If physics init affect block data, then changes must be made to parsing the
      ! SDF (cannot use pointers to point to the central SDF, because a scheme will
      ! be considered as initialized because of its CCPP-internal attribute after a
      ! first call to it), and the is_initialized logic inside the scheme (independent)
      ! of CCPP must be removed for those schemes that need to be initialized multiple
      ! times once per block.

      ! Since physics init can only affect block- (and not thread-) dependent data
      ! structures, it is sufficient to run this over all blocks for one thread only
      nt = 1
      do nb=1,nblks
        !--- Initialize CCPP physics
        call ccpp_physics_init(cdata_block(nb,nt), ierr=ierr)
        if (ierr/=0) then
          write(0,'(2(a,i4))') "An error occurred in ccpp_physics_init for block ", nb, " and thread ", nt
          return
        end if
      end do

      ! Reset number of threads available to physics schemes to default value
      CCPP_shared(:)%nthreads = 1

   else if (trim(step)=="time_vary") then

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP time_vary step'
        ierr = 1
        return
      end if

      ! Loop over blocks, don't use threading on the outside but allowing threading
      ! inside the time_vary routines. This is because the time_vary routines contain
      ! calls to gcycle.f90 and sfcsub.F, which do a lot of I/O and are highly
      !inefficient if executed nthread times.
      CCPP_shared%nthreads = nthrds

      ! Since the time_vary steps only use data structures for all blocks (except the
      ! CCPP-internal variables ccpp_error_flag and ccpp_error_message, which are defined
      ! for all cdata structures independently), we can use cdata_domain here.
      call ccpp_physics_run(cdata_domain, group_name="time_vary", ierr=ierr)
      if (ierr/=0) then
        write(0,'(2(a,i4))') "An error occurred in ccpp_physics_run for group time_vary"
        return
      end if

      ! Reset number of threads available to physics schemes to default value
      CCPP_shared%nthreads = 1

    ! Radiation and stochastic physics
   else if (trim(step)=="radiation" .or. trim(step)=="physics" .or. trim(step)=="stochastics") then

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP ' // trim(step) // ' step'
        ierr = 1
        return
      end if

!$OMP parallel do num_threads (nthrds) &
#ifdef MEMCHECK
!$OMP          schedule (static,nblks),&
#else
!$OMP          schedule (dynamic,1),   &
#endif
!$OMP          default (shared)        &
!$OMP          private (nb,nt,ierr2)   &
!$OMP          reduction (+:ierr)
      do nb = 1,nblks
#ifdef OPENMP
        nt = omp_get_thread_num()+1
#else
        nt = 1
#endif
        !--- Call CCPP radiation group
        call ccpp_physics_run(cdata_block(nb,nt), group_name=trim(step), ierr=ierr2)
        if (ierr2/=0) then
           write(0,'(a,i4,a,i4)') "An error occurred in ccpp_physics_run for group " // trim(step) // ", block ", nb, " and thread ", nt
           ierr = ierr + ierr2
        end if
      end do
!$OMP end parallel do
      if (ierr/=0) return

   ! Finalize
   else if (trim(step)=="finalize") then

      if (.not.present(nblks)) then
        write(0,*) 'Optional argument nblks required for IPD-CCPP finalize step'
        ierr = 1
        return
      end if

      ! Fast physics are finalized in atmosphere_end, loop over
      ! all blocks and threads to finalize all other physics
      do nt=1,nthrds
        do nb=1,nblks
          !--- Finalize CCPP physics
          call ccpp_physics_finalize(cdata_block(nb,nt), ierr=ierr)
          if (ierr/=0) then
            write(0,'(a,i4,a,i4)') "An error occurred in ccpp_physics_finalize for block ", nb, " and thread ", nt
            return
          end if
          !--- Finalize CCPP framework for blocks/threads
          call ccpp_finalize(cdata_block(nb,nt), ierr=ierr)
          if (ierr/=0) then
            write(0,'(a,i4,a,i4)') "An error occurred in ccpp_finalize for block ", nb, " and thread ", nt
            return
          end if
        end do
      end do

      ! Deallocate cdata structure for blocks and threads
      deallocate(cdata_block)

      !--- Finalize CCPP framework for domain
      call ccpp_finalize(cdata_domain, ierr=ierr)
      if (ierr/=0) then
         write(0,'(a)') "An error occurred in ccpp_finalize"
         return
      end if

      !--- Finalize CCPP framework for fast physics last (all other frameworks point to cdata_tile's suite)
      if (ccpp_initialized(cdata_tile)) then
         call ccpp_finalize(cdata_tile, ierr)
         if (ierr/=0) then
            write(0,'(a)') "An error occurred in ccpp_finalize"
            return
         end if
      end if

      ! Deallocate shared CCPP data
      if (allocated(CCPP_shared)) deallocate(CCPP_shared)

    else

      write(0,'(2a)') 'Error, undefined IPD step ', trim(step)
      ierr = 1
      return

    end if

  end subroutine IPD_CCPP_step

end module IPD_CCPP_driver
