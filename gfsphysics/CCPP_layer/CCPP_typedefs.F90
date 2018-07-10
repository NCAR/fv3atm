module CCPP_typedefs

    implicit none

!! \section arg_table_CCPP_typedefs
!! | local_name                                         | standard_name                                                 | long_name                                                                | units   | rank | type                   |    kind   | intent | optional |
!! |----------------------------------------------------|---------------------------------------------------------------|--------------------------------------------------------------------------|---------|------|------------------------|-----------|--------|----------|
!! | CCPP_interstitial                                  | CCPP_Interstitial_type                                        | derived type CCPP_interstitial_type                                      | DDT     |    0 | CCPP_interstitial_type |           | none   | F        |
!!

    private

    public CCPP_shared_type, CCPP_interstitial_type

!! \section arg_table_CCPP_shared_type
!! | local_name                                         | standard_name                                                 | long_name                                                                           | units   | rank | type        |    kind   | intent | optional |
!! |----------------------------------------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------------|---------|------|-------------|-----------|--------|----------|
!! | CCPP_shared(nt)%errflg                             | error_flag                                                    | error flag for error handling in CCPP                                               | flag    |    0 | integer     |           | none   | F        |
!! | CCPP_shared(nt)%errmsg                             | error_message                                                 | error message for error handling in CCPP                                            | none    |    0 | character   | len=512   | none   | F        |
!! | CCPP_shared(nt)%nthreads                           | omp_threads                                                   | number of OpenMP threads available for fast physics schemes                         | count   |    0 | integer     |           | none   | F        |
!!
  type CCPP_shared_type

     integer                             :: errflg
     character(len=512)                  :: errmsg
     integer                             :: nthreads

  contains

    procedure :: create  => shared_create     !<   allocate/set data
    procedure :: reset   => shared_reset      !<   reset data
    procedure :: mprint  => shared_print      !<   print data

  end type CCPP_shared_type

!! \section arg_table_CCPP_interstitial_type
!! | local_name                                         | standard_name                                                 | long_name                                                                           | units   | rank | type        |    kind   | intent | optional |
!! |----------------------------------------------------|---------------------------------------------------------------|-------------------------------------------------------------------------------------|---------|------|-------------|-----------|--------|----------|
!! | CCPP_interstitial%akap                             | kappa_dry_for_fast_physics                                    | modified kappa for fast physics                                                     | none    |    0 | real        |           | none   | F        |
!! | CCPP_interstitial%bdt                              |                                                               | large time step for dynamics                                                        | s       |    0 | real        |           | none   | F        |
!! | CCPP_interstitial%cappa                            | cappa_moist_gas_constant_at_Lagrangian_surface                | cappa(i,j,k) = rdgas / ( rdgas +  cvm(i)/(1.+r_vir*q(i,j,k,sphum)) )                | none    |    3 | real        |           | none   | F        |
!! | CCPP_interstitial%dtdt                             | tendency_of_air_temperature_at_Lagrangian_surface             | air temperature tendency due to fast physics at Lagrangian surface                  | K s-1   |    3 | real        |           | none   | F        |
!! | CCPP_interstitial%do_qa                            | flag_for_inline_cloud_fraction_calculation                    | flag for the inline cloud fraction calculation                                      | flag    |    0 | logical     |           | none   | F        |
!! | CCPP_interstitial%fast_mp_consv                    | flag_for_fast_microphysics_energy_conservation                | flag for fast microphysics energy conservation                                      | flag    |    0 | logical     |           | none   | F        |
!! | CCPP_interstitial%kmp                              | top_layer_index_for_fast_physics                              | top_layer_inder_for_gfdl_mp                                                         | index   |    0 | integer     |           | none   | F        |
!! | CCPP_interstitial%last_step                        | flag_for_the_last_step_of_k_split_remapping                   | flag for the last step of k-split remapping                                         | flag    |    0 | logical     |           | none   | F        |
!! | CCPP_interstitial%mdt                              | time_step_for_remapping_for_fast_physics                      | remapping time step                                                                 | s       |    0 | real        |           | none   | F        |
!! | CCPP_interstitial%npzdelz                          | vertical_dimension_for_thickness_at_Lagrangian_surface        | vertical dimension for thickness at Lagrangian surface                              | count   |    0 | integer     |           | none   | F        |
!! | CCPP_interstitial%out_dt                           | flag_for_tendency_of_air_temperature_at_Lagrangian_surface    | flag for calculating tendency of air temperature due to fast physics                | flag    |    0 | logical     |           | none   | F        |
!! | CCPP_interstitial%te0_2d                           | atmosphere_energy_content_in_column                           | atmosphere total energy in columns                                                  | J m-2   |    2 | real        |           | none   | F        |
!! | CCPP_interstitial%te0                              | atmosphere_energy_content_at_Lagrangian_surface               | atmosphere total energy at Lagrangian surface                                       | J m-2   |    3 | real        |           | none   | F        |
!! | CCPP_interstitial%zvir                             | ratio_of_vapor_to_dry_air_gas_constants_minus_one_default_kind| zvir=rv/rd-1.0                                                                      | none    |    0 | real        |           | none   | F        |
!!
  type CCPP_interstitial_type

     real                                :: akap
     real                                :: bdt
     real, pointer                       :: cappa(:,:,:)
     logical                             :: do_qa
     real, pointer                       :: dtdt(:,:,:)
     logical                             :: fast_mp_consv
     integer                             :: kmp
     logical                             :: last_step
     real                                :: mdt
     integer                             :: npzdelz
     logical                             :: out_dt
     real, pointer                       :: pfull(:)
     real, pointer                       :: te0_2d(:,:) ! called te_2d in fv_dynamics, te0_2d in Lagrangian_to_Eulerian, te0_2d in fv_sat_adj
     real, pointer                       :: te0(:,:,:)  ! called dp1 in fv_dynamics, te in Lagrangian_to_Eulerian, te0 in fv_sat_adj
     real                                :: zvir

  contains

    procedure :: create  => interstitial_create     !<   allocate array data
    procedure :: reset   => interstitial_reset      !<   reset array data
    procedure :: mprint  => interstitial_print      !<   print array data

  end type CCPP_interstitial_type

contains

!-----------------------------
! CCPP_shared_type
!-----------------------------
  subroutine shared_create (Shared)
    !
    implicit none
    !
    class(CCPP_shared_type) :: Shared
    !
    ! Number of OpenMP threads available for schemes, default only one
    Shared%nthreads = 1
    !
    call Shared%reset()
    !
  end subroutine shared_create

  subroutine shared_reset (Shared)
    !
    implicit none
    !
    class(CCPP_shared_type) :: Shared
    !
    Shared%errmsg        = ''
    Shared%errflg        = 0
    !
  end subroutine shared_reset

  subroutine shared_print(Shared)
    !
    implicit none
    !
    class(CCPP_shared_type) :: Shared
    !
    write (0,'(a)') 'Shared_print'
    write (0,*) 'Shared%errflg            = ', Shared%errflg
    write (0,*) 'Shared%errmsg            = ', trim(Shared%errmsg)
    write (0,*) 'Shared%nthreads          = ', Shared%nthreads
    write (0,*) 'Shared_print: end'
    !
  end subroutine shared_print

!-----------------------------
! CCPP_interstitial_type
!-----------------------------
  subroutine interstitial_create (Interstitial, is, ie, isd, ied, js, je, jsd, jed, npz,  &
                                  dt_atmos, p_split, k_split, zvir, p_ref, ak, bk, do_qa, &
                                  kappa, hydrostatic)
    !
    implicit none
    !
    class(CCPP_interstitial_type) :: Interstitial
    integer, intent(in) :: is
    integer, intent(in) :: ie
    integer, intent(in) :: isd
    integer, intent(in) :: ied
    integer, intent(in) :: js
    integer, intent(in) :: je
    integer, intent(in) :: jsd
    integer, intent(in) :: jed
    integer, intent(in) :: npz
    real,    intent(in) :: dt_atmos
    integer, intent(in) :: p_split
    integer, intent(in) :: k_split
    real,    intent(in) :: zvir
    real,    intent(in) :: p_ref
    real,    intent(in) :: ak(:)
    real,    intent(in) :: bk(:)
    logical, intent(in) :: do_qa
    real,    intent(in) :: kappa
    logical, intent(in) :: hydrostatic
    !
#ifdef MOIST_CAPPA
    allocate (Interstitial%cappa  (isd:ied, jsd:jed, 1:npz) )
#else
    allocate (Interstitial%cappa  (isd:isd, jsd:jsd, 1)     )
#endif
    allocate (Interstitial%dtdt   (is:ie, js:je, 1:npz)     )
    allocate (Interstitial%pfull  (1:npz)                   )
    allocate (Interstitial%te0_2d (is:ie, js:je)            )
    allocate (Interstitial%te0    (isd:ied, jsd:jed, 1:npz) )
    !
    ! Initialize variables to default values
#ifdef SW_DYNAMICS
    Interstitial%akap = 1.
#else
    Interstitial%akap = kappa
#endif
    Interstitial%bdt       = dt_atmos/real(abs(p_split))
    Interstitial%do_qa     = do_qa
    Interstitial%mdt       = Interstitial%bdt/real(k_split)
    if (hydrostatic) then
       Interstitial%npzdelz = 1
    else
       Interstitial%npzdelz = npz
    end if
    Interstitial%zvir      = zvir
    !
    ! Calculate vertical pressure levels
    call interstitital_calculate_pressure_levels(Interstitial, npz, p_ref, ak, bk)
    !
    ! Reset all other variables
    call Interstitial%reset()
    !
  end subroutine interstitial_create

  subroutine interstitital_calculate_pressure_levels(Interstitial, npz, p_ref, ak, bk)

     implicit none

     class(CCPP_interstitial_type)       :: Interstitial
     integer, intent(in) :: npz
     real, intent(in) :: p_ref
     real, intent(in) :: ak(:)
     real, intent(in) :: bk(:)

     real :: ph1
     real :: ph2
     integer :: k

#ifdef SW_DYNAMICS
      Interstitial%pfull(1) = 0.5*p_ref
#else
      do k=1,npz
         ph1 = ak(k  ) + bk(k  )*p_ref
         ph2 = ak(k+1) + bk(k+1)*p_ref
         Interstitial%pfull(k) = (ph2 - ph1) / log(ph2/ph1)
      enddo
#endif
      ! DH* This is copied from fv_mapz.F90, does it work with SW_DYNAMICS?
      do k=1,npz
         Interstitial%kmp = k
         if ( Interstitial%pfull(k) > 10.E2 ) exit
      enddo
  end subroutine interstitital_calculate_pressure_levels

  subroutine interstitial_reset (Interstitial)
    !
    implicit none
    !
    class(CCPP_interstitial_type) :: Interstitial
    !
    Interstitial%cappa         = 0.0
    Interstitial%dtdt          = 0.0
    Interstitial%fast_mp_consv = .false.
    Interstitial%last_step     = .false.
    Interstitial%out_dt        = .false.
    Interstitial%te0_2d        = 0.0
    Interstitial%te0           = 0.0
    !
  end subroutine interstitial_reset

  subroutine interstitial_print(Interstitial)
    !
    implicit none
    !
    class(CCPP_interstitial_type) :: Interstitial
    !
    ! Print static variables
    write (0,'(a)') 'Interstitial_print'
    write (0,*) 'Interstitial_print: values that do not change'
    write (0,*) 'Interstitial%akap              = ', Interstitial%akap
    write (0,*) 'Interstitial%bdt               = ', Interstitial%bdt
    write (0,*) 'Interstitial%kmp               = ', Interstitial%kmp
    write (0,*) 'Interstitial%mdt               = ', Interstitial%mdt
    write (0,*) 'sum(Interstitial%pfull       ) = ', sum(Interstitial%pfull      )
    write (0,*) 'Interstitial%zvir              = ', Interstitial%zvir
    ! Print all other variables
    write (0,*) 'Interstitial_print: values that change'
    write (0,*) 'sum(Interstitial%cappa       ) = ', sum(Interstitial%cappa       )
    write (0,*) 'Interstitial%do_qa             = ', Interstitial%do_qa
    write (0,*) 'sum(Interstitial%dtdt        ) = ', sum(Interstitial%dtdt        )
    write (0,*) 'Interstitial%fast_mp_consv     = ', Interstitial%fast_mp_consv
    write (0,*) 'Interstitial%last_step         = ', Interstitial%last_step
    write (0,*) 'Interstitial%out_dt            = ', Interstitial%out_dt
    write (0,*) 'sum(Interstitial%te0_2d      ) = ', sum(Interstitial%te0_2d      )
    write (0,*) 'sum(Interstitial%te0         ) = ', sum(Interstitial%te0         )
    write (0,*) 'Interstitial_print: end'
    !
  end subroutine interstitial_print

end module CCPP_typedefs

