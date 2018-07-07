module CCPP_data

    use ccpp_api,      only: ccpp_t
    use CCPP_typedefs, only: CCPP_interstitial_type, &
                             CCPP_shared_type

    implicit none

    private

    public cdata_tile, cdata_domain, &
           cdata_block, ccpp_suite,  &
           CCPP_interstitial,        &
           CCPP_shared

    !------------------------------------------------------!
    !  CCPP container for domain (entire patch of MPI task !
    !  and individual blocks - dims are nblocks & nthreads !
    !------------------------------------------------------!
    type(ccpp_t), save, target :: cdata_tile
    type(ccpp_t), save, target :: cdata_domain
    type(ccpp_t), dimension(:,:), allocatable, save, target :: cdata_block

    !------------------------------------------------------!
    !  CCPP suite definition file                          !
    !------------------------------------------------------!
    character(len=256)      :: ccpp_suite='undefined.xml'

    type(CCPP_interstitial_type), save, target :: CCPP_interstitial
    type(CCPP_shared_type), dimension(:), allocatable, save, target :: CCPP_shared

end module CCPP_data
