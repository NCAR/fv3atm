module CCPP_data

    use ccpp_api, only: ccpp_t

    implicit none

    private

    public cdata_domain, cdata_block

    !------------------------------------------------------!
    !  CCPP container for domain (entire patch of MPI task !
    !  and individual blocks - dims are nblocks & nthreads !
    !------------------------------------------------------!
    type(ccpp_t), save, target :: cdata_domain
    type(ccpp_t), dimension(:,:), allocatable, save, target :: cdata_block

end module CCPP_data
