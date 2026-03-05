
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module parameters 

  implicit none

  ! numerical parameters
  integer  , parameter :: rp   = kind(1.0d0)
  real(rp) , parameter :: tol  = dble(0.00000000001)
  real(rp) , parameter :: zero = dble(0.00000000000)
  real(rp) , parameter :: one  = dble(1.00000000000)
  real(rp) , parameter :: two  = dble(2.00000000000)
  real(rp) , parameter :: cien = dble(100.000000000)
  real(rp) , parameter :: mil  = dble(1000.00000000)
  real(rp) , parameter :: half = dble(0.50000000000)

  ! working directory
  character(len=120)   :: path                      ! working directory   
  character(len=1)     :: bb
  integer              :: solvemode
  integer              :: threads                   ! num of threads use for parallelization
  integer  , parameter :: iterm = 200               ! max iterations in solving value functions
  real(rp) , parameter :: tolm  = dble(0.000001)    ! tolerance level in solving value functions
  character(len=8)     :: date                      ! date   
  character(len=10)    :: time                      ! time
  real(rp) , parameter :: tau   = 0.3d0             ! corporate tax rate
  real(rp) , parameter :: alpha = 0.30d0

  ! ------------------------------------------------------------
  ! statistics parameters

  type :: moments
    real(rp) :: m,d   ! data
  end type moments
  type :: stats
    type(moments) :: mean,std,p10,p25,p50,p75,p90   ! data
  end type stats
  type(stats)   :: s_lab,s_lab0
  type(stats)   :: s_kap,s_kap0
  type(stats)   :: s_va,s_va0
  type(stats)   :: s_pi,s_pi0
  type(stats)   :: s_z,s_z0
  type(moments) :: shr_firms
  type(moments) :: shr_labor
  type(moments) :: rho_n

  contains

  subroutine pareto(pgrid,pp,uz,xi,sigma,np)
    use toolkit , only : grid
    implicit none
    integer  , intent(in)    :: np
    real(rp) , intent(in)    :: uz,xi,sigma
    real(rp) , intent(inout) :: pgrid(np),pp(np)
    integer                  :: iz
    !pgrid = grid( uz/((one-0.999d0)**(one/xi)) , uz + tol , np , s=one )
    pgrid = grid( uz + ((one/(one-0.999d0))**(one/xi) - one )*sigma*xi  , uz + tol , np , s=one )
    pp    = zero
    pp(1) = paretocdf( uz,xi,sigma, 0.5d0*(pgrid(2)+pgrid(1)) )
    do iz=2,np-1
      pp(iz) = paretocdf( uz,xi,sigma, 0.5d0*(pgrid(iz)+pgrid(iz+1)) ) - paretocdf( uz,xi,sigma, 0.5d0*(pgrid(iz-1)+pgrid(iz)) )
    end do
    pp(np) = one - paretocdf( uz,xi,sigma, 0.5d0*(pgrid(np-1)+pgrid(np)) )
    pp     = pp/sum(pp)
    pgrid  = log(pgrid)
    return
    contains
    function paretocdf(uz0,xi0,sig0,z) result(pcdf)
      implicit none
      real(rp) :: z,xi0,uz0,sig0,pcdf
      !pcdf = 1.0d0 - (uz0/z)**(xi0)
      pcdf = 1.0d0 - (one + (z-uz0)/(xi0*sig0) )**(-xi0)
      return
    end function paretocdf
  end subroutine pareto

  ! subroutine pareto(pgrid,pp,uz,xi,sigma)
  !   use toolkit , only : grid
  !   implicit none
  !   real(rp) , intent(in)             :: uz,xi
  !   real(rp) , intent(in) , optional  :: sigma
  !   real(rp) , intent(inout)          :: pgrid(:),pp(:)
  !   real(rp)                          :: maxz,minz,sig
  !   integer                           :: iz,np
  !   if (present(sigma)) then
  !     sig = sigma
  !   else
  !     sig = uz/xi
  !   end if
  !   np    = size(pgrid)
  !   maxz  = uz + sig*xi*( ((one/(one-0.999d0))**(one/xi)) - one)
  !   minz  = uz + tol
  !   pgrid = grid( maxz , minz , np , s=one )
  !   pp    = zero
  !   maxz  = 0.5d0*(pgrid(2)+pgrid(1))
  !   pp(1) = paretocdf(maxz)
  !   do iz=2,np-1
  !     maxz = 0.5d0*(pgrid(iz)+pgrid(iz+1))
  !     minz = 0.5d0*(pgrid(iz-1)+pgrid(iz))
  !     pp(iz) = paretocdf(maxz) - paretocdf(minz)
  !   end do
  !   minz   = 0.5d0*(pgrid(np-1)+pgrid(np))
  !   pp(np) = one - paretocdf(minz)
  !   maxz   = sum(pp)
  !   pp     = pp/maxz
  !   return
  !   contains
  !   function paretocdf(z) result(pcdf)
  !     implicit none
  !     real(rp) :: z,pcdf
  !     pcdf = 1.0d0 - ( 1.0d0 + (z-uz)/(xi*sig) )**(-xi)
  !     return
  !   end function paretocdf
  ! end subroutine pareto

  subroutine set_data()
    implicit none

    s_lab%mean%d = 14.3747
    s_lab%std%d  = 1.171147
    s_lab%p10%d  = 1.000d0
    s_lab%p25%d  = 2.000d0
    s_lab%p50%d  = 3.630d0
    s_lab%p75%d  = 8.110d0
    s_lab%p90%d  = 19.63d0

    s_kap%mean%d = 4.98197
    s_kap%std%d  = 871.459
    s_kap%p10%d  = 0.04497
    s_kap%p25%d  = 0.10985
    s_kap%p50%d  = 0.31594
    s_kap%p75%d  = 1.03291
    s_kap%p90%d  = 3.41681

    s_va%mean%d = 3.00759
    s_va%std%d  = 82.6978
    s_va%p10%d  = 0.05467
    s_va%p25%d  = 0.12341
    s_va%p50%d  = 0.31122
    s_va%p75%d  = 0.89026
    s_va%p90%d  = 2.72623

    s_pi%mean%d = 0.161466
    s_pi%std%d  = 14.04768
    s_pi%p10%d  = -0.03205
    s_pi%p25%d  = -0.00243
    s_pi%p50%d  =  0.00644
    s_pi%p75%d  =  0.03665
    s_pi%p90%d  =  0.14995

    s_lab0%mean%d = 5.12265
    s_lab0%std%d  = 105.243
    s_lab0%p10%d  = 1.000d0
    s_lab0%p25%d  = 1.250d0
    s_lab0%p50%d  = 2.010d0
    s_lab0%p75%d  = 4.000d0
    s_lab0%p90%d  = 7.250d0

    s_kap0%mean%d = 0.97300
    s_kap0%std%d  = 43.5121
    s_kap0%p10%d  = 0.01507
    s_kap0%p25%d  = 0.03283
    s_kap0%p50%d  = 0.07448
    s_kap0%p75%d  = 0.17883
    s_kap0%p90%d  = 0.46214

    s_va0%mean%d = 0.66770
    s_va0%std%d  = 40.5176
    s_va0%p10%d  = 0.02322
    s_va0%p25%d  = 0.05192
    s_va0%p50%d  = 0.11405
    s_va0%p75%d  = 0.24653
    s_va0%p90%d  = 0.52190

    s_pi0%mean%d = 0.023391
    s_pi0%std%d  = 2.478409
    s_pi0%p10%d  = -0.03578
    s_pi0%p25%d  = -0.01105
    s_pi0%p50%d  =  0.00047
    s_pi0%p75%d  =  0.01149
    s_pi0%p90%d  =  0.04231

    shr_firms%d = 3.26
    shr_labor%d = 54.8312

    rho_n%d = 0.9631323
    
    return
  end subroutine set_data

end module parameters

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%