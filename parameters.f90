
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module parameters 

  use toolkit , only : rp => dp , zero => cero , one => uno , two => dos , cien , mil , half => medio , tol => tolvl

  implicit none

  character(len=120)   :: path                      ! working directory   
  integer              :: threads                   ! num of threads use for parallelization

  integer  , parameter :: iterm = 5000              ! max iterations in solving value functions
  real(rp) , parameter :: tolm  = 0.000001_rp       ! tolerance level in solving value functions
  real(rp) , parameter :: ptrunc = 0.999_rp         ! quantile at which the Pareto distribution is truncated
  
  character(len=8)     :: date                      ! date   
  character(len=10)    :: time                      ! time
  
  real(rp)             :: tau                       ! corporate tax rate
  real(rp)             :: taubase = 0.30_rp         ! profit tax rate of the baseline economy (chosen in the menu)
  real(rp)             :: deltabase = 0.02_rp       ! exogenous exit probability, Hopenhayn economy (chosen in the menu)
  integer              :: noexit = 0                ! 1: no endogenous exit (Hopenhayn economy, fixed cost = 0)

  ! names of the files (see infile, expfile and outfile below)
  character(len=40)    :: ecotag                    ! economy: hopenhayn_tau30_delta02, lucas_tau30, ...
  character(len=20)    :: vname                     ! version of the model: baseline, nocapital, costs_labor
  character(len=20)    :: stamp                     ! date and time at which the program started
  real(rp) , parameter :: alpha = 0.30_rp           ! capital elasticity

  real(rp)             :: minermodel                ! min error for calibration
  real(rp)             :: wcondi                    ! equilibrium condition for wages

  ! model execution parameters
  integer :: withk
  integer :: ecostinw
  integer :: fpop

  ! technology & others
  real(rp) :: gamma
  real(rp) :: beta
  real(rp) :: eprob
  real(rp) :: fcost
  real(rp) :: ecost

  ! fixed productivity distribution
  real(rp) :: xi
  real(rp) :: uz
  real(rp) :: sigma

  ! ar(1) productivity parameters
  real(rp) :: rhoz
  real(rp) :: stdz
  real(rp) :: muz0
  real(rp) :: stdz0

  ! distortions
  real(rp) :: atau  ! average distortion
  real(rp) :: taus  ! distortion for small firms (subsidy rate)
  real(rp) :: taul  ! distortion for large firms (tax rate)
  real(rp) :: zphi  ! level of xz that splits firm distribution
  real(rp) :: phi0  ! distortion as func of z: level
  real(rp) :: phi1  ! distortion as func of z: curvature

  ! disutility of labor supply: level and curvature
  real(rp) :: eta0
  real(rp) :: eta1

  ! results
  real(rp) :: wage
  real(rp) :: rrate
  real(rp) :: M
  real(rp) :: Me
  real(rp) :: aK
  real(rp) :: aN
  real(rp) :: aL
  real(rp) :: aY
  real(rp) :: aC
  real(rp) :: aP
  real(rp) :: aT
  real(rp) :: aZ
  real(rp) :: aZ1
  real(rp) :: aZ2
  real(rp) :: TFP0    ! our definition
  real(rp) :: TFP1    ! nezih's definitoon
  real(rp) :: TFP2    ! econometrician
  real(rp) :: TFP3    ! using total population

  type :: results
    real(rp) :: wage
    real(rp) :: M
    real(rp) :: Me
    real(rp) :: aK
    real(rp) :: aN
    real(rp) :: aL
    real(rp) :: aY
    real(rp) :: aC
    real(rp) :: aP
    real(rp) :: aT
    real(rp) :: aZ
    real(rp) :: aZ1
    real(rp) :: aZ2
    real(rp) :: TFP0    ! our definition
    real(rp) :: TFP1    ! nezih's definitoon
    real(rp) :: TFP2    ! econometrician
    real(rp) :: TFP3    ! using total population
    real(rp) , allocatable :: k(:)
    real(rp) , allocatable :: n(:)
    real(rp) , allocatable :: y(:)
    real(rp) , allocatable :: pi(:)
    real(rp) , allocatable :: e(:)
    real(rp) , allocatable :: t(:)
    real(rp) , allocatable :: tau(:)
    real(rp) , allocatable :: z(:)
    real(rp) , allocatable :: G(:)
    real(rp) , allocatable :: Ge(:)
    logical  , allocatable :: small(:)
    logical  , allocatable :: large(:)
    real(rp) , allocatable :: taxed(:)
  end type results

  ! type(results) :: hugo
  ! type(results) :: lucas

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
  type(moments) :: s_e
  type(moments) :: shr_firms
  type(moments) :: shr_labor
  type(moments) :: rho_n

  abstract interface
    function func1d(x) result(f)
      import :: rp
      real(rp) , intent(in) :: x
      real(rp)              :: f
    end function func1d
  end interface

  contains

  ! ------------------------------------------------------------------------------------------------
  ! names of the files. The economy is set once in main.f90 (set_economy):
  !   ecotag = hopenhayn_tau<tau>_delta<delta> or lucas_tau<tau>, with tau and delta in %, e.g.
  !            hopenhayn_tau30_delta02, hopenhayn_tau0_delta08, lucas_tau30
  !   infile(v)     txtfiles/calib_<ecotag>_<v>.txt           parameters read by the model
  !   expfile(x)    results/<ecotag>_<x>.txt                  symmetric experiment (x = version or
  !                                                           sensitivity exercise)
  !   outfile(k,v)  results/<stamp>_<k>_<ecotag>_<v>.txt      k = calib (calibration) or solve (solution)

  subroutine set_economy(model)
    implicit none
    character(len=*) , intent(in) :: model    ! "hopenhayn" or "lucas"
    character(len=8)              :: ttxt,dtxt
    write(ttxt,'(i0)')   nint(cien*taubase)
    write(dtxt,'(i2.2)') nint(cien*deltabase)
    if (trim(model).eq."hopenhayn") then
      ecotag = "hopenhayn_tau"//trim(ttxt)//"_delta"//trim(dtxt)
    else
      ecotag = trim(model)//"_tau"//trim(ttxt)
    end if
    stamp = date(3:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//time(1:2)//'-'//time(3:4)
    return
  end subroutine set_economy

  function infile(v) result(f)
    implicit none
    character(len=*) , intent(in) :: v
    character(len=120)            :: f
    f = "txtfiles/calib_"//trim(ecotag)//"_"//trim(v)//".txt"
    return
  end function infile

  function expfile(x) result(f)
    implicit none
    character(len=*) , intent(in) :: x
    character(len=120)            :: f
    f = "results/"//trim(ecotag)//"_"//trim(x)//".txt"
    return
  end function expfile

  function outfile(k,v) result(f)
    implicit none
    character(len=*) , intent(in) :: k,v
    character(len=120)            :: f
    f = "results/"//trim(stamp)//"_"//trim(k)//"_"//trim(ecotag)//"_"//trim(v)//".txt"
    return
  end function outfile

  ! ------------------------------------------------------------------------------------------------
  ! findroot: finds x such that f(x)=0 with Brent's method (bracketing + inverse quadratic
  ! interpolation / secant steps, falling back to bisection when needed). The bracket is built
  ! around the initial guess x0 > 0 by multiplicative expansion, so x stays positive.
  ! On exit, the last call to f is always f(x): any side effects of f (policy functions,
  ! distributions, ...) correspond to the returned root.
  subroutine findroot(f,x0,xmax,x,fx,xtol,ftol)
    implicit none
    procedure(func1d)       :: f
    real(rp) , intent(in)   :: x0,xmax,xtol,ftol
    real(rp) , intent(out)  :: x,fx
    real(rp)                :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm,k,xlast
    integer                 :: it
    ! 1. bracket the root around the initial guess
    k  = 1.02_rp
    a  = x0/k           ; fa = f(a)
    b  = min(x0*k,xmax) ; fb = f(b) ; xlast = b
    it = 0
    do while (fa*fb.gt.zero .and. it.lt.40)
      it = it + 1 ; k = k*k
      if (abs(fa).lt.abs(fb)) then
        a = a/k           ; fa = f(a) ; xlast = a
      else
        b = min(b*k,xmax) ; fb = f(b) ; xlast = b
      end if
    end do
    if (fa*fb.gt.zero) then
      print * , ' warning in findroot: root not bracketed, returning the point with the smallest residual'
      if (abs(fa).lt.abs(fb)) then
        x = a
      else
        x = b
      end if
      fx = f(x) ; return
    end if
    ! 2. Brent's method
    c = b ; fc = fb ; d = b - a ; e = d
    do it = 1,200
      if ((fb.gt.zero .and. fc.gt.zero) .or. (fb.lt.zero .and. fc.lt.zero)) then
        c = a ; fc = fa ; d = b - a ; e = d
      end if
      if (abs(fc).lt.abs(fb)) then
        a = b ; b = c ; c = a ; fa = fb ; fb = fc ; fc = fa
      end if
      tol1 = two*epsilon(one)*abs(b) + half*xtol
      xm   = half*(c - b)
      if (abs(xm).le.tol1 .or. abs(fb).le.ftol) exit
      if (abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s = fb/fa
        if (a.eq.c) then
          p = two*xm*s ; q = one - s
        else
          q = fa/fc ; r = fb/fc
          p = s*(two*xm*q*(q - r) - (b - a)*(r - one))
          q = (q - one)*(r - one)*(s - one)
        end if
        if (p.gt.zero) q = -q
        p = abs(p)
        if (two*p.lt.min(3.0_rp*xm*q - abs(tol1*q),abs(e*q))) then
          e = d ; d = p/q
        else
          d = xm ; e = d
        end if
      else
        d = xm ; e = d
      end if
      a = b ; fa = fb
      if (abs(d).gt.tol1) then
        b = b + d
      else
        b = b + sign(tol1,xm)
      end if
      fb = f(b) ; xlast = b
    end do
    x = b ; fx = fb
    ! make sure the last evaluation corresponds to the returned root
    if (x.ne.xlast) fx = f(x)
    return
  end subroutine findroot

  elemental function kstar(z,taui,rrate0,wage0) result(k)
    implicit none
    real(rp) :: k
    real(rp) , intent(in) :: taui,z,rrate0,wage0
    if (withk.eq.1) then
      k = (((one - taui)*exp(z))**(one/(one-gamma)))*&
          gamma*(alpha/rrate0)*&
          ((gamma*((alpha/rrate0)**alpha)*(((one-alpha)/wage0)**(one-alpha)))**(gamma/(one-gamma)))
    else
      k = zero
    end if
    return
  end function kstar
  elemental function nstar(z,taui,rrate0,wage0) result(n)
    implicit none
    real(rp) :: n
    real(rp) , intent(in) :: taui,z,rrate0,wage0
    if (withk.eq.1) then
      n = (((one - taui)*exp(z))**(one/(one-gamma)))*&
          gamma*((one-alpha)/wage0)*&
          ((gamma*((alpha/rrate0)**alpha)*(((one-alpha)/wage0)**(one-alpha)))**(gamma/(one-gamma)))
    else
      n = (((one - taui)*exp(z))**(one/(one-gamma)))*&
          ((gamma/wage0)**(one/(one-gamma)))
    end if
    return
  end function nstar
  elemental function output(z,k,n) result(y)
    implicit none
    real(rp) :: y
    real(rp) , intent(in) :: z,k,n
    if (withk.eq.1) then
      y = exp(z)*( k**(gamma*alpha) )*( n**(gamma*(one-alpha)) )
    else
      y = exp(z)*(n**gamma)
    end if
    return
  end function output
  elemental function profits(z,taui,rrate0,wage0) result(pi)
    implicit none
    real(rp) :: y,k,n,pi
    real(rp) , intent(in) :: z,taui,rrate0,wage0
    if (withk.eq.1) then
      n  = nstar(z,taui,rrate0,wage0)
      k  = kstar(z,taui,rrate0,wage0)
      pi = (one-tau)*( (one-taui)*output(z,k,n) - wage0*n - rrate0*k )
    else
      n  = nstar(z,taui,rrate0,wage0)
      pi = (one-tau)*( (one-taui)*output(z,zero,n) - wage0*n )
    end if
    return
  end function profits

  function aggTFP(mass,distort,dist,zgrids) result (agtfp)
    use toolkit , only : varmean
    implicit none
    real(rp) :: mass,distort(:),dist(:),zgrids(:),agtfp
      agtfp = ( mass**(one-gamma) )* &
                ( varmean( (((one - distort)**gamma)*exp(zgrids))**(one/(one-gamma)) , w = dist )) / &
                ( varmean( (((one - distort)       )*exp(zgrids))**(one/(one-gamma)) , w = dist )**gamma )
    return
  end function aggTFP

  ! this subroutine prints the results of each experiment
  subroutine printexer(numfile,eco0,eco1)
    use toolkit , only : varmean
    implicit none
    type(results) , intent(in)            :: eco0
    type(results) , intent(in) , optional :: eco1
    integer       , intent(in)            :: numfile
    real(rp)                              :: T000,T100
    real(rp)                              :: T010,T001
    real(rp)                              :: T110,T101
    real(rp)                              :: T011,T111
    real(rp)                              :: dM,dphi,dF

    if (.not.present(eco1)) then
      write(numfile,'(30(f9.3))') zero,eco0%wage,eco0%TFP0,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
        eco0%M,eco0%aY,eco0%aK,eco0%aN,eco0%aL+eco0%M,eco0%aC,eco0%TFP1,eco0%TFP2,eco0%TFP3,wcondi,zero
      write(numfile,'(30(a9))') 'atau','wage','TFP','D.TFP','D.M','D.Phi','D.F',&
        'E(t)','E(t|e)','E(t|s)','E(t|l)','taxed','D.E','D.Y','D.K','D.N','D.L','D.C','D.TFP1','D.TFP2','D.TFP3','wcondi',&
        'tax_l'
      write(numfile,'(30(f9.3))') zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,&
        zero,zero,zero,zero,zero,zero,zero,zero,zero,zero,zero
      return
    end if
    
    write(numfile,'(f9.3)',advance="no") cien*atau
    write(numfile,'(f9.3)',advance="no") cien*(eco1%wage/eco0%wage - one)
    write(numfile,'(f9.3)',advance="no") eco1%TFP0
    write(numfile,'(f9.3)',advance="no") cien*( eco1%TFP0/eco0%TFP0 - one)

    T000 = aggTFP( eco0%M , eco0%tau , eco0%G , eco0%z )
    T100 = aggTFP( eco1%M , eco0%tau , eco0%G , eco0%z )
    T010 = aggTFP( eco0%M , eco1%tau , eco0%G , eco0%z )
    T001 = aggTFP( eco0%M , eco0%tau , eco1%G , eco0%z )
    T110 = aggTFP( eco1%M , eco1%tau , eco0%G , eco0%z )
    T101 = aggTFP( eco1%M , eco0%tau , eco1%G , eco0%z )
    T011 = aggTFP( eco0%M , eco1%tau , eco1%G , eco0%z )
    T111 = aggTFP( eco1%M , eco1%tau , eco1%G , eco0%z )

    dM   = (1._rp/3._rp)*(T100 - T000) + (1._rp/6._rp)*((T110 - T010) + (T101 - T001)) + (1._rp/3._rp)*(T111 - T011)
    dphi = (1._rp/3._rp)*(T010 - T000) + (1._rp/6._rp)*((T110 - T100) + (T011 - T001)) + (1._rp/3._rp)*(T111 - T101)
    dF   = (1._rp/3._rp)*(T001 - T000) + (1._rp/6._rp)*((T101 - T100) + (T011 - T010)) + (1._rp/3._rp)*(T111 - T110)

    write(numfile,'(f9.3)',advance="no") cien*dM/T000
    write(numfile,'(f9.3)',advance="no") cien*dphi/T000
    write(numfile,'(f9.3)',advance="no") cien*dF/T000

    ! E(t) is the revenue-neutrality residual: output-weighted average wedge, in % (should be ~0)
    write(numfile,'(f9.3)',advance="no") cien*varmean( eco1%tau   , w = eco1%y*eco1%G                     )
    write(numfile,'(f9.3)',advance="no") cien*varmean( eco1%tau   , w = eco1%Ge                           )
    write(numfile,'(f9.3)',advance="no") cien*varmean( eco1%tau   , w = eco1%y*eco1%G , mask = eco1%small )
    write(numfile,'(f9.3)',advance="no") cien*varmean( eco1%tau   , w = eco1%y*eco1%G , mask = eco1%large )
    write(numfile,'(f9.3)',advance="no") varmean( eco1%taxed , w = eco1%G                            )*cien

    write(numfile,'(f9.3)',advance="no") cien*(eco1%M/eco0%M - one)

    write(numfile,'(f9.3)',advance="no") cien*(eco1%aY/eco0%aY - one)
    write(numfile,'(f9.3)',advance="no") cien*(eco1%aK/eco0%aK - one)
    write(numfile,'(f9.3)',advance="no") cien*(eco1%aL/eco0%aL - one)
    write(numfile,'(f9.3)',advance="no") cien*((eco1%aL+eco1%M)/(eco0%aL+eco0%M) - one)
    write(numfile,'(f9.3)',advance="no") cien*(eco1%aC/eco0%aC - one)

    write(numfile,'(f9.3)',advance="no") cien*(eco1%TFP1/eco0%TFP1 - one)
    write(numfile,'(f9.3)',advance="no") cien*(eco1%TFP2/eco0%TFP2 - one)
    write(numfile,'(f9.3)',advance="no") cien*(eco1%TFP3/eco0%TFP3 - one)

    write(numfile,'(f9.3)',advance="no") wcondi

    ! tax rate paid by firms above the cutoff (in %): equals atau in the symmetric experiment
    write(numfile,'(f9.3)',advance="no") cien*max(zero,maxval(eco1%tau))

    write(numfile,'(a1 )',advance="yes") ' '

    return
  end subroutine printexer

  subroutine pareto(xgrid,xdist,uz0,xi0,sigma0,jj)
    use toolkit , only : grid
    implicit none
    integer  , intent(in)    :: jj
    real(rp) , intent(in)    :: uz0,xi0,sigma0
    real(rp) , intent(inout) :: xgrid(jj),xdist(jj)
    integer                  :: ii
    ! grid up to the ptrunc-quantile of the Pareto, with more points near the minimum (s=2)
    xgrid    = grid( uz0 + ((one/(one-ptrunc))**(one/xi0) - one )*sigma0*xi0  , uz0 + tol , jj , s=one )
    xdist    = zero
    xdist(1) = paretocdf( uz0,xi0,sigma0, 0.5_rp*(xgrid(2)+xgrid(1)) )
    do ii=2,jj-1
      xdist(ii) = paretocdf( uz0,xi0,sigma0, 0.5_rp*(xgrid(ii)+xgrid(ii+1)) ) - &
                  paretocdf( uz0,xi0,sigma0, 0.5_rp*(xgrid(ii-1)+xgrid(ii)) )
    end do
    ! truncated (not winsorized) Pareto: the mass above the ptrunc-quantile is dropped, not piled on the last node
    !xdist(jj) = ptrunc - paretocdf( uz0,xi0,sigma0, 0.5_rp*(xgrid(jj-1)+xgrid(jj)) )
    xdist(jj) = one - paretocdf( uz0,xi0,sigma0, 0.5_rp*(xgrid(jj-1)+xgrid(jj)) )
    xdist = xdist/sum(xdist)
    xgrid = log(xgrid)
    return
    contains
    function paretocdf(uz1,xi1,sig1,z) result(pcdf)
      implicit none
      real(rp) :: z,xi1,uz1,sig1,pcdf
      pcdf = 1.0_rp - (one + (z-uz1)/(xi1*sig1) )**(-xi1)
      return
    end function paretocdf
  end subroutine pareto

  subroutine set_data()
    implicit none

    s_lab%mean%d = 14.3747
    s_lab%std%d  = 1.171147
    s_lab%p10%d  = 1.000_rp
    s_lab%p25%d  = 2.000_rp
    s_lab%p50%d  = 3.630_rp
    s_lab%p75%d  = 8.110_rp
    s_lab%p90%d  = 19.63_rp

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
    s_lab0%std%d  = s_lab%std%d/1.0625574
    s_lab0%p10%d  = 1.000_rp
    s_lab0%p25%d  = 1.250_rp
    s_lab0%p50%d  = 2.010_rp
    s_lab0%p75%d  = 4.000_rp
    s_lab0%p90%d  = 7.250_rp

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

    s_e%d = 8.0_rp
    
    return
  end subroutine set_data

end module parameters

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%