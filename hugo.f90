
! ======================================================================================================
! Hopenhayn (1992) economy with AR(1) productivity
!
!   log z' = -0.5*stdz^2 + rhoz*log z + stdz*eps ,  eps ~ N(0,1)
!
! Entrants draw their initial log productivity from N(muz0, stdz0^2). Firms exit exogenously with
! probability delta (eprob) and endogenously when the expected value of staying is negative.
!
! Calibration (delta is not calibrated: it is chosen in the menu of main.f90):
!   gamma  -> average firm size
!   stdz0  -> share of employment in firms with 50+ workers
!   stdz   -> std. dev. of log employment
!   muz0   -> average employment of entrants
!   wage   -> P25 of employment (the entry cost, ecost, is implied by free entry)
!   fcost  -> exit rate (not calibrated, and set to zero, when there is no endogenous exit)
!   rrate  -> capital-to-output ratio (taken from the baseline calibration when withk = 0)
!   rhoz   -> persistence of employment (outer loop)
!
! Files (names built in parameters.f90 from the economy chosen in the menu, e.g. hopenhayn_tau30_delta02):
!   parameters:   txtfiles/calib_<economy>_<version>.txt
!   calibration:  results/<date>_calib_<economy>_<version>.txt
!   solution:     results/<date>_solve_<economy>_<version>.txt
!   experiments:  results/<economy>_<version>.txt and results/<economy>_<sensitivity exercise>.txt
! ======================================================================================================

module hugo

  use parameters
  use toolkit , only : lmmin,simplex,vect,varmean,varstd,percentile,grid,&
  interpolate,print_err,bisection,olsreg,num2text,error

  implicit none

  private
  public :: solve_hugo

  ! state space
  integer  , parameter :: nz = 1000

  ! productivity variables
  real(rp) , dimension(nz)    :: zgrid    ! log productivity
  real(rp) , dimension(nz)    :: zlo,zhi  ! lower/upper bounds of the cell around each node (midpoints)
  real(rp) , dimension(nz,nz) :: pz       ! transition matrix
  real(rp) , dimension(nz)    :: pz0      ! distribution of entrants

  ! distributions
  real(rp) , dimension(nz)    :: Gz       ! distribution over firms' productivity
  real(rp) , dimension(nz)    :: Gze      ! distribution over firms' productivity: entrants

  ! value and policy functions
  real(rp) , dimension(nz)    :: vfunc
  real(rp) , dimension(nz)    :: vgues
  real(rp) , dimension(nz)    :: pol_k
  real(rp) , dimension(nz)    :: pol_n
  real(rp) , dimension(nz)    :: pol_y
  real(rp) , dimension(nz)    :: pol_pi
  real(rp) , dimension(nz)    :: pol_e
  real(rp) , dimension(nz)    :: pol_t
  real(rp) , dimension(nz)    :: pol_tau

  real(rp)      :: wprev = one      ! wage of the last equilibrium solved (initial guess for findroot)
  real(rp)      :: wfix             ! wage fixed during the calibration (ecost is then implied by free entry)
  logical       :: fixcf = .false.  ! fixed cost set to zero and not calibrated (no endogenous exit)

  type(results) :: base
  type(results) :: eff
  type(results) :: distort


  contains

  subroutine solve_hugo( mode )
    !
    implicit none
    !
    integer , intent(in) :: mode
    integer              :: ifile
    real(rp)             :: par0

    ! --------------------------------------------------------------
    ! initialize some parameter values

    zphi = -99.9_rp ; atau = zero ; ifile = 879 ; minermodel = huge(one) ; vgues = zero

    ! --------------------------------------------------------------
    ! options (menu of main.f90):
    !    1  solve the baseline model                 2  calibrate the baseline model
    !    3  symmetric experiment (baseline)          4  sensitivity of the symmetric experiment
    !    5  solve the model without capital          6  calibrate the model without capital
    !    7  solve the model with costs in labor units
    !    8  calibrate the model with costs in labor units
    ! Each version reads its parameters from txtfiles/calib_<economy>_<version>.txt (see setversion).
    ! When the model without capital is calibrated, the interest rate (and so the discount factor) is
    ! taken from the baseline file of the same economy.

    if (mode.ne.1 .and. mode.ne.5 .and. mode.ne.7) print_err = .false.

    select case (mode)
    case (1) ; call setversion(1) ; call equilibrium( ) ; call writehugo(1)
    case (2) ; call setversion(1) ; call calibmodel( )  ; call writehugo(1)
    case (3) ; call setversion(1) ; call runexperiment( )
    case (4)
      call setversion(1)
      write(*,'(/,a,/)') ' Higher gamma'
      par0 = gamma ; gamma = par0*1.05_rp ; call runexperiment("gamma_up5") ; gamma = par0
      write(*,'(/,a,/)') ' Lower gamma'
      par0 = gamma ; gamma = par0*0.95_rp ; call runexperiment("gamma_down5") ; gamma = par0
      write(*,'(/,a,/)') ' Higher stdz0'
      par0 = stdz0 ; stdz0 = par0*1.10_rp ; call runexperiment("sigma0_up10") ; stdz0 = par0
      write(*,'(/,a,/)') ' Lower stdz0'
      par0 = stdz0 ; stdz0 = par0*0.90_rp ; call runexperiment("sigma0_down10") ; stdz0 = par0
      write(*,'(/,a,/)') ' Frisch elasticity = 0.01'
      par0 = eta1 ; eta1 = one/0.01_rp ; call runexperiment("frisch_0-01") ; eta1 = par0
      write(*,'(/,a,/)') ' Frisch elasticity = 3.0'
      par0 = eta1 ; eta1 = one/3.0_rp ; call runexperiment("frisch_3") ; eta1 = par0
      write(*,'(/,a,/)') ' No capital'
      call setversion(2) ; call runexperiment( )
      write(*,'(/,a,/)') ' Entry/fixed costs in labor units'
      call setversion(3) ; call runexperiment( )
      write(*,'(/,a,/)') ' Interest rate of 3%'
      call setversion(1) ; rrate = 0.03_rp ; vname = "rrate_3" ; call calibmodel( ) ; call runexperiment("rrate_3")
      write(*,'(/,a,/)') ' Interest rate of 5%'
      call setversion(1) ; rrate = 0.05_rp ; vname = "rrate_5" ; call calibmodel( ) ; call runexperiment("rrate_5")
    case (5) ; call setversion(2) ; call equilibrium( ) ; call writehugo(1)
    case (6) ; call setversion(2) ; rrate = baserate( ) ; call calibmodel( ) ; call writehugo(1)
    case (7) ; call setversion(3) ; call equilibrium( ) ; call writehugo(1)
    case (8) ; call setversion(3) ; call calibmodel( )  ; call writehugo(1)
    case default
      write(*,'(a)') ' invalid option for this model'
    end select

    print_err = .true.

    return

    contains

    ! this subroutine runs the symmetric experiment. Results: results/<economy>_<exercise>.txt, where
    ! the exercise is the version of the model (baseline, nocapital, costs_labor) or, if given, the
    ! name of the sensitivity exercise
    subroutine runexperiment(exercise)
      implicit none
      character (len=*) , intent(in) , optional :: exercise

      ! setup file for storing results
      if (present(exercise)) then
        open(unit=ifile,file=expfile(exercise),action='write',status="replace")
        write(*,'(a)') ' results: '//trim(expfile(exercise))
      else
        open(unit=ifile,file=expfile(vname),action='write',status="replace")
        write(*,'(a)') ' results: '//trim(expfile(vname))
      end if

      ! initialize the economy
      atau = zero ; zphi = -99.9_rp ; print * , '   '

      ! baseline economy
      tau = taubase ; call equilibrium(base)

      ! efficient allocation
      tau = 0.0_rp ; call equilibrium(eff)

      ! print initial results
      call printexer(ifile,base,eff) ; call printexer(0,base,eff)
      call printexer(ifile,base)     ; call printexer(0,base)

      ! distorted economies
      tau = taubase ; atau = zero ;
      exper : do ; atau = atau + 0.01_rp
       call equilibrium(distort) ; call printexer(ifile,base,distort) ; call printexer(0,base,distort)
       if (atau.gt.0.20 .and. distort%TFP/base%TFP.lt.0.95_rp) exit exper
       if (atau.gt.0.50) exit exper
      end do exper

      close(ifile)

      return
    end subroutine runexperiment

    !
  end subroutine solve_hugo

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! version of the model: 1 baseline, 2 no capital, 3 entry/fixed costs in labor units.
  ! Sets withk, ecostinw and vname, reads the parameters from txtfiles/calib_<economy>_<version>.txt
  ! and sets delta (and the fixed cost when there is no endogenous exit) as chosen in the menu

  subroutine setversion(iv)
    implicit none
    integer , intent(in) :: iv
    logical              :: exists
    withk = 1 ; ecostinw = 0 ; fpop = 0 ; vname = "baseline"
    if (iv.eq.2) then ; withk    = 0 ; vname = "nocapital"   ; end if
    if (iv.eq.3) then ; ecostinw = 1 ; vname = "costs_labor" ; end if
    inquire(file=infile(vname),exist=exists)
    if (.not.exists) then
      write(*,'(a)') ' parameter file not found: '//trim(infile(vname)) ; stop
    end if
    open(unit=9,file=infile(vname),action="read")
      read(9,*) rrate
      read(9,*) depre
      read(9,*) eprob
      read(9,*) fcost
      read(9,*) ecost
      read(9,*) xi
      read(9,*) uz
      read(9,*) rhoz
      read(9,*) stdz
      read(9,*) muz0
      read(9,*) stdz0
      read(9,*) gamma
      read(9,*) eta0
      read(9,*) eta1
    close(9)
    vgues = zero ; zphi = -99.9_rp ; atau = zero
    eprob = deltabase ; fixcf = (noexit.eq.1)
    if (fixcf) fcost = zero
    return
  end subroutine setversion

  ! interest rate of the baseline model (used to calibrate the model without capital)
  function baserate() result(r)
    implicit none
    real(rp) :: r
    open(unit=9,file=infile("baseline"),action="read")
      read(9,*) r
    close(9)
    return
  end function baserate

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! solve the model

  subroutine equilibrium(res,fixw)
    use omp_lib
    use toolkit , only : grid,tauchen,varmean,varstd,percentile
    implicit none
    type(results) , intent(out) , optional :: res
    integer                                :: iz,jz,iter
    real(rp)                               :: equilw
    real(rp)                               :: maxz
    real(rp)                               :: minz
    real(rp)                               :: equilz
    integer                                :: iterz
    real(rp)                               :: test
    real(rp)                               :: zmu_u,zsd_u,zlo0,zhi0
    logical                                :: fullbr,fixw0
    real(rp)                               :: sw,sx,sy,sxx,sxy,wgt,xx,yy,ww
    real(rp) , allocatable                 :: aux(:,:)
    real(rp) , dimension(nz)               :: Gpz1,Gpz0
    logical , intent(in) , optional :: fixw

    fixw0 = .false. ; if (present(fixw)) fixw0 = fixw

    iterz = 0 ; minz = zero ; maxz = zero ; equilz = zero ; wcondi = 99.9_rp

    ! ------------------------------------------------------------
    ! set productivity grid and process

    ! grid: +-5 unconditional std of the AR(1) and +-5 std of the entrants' distribution
    ! (all heterogeneity comes from this process, so the grid also has to cover the right tail)
    zmu_u = -half*stdz*stdz/(one-rhoz)
    zsd_u = stdz/sqrt(one-rhoz*rhoz)
    zgrid = grid( max(zmu_u+5.0_rp*zsd_u,muz0+5.0_rp*stdz0) , min(zmu_u-5.0_rp*zsd_u,muz0-5.0_rp*stdz0) , nz , s=one )

    ! distribution of entrants: N(muz0,stdz0^2) discretized on the grid
    allocate(aux(nz,nz))
    call tauchen( zgrid , zero , muz0 , max(stdz0,1.0e-6_rp) , nz , aux ) ; pz0 = aux(1,:)
    deallocate(aux)

    ! transition matrix for productivity
    call tauchen( zgrid , rhoz , -half*stdz*stdz , stdz , nz , pz )

    ! cells around each node (the grid is increasing): used to make the distortion cutoff and the
    ! exit threshold move continuously between nodes
    zlo(1)  = zgrid(1)  - half*(zgrid(2)-zgrid(1))
    zhi(nz) = zgrid(nz) + half*(zgrid(nz)-zgrid(nz-1))
    zlo(2:nz)   = half*(zgrid(1:nz-1)+zgrid(2:nz))
    zhi(1:nz-1) = zlo(2:nz)

    Gze = pz0

    ! ------------------------------------------------------------
    ! find the cutoff productivity for distortions

    fullbr = .true.
    if (abs(atau).gt.tol) then
      if (zphi.gt.minval(zgrid)) then
        minz = log(exp(zphi)*0.9_rp)
        maxz = log(exp(zphi)*1.1_rp)
        fullbr = .false.
      else
        minz = minval(zgrid)
        maxz = maxval(zgrid)
      end if
    else
      minz = maxval(zgrid)
      maxz = maxval(zgrid)
    end if
    equilz = huge(one) ; zlo0 = minz ; zhi0 = maxz

    10 continue

    ! ------------------------------------------------------------
    ! distortions

    if (abs(atau).gt.tol) then
      iterz = iterz + 1
      zphi  = half*( maxz + minz )
      ! nodes whose cell lies above zphi are taxed and those below are subsidized. The node whose
      ! cell contains zphi gets the share of its cell above zphi taxed (omega) and the rest subsidized,
      ! applied as the average rate atau*(2*omega-1). Revenue is then continuous in zphi and the
      ! bisection reaches exact neutrality instead of stopping between two nodes.
      do iz=1,nz
        if (zphi.le.zlo(iz)) then
          pol_tau(iz) =   atau
        else if (zphi.ge.zhi(iz)) then
          pol_tau(iz) = - atau
        else
          pol_tau(iz) = atau*( two*(zhi(iz)-zphi)/(zhi(iz)-zlo(iz)) - one )
        end if
      end do
    else
      pol_tau = zero
    end if

    ! ------------------------------------------------------------------------------------------------------------------------
    ! find the equilibrium wage: Brent's method, starting from the wage of the previous equilibrium (wprev)

    if (fixw0) then
      wage   = wfix
      equilw = wres_hugo(wage)
      if (ecostinw.eq.0) ecost = sum(vfunc*Gze)
      if (ecostinw.eq.1) ecost = sum(vfunc*Gze)/wage
      wprev  = wage ; wcondi = zero
    else
      call findroot(wres_hugo,wprev,1.0e4_rp,wage,equilw,tol/cien,tol/cien) ; wprev = wage ; wcondi = equilw
    end if

    ! ------------------------------------------------------------------------------------------------------------------------
    ! find the equilibrium mass of firms

    ! stationary distribution for unit mass of entrants
    M = one ; Gpz1 = zero
    do iter = 1,10000 ; Gpz0 = Gpz1
      !Gpz1 = (one-pol_e)*matmul(Gpz0,pz) + Gze       ! = transpose(pz)*Gpz0, without an nz x nz temporary
      Gpz1 = matmul((one-pol_e)*Gpz0,pz) + Gze
      test = maxval(abs(Gpz1 - Gpz0))
      if (test.lt.tol) exit
    end do

    ! assume unit mass of firms (M=1)
    Me = M/sum(Gpz0) ; Gz = M*Gpz0/sum(Gpz0)

    ! find aggregate labor demand and aggregate consumption under M=1
    if (ecostinw.eq.0) then
      aL = sum(Gz*pol_n)
      aC = wage*aL + sum(Gz*pol_pi) + sum(Gz*pol_t) - Me*ecost
    else
      aL = sum(Gz*pol_n) + Me*ecost + fcost
      aC = wage*aL + sum(Gz*pol_pi) + sum(Gz*pol_t)
    end if

    ! mass of firms that clears labor market
    if (fpop.eq.0) then
      M = (((wage/(eta0*aC))**(one/eta1))/aL)**(eta1/(1._rp+eta1))
    else
      M = one / ( one + aL )
    end if

    ! corresponding mass of entrants
    Me = M/sum(Gpz0) ; Gz = M*Gpz0/sum(Gpz0)

    ! ------------------------------------------------------------
    ! update distortions

    if (abs(atau).gt.tol) then
      equilz = sum(pol_tau*pol_y*Gz)*cien/sum(Gz)
      if (iterz.lt.1000 .and. abs(maxz-minz).gt.tol .and. abs(equilz).gt.tol) then
        if (equilz.gt.zero) minz = zphi
        if (equilz.lt.zero) maxz = zphi
        goto 10
      end if
      if (.not.fullbr .and. (abs(zphi-zlo0).lt.1.0e-6_rp .or. abs(zphi-zhi0).lt.1.0e-6_rp)) then
        fullbr = .true. ; minz = minval(zgrid) ; maxz = maxval(zgrid) ; zlo0 = minz ; zhi0 = maxz ; iterz = 0
        goto 10
      end if
    end if

    ! ------------------------------------------------------------------------------------------------------------------------
    ! aggregates

    M  = sum(Gz)
    aK = sum(Gz*pol_k)
    aN = sum(Gz*pol_n)
    aY = sum(Gz*pol_y)
    aP = sum(Gz*pol_pi)
    aT = sum(Gz*pol_t)

    if (ecostinw.eq.0) then
      aL = aN
      aC = wage*aL + aP + aT - Me*ecost
    else
      aL = aN + Me*ecost + M*fcost
      aC = wage*aL + aP + aT
    end if

    ! ------------------------------------------------------------------------------------------------------------------------
    ! productivity measures

    ! aZ1  = varmean( (((one - pol_tau)**gamma)*exp(zgrid))**(one/(one-gamma)) , w = Gz )
    ! aZ2  = varmean( (((one - pol_tau)       )*exp(zgrid))**(one/(one-gamma)) , w = Gz )**gamma
    ! TFP  = ( M**(one-gamma) )*aZ1/aZ2

    ! main TFP definition (eq. 8)
    TFP   = aY / ( ((aK**alpha) * ( aL**(one-alpha) )) ** gamma ) 

    ! workers + entrepreneurs as labor input
    TFPL  = aY / ( ((aK**alpha) * ( (aL+M)**(one-alpha) )) ** gamma )

    ! output net of entre and fixed costs
    if (ecostinw.eq.0) TFPC  = (aY - Me*ecost - M*fcost ) / ( ((aK**alpha) * ( aL**(one-alpha) )) ** gamma )
    if (ecostinw.eq.1) TFPC  = (aY - Me*ecost*wage - M*fcost*wage ) / ( ((aK**alpha) * ( aL**(one-alpha) )) ** gamma )

    ! output net of entry and fixed costs and with workers+entrepreneurs as input of production
    if (ecostinw.eq.0) TFPLC = (aY - Me*ecost - M*fcost ) / ( ((aK**alpha) * ( (aL+M)**(one-alpha) )) ** gamma ) 
    if (ecostinw.eq.1) TFPLC = (aY - Me*ecost*wage - M*fcost*wage ) / ( ((aK**alpha) * ( (aL+M)**(one-alpha) )) ** gamma )

    ! ------------------------------------------------------------------------------------------------------------------------
    ! store results

    if (present(res)) then

      res%wage  = wage
      res%M     = M
      res%Me    = Me
      res%aK    = aK
      res%aN    = aN
      res%aL    = aL
      res%aY    = aY
      res%aC    = aC
      res%aP    = aP
      res%aT    = aT
      res%TFP   = TFP
      res%TFPL  = TFPL
      res%TFPC  = TFPC
      res%TFPLC = TFPLC

      if (allocated(res%k)   ) deallocate(res%k)   ; allocate (res%k(nz))
      if (allocated(res%n)   ) deallocate(res%n)   ; allocate (res%n(nz))
      if (allocated(res%y)   ) deallocate(res%y)   ; allocate (res%y(nz))
      if (allocated(res%pi)  ) deallocate(res%pi)  ; allocate (res%pi(nz))
      if (allocated(res%e)   ) deallocate(res%e)   ; allocate (res%e(nz))
      if (allocated(res%t)   ) deallocate(res%t)   ; allocate (res%t(nz))
      if (allocated(res%tau) ) deallocate(res%tau) ; allocate (res%tau(nz))
      if (allocated(res%z)   ) deallocate(res%z)   ; allocate (res%z(nz))
      if (allocated(res%G)   ) deallocate(res%G)   ; allocate (res%G(nz))
      if (allocated(res%Ge)  ) deallocate(res%Ge)  ; allocate (res%Ge(nz))

      if (allocated(res%large) ) deallocate(res%large)  ; allocate (res%large(nz))
      if (allocated(res%small) ) deallocate(res%small)  ; allocate (res%small(nz))
      if (allocated(res%taxed) ) deallocate(res%taxed)  ; allocate (res%taxed(nz))

      res%k   = pol_k
      res%n   = pol_n
      res%y   = pol_y
      res%pi  = pol_pi
      res%e   = pol_e
      res%t   = pol_t
      res%tau = pol_tau
      res%z   = zgrid
      res%G   = Gz
      res%Ge  = Gze

      res%small = ( pol_n.lt.percentile(pol_n,0.10_rp,w=Gz) )
      res%large = ( pol_n.gt.percentile(pol_n,0.90_rp,w=Gz) )

      res%taxed = zero ; where (pol_tau.gt.tol) res%taxed = one

    end if

    ! ------------------------------------------------------------------------------------------------------------------------
    ! statistics

    s_lab%mean%m = varmean(    pol_n           , w = Gz )
    s_lab%std%m  = varstd(     log(pol_n)      , w = Gz )
    s_lab%p10%m  = percentile( pol_n , 0.10_rp , w = Gz )
    s_lab%p25%m  = percentile( pol_n , 0.25_rp , w = Gz )
    s_lab%p50%m  = percentile( pol_n , 0.50_rp , w = Gz )
    s_lab%p75%m  = percentile( pol_n , 0.75_rp , w = Gz )
    s_lab%p90%m  = percentile( pol_n , 0.90_rp , w = Gz )

    shr_firms%m  = cien*sum(Gz      ,mask=pol_n.ge.50._rp)/sum(Gz)
    shr_labor%m  = cien*sum(pol_n*Gz,mask=pol_n.ge.50._rp)/sum(pol_n*Gz)

    ! entrants: firms in their first period that survive to the next one
    s_lab0%mean%m = varmean(    pol_n            , w = Gze*(one-pol_e) )
    s_lab0%std%m  = varstd(     log(pol_n)       , w = Gze*(one-pol_e) )
    s_lab0%p10%m  = percentile( pol_n  , 0.10_rp , w = Gze*(one-pol_e) )
    s_lab0%p25%m  = percentile( pol_n  , 0.25_rp , w = Gze*(one-pol_e) )
    s_lab0%p50%m  = percentile( pol_n  , 0.50_rp , w = Gze*(one-pol_e) )
    s_lab0%p75%m  = percentile( pol_n  , 0.75_rp , w = Gze*(one-pol_e) )
    s_lab0%p90%m  = percentile( pol_n  , 0.90_rp , w = Gze*(one-pol_e) )

    ! persistence: population OLS slope of log n(t+1) on log n(t) among continuing firms
    sw = zero ; sx = zero ; sy = zero ; sxx = zero ; sxy = zero
    do iz=1,nz
      wgt = Gz(iz)*(one-pol_e(iz)) ; if (wgt.le.zero) cycle
      xx  = log(pol_n(iz))
      do jz=1,nz
        if (pz(iz,jz).le.zero) cycle
        yy  = log(pol_n(jz)) ; ww = wgt*pz(iz,jz)
        sw  = sw + ww ; sx = sx + ww*xx ; sy = sy + ww*yy ; sxx = sxx + ww*xx*xx ; sxy = sxy + ww*xx*yy
      end do
    end do
    rho_n%m = (sxy - sx*sy/sw)/(sxx - sx*sx/sw)

    s_e%m = cien*sum(pol_e*Gz)/sum(Gz)

    return
  end subroutine equilibrium

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! residual of the free-entry condition at wage w (updates policies, value function and exit)

  function wres_hugo(w) result(f)
    implicit none
    real(rp) , intent(in) :: w
    real(rp)              :: f
    integer               :: iter,j,k
    real(rp)              :: test,zstar
    real(rp)              :: ev(nz)

    wage = w

    ! policy functions
    pol_n = nstar(zgrid,pol_tau,rrate+depre,wage)
    pol_k = kstar(zgrid,pol_tau,rrate+depre,wage)
    pol_y = output(zgrid,pol_k,pol_n)

    if (ecostinw.eq.0) pol_pi = pol_y - wage*pol_n - (rrate+depre)*pol_k - pol_tau*pol_y - fcost
    if (ecostinw.eq.1) pol_pi = pol_y - wage*pol_n - (rrate+depre)*pol_k - pol_tau*pol_y - fcost*wage

    ! pol_t  = tau*(max(0.0d0,pol_pi+pol_tau*pol_y)) + pol_tau*pol_y
    ! pol_pi = pol_pi - tau*(max(0.0d0,pol_pi+pol_tau*pol_y))

    pol_t  = tau*pol_pi + pol_tau*pol_y
    pol_pi = (one-tau)*pol_pi

    ! converge the value function of firms
    iter = 0 ; test = huge(one)
    do while ( iter.lt.iterm .and. test.gt.tol ) ; iter = iter + 1
      vfunc = vgues
      vgues = pol_pi + (one-eprob)*(one/(one+rrate))*max(matmul(pz,vfunc),zero)
      test  = maxval(abs(vgues-vfunc))
    end do
    if (test.gt.tol) print * , ' warning: firm value function did not converge, error = ' , test

    vfunc = vgues ; ev = matmul(pz,vfunc) ; pol_e = eprob ; where(ev.le.zero) pol_e = one

    j = 1
    do while (j.le.nz)
      if (ev(j).gt.zero) exit
      j = j + 1
    end do
    if (j.gt.1 .and. j.le.nz) then
      zstar = zgrid(j-1) + (zgrid(j)-zgrid(j-1))*(-ev(j-1))/(ev(j)-ev(j-1))
      do k=j-1,j
        pol_e(k) = eprob + (one-eprob)*min(one,max(zero,(zstar-zlo(k))/(zhi(k)-zlo(k))))
      end do
    end if

    ! free-entry condition
    f = ecost - sum(vfunc*Gze) ; if (ecostinw.eq.1) f = wage*ecost - sum(vfunc*Gze)

    return
  end function wres_hugo

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! unconditional mean of log productivity (entrants' mean is set relative to it)
  function zmean_u() result(zm)
    implicit none
    real(rp) :: zm
    zm = -half*stdz*stdz/(one-rhoz)
    return
  end function zmean_u

  subroutine calibmodel( )
    implicit none
    integer                    :: npar
    real(rp) , allocatable     :: x0(:),x1(:),y1(:),lb(:),ub(:)
    integer                    :: numiter,exitcode,itrho
    real(rp)                   :: drho,sgn,lastsgn,gap0
    
    ! initial wage: solve the model once with free entry at the current parameters
    if (withk.eq.1) depre = gamma*alpha/3.0_rp - rrate

    call equilibrium( ) ; wfix = wage

    gap0 = max(zmean_u() - muz0, 0.01_rp)
    
    drho = 0.002_rp ; lastsgn = zero ; itrho = 0

    ! delta is not calibrated (it is chosen in the menu). The fixed cost is calibrated to the exit rate,
    ! unless there is no endogenous exit: then it is zero and the exit rate is not a target.
    npar = 6 ; if (fixcf) npar = 5
    allocate(x0(npar),x1(npar),y1(npar),lb(npar),ub(npar))
    x0(1:5) = (/ gamma  , stdz0  , stdz   , gap0   , wfix         /)
    lb(1:5) = (/ 0.30_rp, 0.01_rp, 0.01_rp, 0.00_rp, 0.20_rp*wfix /)
    ub(1:5) = (/ 0.95_rp, 2.00_rp, 1.00_rp, 3.00_rp, 5.00_rp*wfix /)
    if (.not.fixcf) then
      x0(6) = fcost ; lb(6) = 0.00_rp ; ub(6) = 20.0_rp
    end if

    19 continue

    call simplex(funcmodel,x1,y1(1),numiter,exitcode,x0,iprint=2,lb=lb,ub=ub) ; x0 = x1
    call lmmin(funcmodel,x1,y1,numiter,exitcode,x0,shock=0.05_rp,iprint=2,usebro=1,lb=lb,ub=ub) ; x0 = x1
    y1 = funcmodel(x1)

    itrho = itrho + 1
    if ( abs(rho_n%m-rho_n%d).gt.0.002_rp .and. itrho.lt.30 ) then
      sgn = sign(one,rho_n%d-rho_n%m)
      if (sgn*lastsgn.lt.zero) drho = half*drho
      rhoz = min(0.995_rp,max(zero,rhoz + sgn*drho)) ; lastsgn = sgn
      goto 19
    end if

    ! back to the usual solution (free entry) with the implied ecost: it delivers wage = wfix
    wprev = wfix ; call equilibrium( )

    return
  end subroutine calibmodel

  function funcmodel(xp) result(resid)
    implicit none
    real(rp)               :: xp(:)
    integer                :: i
    real(rp) , allocatable :: resid(:) ; i = 1 ; allocate(resid(size(xp)))

    ! gamma =  exp(xp(i)) ; i=i+1
    ! stdz0 =  exp(xp(i)) ; i=i+1
    ! stdz  =  exp(xp(i)) ; i=i+1
    ! muz0  =  zmean_u() - exp(xp(i)) ; i=i+1     ! entrants start below the unconditional mean
    ! ecost =  exp(xp(i)) ; i=i+1
    ! fcost =  exp(xp(i)) ; i=i+1
    ! eprob =  one/(one+exp(-xp(i))) ; i=i+1      ! exogenous exit probability (delta)

    gamma =  xp(i) ; i=i+1
    stdz0 =  xp(i) ; i=i+1
    stdz  =  xp(i) ; i=i+1
    muz0  =  zmean_u() - xp(i) ; i=i+1     ! entrants start below the unconditional mean
    wfix  =  xp(i) ; i=i+1      ! wage (ecost is implied by free entry)
    if (size(xp).ge.6) then
      fcost =  xp(i) ; i=i+1    ! fixed cost (not calibrated without endogenous exit)
    end if
    !eprob =  xp(i) ; i=i+1      ! exogenous exit probability (delta)

    if (withk.eq.1) depre = gamma*alpha/3.0_rp - rrate

    call equilibrium(fixw=.true.) ; i = 1

    resid(i) = 2.0_rp*( one - ( s_lab%mean%m            ) / ( s_lab%mean%d            ) ) ; i=i+1 ! gamma
    resid(i) = 1.0_rp*( one - ( shr_labor%m             ) / ( shr_labor%d             ) ) ; i=i+1 ! stdz0
    resid(i) = 1.0_rp*( one - ( s_lab%std%m             ) / ( s_lab%std%d             ) ) ; i=i+1 ! stdz
    resid(i) = 1.0_rp*( one - ( s_lab0%mean%m           ) / ( s_lab0%mean%d           ) ) ; i=i+1 ! muz0
    resid(i) = 1.0_rp*( one - ( s_lab%p25%m             ) / ( s_lab%p25%d             ) ) ; i=i+1 ! wage
    if (size(xp).ge.6) then
    resid(i) = 2.0_rp*( one - ( s_e%m                   ) / ( s_e%d                   ) ) ; i=i+1 ! fcost
    end if
    !resid(i) = 1.0_rp*( one - ( s_lab%p90%m             ) / ( s_lab%p90%d             ) ) ; i=i+1 ! delta
    !resid(i) = 1.0_rp*( one - ( shr_firms%m             ) / ( shr_firms%d             ) ) ; i=i+1 ! eprob

    if (sum(resid(:)*resid(:)).lt.minermodel) then
      minermodel = sum(resid(:)*resid(:))
      call writehugo(2)
      return
    end if

    return

  end function funcmodel

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! averages by age, from entrants (age 0) to age 20, in the equilibrium just solved.
  ! A cohort starts with the distribution of entrants, Gze, and evolves with the same law of motion
  ! as the stationary distribution in equilibrium: h(a+1) = ((1-pol_e)*h(a))*pz. Summing h(a) over
  ! all ages gives the stationary distribution, Gz (up to the scale M/Me).
  ! Columns: age; survival (% of entrants still active); share of active firms (%); averages of
  ! productivity (z, in levels), output (y) and labor (n); and exit rate (%): share of the firms of
  ! age a that are not active at age a+1, 1 - mass(a+1)/mass(a). The last row is the average over all
  ! active firms; its exit rate is exits over active firms, which in the steady state equals Me/M.
  ! Written to unit un (called by writehugo, after the moments).

  subroutine agestats(un)
    implicit none
    integer  , intent(in) :: un
    integer  , parameter  :: amax = 20
    real(rp)              :: h(nz),hn(nz),mass,mass0,cmass
    integer               :: a
      write(un,'(a)') '  '
      write(un,'(a)') ' Averages by age (age 0 = entrants)'
      write(un,'(a)') ' ------------------------------------------------------------------------------'
      write(un,'(a)') '   age  survival(%)  share(%)           z           y           n   exit(%)'
      write(un,'(a)') ' ------------------------------------------------------------------------------'
      h = Gze ; mass0 = sum(h)
      do a = 0,amax
        mass  = sum(h) ; if (mass.le.zero) exit
        hn    = matmul((one-pol_e)*h,pz)
        cmass = sum(hn)                      ! mass of the cohort at age a+1
        write(un,'(i6,2f13.2,3f12.4,f10.2)') a, cien*mass/mass0, cien*mass*Me/M, &
          sum(exp(zgrid)*h)/mass, sum(pol_y*h)/mass, sum(pol_n*h)/mass, cien*(one-cmass/mass)
        h = hn
      end do
      write(un,'(a)') ' ------------------------------------------------------------------------------'
      write(un,'(a6,2a13,3f12.4,f10.2)') '   all','','', &
        sum(exp(zgrid)*Gz)/sum(Gz), sum(pol_y*Gz)/sum(Gz), sum(pol_n*Gz)/sum(Gz), cien*Me/M
      write(un,'(a)') ' ------------------------------------------------------------------------------'
    return
  end subroutine agestats

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! lambda of the paper (Section 2.6 and Appendix C.3), in the equilibrium just solved (phihat = 0):
  !
  !   lambda = E sum_t d_t iota_t x_t / E sum_t d_t x_t ,   x = z^(1/(1-gamma)) ,
  !
  ! the share of an entrant's expected discounted output that is subsidized, minus the share that is
  ! taxed. d_t = beta^t times the probability of being active at age t, with beta = 1/(1+rrate) as in
  ! the value function, and iota = 1 below the revenue-neutral cutoff z* and -1 above it. As
  ! phihat -> 0, z* splits the output of active firms (sum of x*Gz) in two halves; in the node whose
  ! cell contains z*, iota = 2f-1, with f the share of the cell below z* (the cell rule of the
  ! symmetric experiment). The cohort evolves as in agestats: h(t+1) = ((1-pol_e)*h(t))*pz.
  !
  ! Also reported: the wage response, dlog w/dphihat = lambda/((1-alpha)*gamma), and the wage term of
  ! B(0), (1-alpha*gamma)/((1-alpha)*gamma)*lambda (alpha = 0 without capital). As a check, lambda
  ! without discounting (beta = 1) weights each age as the cross-section does and must be ~0 by
  ! revenue neutrality. Written to unit un (called by writehugo, after the moments).

  subroutine lambdastat(un)
    implicit none
    integer  , intent(in) :: un
    integer  , parameter  :: tmax = 5000
    real(rp)              :: x(nz),iota(nz),h(nz),cum,halfx,f,zst,disc,dd,num,den,lam(2),aeff,cs
    integer               :: k,t,j

    ! revenue-neutral cutoff and iota
    x = exp(zgrid/(one-gamma)) ; halfx = half*sum(x*Gz) ; cum = zero ; iota = -one ; zst = zgrid(nz)
    do k = 1,nz
      if (cum + x(k)*Gz(k) .lt. halfx) then
        iota(k) = one ; cum = cum + x(k)*Gz(k)
      else
        f = (halfx-cum)/(x(k)*Gz(k)) ; iota(k) = two*f - one ; zst = zlo(k) + f*(zhi(k)-zlo(k))
        exit
      end if
    end do
    cs = sum(iota*x*Gz)/sum(x*Gz)                 ! cross-section: zero by construction

    ! lambda with discounting (j = 1) and without (j = 2, check)
    do j = 1,2
      disc = one/(one+rrate) ; if (j.eq.2) disc = one
      h = Gze ; dd = one ; num = zero ; den = zero
      do t = 0,tmax
        num = num + dd*sum(iota*x*h)
        den = den + dd*sum(x*h)
        h   = matmul((one-pol_e)*h,pz) ; dd = dd*disc
        if (sum(h).lt.1.0e-12_rp) exit
      end do
      lam(j) = num/den
    end do

    aeff = alpha*dble(withk)
      write(un,'(a)') '  '
      write(un,'(a)') ' Value of entry and distortions (Section 2.6)'
      write(un,'(a)') ' ------------------------------------------------------------------------------'
      write(un,'(a,f12.6)') ' lambda                                       ',lam(1)
      write(un,'(a,f12.6)') ' dlog w / dphihat = lambda/((1-alpha)*gamma)  ',lam(1)/((one-aeff)*gamma)
      write(un,'(a,f12.6)') ' wage term of B(0)                            ',(one-aeff*gamma)/((one-aeff)*gamma)*lam(1)
      write(un,'(a,f12.6)') ' cutoff z* (in levels)                        ',exp(zst)
      write(un,'(a)') ' checks (should be ~0):'
      write(un,'(a,f12.6)') '   lambda without discounting                 ',lam(2)
      write(un,'(a,f12.6)') '   cross-section, sum(iota*x*Gz)/sum(x*Gz)    ',cs
      write(un,'(a)') ' ------------------------------------------------------------------------------'
    return
  end subroutine lambdastat

  ! output file: the first lines are the parameters, in the order setversion reads them, so the
  ! calibration file can be copied into txtfiles/ (removing the date) to be used as input
  ! (xi and uz are written but not used). unitnum = 1: solution file, 2: calibration file

  subroutine writehugo(unitnum)
    implicit none
    integer , intent(in) , optional :: unitnum
    integer                         :: unitt ; unitt = 0
    !
    if (present(unitnum)) unitt = unitnum
    !
    99 continue
    !
    if (unitt.eq.1) open(unit=unitt,file=outfile("solve",vname),action='write',status="replace")
    if (unitt.eq.2) open(unit=unitt,file=outfile("calib",vname),action='write',status="replace")
    !
    write(unitt,33) rrate,'rrate'
    write(unitt,33) depre,'depre'
    write(unitt,33) eprob,'eprob'
    write(unitt,33) fcost,'fcost'
    write(unitt,33) ecost,'ecost'
    write(unitt,33) xi,'xi (not used)'
    write(unitt,33) uz,'uz (not used)'
    write(unitt,33) rhoz,'rhoz'
    write(unitt,33) stdz,'stdz'
    write(unitt,33) muz0,'muz0'
    write(unitt,33) stdz0,'stdz0'
    write(unitt,33) gamma,'gamma'
    write(unitt,33) eta0,'eta0'
    write(unitt,33) eta1,'eta1'
    write(unitt,33) dble(ecostinw),'ecostinw'
    write(unitt,33) dble(withk),'withk'
    write(unitt,90) '  '
    write(unitt,90) '  ' , 'Equilibrium:     '
    write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
    write(unitt,20) '  ' , 'Wage rate        ' , wage
    write(unitt,20) '  ' , 'Mass of firms    ' , M , cien*M/(aN+M)
    write(unitt,90) '  '
    write(unitt,90) '  ' , 'Productivity:    '
    write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
    write(unitt,20) '  ' , 'TFP              ' , TFP
    write(unitt,20) '  ' , 'TFP w. N+M       ' , TFPL
    write(unitt,20) '  ' , 'TFP w. net Y     ' , TFPC
    write(unitt,20) '  ' , 'TFP w. N+M & net ' , TFPLC
    write(unitt,20) '  ' , 'Output per worker' , aY/aL
    write(unitt,90) '  '
    write(unitt,90) '  ' , 'Aggregates: '
    write(unitt,90) '  ' , '^^^^^^^^^^^ '
    write(unitt,20) '  ' , 'Output      ' , aY
    write(unitt,20) '  ' , 'Capital     ' , aK , aK/aY
    write(unitt,20) '  ' , 'Labor       ' , aL , aN + M
    write(unitt,20) '  ' , 'Consumption ' , aC
    write(unitt,20) '  ' , 'Government  ' , aT
    write(unitt,20) '  ' , 'Exit rate   ' , s_e%m
    write(unitt,20) '  ' , 'Grid: mass at top node (should be ~0) ' , Gz(nz)/sum(Gz)
    write(unitt,20) '  ' , 'Grid: mass at top node (should be ~0) ' , Gz(size(Gz))/sum(Gz)
    write(unitt,90) '  '
    ! all the moments the model computes; * = used in the calibration
    write(unitt,90) '  ' , repeat('*',44)
    write(unitt,90) '  ' , 'MOMENTS             Model    Data    Diff   '
    write(unitt,90) '  ' ,  repeat('*',44)
    write(unitt,10) '  ' , 'Labor: average   ' , s_lab%mean%m , s_lab%mean%d , s_lab%mean%m/s_lab%mean%d - one , ' *'
    write(unitt,10) '  ' , 'Labor: std logs  ' , s_lab%std%m  , s_lab%std%d  , s_lab%std%m/s_lab%std%d - one   , ' *'
    write(unitt,10) '  ' , 'Labor: persisten ' , rho_n%m      , rho_n%d      , rho_n%m/rho_n%d - one             , ' *'
    write(unitt,10) '  ' , 'Labor: p10       ' , s_lab%p10%m  , s_lab%p10%d  , s_lab%p10%m/s_lab%p10%d - one   , '  '
    write(unitt,10) '  ' , 'Labor: p25       ' , s_lab%p25%m  , s_lab%p25%d  , s_lab%p25%m/s_lab%p25%d - one   , ' *'
    write(unitt,10) '  ' , 'Labor: p50       ' , s_lab%p50%m  , s_lab%p50%d  , s_lab%p50%m/s_lab%p50%d - one   , '  '
    write(unitt,10) '  ' , 'Labor: p75       ' , s_lab%p75%m  , s_lab%p75%d  , s_lab%p75%m/s_lab%p75%d - one   , '  '
    write(unitt,10) '  ' , 'Labor: p90       ' , s_lab%p90%m  , s_lab%p90%d  , s_lab%p90%m/s_lab%p90%d - one   , '  '
    write(unitt,90) '  ' ,  repeat('-',44)
    write(unitt,10) '  ' , 'E Labor: average ' , s_lab0%mean%m , s_lab0%mean%d , s_lab0%mean%m/s_lab0%mean%d - one , ' *'
    write(unitt,10) '  ' , 'E Labor: std logs' , s_lab0%std%m  , s_lab0%std%d  , s_lab0%std%m/s_lab0%std%d - one   , '  '
    write(unitt,10) '  ' , 'E Labor: p10     ' , s_lab0%p10%m  , s_lab0%p10%d  , s_lab0%p10%m/s_lab0%p10%d - one   , '  '
    write(unitt,10) '  ' , 'E Labor: p25     ' , s_lab0%p25%m  , s_lab0%p25%d  , s_lab0%p25%m/s_lab0%p25%d - one   , '  '
    write(unitt,10) '  ' , 'E Labor: p50     ' , s_lab0%p50%m  , s_lab0%p50%d  , s_lab0%p50%m/s_lab0%p50%d - one   , '  '
    write(unitt,10) '  ' , 'E Labor: p75     ' , s_lab0%p75%m  , s_lab0%p75%d  , s_lab0%p75%m/s_lab0%p75%d - one   , '  '
    write(unitt,10) '  ' , 'E Labor: p90     ' , s_lab0%p90%m  , s_lab0%p90%d  , s_lab0%p90%m/s_lab0%p90%d - one   , '  '
    write(unitt,90) '  ' ,  repeat('-',44)
    write(unitt,10) '  ' , 'Share firms n>50 ' , shr_firms%m   , shr_firms%d   , shr_firms%m/shr_firms%d - one     , '  '
    write(unitt,10) '  ' , 'Labor in n>50    ' , shr_labor%m   , shr_labor%d   , shr_labor%m/shr_labor%d - one     , ' *'
    write(unitt,90) '  ' ,  repeat('-',44)
    if (fixcf) then
    write(unitt,10) '  ' , 'Exit rate        ' , s_e%m         , s_e%d         , s_e%m/s_e%d - one                 , '  '
    else
    write(unitt,10) '  ' , 'Exit rate        ' , s_e%m         , s_e%d         , s_e%m/s_e%d - one                 , ' *'
    end if
    if (withk.eq.1) then
    write(unitt,10) '  ' , 'Capital/output   ' , aK/aY         , 3.0_rp        , (aK/aY)/3.0_rp - one              , ' *'
    end if
    write(unitt,90) '  ' ,  repeat('*',44)
    write(unitt,90) '  '
    !
    ! age profiles and lambda (Section 2.6)
    call agestats(unitt)
    call lambdastat(unitt)
    write(unitt,90) '  '
    !
    if (unitt.gt.0) close(unitt)
    !
    if (unitt.eq.1) then
      unitt = 0 ; goto 99
    end if
    !
    return
    !
    90 format (a,a)
    10 format (a,a,3(f8.2),a)
    20 format (a,a,3(f12.6))
    33 format (f12.8,' ! ',a)
  end subroutine writehugo

end module hugo
