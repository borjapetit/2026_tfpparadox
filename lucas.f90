
! ======================================================================================================
! Lucas (1978) economy: occupational choice with Pareto-distributed managerial ability
!
! Calibration: gamma -> average firm size, xi -> share of employment in firms with 50+ workers,
! rrate -> capital-to-output ratio (taken from the baseline calibration when withk = 0)
!
! Files (names built in parameters.f90 from the economy chosen in the menu, e.g. lucas_tau30):
!   parameters:   txtfiles/calib_<economy>_<version>.txt
!   calibration:  results/<date>_calib_<economy>_<version>.txt
!   solution:     results/<date>_solve_<economy>_<version>.txt
!   experiments:  results/<economy>_<version>.txt and results/<economy>_<sensitivity exercise>.txt
! ======================================================================================================

module lucas

  use parameters
  use omp_lib
  use toolkit , only : lmmin,simplex,vect,varmean,varstd,percentile,grid,interpolate,&
  correlation,print_err,bisection,brent,olsreg,num2text,tauchen,error
  
  implicit none

  ! state space
  integer  , parameter :: nz = 20000

  ! productivity grid
  real(rp) , dimension(nz) :: zgrid
  real(rp) , dimension(nz) :: pz

  ! distribution of firms
  real(rp) , dimension(nz) :: Gz
  real(rp) , dimension(nz) :: Gze

  ! value and policy function
  real(rp) , dimension(nz) :: vfunc
  real(rp) , dimension(nz) :: vgues
  real(rp) , dimension(nz) :: pol_k
  real(rp) , dimension(nz) :: pol_n
  real(rp) , dimension(nz) :: pol_y
  real(rp) , dimension(nz) :: pol_pi
  real(rp) , dimension(nz) :: pol_e
  real(rp) , dimension(nz) :: pol_t
  real(rp) , dimension(nz) :: pol_tau 

  real(rp) :: bestx(4)
  real(rp) :: wprev = one   ! wage of the last equilibrium solved (initial guess for findroot)

  type(results) :: base
  type(results) :: eff
  type(results) :: distort
  

  contains

  subroutine solve_lucas( mode )
    !
    implicit none
    !
    integer , intent(in) :: mode
    integer              :: ifile
    real(rp)             :: par0

    ! --------------------------------------------------------------
    ! initialize some parameter values

    atau = zero ; ifile = 878 ; minermodel = huge(one)

    ! --------------------------------------------------------------
    ! options (menu of main.f90):
    !    1  solve the baseline model                 2  calibrate the baseline model
    !    3  symmetric experiment (baseline)          4  sensitivity of the symmetric experiment
    !    5  solve the model without capital          6  calibrate the model without capital
    ! Each version reads its parameters from txtfiles/calib_<economy>_<version>.txt (see setversion).
    ! When the model without capital is calibrated, the interest rate is taken from the baseline file
    ! of the same economy.

    if (mode.ne.1 .and. mode.ne.5) print_err = .false.

    select case (mode)
    case (1) ; call setversion(1) ; call equilibrium( ) ; call writelucas(1)
    case (2) ; call setversion(1) ; call calibmodel( )  ; call writelucas(1)
    case (3) ; call setversion(1) ; call runexperiment( )
    case (4)
      call setversion(1)
      write(*,'(/,a,/)') ' Higher gamma'
      par0 = gamma ; gamma = par0*1.05_rp ; call runexperiment("gamma_up5") ; gamma = par0
      write(*,'(/,a,/)') ' Lower gamma'
      par0 = gamma ; gamma = par0*0.95_rp ; call runexperiment("gamma_down5") ; gamma = par0
      write(*,'(/,a,/)') ' Higher xi'
      par0 = xi ; xi = par0*1.10_rp ; call runexperiment("xi_up10") ; xi = par0
      write(*,'(/,a,/)') ' Lower xi'
      par0 = xi ; xi = par0*0.90_rp ; call runexperiment("xi_down10") ; xi = par0
      write(*,'(/,a,/)') ' No capital'
      call setversion(2) ; call runexperiment( )
    case (5) ; call setversion(2) ; call equilibrium( ) ; call writelucas(1)
    case (6) ; call setversion(2) ; rrate = baserate( ) ; call calibmodel( ) ; call writelucas(1)
    case default
      write(*,'(a)') ' invalid option for this model'
    end select

    print_err = .true.

    return
    
    contains

    ! this subroutine runs the symmetric experiment. Results: results/<economy>_<exercise>.txt, where
    ! the exercise is the version of the model (baseline, nocapital) or, if given, the name of the
    ! sensitivity exercise
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
      atau = zero ; print * , '   '

      ! baseline economy
      tau = taubase ; call equilibrium(base)

      ! efficient allocation
      tau = 0.0_rp ; call equilibrium(eff) 

      ! print initial results
      call printexer(ifile,base,eff) ; call printexer(0,base,eff)
      call printexer(ifile,base)     ; call printexer(0,base)

      ! distorted economies
      tau = taubase ; atau = zero ; exper : do ; atau = atau + 0.01_rp
        call equilibrium(distort) ; call printexer(ifile,base,distort) ; call printexer(0,base,distort)
        if (atau.gt.0.20 .and. distort%TFP0/base%TFP0.lt.0.95_rp) exit exper
        if (atau.gt.0.50) exit exper
      end do exper

      close(ifile)
      
      return
    end subroutine runexperiment

    !
  end subroutine solve_lucas

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! version of the model: 1 baseline, 2 no capital.
  ! Sets withk and vname, and reads the parameters from txtfiles/calib_<economy>_<version>.txt

  subroutine setversion(iv)
    implicit none
    integer , intent(in) :: iv
    logical              :: exists
    withk = 1 ; vname = "baseline"
    if (iv.eq.2) then ; withk = 0 ; vname = "nocapital" ; end if
    inquire(file=infile(vname),exist=exists)
    if (.not.exists) then
      write(*,'(a)') ' parameter file not found: '//trim(infile(vname)) ; stop
    end if
    open(unit=9,file=infile(vname),action="read")
      read(9,*) rrate
      read(9,*) xi
      read(9,*) uz
      read(9,*) gamma
    close(9)
    atau = zero
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
  ! solve the model given some equilibrium prices (wage and mass of entrants)

  subroutine equilibrium(res)
    implicit none
    type(results) , intent(out) , optional :: res
    real(rp)                               :: maxw
    real(rp)                               :: minw
    real(rp)                               :: equilw
    real(rp)                               :: maxz
    real(rp)                               :: minz
    real(rp)                               :: equilz
    integer                                :: iterz

    iterz = 0 ; minz = zero ; maxz = zero ; equilz = zero ; wcondi = 99.9_rp

    ! ------------------------------------------------------------
    ! distribution over permanent productivity

    call pareto( zgrid , pz , uz , xi , uz/xi , nz )

    ! ------------------------------------------------------------
    ! find the cutoff productivity for distortions

    if (abs(atau).gt.tol) then
      minz   = minval(zgrid)
      maxz   = maxval(zgrid)
      equilz = huge(one)
    else
      minz   = maxval(zgrid) 
      maxz   = maxval(zgrid)
      equilz = zero
    end if

    10 continue

    ! ------------------------------------------------------------
    ! distortions

    if (abs(atau).gt.tol) then
      iterz = iterz + 1
      zphi  = half*( maxz + minz )
      where (zgrid.gt.zphi) pol_tau =   atau
      where (zgrid.le.zphi) pol_tau = - atau 
    else
      pol_tau = zero 
    end if

    ! ------------------------------------------------------------
    ! find the equilibrium wage

    ! Brent's method, starting from the wage of the previous equilibrium (wprev)
    call findroot(wres_lucas,wprev,1.0e5_rp,wage,equilw,tol/cien,tol/cien)
    wprev = wage

    wcondi = equilw

    ! ------------------------------------------------------------
    ! update distortions

    if (abs(atau).gt.tol) then
      equilz = sum(pol_tau*pol_y*Gz)*cien/sum(Gz)
      if (iterz.lt.1000 .and. abs(maxz-minz).gt.tol .and. abs(equilz).gt.tol) then
        if (equilz.gt.zero) minz = zphi
        if (equilz.lt.zero) maxz = zphi
        goto 10
      end if
    end if

    ! ------------------------------------------------------------
    ! equilibrium objects

    M  = sum(Gz)                ! mass of firms
    pol_e = zero                ! no exit in the static economy
    Me = M                      ! mass of entrants
    aK = sum(Gz*pol_k)          ! aggregate capital
    aN = sum(Gz*pol_n)          ! aggregate labor in production
    aL = aN                     ! aggregate labor supply
    aP = sum(Gz*pol_pi)         ! aggregate profits
    aT = sum(Gz*pol_t)          ! aggregate tax revenues
    aY = sum(Gz*pol_y)          ! aggregate output
    aC = wage*aN + aP + aT      ! aggregate consumption

    ! aggregate TFP
    aZ1 = varmean( (((one - pol_tau)**gamma)*exp(zgrid))**(one/(one-gamma)) , w = Gz )
    aZ2 = varmean( (((one - pol_tau)       )*exp(zgrid))**(one/(one-gamma)) , w = Gz )**gamma
    TFP0  = ( M**(one-gamma) )*aZ1/aZ2
    aZ    = ( TFP0 / ( M**(one-gamma) )  )**(one/(one-gamma))
    
    ! TFP as in Guner et al
    TFP1 = aY / ( ( aK**(gamma*alpha) ) * ( aL**(one-(gamma*alpha)) ) )

    ! measured TFP as an econometrician
    TFP2 = aY / ( ( aK**(one-(wage*aL/aY)) ) * ( aL**(wage*aL/aY) ) )

    ! TFP using total population
    TFP3 = aY / ( ( (aK**alpha) * ( (aL+M)**(one-alpha) ) ) ** gamma )

    ! ------------------------------------------------------------
    ! store results

    if (present(res)) then

      res%wage = wage
      res%M    = M
      res%Me   = M
      res%aK   = aK
      res%aN   = aN
      res%aL   = aL
      res%aY   = aY
      res%aC   = aC
      res%aP   = aP
      res%aT   = aT
      res%aZ   = aZ
      res%aZ1  = aZ1
      res%aZ2  = aZ2
      res%TFP0 = TFP0
      res%TFP1 = TFP1
      res%TFP2 = TFP2
      res%TFP3 = TFP3

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

      if (allocated(res%large) ) deallocate(res%large) ; allocate (res%large(nz))
      if (allocated(res%small) ) deallocate(res%small) ; allocate (res%small(nz))
      if (allocated(res%taxed) ) deallocate(res%taxed) ; allocate (res%taxed(nz))

      res%k   = pol_k
      res%n   = pol_n
      res%y   = pol_y
      res%pi  = pol_pi
      res%e   = pol_e
      res%t   = pol_t
      res%tau = pol_tau
      res%z   = zgrid
      res%G   = Gz
      res%Ge  = Gz

      res%small = ( vect(pol_n).lt.percentile(vect(pol_n),0.10_rp,w=vect(Gz)) )
      res%large = ( vect(pol_n).gt.percentile(vect(pol_n),0.90_rp,w=vect(Gz)) )

      res%taxed = zero ; where (vect(pol_tau).gt.tol) res%taxed = one

    end if

    ! ------------------------------------------------------------
    ! statistics and aggregate variables

    s_lab%mean%m = varmean(    pol_n           , w = Gz )
    s_lab%std%m  = varstd(     log(pol_n)      , w = Gz )
    s_lab%p10%m  = percentile( pol_n  , 0.10_rp , w = Gz )
    s_lab%p25%m  = percentile( pol_n  , 0.25_rp , w = Gz )
    s_lab%p50%m  = percentile( pol_n  , 0.50_rp , w = Gz )
    s_lab%p75%m  = percentile( pol_n  , 0.75_rp , w = Gz )
    s_lab%p90%m  = percentile( pol_n  , 0.90_rp , w = Gz )

    shr_firms%m  = cien*sum(Gz      ,mask=pol_n.ge.50.0_rp)/sum(Gz)
    shr_labor%m  = cien*sum(pol_n*Gz,mask=pol_n.ge.50.0_rp)/sum(pol_n*Gz)

    s_lab0%mean%m = zero
    s_lab0%std%m  = zero
    s_lab0%p10%m  = zero
    s_lab0%p25%m  = zero
    s_lab0%p50%m  = zero
    s_lab0%p75%m  = zero
    s_lab0%p90%m  = zero
    
    rho_n%m = zero

    return
  end subroutine equilibrium

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! excess labor demand at wage w. It also updates policies, occupational choice and the
  ! distribution of active firms (module variables), used by equilibrium() afterwards.

  function wres_lucas(w) result(f)
    implicit none
    real(rp) , intent(in) :: w
    real(rp)              :: f

    wage = w

    ! policy functions
    pol_n  = nstar(zgrid,pol_tau,rrate,wage)
    pol_k  = kstar(zgrid,pol_tau,rrate,wage)
    pol_y  = output(zgrid,pol_k,pol_n)
    pol_pi = (one-pol_tau)*pol_y - wage*pol_n - rrate*pol_k
    pol_t  = tau*pol_pi + pol_tau*pol_y
    pol_pi = (one-tau)*pol_pi

    ! mass of firms and equilibrium
    Gz  = zero ; where (pol_pi.gt.wage) Gz = pz
    Gze = Gz
    M   = sum(Gz)
    aN  = one - M

    ! labor market clearing condition
    f = sum(pol_n*Gz) - aN

    return
  end function wres_lucas

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine calibmodel( )
    use toolkit , only : simplex,lmmin,normalize,denormalize
    implicit none
    real(rp) , dimension(3) :: x0,x1
    real(rp) , dimension(3) :: y1
    real(rp)                :: toty
    integer                 :: numiter,exitcode,ipp
    x0 = (/ xi,gamma*xi,rrate /) ; x0 = log(x0)
    do ipp = 0,size(x0)
      if (ipp.gt.0) then
        call random_number(x1) ; x0 = x0 + 0.1_rp*(x1 - 0.5_rp)*abs(x0)
      end if
      call lmmin(funcmodel,x1,y1,numiter,exitcode,x0,usebro=0,iprint=2) ; x0 = x1
      call simplex(sumfuncmodel,x1,toty,numiter,exitcode,x0,iprint=2)
      x0 = x1
    end do
    do ipp = 1,500
      call random_number(x1)
      x0   = bestx(1:size(x0)) + 0.1_rp*(x1 - 0.5_rp)*abs(bestx(1:size(x0)))
      toty = sumfuncmodel(x0)
      print * , ipp , toty , minermodel
    end do
    return
  end subroutine calibmodel

  function funcmodel(xp) result(resid)
    use toolkit , only : denormalize
    implicit none
    real(rp)               :: xp(:)
    real(rp) , allocatable :: resid(:) ; allocate(resid(size(xp))) ; resid = zero

    xi    = exp(xp(1))
    gamma = exp(xp(2))/xi
    
    if (withk.eq.1) rrate = gamma*alpha/3.0_rp

    call equilibrium( )

    !if (size(xp).le.3)
    resid(1) = 5.0_rp*( one - s_lab%mean%m / s_lab%mean%d )
    !if (size(xp).le.3) 
    resid(2) = 1.0_rp*( one - shr_labor%m  / shr_labor%d  ) !
    !if (size(xp).eq.3) resid(3) = 1.0_rp*( one - (aK/aY) / 3.0_rp )

    !if (size(xp).le.3) resid(2) = 1.0_rp*( one - (s_lab%p10%m/s_lab%mean%m) / (s_lab%p10%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(3) = 1.0_rp*( one - (s_lab%p25%m/s_lab%mean%m) / (s_lab%p25%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(4) = 1.0_rp*( one - (s_lab%p50%m/s_lab%mean%m) / (s_lab%p50%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(5) = 1.0_rp*( one - (s_lab%p75%m/s_lab%mean%m) / (s_lab%p75%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(6) = 1.0_rp*( one - (s_lab%p90%m/s_lab%mean%m) / (s_lab%p90%d/s_lab%mean%d)  )

    if (sum(resid(:)*resid(:)).lt.minermodel) then
      minermodel = sum(resid(:)*resid(:))
      bestx      = xp
      call writelucas(2)
      return
    end if

    return
  end function funcmodel
  function sumfuncmodel(xp) result(resid)
    use toolkit , only : denormalize
    implicit none
    real(rp) :: xp(:),res(size(xp)+0),resid
    res = funcmodel(xp) ; resid = sum(res(:)*res(:))
    return
  end function sumfuncmodel

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine writelucas(unitnum)
    use toolkit , only : varmean,percentile,vect
    implicit none
    integer , intent(in) , optional :: unitnum
    integer                         :: unitt ; unitt = 0
    !
    if (present(unitnum)) unitt = unitnum
    !
    ! unitnum = 1: solution file, 2: calibration file
    if (unitt.eq.1) open(unit=unitt,file=outfile("solve",vname),action='write',status="replace")
    if (unitt.eq.2) open(unit=unitt,file=outfile("calib",vname),action='write',status="replace")
    !
    99 continue
    !
    ! parameters first, in the order setversion reads them (the calibration file can be copied into
    ! txtfiles/, removing the date, to be used as input)
    write(unitt,33) rrate,'rrate'
    write(unitt,33) xi,'xi'
    write(unitt,33) uz,'uz'
    write(unitt,33) gamma,'gamma'
    write(unitt,33) dble(withk),'withk (not read)'
    write(unitt,90) '  '
    write(unitt,90) '  '
    write(unitt,90) '  ' , 'Equilibrium:     '
    write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
    write(unitt,20) '  ' , 'Wage rate        ' , wage
    write(unitt,20) '  ' , 'Share managers   ' , M*cien
    write(unitt,90) '  '
    write(unitt,90) '  ' , 'Productivity:    '
    write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
    write(unitt,20) '  ' , 'Z                ' , aZ   
    write(unitt,20) '  ' , 'TFP              ' , TFP0 
    write(unitt,20) '  ' , 'TFP as in Guner  ' , TFP1
    write(unitt,20) '  ' , 'TFP as in Solow  ' , TFP2
    write(unitt,20) '  ' , 'TFP with N+M     ' , TFP3 
    write(unitt,20) '  ' , 'Output per worker' , aY/aL
    write(unitt,90) '  '
    write(unitt,90) '  ' , 'Aggregates: '
    write(unitt,90) '  ' , '^^^^^^^^^^^ '
    write(unitt,20) '  ' , 'Output      ' , aY
    write(unitt,20) '  ' , 'Capital     ' , aK , aK/aY
    write(unitt,20) '  ' , 'Labor       ' , aL 
    write(unitt,20) '  ' , 'Consumption ' , aC
    write(unitt,20) '  ' , 'Government  ' , aT
    write(unitt,20) '  ' , 'Labor sh.   ' , wage*aL/(rrate*aK + aP + wage*aL + aT)
    write(unitt,20) '  ' , 'Capital sh. ' , (rrate*aK + aP)/(rrate*aK + aP + wage*aL + aT)
    write(unitt,90) '  '
    write(unitt,90) '  '
    write(unitt,90) '  ' , repeat('*',41)
    write(unitt,90) '  ' , 'MOMENTS             Model    Data    Diff '
    write(unitt,90) '  ' ,  repeat('*',41)
    write(unitt,10) '  ' , 'Labor: average   ' , s_lab%mean%m , s_lab%mean%d , s_lab%mean%m/s_lab%mean%d - one
    write(unitt,10) '  ' , 'Labor: std logs  ' , s_lab%std%m  , s_lab%std%d  , s_lab%std%m/s_lab%std%d - one
    write(unitt,90) '  ' ,  repeat('-',41)
    write(unitt,10) '  ' , 'Labor: p10       ' , s_lab%p10%m  , s_lab%p10%d  , s_lab%p10%m/s_lab%p10%d - one
    write(unitt,10) '  ' , 'Labor: p25       ' , s_lab%p25%m  , s_lab%p25%d  , s_lab%p25%m/s_lab%p25%d - one
    write(unitt,10) '  ' , 'Labor: p50       ' , s_lab%p50%m  , s_lab%p50%d  , s_lab%p50%m/s_lab%p50%d - one
    write(unitt,10) '  ' , 'Labor: p75       ' , s_lab%p75%m  , s_lab%p75%d  , s_lab%p75%m/s_lab%p75%d - one
    write(unitt,10) '  ' , 'Labor: p90       ' , s_lab%p90%m  , s_lab%p90%d  , s_lab%p90%m/s_lab%p90%d - one
    write(unitt,90) '  ' ,  repeat('-',41)
    write(unitt,10) '  ' , 'Share firms n>50 ' , shr_firms%m  , shr_firms%d  , shr_firms%m/shr_firms%d - one
    write(unitt,10) '  ' , 'Labor in n>50    ' , shr_labor%m  , shr_labor%d  , shr_labor%m/shr_labor%d - one
    write(unitt,90) '  ' ,  repeat('*',41)
    write(unitt,90) '  '
    write(unitt,90) '  '
    !
    if (unitt.gt.0) close(unitt)
    !
    if (unitt.gt.0) then
      unitt = 0 ; goto 99
    end if
    !
    return
    90 format (a,a)
    10 format (a,a,10(f8.2))
    20 format (a,a,3(f12.6))
    33 format (f12.8,' ! ',a)
  end subroutine writelucas

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module lucas