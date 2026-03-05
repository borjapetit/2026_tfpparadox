module lucas

  use parameters
  use toolkit , only : lmmin,simplex,vect,varmean,varstd,percentile,grid,interpolate,&
  correlation,print_err,bisection,brent,olsreg,num2text
  
  implicit none

  ! state space
  integer  , parameter     :: nz = 5000
  integer                  :: fixedN
  real(rp) , dimension(nz) :: zgrid
  real(rp) , dimension(nz) :: pz
  real(rp) , dimension(nz) :: Gz
  real(rp)                 :: taus,taul,zphi

  ! ------------------------------------------------------------
  ! model parameters
  ! ------------------------------------------------------------ 

  real(rp) :: minermodel,bestx(4)
  real(rp) :: condi
  integer  :: withk
  ! equilibrium objects: wage, mass of firms, mass of entrants, totallabor supply
  real(rp) :: wage,M,Lbar,zhat
  ! caPtal cost
  real(rp) :: rrate
  ! technology
  real(rp) :: gamma,eprob
  ! productivity distribution
  real(rp) :: xi,uz,sigma
  ! distortions: level and curvature
  real(rp) :: atau,phi1,phi0
  ! aggregates
  real(rp) :: aK,aN,aL,aY,aC,aP,aT,aZ,aTFP,nTFP,mTFP,eZ,eTFP  ! distorted economy
  ! value and policy function
  real(rp) , dimension(nz) :: vfunc,vgues,pol_k,pol_n,pol_y,pol_pi,pol_e,pol_t,pol_tau 

  contains

  subroutine solve_lucas( mode )
    implicit none
    integer , intent(in) :: mode
    real(rp) :: aux,vecr(14),vecr0(14)
    integer  :: ia,iy,ind,ifile,modsol
    real(rp) :: M0,dist0(nz),distort0(nz),wage0,aK0,TFP0,par0
    character(len=100) :: filename

    modsol = mode

    ! load model settings
    open(unit=9,file="txtfiles/in_lucas_0.txt",action="read") 
      read(9,*) aux ; withk = int(aux)
    close(9)

    85 continue

    ! baseline economy
    if (withk.eq.1) then
      open(unit=9,file="txtfiles/in_lucas_base.txt",action="read")
    ! economy without capital
    else if (withk.eq.0) then
      open(unit=9,file="txtfiles/in_lucas_nok.txt",action="read")
    ! bad combination of model settings
    else
      write(*,*) ' Error: combination of parameters not implemented '
      read * , aux
    end if

    ! read parameters of the model
    read(9,*) rrate
    read(9,*) xi
    read(9,*) uz
    read(9,*) sigma
    read(9,*) gamma
    close(9)

    ! initialize some parameter values
    atau = zero ; phi0 = zero ; phi1 = zero

    if (modsol.eq.1) then
      
      call equilibrium( )
      call writelucas(1)
      call writelucas( )

    elseif (modsol.eq.2) then

      print_err = .false. ; ifile = 878

      write(*,'(/,a,/)') ' Revenue neutral distortions'
      filename = "base" ; call runexperiment( )

      print_err = .true.

    elseif (modsol.eq.3) then
    
      print_err  = .false.
      minermodel = huge(one)

      call calibmodel( )
      call writelucas( )

      print_err = .true.

    elseif (modsol.eq.4) then

      print_err = .false. ; ifile = 878

      write(*,'(/,a,/)') ' Baseline'
      filename = "s_base" ; call runexperiment( )

      write(*,'(/,a,/)') ' Higher gamma'
      par0 = gamma ; gamma = par0*1.05d0
      filename = "s_gamma_1" ; call runexperiment( )
      gamma = par0

      write(*,'(/,a,/)') ' Lower gamma'
      par0 = gamma ; gamma = par0*0.95d0
      filename = "s_gamma_0" ; call runexperiment( )
      gamma = par0

      write(*,'(/,a,/)') ' Higher xi'
      par0 = xi ; xi = par0*1.10d0
      filename = "s_xi_1" ; call runexperiment( )
      xi = par0

      write(*,'(/,a,/)') ' Lower xi'
      par0 = xi ; xi = par0*0.90d0
      filename = "s_xi_0" ; call runexperiment( )
      xi = par0

      print_err = .true.
    
    end if

    if (modsol.eq.2) then
      if (withk.eq.1) then
        withk = 0
        goto 85
      else if (withk.eq.0) then
        withk  = 1
        modsol = 4
        goto 85
      end if
    end if
    
    return
    contains
    ! this subroutine runs each experiment
    subroutine runexperiment( )
      implicit none
      character (len=100) :: ecotype
      real(rp)            :: tfpp

      ! initialize the economy
      atau = zero ; taus = zero ; taul = zero

      if (withk.eq.1) then
        ecotype = 'base'
      else if (withk.eq.0) then
        ecotype = 'nok'
      end if

      ! setup file for storing results
      open(unit=ifile,file="results/l_"//trim(adjustl(ecotype))//"_"//&
      trim(adjustl(filename))//".txt",action='write',status="replace")

      ! undistorted economy
      call equilibrium( ) ; call printexer(1) ; tfpp = aggTFP(M,pol_tau,Gz)

      ! distorted economies
      do ia = 1,20 ; atau = dble(ia)/cien
        call bisection(findzphi,zphi,iy,ind,maxval(zgrid),minval(zgrid),iprint=0)
        condi = findzphi(zphi) ; call printexer( ) ; tfpp = aggTFP(M,pol_tau,Gz)
      end do

      ! if tfp increases with 20% distortion, keep increasing it until TFP falls
      if (tfpp/TFP0.gt.0.95d0) then
        84 atau = atau + 0.01d0
        call bisection(findzphi,zphi,iy,ind,maxval(zgrid),minval(zgrid),iprint=0)
        condi = findzphi(zphi) ; call printexer(0) ; tfpp = aggTFP(M,pol_tau,Gz)
        if (tfpp/TFP0.gt.0.95d0 .and. atau.lt.0.5d0) goto 84
      end if

      close(ifile)
      
      return
    end subroutine runexperiment
    ! this subroutine prints the results of each experiment
    subroutine printexer(indic)
      implicit none
      integer , optional , intent(in) :: indic
      real(rp) :: Expz,Etau,Etau0,Etaulow,Etauhigh
      real(rp) :: T000,T100,T010,T001,T110,T101,T011,T111
      real(rp) :: dM,dphi,dF,dwage

      Expz     = sum( exp(zgrid)*Gz )/sum( Gz )
      Etau     = sum( pol_y*pol_tau*Gz                             )/sum( pol_y*Gz                             )
      Etau0    = sum( pol_y*pol_tau*Gz                             )/sum( pol_y*Gz                             )
      Etaulow  = sum( pol_y*pol_tau*Gz , mask=pol_n.lt.s_lab%p10%m )/sum( pol_y*Gz , mask=pol_n.lt.s_lab%p10%m )
      Etauhigh = sum( pol_y*pol_tau*Gz , mask=pol_n.gt.s_lab%p90%m )/sum( pol_y*Gz , mask=pol_n.gt.s_lab%p90%m )
      
      vecr  = (/ Etau,Etau0,Etaulow,Etauhigh,Expz,M,aY,aN,aN+M,aK,aY/aN,aC,mTFP,nTFP /)

      if (present(indic) .and. indic.eq.1) then  
        M0       = M
        dist0    = Gz
        distort0 = zero
        aK0      = aK
        wage0    = wage
        vecr0    = vecr
        TFP0     = aggTFP(M,pol_tau,Gz)
        write(ifile,'(a1,30(f9.3))') 'X',atau*cien,wage,TFP0,TFP0,zero,zero,zero,zero,vecr,condi,zphi 
        write(ifile,'(50(a9))') 'atau','wage','TFP0','TFP1','DTFP','DM','DPhi','DF',&
          'Etau','Etau0','Etaul','Etauh','Expz','M','aY','aN','aL','aK','YN','aC','mTFP','nTFP','check','zphi'
        write(*,'(50(a9))') 'atau','wage','TFP0','TFP1','DTFP','DM','DPhi','DF',&
          'Etau','Etau0','Etaul','Etauh','Expz','M','aY','aN','aL','aK','YN','aC','mTFP','nTFP','check','zphi'
      end if
      
      dwage      = cien*( wage/wage0             - one )
      vecr(5:14) = cien*( vecr(5:14)/vecr0(5:14) - one )

      T000 = aggTFP( M0 , distort0 , dist0 )
      T100 = aggTFP( M  , distort0 , dist0 )
      T010 = aggTFP( M0 , pol_tau  , dist0 )
      T001 = aggTFP( M0 , distort0 , Gz    )
      T110 = aggTFP( M  , pol_tau  , dist0 )
      T101 = aggTFP( M  , distort0 , Gz    )
      T011 = aggTFP( M0 , pol_tau  , Gz    )
      T111 = aggTFP( M  , pol_tau  , Gz    )

      dM   = (1.d0/3.d0)*(T100 - T000) + (1.d0/6.d0)*((T110 - T010) + (T101 - T001)) + (1.d0/3.d0)*(T111 - T011)
      dphi = (1.d0/3.d0)*(T010 - T000) + (1.d0/6.d0)*((T110 - T100) + (T011 - T001)) + (1.d0/3.d0)*(T111 - T101)
      dF   = (1.d0/3.d0)*(T001 - T000) + (1.d0/6.d0)*((T101 - T100) + (T011 - T010)) + (1.d0/3.d0)*(T111 - T110)

      write(ifile,'(30(f9.3))') atau*cien,dwage,T000,T111,cien*(T111-T000)/T000,cien*dM/T000,cien*dphi/T000,cien*dF/T000,&
        vecr,condi,zphi
      write(*,'(30(f9.3))') atau*cien,dwage,T000,T111,cien*(T111-T000)/T000,cien*dM/T000,cien*dphi/T000,cien*dF/T000,&
        vecr,condi,zphi

      return
    end subroutine printexer
    ! this function computes TFP given mass of firms, distortions and distribution of firms
    function aggTFP(mass,distort,dist) result (agtfp)
      implicit none
      real(rp) :: mass,distort(:),dist(:),agtfp
        agtfp = ( mass**(one-gamma) )* &
                 ( varmean( (((one - distort)**gamma)*exp(zgrid))**(one/(one-gamma)) , w = dist )) / &
                 ( varmean( (((one - distort)       )*exp(zgrid))**(one/(one-gamma)) , w = dist )**gamma )
      return
    end function aggTFP
    ! this function is used to find the threshold productivity for each experiment
    function findzphi(xp) result(resid)
      implicit none
      real(rp) :: xp
      real(rp) :: resid
      zphi = xp ; taul = atau ; taus = -atau ; call equilibrium( )
      resid = sum(pol_tau*pol_y*Gz)*cien
      return
    end function findzphi
    ! this function is used to find the level of distortion that makes TFP equal to the baseline
    function findatau0(xp) result(resid)
      implicit none
      real(rp) :: xp
      real(rp) :: resid,tfpa0,tfpa1
      atau = xp
      call bisection(findzphi,zphi,iy,ind,maxval(zgrid),minval(zgrid),iprint=0)
      tfpa0 = aggTFP( M0 , distort0 , dist0 )
      tfpa1 = aggTFP( M , pol_tau , Gz )
      resid = cien*( one - tfpa1/tfpa0 )
      return
    end function findatau0
  end subroutine solve_lucas

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! solve the model given some equilibrium prices (wage and mass of entrants)

  subroutine equilibrium( )
    use omp_lib
    use toolkit , only : grid,tauchen,olsreg,varmean,varstd,percentile,vect,interpolate
    implicit none
    real(rp) :: maxw,minw,equilw
    real(rp) :: lshare,zaux1,zaux2

    ! ------------------------------------------------------------
    ! distribution over permanent productivity
    call pareto(zgrid,pz,uz,xi,uz/xi,nz)

    ! ------------------------------------------------------------
    ! distortions

    if (abs(taus)+abs(taul).lt.tol) then
      pol_tau = one - exp( phi0 + phi1*(zgrid(nz)-zgrid) )
    else
      where (zgrid.gt.zphi) pol_tau = taul
      where (zgrid.le.zphi) pol_tau = taus
    end if

    ! ------------------------------------------------------------
    ! find the equilibrium wage

    maxw = 3000.0d0 ; minw = 0.0d0 ; equilw = huge(one)
    do while (abs(maxw-minw)*10.0d0.gt.tol .and. abs(equilw)*10.0d0.gt.tol)

      ! set the wage rate 
      wage  = 0.5d0*(maxw + minw)

      ! policy functions
      pol_n  = nstar(zgrid,pol_tau)
      pol_k  = kstar(zgrid,pol_tau)
      pol_y  = output(zgrid,pol_k,pol_n,pol_tau)
      pol_pi = (one-pol_tau)*pol_y - wage*pol_n - rrate*pol_k
      pol_t  = tau*pol_pi + pol_tau*pol_y
      pol_pi = (one-tau)*pol_pi

      ! mass of firms and equilibrium
      Gz = zero ; where (pol_pi.gt.wage) Gz = pz
      M  = sum(Gz)
      aL = one - M

      ! labor market clearing condition
      equilw = sum(pol_n*Gz) - aL

      ! update wage rate
      if (equilw.lt.zero) then
        maxw = 0.5d0*(wage+maxw)
      else
        minw = 0.5d0*(wage+minw)
      end if

    end do

    ! ------------------------------------------------------------
    ! find the equilibrium mass of firms

    M  = sum(Gz)
    aL = one - M

    ! ------------------------------------------------------------
    ! statistics and aggregate variables

    s_lab%mean%m = varmean(    pol_n           , w = Gz )
    s_lab%std%m  = varstd(     log(pol_n)      , w = Gz )
    s_lab%p10%m  = percentile( pol_n  , 0.10d0 , w = Gz )
    s_lab%p25%m  = percentile( pol_n  , 0.25d0 , w = Gz )
    s_lab%p50%m  = percentile( pol_n  , 0.50d0 , w = Gz )
    s_lab%p75%m  = percentile( pol_n  , 0.75d0 , w = Gz )
    s_lab%p90%m  = percentile( pol_n  , 0.90d0 , w = Gz )
    shr_firms%m  = cien*sum(Gz      ,mask=pol_n.ge.50.d0)/sum(Gz)
    shr_labor%m  = cien*sum(pol_n*Gz,mask=pol_n.ge.50.d0)/sum(pol_n*Gz)

    ! ------------------------------------------------------------
    ! aggregate variables

    M  = sum(Gz)                ! mass of firms
    aK = sum(Gz*pol_k)          ! aggregate capital
    aN = sum(Gz*pol_n)          ! aggregate labor
    aP = sum(Gz*pol_pi)         ! aggregate profits
    aT = sum(Gz*pol_t)          ! aggregate tax revenues
    aY = sum(Gz*pol_y)          ! aggregate output
    aC = wage*aN + aP + aT      ! aggregate consumption
    aL = aN                     ! aggregate labor supply

    ! aggregate TFP
    zaux1 = varmean( (((one - pol_tau)**gamma)*exp(zgrid))**(one/(one-gamma)) , w = Gz )
    zaux2 = varmean( (((one - pol_tau)       )*exp(zgrid))**(one/(one-gamma)) , w = Gz )**gamma
    aTFP  = ( M**(one-gamma) )*zaux1/zaux2
    aZ    = ( aTFP / ( M**(one-gamma) )  )**(one/(one-gamma))
    
    ! efficient TFP – given aggregate inputs and number of firms
    eZ   = varmean(exp(zgrid/(one-gamma)),w=Gz)
    eTFP = ( eZ**(one-gamma) )*( M**(one-gamma) ) 

    ! measured TFP as an econometrician
    lshare = wage*aN/aY
    mTFP   = aY/( (aK**(one-lshare))*(aN**lshare) )

    ! TFP as in Guner et al
    nTFP = aY/( (aK**(gamma*alpha))*(aN**(one-(gamma*alpha))) )

    return
  end subroutine equilibrium

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  elemental function kstar(z,taui) result(k)
    implicit none
    real(rp) :: k
    real(rp) , intent(in) :: taui,z
    if (withk.eq.1) then
      k = (((one - taui)*exp(z))**(one/(one-gamma)))*&
          gamma*(alpha/rrate)*&
          ((gamma*((alpha/rrate)**alpha)*(((one-alpha)/wage)**(one-alpha)))**(gamma/(one-gamma)))
    else
      k = zero
    end if
    return
  end function kstar
  elemental function nstar(z,taui) result(n)
    implicit none
    real(rp) :: n
    real(rp) , intent(in) :: taui,z
    if (withk.eq.1) then
      n = (((one - taui)*exp(z))**(one/(one-gamma)))*&
          gamma*((one-alpha)/wage)*&
          ((gamma*((alpha/rrate)**alpha)*(((one-alpha)/wage)**(one-alpha)))**(gamma/(one-gamma)))
    else
      n = (((one - taui)*exp(z))**(one/(one-gamma)))*&
          ((gamma/wage)**(one/(one-gamma)))
    end if
    return
  end function nstar
  elemental function output(z,k,n,taui) result(y)
    implicit none
    real(rp) :: y
    real(rp) , intent(in) :: z,k,n,taui
    if (withk.eq.1) then
      y = (one-taui)*exp(z)*( k**(gamma*alpha) )*( n**(gamma*(one-alpha)) )
    else
      y = exp(z)*(n**gamma)
    end if
    return
  end function output
  elemental function profits(z,taui) result(pi)
    implicit none
    real(rp) :: y,k,n,pi
    real(rp) , intent(in) :: z,taui
    if (withk.eq.1) then
      n  = nstar(z,taui)
      k  = kstar(z,taui)
      pi = (one-tau)*( (one-taui)*output(z,k,n,taui) - wage*n - rrate*k )
    else
      n  = nstar(z,taui)
      pi = (one-tau)*( (one-taui)*output(z,zero,n,taui) - wage*n )
    end if
    return
  end function profits

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
        call random_number(x1) ; x0 = x0 + 0.1d0*(x1 - 0.5d0)*abs(x0)
      end if
      if (withk.eq.0) then
        call lmmin(funcmodel,x1(1:size(x0)-1),y1(1:size(x0)-1),numiter,exitcode,x0(1:size(x0)-1),iprint=2) 
        !call simplex(sumfuncmodel,x1(1:size(x0)-1),toty,numiter,exitcode,x0(1:size(x0)-1),iprint=2)
      else
        call lmmin(funcmodel,x1,y1,numiter,exitcode,x0,usebro=0,iprint=2) ; x0 = x1
        call simplex(sumfuncmodel,x1,toty,numiter,exitcode,x0,iprint=2)
      end if
      x0 = x1
    end do
    do ipp = 1,500
      call random_number(x1)
      x0   = bestx(1:size(x0)) + 0.1d0*(x1 - 0.5d0)*abs(bestx(1:size(x0)))
      toty = sumfuncmodel(x0)
      print * , ipp , toty , minermodel
    end do
    return
  end subroutine calibmodel
  function funcmodel(xp) result(resid)
    use toolkit , only : denormalize
    implicit none
    real(rp)               :: xp(:)
    real(rp) , allocatable :: resid(:) ; allocate(resid(size(xp)+0))

    if (size(xp).le.3) xi    = exp(xp(1))
    if (size(xp).le.3) gamma = exp(xp(2))/xi
    if (size(xp).eq.3) rrate = min(0.10d0,exp(xp(3)))

    call equilibrium( )

    if (size(xp).le.3) resid(1) = 5.0d0*( one - s_lab%mean%m / s_lab%mean%d )
    if (size(xp).le.3) resid(2) = 1.0d0*( one - shr_labor%m  / shr_labor%d  ) !
    if (size(xp).eq.3) resid(3) = 1.0d0*( one - (aK/aY) / 3.0d0 )

    !if (size(xp).le.3) resid(2) = 1.0d0*( one - (s_lab%p10%m/s_lab%mean%m) / (s_lab%p10%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(3) = 1.0d0*( one - (s_lab%p25%m/s_lab%mean%m) / (s_lab%p25%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(4) = 1.0d0*( one - (s_lab%p50%m/s_lab%mean%m) / (s_lab%p50%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(5) = 1.0d0*( one - (s_lab%p75%m/s_lab%mean%m) / (s_lab%p75%d/s_lab%mean%d)  )
    !if (size(xp).le.3) resid(6) = 1.0d0*( one - (s_lab%p90%m/s_lab%mean%m) / (s_lab%p90%d/s_lab%mean%d)  )

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
    character(len=100)              :: filename
    integer                         :: i,unitt ; unitt = 0
    filename="results/"//date(3:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//time(1:2)//'-'//time(3:4)
    if (present(unitnum)) unitt = unitnum
    if (unitt.eq.1) open(unit=unitt,file=trim(adjustl(filename))//"_lucas.txt",action='write',status="replace")
    if (unitt.eq.2) open(unit=unitt,file=trim(adjustl(filename))//"_lucas_calib.txt",action='write',status="replace")
    if (unitt.eq.3) open(unit=unitt,file=trim(adjustl(filename))//"_lucas_calib_0.txt",action='write',status="replace")
      write(unitt,90) '                      '
      write(unitt,90) '  parameters:         '
      write(unitt,90) '  ----------          '
      write(unitt,33) dble(withk),'withk'
      write(unitt,33) rrate,'rrate'
      write(unitt,33) xi,'xi'
      write(unitt,33) uz,'uz'
      write(unitt,33) sigma,'sigma'
      write(unitt,33) gamma,'gamma'
      write(unitt,90) '  '
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Equilibrium:     '
      write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
      write(unitt,20) '  ' , 'Wage rate        ' , wage
      write(unitt,20) '  ' , 'Share managers   ' , M*cien
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Productivity:    '
      write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
      write(unitt,20) '  ' , 'Z                ' , aZ   , cien*aZ/eZ - cien
      write(unitt,20) '  ' , 'TFP              ' , aTFP , cien*aTFP/eTFP - cien
      write(unitt,20) '  ' , 'TFP as in Guner  ' , nTFP , cien*nTFP/aTFP - cien
      write(unitt,20) '  ' , 'TFP as in Solow  ' , mTFP , cien*mTFP/aTFP - cien
      write(unitt,20) '  ' , 'Output per worker' , aY/aL
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Distortions: '
      write(unitt,90) '  ' , '^^^^^^^^^^^^ '
      write(unitt,20) '  ' , 'Average      ' , varmean( pol_tau , w = Gz )
      write(unitt,20) '  ' , 'Small firms  ' , varmean( pol_tau , w = Gz , mask=pol_n.lt.s_lab%p10%m )
      write(unitt,20) '  ' , 'Large firms  ' , varmean( pol_tau , w = Gz , mask=pol_n.gt.s_lab%p90%m )
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Aggregates: '
      write(unitt,90) '  ' , '^^^^^^^^^^^ '
      write(unitt,20) '  ' , 'Output      ' , aY
      write(unitt,20) '  ' , 'Capital     ' , aK , aK/aY
      write(unitt,20) '  ' , 'Labor       ' , aL 
      write(unitt,20) '  ' , 'Consumption ' , aC
      write(unitt,20) '  ' , 'Government  ' , aT
      write(unitt,90) '  '
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
      write(unitt,90) '  ' ,  repeat('-',41)
      write(unitt,10) '  ' , 'Capital: average ' , s_kap%mean%m , s_kap%mean%d , s_kap%mean%m/s_kap%mean%d - one
      write(unitt,90) '  ' ,  repeat('*',41)
      write(unitt,90) '  '
      write(unitt,90) '  '
      if (unitt.gt.0) close(unitt)

    return
    90 format (a,a)
    10 format (a,a,10(f8.2))
    20 format (a,a,3(f12.6))
    33 format (f12.8,' ! ',a)
  end subroutine writelucas

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module lucas