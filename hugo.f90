
module hugo

  use parameters
  use toolkit , only : lmmin,simplex,vect,varmean,varstd,percentile,grid,&
  interpolate,print_err,bisection,olsreg,num2text
  
  implicit none

  ! state space
  integer  , parameter        :: nz  = 10
  integer  , parameter        :: np  = 1000
  real(rp) , dimension(nz)    :: zgrid
  real(rp) , dimension(np)    :: pgrid
  real(rp) , dimension(np,nz) :: zpgrid,pe
  real(rp) , dimension(nz)    :: pz0
  real(rp) , dimension(nz,nz) :: pz
  real(rp) , dimension(np)    :: pp
  real(rp) , dimension(np,nz) :: Gpz
  real(rp) , dimension(np,nz) :: pzp0
  real(rp)                    :: taus,taul,zphi

  ! ------------------------------------------------------------
  ! model parameters
  ! ------------------------------------------------------------ 

  real(rp) :: minermodel
  real(rp) :: condi
  integer  :: ecostinw,elasticL,withk
  ! equilibrium objects: wage, mass of firms, mass of entrants, totallabor supply
  real(rp) :: wage,M,Me,Lbar
  ! costs of production: fixed and entry
  real(rp) :: fcost,ecost
  ! caPtal cost
  real(rp) :: beta,rrate
  ! technology
  real(rp) :: gamma,eprob
  ! productivity distribution
  real(rp) :: rhoz,stdz,muz0,stdz0,xi,uz
  ! distortions: level and curvature
  real(rp) :: phi0,phi1,atau
  ! disutility of labor supply: level and curvature
  real(rp) :: eta0,eta1
  ! aggregates
  real(rp) :: aK,aN,aL,aY,aC,aP,aT,aZ,aTFP,nTFP,mTFP,eZ,eTFP  ! distorted economy
  ! value and policy function
  real(rp) , dimension(np,nz) :: vfunc,vgues,pol_k,pol_n,pol_y,pol_pi,pol_e,pol_t,pol_tau 

  contains

  subroutine solve_hugo( mode )
    implicit none
    integer , intent(in)   :: mode
    character(len=100)     :: filename
    real(rp)               :: aux,vecr(14),vecr0(14)
    integer                :: ia,iy,ind,ifile,modsol
    real(rp)               :: M0,wage0,aK0,TFP0,par0
    real(rp) , allocatable :: dist0(:,:),distort0(:,:)

    allocate (dist0(np,nz),distort0(np,nz))
    
    modsol = mode

    ! load model settings
    open(unit=9,file="txtfiles/in_hugo_0.txt",action="read") 
      read(9,*) aux ; ecostinw = int(aux)
      read(9,*) aux ; elasticL = int(aux)
      read(9,*) aux ; withk    = int(aux)
    close(9)

    85 continue

    ! baseline economy
    if (withk.eq.1 .and. ecostinw.eq.0) then
      open(unit=9,file="txtfiles/in_hugo_base.txt",action="read")
    ! economy with entry cost in labor units
    else if (withk.eq.1 .and. ecostinw.eq.1) then
      open(unit=9,file="txtfiles/in_hugo_ew.txt",action="read")
    ! economy without capital
    else if (withk.eq.0 .and. ecostinw.eq.0) then
      open(unit=9,file="txtfiles/in_hugo_nok.txt",action="read")
    ! bad combination of model settings
    else
      write(*,*) ' Error: combination of parameters not implemented '
      read * , aux
    end if

    ! load parameter values
    read(9,*) rrate ; beta = one/(one + rrate)
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

    ! initialize some parameter values
    Lbar = one ; atau = zero ; phi0 = zero ; phi1 = zero

    if (modsol.eq.1) then

      call equilibrium( )
      call writehugo(1)
      call writehugo( )

    elseif (modsol.eq.2) then

      print_err = .false. ; ifile = 878

      write(*,'(/,a,/)') ' Revenue neutral distortions'
      filename = "base" ; call runexperiment( )

      print_err = .true.

    elseif (modsol.eq.3) then
    
      print_err = .false.
      minermodel = huge(one)
      call calibmodel( )
      call writehugo( )
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

      write(*,'(/,a,/)') ' Frisch elasticity = 0.01'
      par0 = eta1 ; eta1 = 100.0d0
      filename = "s_eta1_0-01" ; call runexperiment( )
      eta1 = par0

      write(*,'(/,a,/)') ' Frisch elasticity = 3.0'
      par0 = eta1 ; eta1 = one/3.0d0
      filename = "s_eta1_3-00" ; call runexperiment( )
      eta1 = par0

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

      if (withk.eq.1 .and. ecostinw.eq.0) then
        ecotype = 'base'
      else if (withk.eq.1 .and. ecostinw.eq.1) then
        ecotype = 'einw'
      else if (withk.eq.0 .and. ecostinw.eq.0) then
        ecotype = 'nok'
      end if

     ! setup file for storing results
     open(unit=ifile,&
     file="results/h_"//trim(adjustl(ecotype))//"_"//trim(adjustl(filename))//".txt",&
     action='write',status="replace")

      ! undistorted economy
      call equilibrium( ) ; call printexer(1) ; tfpp = aggTFP(M,pol_tau,Gpz)

      ! distorted economies
      do ia = 1,20 ; atau = dble(ia)/cien
        call bisection(findzphi,zphi,iy,ind,maxval(zpgrid),minval(zpgrid),iprint=0)
        condi = findzphi(zphi) ; call printexer( ) ; tfpp = aggTFP(M,pol_tau,Gpz)
      end do

      ! if tfp increases with 20% distortion, keep increasing it until TFP falls
      if (tfpp/TFP0.gt.0.95d0) then
        84 atau = atau + 0.01d0
        call bisection(findzphi,zphi,iy,ind,maxval(zpgrid),minval(zpgrid),iprint=0)
        condi = findzphi(zphi) ; call printexer(0) ; tfpp = aggTFP(M,pol_tau,Gpz)
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

      Expz     = sum( exp(zpgrid)*Gpz                               )/sum( Gpz )
      Etau     = sum( pol_y*pol_tau*Gpz                             )/sum( pol_y*Gpz                             )
      Etau0    = sum( pol_y*pol_tau*pzp0                            )/sum( pol_y*pzp0                            )
      Etaulow  = sum( pol_y*pol_tau*Gpz , mask=pol_n.lt.s_lab%p10%m )/sum( pol_y*Gpz , mask=pol_n.lt.s_lab%p10%m )
      Etauhigh = sum( pol_y*pol_tau*Gpz , mask=pol_n.gt.s_lab%p90%m )/sum( pol_y*Gpz , mask=pol_n.gt.s_lab%p90%m )
      
      vecr = (/ Etau,Etau0,Etaulow,Etauhigh,Expz,M,aY,aN,aN+M,aK,aY/aN,aC,mTFP,nTFP /)

      if (present(indic) .and. indic.eq.1) then  
        M0       = M
        dist0    = Gpz
        distort0 = zero
        wage0    = wage
        aK0      = aK
        vecr0    = vecr
        TFP0     = aggTFP(M,pol_tau,Gpz)
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
      T001 = aggTFP( M0 , distort0 , Gpz   )
      T110 = aggTFP( M  , pol_tau  , dist0 )
      T101 = aggTFP( M  , distort0 , Gpz   )
      T011 = aggTFP( M0 , pol_tau  , Gpz   )
      T111 = aggTFP( M  , pol_tau  , Gpz   )

      dM   = (1.d0/3.d0)*(T100 - T000) + (1.d0/6.d0)*((T110 - T010) + (T101 - T001)) + (1.d0/3.d0)*(T111 - T011)
      dphi = (1.d0/3.d0)*(T010 - T000) + (1.d0/6.d0)*((T110 - T100) + (T011 - T001)) + (1.d0/3.d0)*(T111 - T101)
      dF   = (1.d0/3.d0)*(T001 - T000) + (1.d0/6.d0)*((T101 - T100) + (T011 - T010)) + (1.d0/3.d0)*(T111 - T110)

      write(ifile,'(30(f9.3))') atau*cien,dwage,T000,T111,cien*(T111-T000)/T000,cien*dM/T000,cien*dphi/T000,cien*dF/T000,&
        vecr,condi,zphi
      write(*,'(30(f9.3))') atau*cien,dwage,T000,T111,cien*(T111-T000)/T000,cien*dM/T000,cien*dphi/T000,cien*dF/T000,&
        vecr,condi,zphi,aT

      return
    end subroutine printexer
    ! this function computes TFP given mass of firms, distortions and distribution of firms
    function aggTFP(mass,distort,dist) result (agtfp)
      implicit none
      real(rp) :: mass,distort(:,:),dist(:,:),agtfp
        agtfp = ( mass**(one-gamma) )* &
                ( varmean( vect((((one - distort)**gamma)*exp(zpgrid))**(one/(one-gamma))) , w = vect(dist) )) / &
                ( varmean( vect((((one - distort)       )*exp(zpgrid))**(one/(one-gamma))) , w = vect(dist) )**gamma )
      return
    end function aggTFP
    function findzphi(xp) result(resid)
      implicit none
      real(rp) :: xp
      real(rp) :: resid
      zphi = xp ; taul = atau ; taus = -atau ; call equilibrium( )
      resid = sum(pol_tau*pol_y*Gpz)*cien
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
      tfpa1 = aggTFP( M , pol_tau , Gpz )
      resid = one - tfpa1/tfpa0
      resid = resid*cien
      return
    end function findatau0
  end subroutine solve_hugo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! solve the model given some equilibrium prices (wage and mass of entrants)

  subroutine equilibrium( )
    use omp_lib
    use toolkit , only : grid,tauchen,olsreg,varmean,varstd,percentile,vect
    implicit none
    integer                :: iz,ip,iter
    real(rp)               :: equilw,maxw,minw
    real(rp)               :: test,coefs(2),zhat,ztil
    real(rp) , allocatable :: Gpz1(:,:),Gpz0(:,:),Gpze(:,:),vnp(:,:)

    allocate(Gpz1(np,nz),Gpz0(np,nz),Gpze(np,nz),vnp(np,nz))

    ! ------------------------------------------------------------
    ! set productivity grid and process
    
    ! grid for productivity
    zgrid = grid( stdz*6.0d0 , muz0 - max(stdz,stdz0)*6.0d0 , nz  , s=one )
    ! distribitopn of entrants
    call tauchen(zgrid , zero , muz0 - half*stdz0*stdz0 , stdz0 , nz , pz ) ; pz0 = pz(1,:)
    ! transition matrix for productivity
    call tauchen(zgrid,rhoz,-half*stdz*stdz,stdz,nz,pz) 
    ! distribution over permanent productivity
    call pareto(pgrid,pp,uz,xi,uz/xi,np)
    ! total firm-level productivity
    do iz = 1,nz
      zpgrid(:,iz) = zgrid(iz) + pgrid(:)
      pzp0(:,iz)   = pz0(iz)*pp(:)
    end do
  
    ! ------------------------------------------------------------
    ! distortions

    if (abs(taus)+abs(taul).lt.tol) then
      pol_tau = one - exp( phi0 + phi1*(zpgrid(ip,nz)-zpgrid) )
    else 
      where (zpgrid.gt.zphi) pol_tau = taul
      where (zpgrid.le.zphi) pol_tau = taus
    end if

    ! ------------------------------------------------------------
    ! find the equilibrium wage

    maxw = 100.0d0 ; minw = 0.0d0 ; equilw = huge(one)
    do while (abs(maxw-minw)*cien.gt.tol .and. abs(equilw)*cien.gt.tol)

      ! set the wage rate 
      wage  = 0.5d0*(maxw + minw)

      ! policy functions
      pol_n  = nstar(zpgrid,pol_tau)
      pol_k  = kstar(zpgrid,pol_tau)
      pol_y  = output(zpgrid,pol_k,pol_n,pol_tau)
      pol_pi = (one-pol_tau)*pol_y - wage*pol_n - rrate*pol_k 
      pol_t  = tau*pol_pi + pol_tau*pol_y
      pol_pi = (one-tau)*pol_pi - fcost

      ! converge the value function of firms
      iter = 0 ; test  = huge(one) ; vgues = zero ; vfunc = vgues ; ip = 1
      do while ( iter.lt.iterm .and. test.gt.tol ) ; iter  = iter + 1
        vfunc = vgues
        !$omp parallel num_threads(threads) default(shared)
        !$omp do schedule(dynamic)
        do ip=1,np
          vgues(ip,:) = pol_pi(ip,:) + (one-eprob)*beta*max(matmul(pz,vfunc(ip,:)),zero)
          !vgues(ip,:) = max(zero,pol_pi(ip,:) + (one-eprob)*beta*matmul(pz,vfunc(ip,:)))
        end do
        !$omp end do
        !$omp end parallel
        test  = maxval(abs(vgues-vfunc))
      end do
      vfunc = vgues
      do ip=1,np
        pol_e(ip,:) = eprob ; where(matmul(pz,vfunc(ip,:)).le.zero) pol_e(ip,:) = one
        !pol_e(ip,:) = zero ; where(pol_pi(ip,:)+(one-eprob)*beta*matmul(pz,vfunc(ip,:)).le.zero) pol_e(ip,:) = one
      end do
      pe = zero

      ! free-entry condition
      if (ecostinw.eq.1) then
        equilw = wage*ecost - sum(vfunc*pzp0)
      else
        equilw = ecost - sum(vfunc*pzp0)
      end if

      ! update wage rate
      if (equilw.gt.zero) then
        maxw = wage
      else
        minw = wage
      end if

    end do

    ! ------------------------------------------------------------
    ! find the equilibrium mass of firms

    ! stationary distribution for unit mass of entrants
    Gpz0 = pzp0 ; Gpz1 = zero
    do iter = 1,10000 ; Gpz0 = Gpz1
      Gpze = pzp0*(one-pe)
      !$omp parallel num_threads(threads) default(shared)
      !$omp do schedule(dynamic)
      do ip=1,np
        Gpz1(ip,:) = (one-pol_e(ip,:))*matmul(transpose(pz),Gpz0(ip,:))
      end do
      !$omp end do
      !$omp end parallel
      Gpz1 = Gpz1 + Gpze
      test = maxval(abs(vect(Gpz1 - Gpz0)))
      if (test.lt.tol) exit
    end do

    ! assume unit mass of firms (M=1)
    Me  = one/sum(Gpz0)
    Gpz = Gpz0*Me

    ! find aggregate labor demand and aggregate consumption under M=1
    if (ecostinw.eq.1) then
      aN = sum(Gpz*pol_n) + Me*ecost
      aC = wage*aN + sum(Gpz*pol_pi) + sum(Gpz*pol_t)
    else
      aN = sum(Gpz*pol_n)
      aC = wage*aN + sum(Gpz*pol_pi) + sum(Gpz*pol_t) - Me*ecost
    end if

    ! mass of firms that clears labor market
    if (elasticL.eq.1) then
      M = ( ((wage/(eta0*aC))**(one/eta1)) / aN )**(eta1/(1.d0+eta1))
    else
      M = Lbar/aN
    end if 

    ! corresponding mass of entrants
    Me  = M/sum(Gpz0)

    ! distribution of firms
    Gpz = Gpz0*Me

    ! ------------------------------------------------------------
    ! statistics and aggregate variables    

    s_lab%mean%m = varmean(    vect(pol_n)           , w = vect(Gpz) )
    s_lab%std%m  = varstd(     log(vect(pol_n))      , w = vect(Gpz) )
    s_lab%p10%m  = percentile( vect(pol_n)  , 0.10d0 , w = vect(Gpz) )
    s_lab%p25%m  = percentile( vect(pol_n)  , 0.25d0 , w = vect(Gpz) )
    s_lab%p50%m  = percentile( vect(pol_n)  , 0.50d0 , w = vect(Gpz) )
    s_lab%p75%m  = percentile( vect(pol_n)  , 0.75d0 , w = vect(Gpz) )
    s_lab%p90%m  = percentile( vect(pol_n)  , 0.90d0 , w = vect(Gpz) )

    s_lab0%mean%m = varmean(    vect(pol_n)           , w = vect(pzp0*(one-pol_e)) )
    s_lab0%std%m  = varstd(     log(vect(pol_n))      , w = vect(pzp0*(one-pol_e)) )
    s_lab0%p10%m  = percentile( vect(pol_n)  , 0.10d0 , w = vect(pzp0*(one-pol_e)) )
    s_lab0%p25%m  = percentile( vect(pol_n)  , 0.25d0 , w = vect(pzp0*(one-pol_e)) )
    s_lab0%p50%m  = percentile( vect(pol_n)  , 0.50d0 , w = vect(pzp0*(one-pol_e)) )
    s_lab0%p75%m  = percentile( vect(pol_n)  , 0.75d0 , w = vect(pzp0*(one-pol_e)) )
    s_lab0%p90%m  = percentile( vect(pol_n)  , 0.90d0 , w = vect(pzp0*(one-pol_e)) )
    
    shr_firms%m  = cien*sum(Gpz      ,mask=pol_n.ge.50.d0)/sum(Gpz)
    shr_labor%m  = cien*sum(pol_n*Gpz,mask=pol_n.ge.50.d0)/sum(pol_n*Gpz)

    vnp = zero
    do ip=1,np ; do iz=1,nz
      vnp(ip,iz) = sum(pz(iz,:)*pol_n(ip,:))
    end do ; end do

    call olsreg(coefs,log(vect(vnp)),log(vect(pol_n)),w=vect((one-pol_e)*Gpz),mask=vect(vnp.gt.0.10d0)) ; rho_n%m = coefs(2)
    
    ! ------------------------------------------------------------
    ! aggregate variables

    M  = sum(Gpz)          ! mass of firms
    aK = sum(Gpz*pol_k)    ! aggregate capital
    aN = sum(Gpz*pol_n)    ! aggregate labor
    aP = sum(Gpz*pol_pi)   ! aggregate profits
    aT = sum(Gpz*pol_t)    ! aggregate tax revenues
    aY = sum(Gpz*pol_y)

    ! consumption and aggregat elabro
    if (ecostinw.eq.1) then
      aL = aN + Me*ecost 
      aC = wage*aL + aP + aT
    else
      aL = aN
      aC = wage*aL + aP + aT - Me*ecost
    end if

    ! aggregate TFP
    zhat = varmean( (((one - vect(pol_tau))**gamma)*exp(vect(zpgrid)))**(one/(one-gamma)) , w = vect(Gpz) )
    ztil = varmean( (((one - vect(pol_tau))       )*exp(vect(zpgrid)))**(one/(one-gamma)) , w = vect(Gpz) )**gamma
    aTFP = ( M**(one-gamma) )*zhat/ztil
    aZ   = ( aTFP / ( M**(one-gamma) )  )**(one/(one-gamma))
    
    ! measured TFP as an econometrician
    mTFP = aY / ( ( aK**(one-(wage*aN/aY)) ) * ( aN**(wage*aN/aY) ) )

    ! TFP as in Guner et al
    nTFP = aY / ( ( aK**(gamma*alpha) ) * ( aN**(one-(gamma*alpha)) ) )

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

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine calibmodel( )
    use toolkit , only : simplex,lmmin,normalize,denormalize
    implicit none
    real(rp) , dimension(8)  :: x0,x1
    real(rp) , dimension(8)  :: y1
    real(rp)                 :: toty
    integer                  :: numiter,exitcode,ipp
    x0 = (/ gamma,xi,stdz,-muz0,stdz0,ecost,fcost,rrate /)
    x0 = log(x0) 
    do ipp = 0,size(x0)
      if (ipp.gt.0) then
        call random_number(x1)
        x0 = x0 + 0.1d0*(x1 - 0.5d0)*x0
      end if
      if (withk.eq.0) then
        ! call lmmin(funcmodel,x1(1:size(x0)-1),y1(1:size(x0)-1),numiter,exitcode,x0(1:size(x0)-1),shock=0.05d0,iprint=2,usebro=1)
        ! y1(1:8) = funcmodel(x1(1:size(x0)-1))
        call simplex(sumfuncmodel,x1(1:size(x0)-1),toty,numiter,exitcode,x0(1:size(x0)-1),iprint=2)
      else
        call lmmin(funcmodel,x1,y1,numiter,exitcode,x0,shock=0.05d0,iprint=2,usebro=1) ; x0 = x1
        call simplex(sumfuncmodel,x1,toty,numiter,exitcode,x0,iprint=2)
      end if
      x0 = x1
    end do
    return
  end subroutine calibmodel
  function funcmodel(xp) result(resid)
    use toolkit , only : denormalize
    implicit none
    real(rp)               :: xp(:)
    real(rp) , allocatable :: resid(:) ; allocate(resid(size(xp)))

    if (size(xp).le.8) gamma =  exp(xp(1)) 
    if (size(xp).le.8) xi    =  exp(xp(2)) 
    if (size(xp).le.8) stdz  =  exp(xp(3))
    if (size(xp).le.8) muz0  = -exp(xp(4))
    if (size(xp).le.8) stdz0 =  exp(xp(5))
    if (size(xp).le.8) ecost =  exp(xp(6))
    if (size(xp).le.8) fcost =  exp(xp(7))
    if (size(xp).eq.8) rrate =  exp(xp(8))

    call equilibrium( )

    if (size(xp).le.8) resid(1) = 5.0d0*( one - ( s_lab%mean%m     ) / ( s_lab%mean%d  ) ) ! gamma
    if (size(xp).le.8) resid(2) = 1.0d0*( one - ( shr_labor%m      ) / ( shr_labor%d   ) ) ! xi 

    if (size(xp).le.8) resid(3) = 1.0d0*( one - ( s_lab%std%m      ) / ( s_lab%std%d   ) ) ! stdz
    if (size(xp).le.8) resid(4) = 1.0d0*( one - ( s_lab0%mean%m    ) / ( s_lab0%mean%d ) ) ! muz0
    if (size(xp).le.8) resid(5) = 1.0d0*( one - ( s_lab0%p90%m     ) / ( s_lab0%p90%d  ) ) ! stdz0
    if (size(xp).le.8) resid(6) = 1.0d0*( one - ( s_lab0%p10%m     ) / ( s_lab0%p10%d  ) ) ! ecost
    if (size(xp).le.8) resid(7) = 4.0d0*( one - ( sum(pol_e*Gpz)/M ) / ( 0.08d0        ) ) ! fcost
    if (size(xp).eq.8) resid(8) = 1.0d0*( one - ( aK/aY            ) / ( 3.0d0         ) ) ! rrate


    if (sum(resid(:)*resid(:)).lt.minermodel) then
      minermodel = sum(resid(:)*resid(:))
      call writehugo(2)
      return
    end if

    return
  end function funcmodel
  function sumfuncmodel(xp) result(resid)
    implicit none
    real(rp) :: xp(:),res(size(xp)),resid
    res  = funcmodel(xp) ; resid = sum(res(:)*res(:))
    return
  end function sumfuncmodel

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine writehugo(unitnum)
    use toolkit , only : varmean,percentile,vect
    implicit none
    integer , intent(in) , optional :: unitnum
    character(len=100)              :: filename
    integer                         :: i,unitt ; unitt = 0
    filename="results/"//date(3:4)//'-'//date(5:6)//'-'//date(7:8)//'_'//time(1:2)//'-'//time(3:4)
    if (present(unitnum)) unitt = unitnum
    if (unitt.eq.1) open(unit=unitt,file=trim(adjustl(filename))//"_hugo.txt",action='write',status="replace")
    if (unitt.eq.2) open(unit=unitt,file=trim(adjustl(filename))//"_hugo_calib.txt",action='write',status="replace")
    if (unitt.eq.3) open(unit=unitt,file=trim(adjustl(filename))//"_hugo_calib_0.txt",action='write',status="replace")
      write(unitt,90) '                      '
      write(unitt,90) '  parameters:         '
      write(unitt,90) '  ----------          '
      write(unitt,33) dble(ecostinw),'ecostinw'
      write(unitt,33) dble(elasticL),'elasticL'
      write(unitt,33) dble(withk),'withk'
      write(unitt,33) rrate,'rrate'
      write(unitt,33) eprob,'eprob'
      write(unitt,33) fcost,'fcost'
      write(unitt,33) ecost,'ecost'
      write(unitt,33) xi,'xi'
      write(unitt,33) uz,'uz'
      write(unitt,33) rhoz,'rhoz'
      write(unitt,33) stdz,'stdz'
      write(unitt,33) muz0,'muz0'
      write(unitt,33) stdz0,'stdz0'
      write(unitt,33) gamma,'gamma'
      write(unitt,33) eta0,'eta0'
      write(unitt,33) eta1,'eta1'
      write(unitt,90) '  '
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Equilibrium:     '
      write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
      write(unitt,20) '  ' , 'Wage rate        ' , wage
      write(unitt,20) '  ' , 'Mass of entrants ' , M , cien*M/(aN+M)
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Productivity:    '
      write(unitt,90) '  ' , '^^^^^^^^^^^^^    '
      write(unitt,20) '  ' , 'Z                ' , aZ   
      write(unitt,20) '  ' , 'TFP              ' , aTFP 
      write(unitt,20) '  ' , 'TFP as in Guner  ' , nTFP 
      write(unitt,20) '  ' , 'TFP as in Solow  ' , mTFP 
      write(unitt,20) '  ' , 'Output per worker' , aY/aL
      write(unitt,90) '  '
      write(unitt,90) '  ' , 'Aggregates: '
      write(unitt,90) '  ' , '^^^^^^^^^^^ '
      write(unitt,20) '  ' , 'Output      ' , aY
      write(unitt,20) '  ' , 'Capital     ' , aK , aK/aY
      write(unitt,20) '  ' , 'Labor       ' , aL , aN + M
      write(unitt,20) '  ' , 'Consumption ' , aC
      write(unitt,20) '  ' , 'Government  ' , aT
      write(unitt,20) '  ' , 'Labor sh.   ' , wage*aL/(rrate*aK + aP + wage*aL + aT)
      write(unitt,20) '  ' , 'Capital sh. ' , (rrate*aK + aP)/(rrate*aK + aP + wage*aL + aT)
      write(unitt,20) '  ' , 'Exit rate   ' , cien*sum(pol_e*Gpz)/sum(Gpz)
      write(unitt,90) '  '
      write(unitt,90) '  '
      write(unitt,90) '  ' , repeat('*',41)
      write(unitt,90) '  ' , 'MOMENTS             Model    Data    Diff '
      write(unitt,90) '  ' ,  repeat('*',41)
      write(unitt,10) '  ' , 'Labor: average   ' , s_lab%mean%m , s_lab%mean%d , s_lab%mean%m/s_lab%mean%d - one
      write(unitt,10) '  ' , 'Labor: std logs  ' , s_lab%std%m  , s_lab%std%d  , s_lab%std%m/s_lab%std%d - one
      write(unitt,10) '  ' , 'Labor: persisten ' , rho_n%m      , rho_n%d      , rho_n%m/rho_n%d - one
      write(unitt,10) '  ' , 'Labor: p10       ' , s_lab%p10%m  , s_lab%p10%d  , s_lab%p10%m/s_lab%p10%d - one
      write(unitt,10) '  ' , 'Labor: p25       ' , s_lab%p25%m  , s_lab%p25%d  , s_lab%p25%m/s_lab%p25%d - one
      write(unitt,10) '  ' , 'Labor: p50       ' , s_lab%p50%m  , s_lab%p50%d  , s_lab%p50%m/s_lab%p50%d - one
      write(unitt,10) '  ' , 'Labor: p75       ' , s_lab%p75%m  , s_lab%p75%d  , s_lab%p75%m/s_lab%p75%d - one
      write(unitt,10) '  ' , 'Labor: p90       ' , s_lab%p90%m  , s_lab%p90%d  , s_lab%p90%m/s_lab%p90%d - one
      write(unitt,90) '  ' ,  repeat('-',41)
      write(unitt,10) '  ' , 'E Labor: average ' , s_lab0%mean%m , s_lab0%mean%d , s_lab0%mean%m/s_lab0%mean%d - one
      write(unitt,10) '  ' , 'E Labor: p10     ' , s_lab0%p10%m  , s_lab0%p10%d  , s_lab0%p10%m/s_lab0%p10%d - one
      write(unitt,10) '  ' , 'E Labor: p25     ' , s_lab0%p25%m  , s_lab0%p25%d  , s_lab0%p25%m/s_lab0%p25%d - one
      write(unitt,10) '  ' , 'E Labor: p50     ' , s_lab0%p50%m  , s_lab0%p50%d  , s_lab0%p50%m/s_lab0%p50%d - one
      write(unitt,10) '  ' , 'E Labor: p75     ' , s_lab0%p75%m  , s_lab0%p75%d  , s_lab0%p75%m/s_lab0%p75%d - one
      write(unitt,10) '  ' , 'E Labor: p90     ' , s_lab0%p90%m  , s_lab0%p90%d  , s_lab0%p90%m/s_lab0%p90%d - one
      write(unitt,90) '  ' ,  repeat('-',41)
      write(unitt,10) '  ' , 'Share firms n>50 ' , shr_firms%m   , shr_firms%d   , shr_firms%m/shr_firms%d - one
      write(unitt,10) '  ' , 'Labor in n>50    ' , shr_labor%m   , shr_labor%d   , shr_labor%m/shr_labor%d - one
      write(unitt,90) '  ' ,  repeat('*',41)
      write(unitt,90) '  '
      write(unitt,90) '  '
      if (unitt.gt.0) close(unitt)

    return
    90 format (a,a)
    10 format (a,a,10(f8.2))
    20 format (a,a,3(f12.6))
    33 format (f12.8,' ! ',a)
  end subroutine writehugo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module hugo
