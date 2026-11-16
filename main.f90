

! ***********************************************************************************
! this program controls the execution of the codes to solve the model in the paper
! "Firm Entry and the TFP Paradox of Productivity-Dependent Distortions", by
! M. Almunia, JF Jimeno, D. Lopez-Rodriguez, and B. Petit.
!
! Menu: (1) economy, (2) exercise, (3) profit tax rate and, in the Hopenhayn economy,
! (4) exogenous exit probability. The answers set the names of all the files:
!
!   parameters (input):  txtfiles/calib_<economy>_<version>.txt
!   calibration:         results/<date>_calib_<economy>_<version>.txt
!   solution:            results/<date>_solve_<economy>_<version>.txt
!   experiments:         results/<economy>_<version>.txt
!                        results/<economy>_<sensitivity exercise>.txt
!
! where <economy> is hopenhayn_tau<tau>_delta<delta> or lucas_tau<tau> (tau and delta in %),
! <version> is baseline, nocapital or costs_labor, and the sensitivity exercises are
! gamma_up5, gamma_down5, sigma0_up10, sigma0_down10, frisch_0-01, frisch_3 (Hopenhayn) and
! gamma_up5, gamma_down5, xi_up10, xi_down10 (Lucas). To use a calibration as input, copy it
! to txtfiles/ and remove the date from its name.
! ***********************************************************************************

program main

  use omp_lib
  use parameters
  use toolkit    , only : num2text
  use hugo       , only : solve_hugo
  use lucas      , only : solve_lucas

  implicit none

  logical  :: exists
  real(rp) :: timing0,timing1
  integer  :: which,what,itau,idelta

  ! ***********************************************************************************
  ! set working directory and number of threads

  ! if the program is run from a folder that contains the code, work there
  threads = omp_get_max_threads( )
  inquire(file="main.f90",exist=exists)
  if (exists) then
    path = "./" ; goto 99
  end if

  ! working on mac laptop
  path = "/Users/borjapetit/Library/CloudStorage/Dropbox/jobmarket/jmp/code/"
  inquire(file=trim(adjustl(path))//"main.f90",exist=exists)
  if (exists) goto 99

  ! working on my office's computer
  path = "/mnt/c/Users/bpetit/Dropbox/jobmarket/jmp/code/"
  inquire(file=trim(adjustl(path))//"main.f90",exist=exists)
  if (exists) goto 99

  99 continue

  call system('rm -f *.mod')
  call system('rm -f *.log')

  call chdir(trim(adjustl(path)))
  call system('mkdir -p results txtfiles')

  ! ***********************************************************************************
  ! print header

  call date_and_time(date=date,time=time)

  write(*,fmt="(a)",advance="yes") '  '
  write(*,fmt="(a)",advance="yes") ' # ************************************************************************ # '
  write(*,fmt="(a)",advance="yes") ' # Firm Entry and the TFP Paradox of Productivity-Dependent Distortions     # '
  write(*,fmt="(a)",advance="yes") ' # M. Almunia, JF Jimeno, D. Lopez-Rodriguez, B. Petit                      # '
  write(*,fmt="(a)",advance="yes") ' # November 2026                                                            # '
  write(*,fmt="(a)",advance="yes") ' # ************************************************************************ # '
  write(*,fmt="(a)",advance="yes") '  '
  write(*,fmt="(a)",advance="yes") ' date: '//date(7:8)//'/'//date(5:6)//'/'//date(1:4)
  write(*,fmt="(a)",advance="yes") ' time: '//time(1:2)//':'//time(3:4)
  write(*,fmt="(a)",advance="yes") '  '
  write(*,fmt="(a)",advance="yes") ' threads: '//num2text(threads)
  write(*,fmt="(a)",advance="yes") '  '

  ! -----------------------------------------------------------
  ! 1. economy

  write(*,fmt="(a)",advance="yes") ' which model do you want to work with?      '
  write(*,fmt="(a)",advance="yes") '                                            '
  write(*,fmt="(a)",advance="yes") '   [ 1 ] Hopenhayn (1992)                   '
  write(*,fmt="(a)",advance="yes") '   [ 2 ] Lucas (1978)                       '
  write(*,fmt="(a)",advance="yes") '                                            '
  write(*,fmt="(a)",advance="no" ) ' your choice: '; read (*,*) which
  write(*,fmt="(a)",advance="yes") '                                    '
  if (which.lt.1 .or. which.gt.2) then
    write(*,fmt="(a)",advance="yes") ' invalid choice' ; stop
  end if

  ! -----------------------------------------------------------
  ! 2. exercise

  write(*,fmt="(a)",advance="yes") ' what do you want to do?                              '
  write(*,fmt="(a)",advance="yes") '                                                      '
  write(*,fmt="(a)",advance="yes") '   [ 1 ] solve the baseline model                     '
  write(*,fmt="(a)",advance="yes") '   [ 2 ] calibrate the baseline model                 '
  write(*,fmt="(a)",advance="yes") '   [ 3 ] symmetric experiment                         '
  write(*,fmt="(a)",advance="yes") '   [ 4 ] sensitivity of the symmetric experiment      '
  write(*,fmt="(a)",advance="yes") '   [ 5 ] solve the model without capital              '
  write(*,fmt="(a)",advance="yes") '   [ 6 ] calibrate the model without capital          '
  if (which.eq.1) then
  write(*,fmt="(a)",advance="yes") '   [ 7 ] solve the model with costs in labor units    '
  write(*,fmt="(a)",advance="yes") '   [ 8 ] calibrate the model with costs in labor units'
  end if
  write(*,fmt="(a)",advance="yes") '                                                      '
  write(*,fmt="(a)",advance="no" ) ' your choice: '; read (*,*) what
  write(*,fmt="(a)",advance="yes") '                                    '
  if (what.lt.1 .or. what.gt.8 .or. (which.eq.2 .and. what.gt.6)) then
    write(*,fmt="(a)",advance="yes") ' invalid option for this model' ; stop
  end if

  ! -----------------------------------------------------------
  ! 3. profit tax rate

  write(*,fmt="(a)",advance="yes") ' profit tax rate of the baseline economy?   '
  write(*,fmt="(a)",advance="yes") '                                            '
  write(*,fmt="(a)",advance="yes") '   [ 1 ] 30%                                '
  write(*,fmt="(a)",advance="yes") '   [ 2 ] 0%                                 '
  write(*,fmt="(a)",advance="yes") '                                            '
  write(*,fmt="(a)",advance="no" ) ' your choice: '; read (*,*) itau
  write(*,fmt="(a)",advance="yes") '                                    '
  select case (itau)
  case (1) ; taubase = 0.30_rp
  case (2) ; taubase = zero
  case default
    write(*,fmt="(a)",advance="yes") ' invalid choice' ; stop
  end select

  ! -----------------------------------------------------------
  ! 4. exogenous exit probability (Hopenhayn economy)

  noexit = 0 ; deltabase = zero
  if (which.eq.1) then
    write(*,fmt="(a)",advance="yes") ' exogenous exit probability (delta)?                                  '
    write(*,fmt="(a)",advance="yes") '                                                                      '
    write(*,fmt="(a)",advance="yes") '   [ 1 ] 0%                                                           '
    write(*,fmt="(a)",advance="yes") '   [ 2 ] 2%                                                           '
    write(*,fmt="(a)",advance="yes") '   [ 3 ] no endogenous exit (delta = exit rate in the data, no fixed cost)'
    write(*,fmt="(a)",advance="yes") '                                                                      '
    write(*,fmt="(a)",advance="no" ) ' your choice: '; read (*,*) idelta
    write(*,fmt="(a)",advance="yes") '                                    '
    select case (idelta)
    case (1) ; deltabase = zero
    case (2) ; deltabase = 0.02_rp
    case (3) ; noexit = 1          ! delta is set below, once the data moments are loaded
    case default
      write(*,fmt="(a)",advance="yes") ' invalid choice' ; stop
    end select
  end if

  ! -----------------------------------------------------------
  ! start timing
  timing0 = omp_get_wtime()

  ! -----------------------------------------------------------
  ! setup the problem

  call set_data( ) ; tau = taubase ; withk = 1 ; ecostinw = 0 ; fpop = 0

  if (noexit.eq.1) deltabase = s_e%d/cien      ! exit rate in the data

  if (which.eq.1) call set_economy("hopenhayn")
  if (which.eq.2) call set_economy("lucas")

  
  write(*,'(a)'                ) ' ------------------------------------------------- '
  if (which.eq.1) write(*,'(a)') ' economy:   Hopenhayn'
  if (which.eq.2) write(*,'(a)') ' economy:   Lucas'
  write(*,'(a)'                ) ' tax rate:  '//trim(num2text(cien*taubase,1))//"%"
  if (which.eq.1) write(*,'(a)') ' exit prob: '//trim(num2text(cien*deltabase,1))//"%"
  write(*,'(a)'                ) ' file name: '//trim(ecotag)
  write(*,'(a)'                ) ' ------------------------------------------------- '
  write(*,'(a)'                )      '  '

  ! -----------------------------------------------------------
  ! each solve_* routine sets the version of the model and reads its parameters (option "what")

  if (which.eq.1) call solve_hugo(what)      ! Hopenhayn  (hugo.f90)
  if (which.eq.2) call solve_lucas(what)     ! Lucas      (lucas.f90)

  ! compute and print the elapsed time
  timing1 = omp_get_wtime()
  write(*,*) '  '
  write(*,'(a,f10.2,a)') ' elapsed time: ',(timing1-timing0)/dble(60.000),' minutes'
  write(*,*) '  '
  write(*,*) ' press any key to finish the program'
  read( *,*)

end program main

! ***********************************************************************************
! columns of the results files of the experiments
!
!  1 atau      subsidy/tax rate (%)
!  2 wage      % change
!  3 TFP       level
!  4 D.TFP     % change in TFP
!  5 D.M       scale effect
!  6 D.Phi     misallocation effect
!  7 D.F       selection effect
!  8 E(t)      revenue-neutrality residual
!  9 E(t|e)
! 10 E(t|s)
! 11 E(t|l)
! 12 taxed
! 13 D.E       % change in the mass of firms
! 14 D.Y
! 15 D.K
! 16 D.N
! 17 D.L
! 18 D.C
! 19 D.TFP1
! 20 D.TFP2
! 21 D.TFP3    TFP with workers and entrepreneurs in the labor input
! 22 wcondi
! 23 tax_l     tax rate paid by firms above the cutoff (in %)
! ***********************************************************************************
