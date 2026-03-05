


! ***********************************************************************************
! this program controls the execution of the codes to solve the model in the paper 
! "Size-Dependent Regulations, Tax Misreporting, and Aggregate Productivity", by 
! M. Almunia, JF Jimeno, D. Lopez-Rodriguez, and B. Petit.
! ***********************************************************************************

program main

  use omp_lib
  use parameters
  use hugo       , only : solve_hugo
  use lucas      , only : solve_lucas
  use toolkit    , only : num2text
  
  implicit none
  
  logical  :: exists
  real(rp) :: timing0,timing1,aux
  
  ! ***********************************************************************************
  ! set working directory and number of threads

  ! working on mac laptop
  path = "/Users/borjapetit/Library/CloudStorage/Dropbox/research/projects/2021_bunching/code_tfp/"
  inquire(file=trim(adjustl(path))//"main.f90",exist=exists)
  threads = omp_get_max_threads( )
  if (exists) goto 99
  
  ! working on my office's computer
  path = "/mnt/c/Users/bpetit/Dropbox/research/projects/2021_bunching/code_tfp/"
  inquire(file=trim(adjustl(path))//"main.f90",exist=exists)
  threads = omp_get_max_threads( )
  if (exists) goto 99
  
  99 continue
  
  call system('rm *.mod')
  call system('rm *.log')
  
  call chdir(trim(adjustl(path)))
  
  ! ***********************************************************************************
  ! print header

  call date_and_time(date=date,time=time)

  write(*,fmt="(a)",advance="yes") '  '
  write(*,fmt="(a)",advance="yes") ' # ************************************************************************ # '
  write(*,fmt="(a)",advance="yes") ' # Size-Dependent Regulations, Tax Misreporting, and Aggregate Productivity # '
  write(*,fmt="(a)",advance="yes") ' # M. Almunia, JF Jimeno, D. Lopez-Rodriguez, B. Petit                      # '
  write(*,fmt="(a)",advance="yes") ' # Mayo 2025                                                                # '
  write(*,fmt="(a)",advance="yes") ' # ************************************************************************ # '
  write(*,fmt="(a)",advance="yes") '  '
  write(*,fmt="(a)",advance="yes") ' date: '//date(7:8)//'/'//date(5:6)//'/'//date(1:4)
  write(*,fmt="(a)",advance="yes") ' time: '//time(1:2)//':'//time(3:4)
  write(*,fmt="(a)",advance="yes") '  '
  write(*,fmt="(a)",advance="yes") ' threads: '//num2text(threads)
  write(*,fmt="(a)",advance="yes") '                                    '
  write(*,fmt="(a)",advance="yes") ' solving the model. choose a mode:  '
  write(*,fmt="(a)",advance="yes") '                                    '
  write(*,fmt="(a)",advance="yes") '  model a la hopenhayn (1992)       '
  write(*,fmt="(a)",advance="yes") '  ^^^^^^^^^^^^^^^^^^^^^^^^^^^       '
  write(*,fmt="(a)",advance="yes") '   [ 11 ] model                     '
  write(*,fmt="(a)",advance="yes") '   [ 12 ] experiment                '
  write(*,fmt="(a)",advance="yes") '   [ 13 ] calibrate                 '
  write(*,fmt="(a)",advance="yes") '   [ 14 ] sensitivity               '
  write(*,fmt="(a)",advance="yes") '                                    '
  write(*,fmt="(a)",advance="yes") '  model a la lucas (1976)           '
  write(*,fmt="(a)",advance="yes") '  ^^^^^^^^^^^^^^^^^^^^^^^           '
  write(*,fmt="(a)",advance="yes") '   [ 21 ] model                     '
  write(*,fmt="(a)",advance="yes") '   [ 22 ] experiment                '
  write(*,fmt="(a)",advance="yes") '   [ 23 ] calibrate                 '
  write(*,fmt="(a)",advance="yes") '   [ 24 ] sensitivity               '
  write(*,fmt="(a)",advance="yes") '                                    '
  write(*,fmt="(a)",advance="no" ) ' your choice: '; read (*,'(i2)') solvemode
  write(*,fmt="(a)",advance="yes") '                                    '

  ! ***********************************************************************************
  ! setup the problem

  ! set parameter values
  call set_data( )
  
  ! start timing
  timing0 = omp_get_wtime()

  ! ***********************************************************************************
  ! solve the model according to user's choice

   ! solve the model with given prices
  if (solvemode.gt.10 .and. solvemode.lt.20) then
    call solve_hugo(solvemode-10)
  end if
  if (solvemode.gt.20 .and. solvemode.lt.30) then
    call solve_lucas(solvemode-20)
  end if

  91 continue

  ! compute and print the elapsed time
  timing1 = omp_get_wtime()
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) (timing1-timing0)/dble(60.000)
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) '  '
  write(*,*) ' press any key to finish the program'
  read( *,*) 

  return
end program main

! ***********************************************************************************










 