program ising
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)
!********************************************************************
!************************* DEFINING OBJECTS *************************
!********************************************************************
  !!! Deines constants !!!
  real (kind=dp), parameter :: kB = 1.38064852E-23
  real (kind=dp), parameter :: h = 0.01

  !!! Defines arrays !!!
  real (kind=dp), dimension(0:4,2) :: probArray
  real (kind=dp), dimension(:,:), allocatable :: isingGrid

  !!! Defines variables
  real (kind=dp) :: temp
  integer :: isingWidth, isingHeight, i, j, spinSi

  !!! Deines the functions names and types !!!
  real (kind=dp) :: Probability, Energy

!********************************************************************
!************************** MAIN PROGRAMME **************************
!********************************************************************

  isingWidth = 3
  isingHeight = 3
  allocate(isingGrid(isingWidth,IsingHeight))

  do i = 1, isingWidth
    do j = 1, isingHeight
      call CreateIsingGrid
      isingGrid(3,3) = -1
      call CheckSpin
    end do
  end do

!********************************************************************
!******************** FUNCTIONS AND SUBROUTINES *********************
!********************************************************************
contains
  !function Energy()
    !!! Function to calculate Energy of state at current lattice point !!!
    !Energy = -Si*sum(Sj)-hSi
  !end function

  !function Probability()
    !!! Function to calculate probability of spin flip using Energy function !!!
    !Probability = exp(-(Energy(new)-Energy(old))/(kB*temp))
  !end function

  subroutine CreateIsingGrid
    !!! Subroutine to create the Ising Grid system of lattice points, each with a magnetic spin of 1 !!!
    isingGrid(i,j) = 1
  end subroutine

  subroutine CheckSpin
    !!! Subroutine to check what the spin of the current lattice point has (Si) !!!
    if (isingGrid(i,j) ==  1) then
      spinSi = 1
    else if (isingGrid(i,j) == -1) then
      spinSi = -1
    else
      print *, 'ERROR: magnetic spin at position', i, j, 'is not equal to 1 or -1'
    end if

    print *, spinSi
  end subroutine

  !subroutine PreCalcProbs
    !!! Subroutine to precalculate the probabilities of flipping spin !!!
   ! integer :: d

  !  do d = 0, 4
  !    if (spinSi == 1) then
  !      probArray(d,1) = Probability()
  !    else if (spinSi == -1) then
 !       probArray(d,2) = Probability()
     ! end if
    !end do
  !end subroutine
end program ising

!!!!!!!!!!!!!
!!! NOTES !!!
!!!!!!!!!!!!!
! – Maybe move allocation of isingGrid to subroutine CreateIsingGrid (call CreateIsingGrid be removed from nested loop in main programme)
