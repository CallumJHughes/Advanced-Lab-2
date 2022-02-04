program ising
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)

!********************************************************************
!************************* DEFINING OBJECTS *************************
!********************************************************************
  !!! Defines constants !!!
  real (kind=dp), parameter :: kB = 1.38064852E-23
  real (kind=dp), parameter :: h = 0.01

  !!! Defines arrays !!!
  real (kind=dp), dimension(0:4,2) :: probArray
  integer, dimension(:,:), allocatable :: isingGrid

  !!! Defines variables
  real (kind=dp) :: temp
  integer :: isingWidth, isingHeight, i, j, spinSi, downSpins
  integer :: neighbour1, neighbour2, neighbour3, neighbour4, neighbourSum

!********************************************************************
!************************** MAIN PROGRAMME **************************
!********************************************************************

  temp = 3

  spinSi = 1
  call PreCalcProbs
  spinSi = -1
  call PreCalcProbs

  print *, probArray

  call CreateIsingGrid
  isingGrid(2,1) = -1
  isingGrid(1,2) = -1
  isingGrid(3,2) = -1
  isingGrid(2,3) = -1
  print *, isingGrid

  !do j = 1, isingHeight
    !do i = 1, isingWidth
      !call CheckSpin
      !print *, 'Current spin is: ' ,spinSi
      !call SumNeighbourSpin
      !call CheckDownSpins
      !print *, 'energy is: ', Energy(spinSi)
      !print *, 'probability is: ', Probability(spinSi)
      !call PreCalcProbs
      !print *, ' '
    !end do
  !xend do

!********************************************************************
!******************** FUNCTIONS AND SUBROUTINES *********************
!********************************************************************
contains
  function Energy(Si)
    !!! Function to calculate Energy of state at current lattice point !!!
    real(kind=dp) :: Energy
    integer :: Si
    real(kind=dp), parameter :: J = 1
    real(kind=dp), parameter :: h = 0.01

    Energy = -(Si*neighbourSum)-(h*Si)
  end function Energy

  function Probability(Si)
    !!! Function to calculate probability of spin flip using Energy function !!!
    real(kind=dp) :: Probability
    integer :: Si

    Probability = exp(-(Energy(Si*(-1))-Energy(Si))/temp)
  end function

  subroutine CreateIsingGrid
    !!! Subroutine to create the Ising Grid system of lattice points, each with a magnetic spin of 1 !!!
    isingWidth = 3
    isingHeight = 3
    allocate(isingGrid(isingWidth,isingHeight))

    do j = 1, isingHeight
      do i = 1, isingWidth
        isingGrid(i,j) = 1
      end do
    end do
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
  end subroutine

  subroutine SumNeighbourSpin
    !!! Subroutine to sum the magnetic spins of the neighbouring lattice points !!!
    neighbourSum = 0

    if ((j-1).LT.1) then
      neighbour1 = isingGrid(i,isingHeight)
    else
      neighbour1 = isingGrid(i,j-1)
    end if

    if ((i+1).GT.isingWidth) then
      neighbour2 = isingGrid(1,j)
    else
      neighbour2 = isingGrid(i+1,j)
    end if

    if ((j+1).GT.isingHeight) then
      neighbour3 = isingGrid(i,1)
    else
      neighbour3 = isingGrid(i,j+1)
    end if

    if ((i-1).LT.1) then
      neighbour4 = isingGrid(isingWidth,j)
    else
      neighbour4 = isingGrid(i-1,j)
    end if

    print *, 'neighbour spins are: ' , neighbour1, neighbour2, neighbour3, neighbour4

    neighbourSum = neighbour1 + neighbour2 + neighbour3 + neighbour4

    print *,'neighbour sum is: ' , neighbourSum
  end subroutine

  subroutine CheckDownSpins
    !!! Subroutine to check the number of down spins in the set of neighbours !!!
    if (downSpins.EQ.0) then
      neighbourSum = 4
    else if (downSpins.EQ.1) then
      neighbourSum = 2
    else if (downSpins.EQ.2) then
      neighbourSum = 0
    else if (downSpins.EQ.3) then
      neighbourSum = -2
    else if (downSpins.EQ.4) then
      neighboursum = -4
    end if

    !print *, 'number of down spins is: ', downSpins
  end subroutine

  subroutine PreCalcProbs
    !!! Subroutine to precalculate the probabilities of flipping spin !!!
    do downSpins = 0, 4
      if (spinSi == 1) then
        call CheckDownSpins
        probArray(downSpins,1) = Probability(spinSi)
        print *, neighbourSum
      else if (spinSi == -1) then
        call CheckDownSpins
        probArray(downSpins,2) = Probability(spinSi)
        print *, neighbourSum
      end if
    end do
  end subroutine
end program ising

!!!!!!!!!!!!!
!!! NOTES !!!
!!!!!!!!!!!!!
! – Maybe move allocation of isingGrid to subroutine CreateIsingGrid (call CreateIsingGrid be removed from nested
!   loop in main programme)
! – Remember to remove all print statements in each subroutine used for testing
! – IOSTAT=err for input and what not
! – Will kB be used in Probability function?
! – When all neighbour flips are up, the probability of changing to flip down is very high
