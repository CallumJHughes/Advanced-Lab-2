program ising
  implicit none
  integer, parameter :: dp=selected_real_kind(15,300)

!********************************************************************
!************************* DEFINING OBJECTS *************************
!********************************************************************
  !!! Defines constants !!!
  !real (kind=dp), parameter :: kB = 1
  !real (kind=dp), parameter :: h = 0.01

  !!! Defines arrays !!!
  real (kind=dp), dimension(0:4,2) :: probArray
  integer, dimension(:,:), allocatable :: isingGrid
  real (kind=dp), dimension(:), allocatable :: magData, energyData

  !!! Defines variables
  real (kind=dp) ::  probValue, magAvg, energyAvg, totalLatticeEnergy
  integer :: isingWidth, isingHeight, i, j, timestep, spinSi, downSpins, numTimeSteps, temp
  integer :: neighbour1, neighbour2, neighbour3, neighbour4, neighbourSum, sumSi
  logical :: flipAllowed

!********************************************************************
!************************** MAIN PROGRAMME **************************
!********************************************************************

  numTimeSteps = 1000
  allocate(magData(numTimeSteps))
  allocate(energyData(numTimeSteps))
  temp = 1

  ! Precalculate the probability array
  spinSi = 1
  call PreCalcProbs
  spinSi = -1
  call PreCalcProbs

  call CreateIsingGrid
  !isingGrid(3,4) = -1
  !isingGrid(1,5) = -1
  !isingGrid(5,3) = -1
  !isingGrid(3,4) = -1
  print *, 'initial ising: ', isingGrid

  timestep=0
  sumSi=0
  call SumLatticeSpins
  print *, sumSi
  call AssignMagData

  do timestep = 1, numTimeSteps
    sumSi = 0 ! Resets the sum of the spins of lattice points of the system
    !call SumLatticeSpins
    !print *, 'sum is: ', sumSi
    !call AssignMagData
    !call AssignEnergyData
    totalLatticeEnergy = 0
    do j = 1, isingHeight
      do i = 1, isingWidth
        !print *, Magnetisation(sumSi,isingWidth,isingHeight)
        call CheckSpin
        call SumNeighbourSpin
        call CheckNeighbourSum
        call FindProbability
        call FlipSpin
        call UpdateIsingGrid
        !call StateEnergy
        !print *, ' ' ! Creates gap between output statements for each step
      end do
    end do
    call SumLatticeSpins
    call AssignMagData
    !print *, 'total lattice energy is: ', totalLatticeEnergy
  end do

  print *, isingGrid

  !print *, magData

  print *, probArray

  call MagDataAvg

  print *, 'Mag Avg is: ', magAvg

  !call EnergyDataAvg

  !print *, 'Energy Avg is: ', energyAvg

!********************************************************************
!******************** FUNCTIONS AND SUBROUTINES *********************
!********************************************************************
contains
  function Energy(Si,neighbours)
    !!! Function to calculate Energy of state at current lattice point !!!
    real (kind=dp) :: Energy
    integer :: Si, neighbours
    real (kind=dp), parameter :: J = 1.0
    real (kind=dp), parameter :: h = -0.01

    Energy = -(J*Si*neighbours)-(h*Si)
  end function Energy

  function Probability(Si,neighbours,temp)
    !!! Function to calculate probability of spin flip using Energy function !!!
    real (kind=dp) :: Probability
    integer :: Si, neighbours, temp

    Probability = exp(-(Energy(Si*(-1),neighbours)-Energy(Si,neighbours))/temp)
  end function

  function Magnetisation(sumSi,width,height)
    !!! Function to calculate the Magnetisation of the system !!!
    real (kind=dp) :: Magnetisation
    integer :: sumSi, width, height

    Magnetisation = (1.0/(width*height))*(sumSi)
  end function

  subroutine CreateIsingGrid
    !!! Subroutine to create the Ising Grid system of lattice points, each with a magnetic spin of 1 !!!
    isingWidth = 6
    isingHeight = 6
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

    !print *, 'Current spin is: ', spinSi
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

    !print *, 'neighbour spins are: ' , neighbour1, neighbour2, neighbour3, neighbour4

    neighbourSum = neighbour1 + neighbour2 + neighbour3 + neighbour4

    !print *,'neighbour sum is: ' , neighbourSum
  end subroutine

  subroutine CheckNeighbourSum
    !!! Subroutine to check the neighbour sum and assign relevant number of down spins !!!
    !!! Used to in main programme !!!
    if (neighbourSum.EQ.4) then
      downSpins = 0
    else if (neighbourSum.EQ.2) then
      downSpins = 1
    else if (neighbourSum.EQ.0) then
      downSpins = 2
    else if (neighbourSum.EQ.-2) then
      downSpins = 3
    else if (neighbourSum.EQ.-4) then
      downSpins = 4
    end if

    !print *, 'number of down spins is: ', downSpins

  end subroutine

  subroutine CheckDownSpins
    !!! Subroutine to check the number of down spin neighbours and assign relavant neighbour sum !!!
    !!! Used to precalculate the probability array !!!
    if (downSpins.EQ.0) then
      neighbourSum = 4
    else if (downSpins.EQ.1) then
      neighbourSum = 2
    else if (downSpins.EQ.2) then
      neighbourSum = 0
    else if (downSpins.EQ.3) then
      neighbourSum = -2
    else if (downSpins.EQ.4) then
      neighbourSum = -4
    end if

    !print *, 'Neighbour sum is : ', neighbourSum
  end subroutine

  subroutine PreCalcProbs
    !!! Subroutine to precalculate the probabilities of flipping spin !!!
    do downSpins = 0, 4
      if (spinSi == 1) then
        call CheckDownSpins
        probArray(downSpins,1) = Probability(spinSi,neighbourSum,temp)
        !print *, neighbourSum
      else if (spinSi == -1) then
        call CheckDownSpins
        probArray(downSpins,2) = Probability(spinSi,neighbourSum,temp)
        !print *, neighbourSum
      end if
    end do
  end subroutine

  subroutine FindProbability
    !!! Subroutine to find the relavant probability value in the probability array !!!
    if (spinSi == 1) then
      probValue = probArray(downSpins,1)
    else if (spinSi == -1) then
      probValue = probArray(downSpins,2)
    end if

    !print *, 'Prob value: ', probValue
  end subroutine

  subroutine FlipSpin
    !!! Subroutine to check if a random number is less than the probability, if so allow flip !!!
    real (kind=dp) :: r

    call Random_Number(r)

    print *, 'Random number is: ', r

    if (r.LT.probValue) then
      spinSi = spinSi * (-1)
      flipAllowed = .True.
      !print *, 'Flip allowed'
    else
      spinSi = spinSi
      flipAllowed = .False.
      !print *, 'Flip not allowed'
    end if

    !print *, 'New Si value: ', spinSi

    !print *, isingGrid
  end subroutine

  subroutine UpdateIsingGrid
    !!! Subroutine to update the Ising grid if the spin flip was allowed !!!
    if (flipAllowed .eqv. .True.) then
      isingGrid(i,j) = spinSi
    end if

    !print *, isingGrid
  end subroutine

  subroutine SumLatticeSpins
    !!! Subroutine to sum the spins of all the lattice points in the system
    do j = 1, isingHeight
      do i = 1, isingWidth
        sumSi = sumSi + isingGrid(i,j)
        !print *, 'current spin is: ', isingGrid(i,j)
        !print *, 'sum is: ', sumSi
      end do
    end do
  end subroutine

  subroutine AssignMagData
    !!! Assigns values of Magnetisation to magData array for each time step !!!
    !!! Writes magnetisation values out to a data file alongside the current timestep !!!
    magData(timestep) = Magnetisation(sumSi,isingWidth,isingHeight)

    open(unit=1, file='magdata.dat')

    write(1,*) timestep, magData(timestep)
  end subroutine

  subroutine MagDataAvg
    !!! Subroutine to calculate the average of Magnetisation data !!!
    integer :: k
    real (kind=dp) :: magTot

    do k = 1, size(magData)
      magTot = magTot + magData(k)
    end do

    magAvg = magTot / size(magData)
  end subroutine

  subroutine AssignEnergyData
    !!! Calculates Energy of each lattice point and sums them to find Total Lattice Energy, then assigns values to energyData array for each time step !!!
    !!! Writes energy values alongside current timestep out to a data file !!!
    energyData(timestep) = totalLatticeEnergy

    print *, 'Energy is: ', energyData(timestep)

    open(unit=2, file='energydata.dat')

    write(2,*) timestep, energyData(timestep)
  end subroutine

  subroutine StateEnergy
    !!! Calculates energy of state at current lattice point !!!
    totalLatticeEnergy = totalLatticeEnergy + Energy(spinSi,neighbourSum)
  end subroutine

  subroutine EnergyDataAvg
    !!! Calculates the average of the energy data !!!
    integer :: k
    real (kind=dp) :: energyTot

    do k = 1, size(energyData)
      energyTot = energyTot + energyData(k)
    end do

    energyAvg = energyTot / size(energyData)
  end subroutine

end program ising

!!!!!!!!!!!!!
!!! NOTES !!!
!!!!!!!!!!!!!
! – Remember to remove all print statements in each subroutine used for testing
! – IOSTAT=err for input and what not
! – When all neighbour flips are up, the probability of changing to flip down is very high
! – Add code for user to choose lattice size and number of time steps?
! – Add code to AssignMagData subroutine to write out to mag data and time steps to a data file
! – Add subroutine for garbage collection
! – Add graph for Magnetisation vs Temp and include reference
! – Calculating total lattice energy???
! – If add modules, functions and subroutines start with module name
