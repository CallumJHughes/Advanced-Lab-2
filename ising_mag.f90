program ising_mag
  implicit none
  integer, parameter :: dp=real_selected_kind(15,300)
  !!! Deines constants !!!
  real (kind=dp), parameter :: kB

  !!! Defines arrays !!!
  real (kind=dp), dimension(0:4,2) :: probArray

  !!! Defines variables
  real (kind=dp) :: temp

  !!! Deines the functions names and types !!!
  real (kind=dp) :: Probability, Energy

  call CheckSpin

contains
  function Energy()
    !!! Function to calculate Energy of state at current lattice point !!!
  end function

  function Probability()
    !!! Function to calculate probability of spin flip using Energy function !!!
    Probability = exp(-(Energy(new)-Energy(old))/(kB*temp))
  end function

  subroutine preCalcProbs
    !!! Subroutine to precalculate the probabilities of flipping spin !!!
    integer :: d

    do d = 0, 4
      if (spinSi == 1) then
        probArray(d,1) = Probability()
      else if (spinSi == -1) then
        probArray(d,2) = Probability()
      end if
    end do
  end subroutine
end program ising_mag
