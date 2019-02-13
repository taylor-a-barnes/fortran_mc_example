FUNCTION calculate_energy_lj(coords,natoms)

  IMPLICIT DOUBLE PRECISION(A-Z)

  INTEGER :: natoms
  DOUBLE PRECISION :: coords(natoms,3)

  ! Boltzmann constant
  DOUBLE PRECISION :: boltzmann = 0.0019872041    ! kcal/(mol*Kelvin)

  ! lennard-jones parameters
  DOUBLE PRECISION :: eps = 10.2                  ! Kelvin
  DOUBLE PRECISION :: sigma = 2.28                ! Angstrom

  COMMON /energy/ boltzmann, eps, sigma

  calculate_energy_lj = 0.0

  ! loop over all unique atom pairs
  DO 10 i=1, natoms
     DO 1234 j=i+1, natoms

        ! calculate the distance between this pair
        dr = ( coords(i,1) - coords(j,1) ) ** 2 + &
             ( coords(i,2) - coords(j,2) ) ** 2 + &
             ( coords(i,3) - coords(j,3) ) ** 2
        dr = dr**(0.5)

        ! calculate the Lennard-Jones interaction between this pair
        s6 = ( sigma / dr )**6
        calculate_energy_lj = calculate_energy_lj + &
             4.0*eps*boltzmann*( s6**2 - s6 )

1234 CONTINUE
10 CONTINUE

END FUNCTION calculate_energy_lj
