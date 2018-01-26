PROGRAM mccode
  USE energy

  INTEGER :: i, j
  INTEGER :: iter
  INTEGER :: natoms
  INTEGER :: iout
  DOUBLE PRECISION, ALLOCATABLE :: coords(:,:)
  INTEGER :: atom
  DOUBLE PRECISION :: mxdplc = 0.1
  DOUBLE PRECISION :: displc(3)
  DOUBLE PRECISION :: e, eold
  DOUBLE PRECISION :: r
  DOUBLE PRECISION :: k
  DOUBLE PRECISION :: temp
  DOUBLE PRECISION :: prob

  ! init
  natoms = 10
  ALLOCATE( coords(natoms,3) )
  DO 10 i=1, natoms
     DO 11 j=1, 3
        coords(i,j) = 1.0d0*i
11 CONTINUE
10 CONTINUE

  ! calculate the energy of this configuration
  e = caleng(COORDS,NATOMS)
  eold = e

  ! init rand
  r = RAND(1)

  ! calculate the energy of this configuration
  DO 20 iter=1, 10

     WRITE(6,*)'Iteration: ',iter, e

     ! randomly select an atom
     r = RAND(0)
     atom = FLOOR((r*DBLE(natoms))+1.d0)

     ! randomly select a direction
     displc(1) = mxdplc*2.0*(RAND(0)-0.5)
     displc(2) = mxdplc*2.0*(RAND(0)-0.5)
     displc(3) = mxdplc*2.0*(RAND(0)-0.5)

     ! move the atom in this direction
     coords(atom,1) = coords(atom,1) + displc(1)
     coords(atom,2) = coords(atom,2) + displc(2)
     coords(atom,3) = coords(atom,3) + displc(3)

     ! calculate the energy of this configuration
     e = caleng(coords,natoms)

     ! accept of reject the trial step
     temp = 1.0 ! K
     k = 0.0019872041 ! kcal/(mol*K)
     prob = EXP(-(e-eold)/(k*temp))
     GOTO 50
40   eold=e
     GOTO 20
45   e = eold
     coords(atom,1) = coords(atom,1) - displc(1)
     coords(atom,2) = coords(atom,2) - displc(2)
     coords(atom,3) = coords(atom,3) - displc(3)
     GOTO 20
50   CONTINUE
     IF(RAND(0).lt.prob) GOTO 40
     GOTO 45

20 CONTINUE

  ! write the output
  iout = 31
  OPEN(newunit=iout, file="output.xyz", status="replace")
  WRITE(iout,*)natoms
  WRITE(iout,*)'Helium Cluster'
  DO 21 i=1, natoms
     WRITE(iout,*)"He ",coords(i,1),coords(i,2),coords(i,3)
21 CONTINUE

END PROGRAM mccode
