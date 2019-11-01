    
!|-------------------------------------------------|
!| /////////////////////////////////////////////// |
!| ///// T R A J E C T O R Y     M O D U L E ///// |
!| /////////////////////////////////////////////// |
!|-------------------------------------------------|
    
MODULE trajectories
IMPLICIT NONE

!$$$$$$$$$$ USER DERIVED TYPES $$$$$$$$$$

TYPE istream
    CHARACTER(50) :: filename
    INTEGER :: atoms
    CHARACTER(30) :: genbyvmd
    INTEGER :: dimen = 3
END TYPE istream

TYPE istrArr
    CHARACTER(3), ALLOCATABLE :: atomid(:)
    REAL, ALLOCATABLE :: coor(:,:)
END TYPE istrArr

TYPE :: ostream
    INTEGER :: molec, size
    CHARACTER(50) :: outname
END TYPE ostream

PUBLIC :: garray, printarray, centermass, ordpar

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!********** FUNCTIONS AND SUBROUTINES **********

    CONTAINS
    
    !#############################################
    !##### G A R R A Y     P R O C E D U R E #####
    !#############################################

    
    SUBROUTINE garray(inpt,atms,mol,sz,frm)
    IMPLICIT NONE

    ! Dummy args
    CHARACTER(50), INTENT(IN) :: inpt
    INTEGER, INTENT(IN) :: atms
    INTEGER, INTENT(OUT) :: mol, sz, frm
    
    ! Constructors
    TYPE(istream) :: garrayIn
    TYPE(istrArr) :: garrayInArr
    TYPE(ostream) :: garrayOut
    
    ! Subroutine vars
    INTEGER :: nxtFrm, duplTot
    REAL :: fstCoor
    INTEGER(KIND = 4) :: sumFrame = 0
    CHARACTER(LEN = 3) :: fstAtom
    INTEGER :: i
    
    garrayIn%filename = inpt
    garrayIn%atoms = atms
    
    !--------------------

    OPEN(20, FILE = garrayIn%filename)
    READ(20,*) garrayIn%atoms
    READ(20,'(A)') garrayIn%genbyvmd

    sumframe = sumframe + garrayIn%atoms

        DO
            ALLOCATE(garrayInArr%atomid(garrayIn%atoms), &
                    & garrayInArr%coor(garrayIn%atoms, garrayIn%dimen))

                DO i = 1, garrayIn%atoms

                    READ(20,*) garrayInArr%atomid(i), &
                    & garrayInArr%coor(i,1:garrayIn%dimen)
        
                    IF (i == 1) THEN
                        fstAtom = garrayInArr%atomid(i)
                        fstCoor = garrayInArr%coor(i,1)
                    ENDIF
                    IF ((fstAtom == garrayInArr%atomid(i)) .AND. (fstCoor /= garrayInArr%coor(i,1))) THEN
                        dupltot = i - 1
                    ENDIF  

                ENDDO
            DEALLOCATE(garrayInArr%atomid, garrayInArr%coor)

            garrayOut%size = garrayIn%atoms - dupltot

            READ(20,*) nxtfrm
            READ(20,*) garrayIn%genbyvmd
            IF (nxtfrm /= garrayIn%atoms) THEN
                EXIT
            ENDIF

            sumframe = sumframe + nxtfrm

        ENDDO

    frm = sumframe / garrayIn%atoms
    garrayOut%molec = garrayIn%atoms / garrayOut%size

    mol = garrayOut%molec
    sz = garrayOut%size



    WRITE(6,*) '========================================================'
    WRITE(6,*) 'CALVIN : The data has the following array dimensions... '
    WRITE(6,*) '        '
    WRITE(6,*) 'Simulation Frames: ', frm
    WRITE(6,*) 'Molecules per Frame: ', mol
    WRITE(6,*) 'Atoms per Molecule: ', sz
    WRITE(6,*) '========================================================'

    WRITE(6,*) '        '

    CLOSE(20)
    
    RETURN

    END SUBROUTINE garray

    !#####################################################
    !##### P R I N T A R R A Y     P R O C E D U R E #####
    !#####################################################
    
    SUBROUTINE printarray(inpt,frm,mol,sz)
    IMPLICIT NONE
    
    ! Dummy args
    CHARACTER(50), INTENT(IN) :: inpt
    INTEGER, INTENT(IN) :: frm, mol, sz
    
    ! Constructors
    TYPE(istream) :: parrayIn
    TYPE(istrArr) :: parrayInArr
    TYPE(ostream) :: parrayOut
    
    ! Subroutine vars
    INTEGER :: nxtfrm, frame = 1
    INTEGER :: i, j
     
    parrayIn%filename = inpt
    parrayOut%molec = mol
    parrayOut%size = sz
    
    !===============================================================================

    WRITE(6,*) '========================================================'
    WRITE(6,*) 'CALVIN : What name would like the output to be? '
    WRITE(6,*) '        '
    WRITE(6,*) '========================================================'
    READ(5,*) parrayOut%outname
    
    OPEN(21, FILE = parrayIn%filename)
    READ(21,*) parrayIn%atoms
    READ(21,'(A)') parrayIn%genbyvmd

    
    OPEN(51, file = parrayOut%outname, status = 'new')
    WRITE(51,*) 'frames', '     ', 'atom', '     ', 'xcoord', '        ', 'ycoord', '        ', 'zcoord'
    
    
        DO
            ALLOCATE(parrayInArr%atomid(parrayIn%atoms), &
            & parrayInArr%coor(parrayIn%atoms, parrayIn%dimen))
    
                DO i = 1, parrayOut%molec
                    DO j = 1, parrayOut%size
                        READ(21,*) parrayInArr%atomid(j), &
                        & parrayInArr%coor(j,1:parrayIn%dimen)

                        WRITE(51,*) frame, parrayInArr%atomid(j), &
                        & parrayInArr%coor(j,1:parrayIn%dimen) 
                    ENDDO
                ENDDO
    
            DEALLOCATE(parrayInArr%atomid, parrayInArr%coor)
    
            READ(21,*) nxtfrm
            READ(21,*) parrayIn%genbyvmd

            IF (nxtfrm /= parrayIn%atoms) THEN
                EXIT
            ENDIF

            frame = frame + 1
    
        ENDDO
    
    CLOSE(21)
    CLOSE(51)

    frame = 1
    
    RETURN
    
    END SUBROUTINE printarray

    !#####################################################
    !##### C E N T E R M A S S     P R O C E D U R E #####
    !#####################################################

    
    SUBROUTINE centermass(inpt, mol, sz)
    IMPLICIT NONE

    ! Dummy args
    CHARACTER(50), INTENT(IN) :: inpt
    INTEGER, INTENT(IN) :: mol, sz
    
    ! Constructors
    TYPE(istream) :: centerIn
    TYPE(istrArr) :: centerInArr
    TYPE(ostream) :: centerOut
    
    ! Subroutine vars
    INTEGER :: nxtfrm, frame = 1
    INTEGER :: i, j
    CHARACTER(2) :: atom1, atom2
    REAL , PARAMETER :: Cmass = 12.0107 , Hmass = 1.00784 , Nmass = 14.0067, Omass  = 15.999
    CHARACTER, PARAMETER :: Cstring = 'C' , Hstring = 'H' , Nstring = 'N' , Ostring = 'O'

    REAL, ALLOCATABLE :: mass(:), momenta(:,:), com(:,:)
    REAL :: masstot = 0, xtot = 0, ytot = 0, ztot = 0    
    
    centerIn%filename = inpt
    centerOut%molec = mol
    centerOut%size = sz

    !===================================================================================================

    WRITE(6,*) '========================================================'
    WRITE(6,*) 'CALVIN : What name would like the output to be? '
    WRITE(6,*) '        '
    WRITE(6,*) '========================================================'
    WRITE(6,*) '        '
    READ(5,*) centerOut%outname

    WRITE(6,*) '==============================================================================================='
    WRITE(6,*) 'CALVIN: find the center of mass between what 2 atoms? (separate your answers by a tab or space)'
    WRITE(6,*) 'Example: The center of mass between C1 and C2 can be entered as C1  C2'
    READ(5,*) atom1, atom2
    WRITE(6,*) 'CALVIN is crunching the numbers....'




    OPEN(22, FILE = centerIn%filename)
    READ(22,*) centerIn%atoms
    READ(22,'(A)') centerIn%genbyvmd


    OPEN(52, file = centerOut%outname, status = 'new')
    WRITE(52,*) 'frames', '     ', 'mol', '     ', 'xcom', '        ', 'ycom', '        ', 'zcom'


        DO
            ALLOCATE(centerInArr%atomid(centerIn%atoms), &
            & centerInArr%coor(centerIn%atoms, centerIn%dimen), &
            & mass(centerIn%atoms), &
            & momenta(centerIn%atoms, centerIn%dimen), &
            & com(centerOut%molec, centerIn%dimen))


                DO i = 1, centerOut%molec

                    DO j = 1, centerOut%size
                        READ(22,*) centerInArr%atomid(j), &
                        & centerInArr%coor(j,1:centerIn%dimen)
            
                        IF (INDEX(centerInArr%atomid(j) , Cstring) == 1) THEN
                            mass(j) = Cmass
                        ELSEIF (INDEX(centerInArr%atomid(j) , Hstring) == 1) THEN
                            mass(j) = Hmass
                        ELSEIF (INDEX(centerInArr%atomid(j) , Nstring) == 1) THEN
                            mass(j) = Nmass
                        ELSE
                            mass(j) = Omass
                        ENDIF


                        IF (atom1 == centerInArr%atomid(j)) THEN
                            masstot = masstot + mass(j)
                            xtot = xtot + (centerInArr%coor(j,1) * mass(j))
                            ytot = ytot + (centerInArr%coor(j,2) * mass(j))
                            ztot = ztot + (centerInArr%coor(j,centerIn%dimen) * mass(j))

                        ELSEIF (atom2 == centerInArr%atomid(j)) THEN
                            masstot = masstot + mass(j)
                            xtot = xtot + (centerInArr%coor(j,1) * mass(j))
                            ytot = ytot + (centerInArr%coor(j,2) * mass(j))
                            ztot = ztot + (centerInArr%coor(j,centerIn%dimen) * mass(j))
                        ENDIF

                        ! masstot = masstot + mass(j)

                        ! momenta(j,1) = mass(j) * centerInArr%coor(j,1)
                        ! momenta(j,2) = mass(j) * centerInArr%coor(j,2)
                        ! momenta(j,3) = mass(j) * centerInArr%coor(j,centerIn%dimen)
            
                        ! xtot = xtot + momenta(j,1)
                        ! ytot = ytot + momenta(j,2)
                        ! ztot = ztot + momenta(j,3) 
        
                    ENDDO

                    com(i,1) = xtot / masstot
                    com(i,2) = ytot / masstot
                    com(i,3) = ztot / masstot
        
                    WRITE(52,*) frame, i, com(i,1:centerIn%dimen)
        
                    masstot = 0.0
                    xtot = 0.0
                    ytot = 0.0
                    ztot = 0.0

                ENDDO

            DEALLOCATE(centerInArr%atomid, centerInArr%coor, mass, momenta, com)

            READ(22,*) nxtfrm
            READ(22,*) centerIn%genbyvmd

            IF (nxtfrm /= centerIn%atoms) THEN
                EXIT
            ENDIF

            frame = frame + 1

        ENDDO
    CLOSE(22)
    CLOSE(52)

    frame = 1
    
    RETURN

    END SUBROUTINE centermass

    !#############################################
    !##### O R D P A R     P R O C E D U R E #####
    !#############################################    
    
    SUBROUTINE ordpar(inpt, mol, sz)
    IMPLICIT NONE

    ! Dummy args
    CHARACTER(50), INTENT(IN) :: inpt
    INTEGER, INTENT(IN) :: mol, sz
    
    ! Constructors
    TYPE(istream) :: ordIn
    TYPE(istrArr) :: ordInArr
    TYPE(ostream) :: ordOut
    
    ! Subroutine vars
    INTEGER :: nxtfrm, frame = 1
    INTEGER :: i, j

    REAL, ALLOCATABLE :: ring(:,:,:), vector(:,:)
    REAL, ALLOCATABLE :: angle(:), ord(:)
    
    CHARACTER(2) :: string1, string2

    REAL :: xvec1, xvec2
    REAL :: yvec1, yvec2
    REAL :: zvec1, zvec2

    
    ordIn%filename = inpt
    ordOut%molec = mol
    ordOut%size = sz


    !========================================================================================================

    WRITE(6,*) '========================================================'
    WRITE(6,*) 'CALVIN : What name would like the output to be? '
    WRITE(6,*) '        '
    WRITE(6,*) '========================================================'
    WRITE(6,*) '        '
    READ(5,*) ordOut%outname

    WRITE(6,*) '==============================================================================================='
    WRITE(6,*) 'CALVIN: what would you like to find the cross product of? (separate your answers by a tab or space)'
    WRITE(6,*) 'Example: Atoms C1<=== center ===>C2 can be entered as C1  C2'
    READ(5,*) string1, string2
    WRITE(6,*) 'CALVIN is crunching the numbers....'


    OPEN(23, FILE = ordIn%filename)
    READ(23,*) ordIn%atoms
    READ(23,'(A)') ordIn%genbyvmd


    OPEN(53, file = ordOut%outname, status = 'new')
    WRITE(53,*) 'frame', '      ', 'mol', '     ', 'rads', '        ', 'degs', '        ', 'ordpar'


        DO
            ALLOCATE(ordInArr%atomid(ordIn%atoms), &
            & ordInArr%coor(ordIn%atoms, ordIn%dimen), &
            & ring(ordOut%molec, 2, ordIn%dimen), &
            & vector(ordOut%molec, ordIn%dimen), &
            & angle(ordOut%molec), &
            & ord(ordOut%molec))


                DO i = 1, ordOut%molec

                    DO j = 1, ordOut%size
                        READ(23,*) ordInArr%atomid(j), ordInArr%coor(j, 1:ordIn%dimen)
            
                        IF (ordInArr%atomid(j) == string1 ) THEN
                            xvec1 = ordInArr%coor(j, 1)
                            yvec1 = ordInArr%coor(j, 2)
                            zvec1 = ordInArr%coor(j, ordIn%dimen)
                        ELSEIF (ordInArr%atomid(j) == string2 ) THEN
                            xvec2 = ordInArr%coor(j, 1)
                            yvec2 = ordInArr%coor(j, 2)
                            zvec2 = ordInArr%coor(j, ordIn%dimen)
                        ELSE
                            CYCLE
                        ENDIF
        
                    ENDDO
                    
                    ring(i , 1 , 1) = xvec1
                    ring(i , 2 , 1) = xvec2
                    ring(i , 1 , 2) = yvec1
                    ring(i , 2 , 2) = yvec2
                    ring(i , 1 , ordIn%dimen) = zvec1
                    ring(i , 2 , ordIn%dimen) = zvec2

                    vector(i , 1) = ring(i , 2 , 1) - ring(i , 1 , 1)    
                    vector(i , 2) = ring(i , 2 , 2) - ring(i , 1 , 2)    
                    vector(i , 3) = ring(i , 2 , ordIn%dimen) - ring(i , 1 , ordIn%dimen)   
                

                    angle(i) = ACOS( (vector(i, ordIn%dimen)) / &
                    & (( ((vector(i, 1)) ** 2) + &
                    & ((vector(i, 2)) ** 2) + &
                    & ((vector(i, ordIn%dimen)) ** 2) ) ** 0.5))
                    ! In order to convert radians to degrees, multiply by 57.296
                                   
                    ord(i) = (((((COS((angle(i)))) ** 2) * 3) - 1) / 2)

                    
                    
                
                    WRITE(53,*) frame, i, angle(i), (angle(i) * 57.296), ord(i)
                
                    Xvec1 = 0.0
                    Xvec2 = 0.0
                    Yvec1 = 0.0
                    Yvec2 = 0.0
                    Zvec1 = 0.0
                    Zvec2 = 0.0

                ENDDO

            DEALLOCATE(ordInArr%atomid, ordInArr%coor, ring, vector, angle, ord)

            READ(23,*) nxtfrm
            READ(23,*) ordIn%genbyvmd

            IF (nxtfrm /= ordIn%atoms) THEN
                EXIT
            ENDIF

            frame = frame + 1

        ENDDO

    CLOSE(23)
    CLOSE(53)

    frame = 1
    
    RETURN
    
    END SUBROUTINE ordpar
    
    
END MODULE trajectories

!|-----------------------------------------|
!| /////////////////////////////////////// |
!| ///// E N E R G Y     M O D U L E ///// |
!| /////////////////////////////////////// |
!|-----------------------------------------|
    
MODULE energies
IMPLICIT NONE

TYPE energyIn
    CHARACTER(50) :: filename
    CHARACTER(100) :: header
    INTEGER :: head, tail
    INTEGER :: frame, time
    REAL :: elec, vdw, nonbond, total
END TYPE energyIn

TYPE energyOut
    REAL :: elecMean, vdwMean, nonbondMean, totalMean
END TYPE energyOut


CONTAINS
    
    
    SUBROUTINE avgenergy(nrgin)
    IMPLICIT NONE

    CHARACTER(50),INTENT(IN) :: nrgin
    
    TYPE(energyIn) :: avg
    TYPE(energyOut) :: means

    INTEGER :: stat, count = 0
    REAL :: elecsum = 0, vdwsum = 0, nonbondsum = 0, totalsum = 0


    avg%filename = nrgin

    !========================================================================================================

    OPEN(30, FILE = avg%filename)
    READ(30, '(A)') avg%header

    DO 

        READ(30,*,IOSTAT = stat) avg%frame, avg%time, avg%elec, avg%vdw, avg%nonbond, avg%total
        IF (stat /= 0) THEN
            EXIT
        ENDIF

        elecsum = elecsum + avg%elec
        vdwsum = vdwsum + avg%vdw
        nonbondsum = nonbondsum + avg%nonbond
        totalsum = totalsum + avg%total

        count = count + 1

    END DO

    means%elecMean = elecsum / count
    means%vdwMean = vdwsum / count
    means%nonbondMean = nonbondsum / count
    means%totalMean = totalsum / count


    WRITE(6,*) '========================================================'
    WRITE(6,*) 'CALVIN : The data has the following average energies... '
    WRITE(6,*) '        '
    WRITE(6,*) 'Total Frames: ', count
    WRITE(6,*) 'Mean Electrostatic Energy: ', means%elecMean
    WRITE(6,*) 'Mean VDW Energy: ', means%vdwMean
    WRITE(6,*) 'Mean Non-bonded Energy: ', means%nonbondMean
    WRITE(6,*) 'Mean Total Energy: ', means%totalMean
    WRITE(6,*) '========================================================'
    
    
    CLOSE(30)

    RETURN
    END SUBROUTINE avgenergy




END MODULE energies
    
    
    
    
    
    
    
    
    
    
    
    
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@ M A I N @@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
   
PROGRAM main
USE trajectories
USE energies
IMPLICIT NONE

!================================================================
INTEGER :: answer1, answer2
INTEGER :: atomCount, molCount, molSize, simFrame
CHARACTER(1) :: choice, choice2
CHARACTER(50) :: trajInput, nrginput
!================================================================

WRITE(6,*) '        '
WRITE(6,*) 'Hello, I am CALVIN. I am an Interface Designed to Process Large-Scale &
            &Simulation Trajectory data for analysis and presentation. The maker &
            &created me with one purpose: speed and efficiency'
WRITE(6,*) '        '



    DO

        WRITE(6,*) '============================================================'
        WRITE(6,*) 'CALVIN : Am I analyzing TRAJECTORIES, PROTEINS, OR ENERGIES?'
        WRITE(6,*) '0 - Atomic Trajectories'
        WRITE(6,*) '1 - Protein Sequence'
        WRITE(6,*) '2 - Energies'
        WRITE(6,*) '============================================================'
        READ(5,*) answer2

    
        !/////////////////////////////////////////////
        !////////// T R A J E C T O R I E S //////////
        !/////////////////////////////////////////////

        IF (answer2 == 0) THEN
            
            WRITE(6,*) '        '
            WRITE(6,*) 'CALVIN : Please enter the name of the trajectory file?'
            READ(5,*) trajInput
            WRITE(6,*) '        '
            WRITE(6,*) 'CALVIN: Before I can analyze the trajectories, I need to measure&
                        & the dimensions of this file so I can better understand it,&
                        & give me a few seconds...'

            
            CALL garray(trajInput, atomCount, molCount, molSize, simFrame)

            DO

                WRITE(6,*) '=========================================================='
                WRITE(6,*) '        '
                WRITE(6,*) 'CALVIN : What method would you like to perform?'
                WRITE(6,*) '0 - Print the atomic trajectories'
                WRITE(6,*) '1 - Center of Mass'
                WRITE(6,*) '2 - Order Parameter'
                WRITE(6,*) '3 - Density'
                WRITE(6,*) '        '
                READ(5,*) answer1

                IF (answer1 == 0) THEN
                    WRITE(6,*) 'CALVIN is crunching the numbers....'
                    CALL printarray(trajInput, simFrame, molCount, molSize)
                    WRITE(6,*) 'CALVIN : An ATOMIC TRAJECTORY output has&
                                & been successfully generated....'
                    WRITE(6,*) '        '
                    WRITE(6,*) '=========================================================='
                    WRITE(6,*) 'CALVIN : Would you like to perform another analysis? (y/n)'
                        READ(5,*) choice
                        IF (choice == 'y') THEN
                            CYCLE
                        ELSEIF (choice == 'n') THEN
                            EXIT
                        ENDIF
                ELSEIF (answer1 == 1) THEN
                    WRITE(6,*) '        '
                    WRITE(6,*) 'CALVIN is crunching the numbers....'
                    CALL centermass(trajInput, molCount, molSize)
                    WRITE(6,*) 'CALVIN : A CENTER OF MASS output file has&
                                & been successfully generated....'
                    WRITE(6,*) '        '
                    WRITE(6,*) '=========================================================='
                    WRITE(6,*) 'CALVIN : Would you like to perform another analysis? (y/n)'
                        READ(5,*) choice
                        IF (choice == 'y') THEN
                            CYCLE
                        ELSEIF (choice == 'n') THEN
                            EXIT
                        ENDIF
                ELSEIF (answer1 == 2) THEN
                    WRITE(6,*) '        '
                    CALL ordpar(trajInput, molCount, molSize)
                    WRITE(6,*) 'CALVIN : An ORDER PARAMETER output has &
                                & been successfully generated....'
                    WRITE(6,*) '        '
                    WRITE(6,*) '=========================================================='
                    WRITE(6,*) 'CALVIN : Would you like to perform another analysis? (y/n)'
                        READ(5,*) choice
                        IF (choice == 'y') THEN
                            CYCLE
                        ELSEIF (choice == 'n') THEN
                            EXIT
                        ENDIF
                ELSEIF (answer1 == 3) THEN
                    WRITE(6,*) 'CALVIN : This analysis is currently being built, pardon our dust...'
                    WRITE(6,*) '        '
                    WRITE(6,*) '=========================================================='
                    WRITE(6,*) 'CALVIN : Would you like to perform another analysis? (y/n)'
                    READ(5,*) choice
                        IF (choice == 'y') THEN
                            CYCLE
                        ELSEIF (choice == 'n') THEN
                            EXIT
                        ENDIF
                ELSE
                    WRITE(6,*) 'CALVIN : That is not a valid answer, please try again'
                ENDIF


            ENDDO


            WRITE(6,*) '        '
            WRITE(6,*) '=========================================================='
            WRITE(6,*) 'CALVIN : Would you like to look at another file? (y/n)'
            READ(5,*) choice2
                IF (choice2 == 'y') THEN
                    CYCLE
                ELSEIF (choice2 == 'n') THEN
                    EXIT
                ENDIF



        !///////////////////////////////////////////////////////
        !////////// P R O T E I N     S E Q U E N C E //////////
        !///////////////////////////////////////////////////////

        ELSEIF (answer2 == 1) THEN
            WRITE(6,*) 'You have chosen to analyze protein structure...'




        !/////////////////////////////////////
        !////////// E N E R G I E S //////////
        !/////////////////////////////////////

        ELSEIF (answer2 == 2) THEN

            WRITE(6,*) '        '
            WRITE(6,*) 'CALVIN : Please enter the name of the Energy File?'
            READ(5,*) nrginput
            WRITE(6,*) '        '


            ! WRITE(6,*) 'CALVIN is crunching the numbers....'
            CALL avgenergy(nrginput)

            ! WRITE(6,*) 'CALVIN : An ENERGY AVERAGE output file has&
            !             & been successfully generated....'
            WRITE(6,*) '        '
            WRITE(6,*) '=========================================================='
            WRITE(6,*) 'CALVIN : Would you like to perform another analysis? (y/n)'
            READ(5,*) choice
                IF (choice == 'y') THEN
                    CYCLE
                ELSEIF (choice == 'n') THEN
                    EXIT
                ENDIF

        ENDIF



    ENDDO


END PROGRAM main
