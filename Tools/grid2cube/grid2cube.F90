!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2013,2015, Ignacio Fdez. Galvan                        *
!***********************************************************************
!===============================================================================
! Fortran 2003 program for converting Molcas "grid" files into
!                                     Gaussian "cube" format
!
! Last modified: 2015 June 18
!            by: Ignacio Fdez. Galván
!===============================================================================

#define STRLEN 256

PROGRAM Grid2Cube
IMPLICIT NONE
CHARACTER(LEN=STRLEN), DIMENSION(:), ALLOCATABLE :: Args
INTEGER :: FileIn, FileOut, Natom, Block_Size, N_of_Grids, Grid, NP
LOGICAL :: Binary
INTEGER, DIMENSION(3) :: Net, NPt
REAL(KIND=8), DIMENSION(3) :: Origin, Axis_1, Axis_2, Axis_3
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: Coor
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PtData
CHARACTER(6), DIMENSION(:), ALLOCATABLE :: Label
CHARACTER(LEN=STRLEN), DIMENSION(:), ALLOCATABLE :: GridName
CHARACTER(LEN=2), DIMENSION(0:112) :: Symbol

! Define the atomic symbols for getting the atomic numbers
Symbol(:) = (/                                          "X ", &
  "H ", "HE", "LI", "BE", "B ", "C ", "N ", "O ", "F ", "NE", &
  "NA", "MG", "AL", "SI", "P ", "S ", "CL", "AR", "K ", "CA", &
  "SC", "TI", "V ", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", &
  "GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y ", "ZR", &
  "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", &
  "SB", "TE", "I ", "XE", "CS", "BA", "LA", "CE", "PR", "ND", &
  "PM", "SM", "EU", "GD", "TB", "DY", "HO", "ER", "TM", "YB", &
  "LU", "HF", "TA", "W ", "RE", "OS", "IR", "PT", "AU", "HG", &
  "TL", "PB", "BI", "PO", "AT", "RN", "FR", "RA", "AC", "TH", &
  "PA", "U ", "NP", "PU", "AM", "CM", "BK", "CF", "ES", "FM", &
  "MD", "NO", "LR", "RF", "DB", "SG", "BH", "HS", "MT", "DS", &
  "RG", "CN"                                                  /)

! Read the filenames, grid header, and let the user select
! the desired grid
CALL Process_Arguments()
CALL Read_Header()
CALL Select_Grid(Grid)

! Total number of points
NP = PRODUCT(NPt(:))
ALLOCATE(PtData(NP))

! Read grid data
CALL Read_Data(Grid)
CLOSE(FileIn)

! Write output file
CALL Write_Data(TRIM(GridName(Grid)))
CLOSE(FileOut)

DEALLOCATE(Args,Coor,Label,GridName,PtData)

!===============================================================================

CONTAINS

!-------------------------------------------------------------------------------
! Read the command-line arguments and perform some checks

SUBROUTINE Process_Arguments()
INTEGER :: NumArgs, i, Error

NumArgs = COMMAND_ARGUMENT_COUNT()
IF (NumArgs < 2) THEN
  CALL Write_Usage()
  STOP
END IF

ALLOCATE(Args(NumArgs))

DO i=1,NumArgs
  CALL GET_COMMAND_ARGUMENT(i,Args(i))
END DO

#ifdef _DEBUG_
WRITE(6,1) 'Input file: ',TRIM(Args(1))
WRITE(6,1) 'Output file:',TRIM(Args(2))
#endif

! Check that the input file can be opened, and whether it is binary or ASCII
CALL Binary_or_ASCII(20,TRIM(Args(1)),Error)
IF (Error /= 0 ) THEN
  WRITE(6,1) '** Error opening file',TRIM(Args(1))
  STOP
 ELSE
  FileIn=20
END IF

! Check that the output file can be opened
OPEN(30,FILE=TRIM(Args(2)),STATUS='REPLACE',ACTION='WRITE',IOSTAT=Error)
IF (Error /= 0 ) THEN
  WRITE(6,1) '** Error opening file',TRIM(Args(2))
  STOP
 ELSE
  FileOut=30
END IF

1 FORMAT (A,1X,A)

END SUBROUTINE Process_Arguments

!-------------------------------------------------------------------------------
! Write usage information

SUBROUTINE Write_Usage()

WRITE(6,1)
WRITE(6,1) 'Convert a Molcas grid file into a Gaussian cube file.'
WRITE(6,1) 'The input file can be ASCII or binary (non-packed),'
WRITE(6,1) 'the generated output file will be ASCII.'
WRITE(6,1)
WRITE(6,1) 'USAGE:'
WRITE(6,1) '  grid2cube input_file output_file'
WRITE(6,1)
WRITE(6,1) 'If there are several grids in the input file,'
WRITE(6,1) 'the user will be asked which one to convert.'
WRITE(6,1)

1 FORMAT (A)

END SUBROUTINE Write_Usage

!-------------------------------------------------------------------------------
! Detect whether the input file is binary or ASCII

SUBROUTINE Binary_or_ASCII(U,FileName,Error)
INTEGER, INTENT(IN) :: U
CHARACTER(LEN=*), INTENT(IN) :: FileName
INTEGER, INTENT(OUT) :: Error
INTEGER :: Test
CHARACTER :: TestChar
REAL(KIND=8) :: TestFloat

! Try to read as ASCII, the first line should be 0
OPEN(U,FILE=FileName,STATUS='OLD',ACTION='READ',IOSTAT=Error,FORM='FORMATTED')
IF (Error /= 0) RETURN
READ(U,*,IOSTAT=Error) Test
IF (Error == 0 .AND. Test == 0) THEN
#ifdef _DEBUG_
  WRITE(6,1) 'Input format: ASCII'
#endif
  Binary = .FALSE.
  RETURN
END IF
CLOSE(U)

! Try to read as binary, the first two records should be "a" and 1999.0
OPEN(U,FILE=FileName,STATUS='OLD',ACTION='READ',IOSTAT=Error,FORM='UNFORMATTED')
IF (Error /= 0) RETURN
READ(U,IOSTAT=Error) TestChar
IF (Error == 0) READ(U,IOSTAT=Error) TestFloat
IF (Error == 0 .AND. TestChar == 'a' .AND. TestFloat == 1999.0) THEN
#ifdef _DEBUG_
  WRITE(6,1) 'Input format: binary'
#endif
  Binary = .TRUE.
  RETURN
END IF
CLOSE(U)

#ifdef _DEBUG_
WRITE(6,1) 'Unknown input format'
#endif
Binary = .FALSE.
Error=1

#ifdef _DEBUG_
1 FORMAT (A)
#endif

END SUBROUTINE Binary_or_ASCII

!-------------------------------------------------------------------------------
! Read a string from the input file

SUBROUTINE Read_String(U,String)
INTEGER, INTENT(IN) :: U
CHARACTER(LEN=*), INTENT(OUT) :: String
INTEGER :: i, Error

! If the file is binary, try to read the string with decreasing lengths until
! it gives no error
IF (Binary) THEN
  DO i=LEN(String),0,-1
    READ(U,IOSTAT=Error) String(:i)
    IF (Error == 0) THEN
      String(i+1:) = ''
      EXIT
    END IF
    BACKSPACE(U)
  END DO
! For ASCII files it is much easier
 ELSE
  READ(U,1) String
END IF

1 FORMAT (A)

END SUBROUTINE Read_String

!-------------------------------------------------------------------------------
! Convert string to uppercase letters

SUBROUTINE To_Upper(String)
CHARACTER(LEN=*), INTENT(INOUT) :: String
INTEGER :: i, Letter

DO i=1,LEN(TRIM(String))
  Letter=ICHAR(String(i:i))
  IF ((Letter > 96) .AND. (Letter < 193)) String(i:i)=CHAR(Letter-32)
END DO

END SUBROUTINE

!-------------------------------------------------------------------------------
! Read the header of the input file

SUBROUTINE Read_Header()
CHARACTER(LEN=STRLEN) :: Line, Word
INTEGER :: i, ngrid

ngrid = 0

! Read line by line until "Title=" is found and store some useful values
DO
  CALL Read_String(FileIn,Line)
  IF (LEN(TRIM(Line)) == 0) CYCLE
  READ(Line,*) Word
  SELECT CASE (TRIM(Word))

   CASE("Natom=")
    READ(Line,*) Word,Natom
    ALLOCATE(Coor(Natom,3),Label(Natom))
    DO i=1,Natom
      CALL Read_String(FileIn,Line)
      READ(Line,*) Label(i),Coor(i,:)
    END DO

   CASE("N_of_Grids=")
    READ(Line,*) Word,N_of_Grids
    ALLOCATE(GridName(N_of_Grids))

   CASE("Block_Size=")
    READ(Line,*) Word,Block_Size

   CASE("Net=")
    READ(Line,*) Word,Net(:)
    NPt = Net(:)+1

   CASE("Origin=")
    READ(Line,*) Word,Origin(:)

   CASE("Axis_1=")
    READ(Line,*) Word,Axis_1(:)

   CASE("Axis_2=")
    READ(Line,*) Word,Axis_2(:)

   CASE("Axis_3=")
    READ(Line,*) Word,Axis_3(:)

   CASE("GridName=")
    ngrid = ngrid+1
    GridName(ngrid) = ADJUSTL(Line(11:))

   CASE("Title=")
    EXIT

  END SELECT
END DO
BACKSPACE(FileIn)

#ifdef _DEBUG_
WRITE(6,2) 'Natom=     ',Natom
DO i=1,Natom
  WRITE(6,4) Label(i),Coor(i,:)
END DO
WRITE(6,2) 'Block_Size=',Block_Size
WRITE(6,2) 'Net=       ',Net
WRITE(6,3) 'Origin=    ',Origin(:)
WRITE(6,3) 'Axis_1=    ',Axis_1(:)
WRITE(6,3) 'Axis_2=    ',Axis_2(:)
WRITE(6,3) 'Axis_3=    ',Axis_3(:)
WRITE(6,2) 'N_of_Grids=',N_of_Grids
DO i=1,N_of_Grids
  WRITE(6,5) i,TRIM(GridName(i))
END DO
2 FORMAT (A,3(1X,I6))
3 FORMAT (A,3(1X,F11.3))
4 FORMAT (2X,A,3(1X,F10.5))
5 FORMAT (2X,I3,' -- ',A)
#endif

IF (N_of_Grids < 1 .OR. ngrid < 1) THEN
  WRITE(6,1) '** No grids found in input file'
  STOP
END IF

1 FORMAT (A)

END SUBROUTINE Read_Header

!-------------------------------------------------------------------------------
! Print the grid names found in the file, and ask the user to choose one

SUBROUTINE Select_Grid(ngrid)
INTEGER, INTENT(OUT) :: ngrid
INTEGER :: i, Error

WRITE(6,1) 'Grids in input file:'
DO i=1,N_of_Grids
  WRITE(6,2) i,TRIM(GridName(i))
END DO

WRITE(6,1,ADVANCE='NO') 'Which one to convert? '
READ(5,*,IOSTAT=Error) ngrid

IF (Error /= 0) THEN
  WRITE(6,1) "** Sorry, I don't understand"
  STOP
END IF
IF (ngrid < 1 .OR. ngrid > N_of_Grids) THEN
  WRITE(6,1) '** Sorry, no such grid'
  STOP
END IF

1 FORMAT (A)
2 FORMAT (I3,': ',A)

END SUBROUTINE Select_Grid

!-------------------------------------------------------------------------------
! Read the data corresponding to the desired grid

SUBROUTINE Read_Data(ngrid)
INTEGER, INTENT(IN) :: ngrid
CHARACTER(LEN=STRLEN) :: Line
INTEGER :: i, Current_Block, Current_Grid, Data_Read, Data_in_Block

Current_Block = 0
Current_Grid = 0
Data_Read = 0

! The data is stored in blocks of at most Block_Size values for each grid,
! separated by "Title=" lines

DO WHILE (Data_Read < NP)

  CALL Read_String(FileIn,Line)
  ! Advance grid number when "Title" is found,
  ! and block number when the first grid is found again
  IF (Line(1:6) == 'Title=') THEN
    IF (MOD(Current_Grid,N_of_Grids) == 0) Current_Block = Current_Block+1
    Current_Grid = MOD(Current_Grid,N_of_Grids)+1
#ifdef _DEBUG_
    WRITE(6,*) TRIM(Line)
    WRITE(6,*) 'Current_Block: ',Current_Block
    WRITE(6,*) 'Current_Grid:  ',Current_Grid
#endif
  END IF
  IF (Current_Grid /= ngrid) CYCLE

  ! Number of values in this block, is Block_Size or the remaining number
  Data_in_Block = MIN(NP-(Current_Block-1)*Block_Size,Block_Size)
#ifdef _DEBUG_
  WRITE(6,*) 'Data_in_Block: ',Data_in_Block
#endif
  IF (Binary) THEN
    READ(FileIn) PtData(Data_Read+1:Data_Read+Data_in_Block)
   ELSE
    DO i=Data_Read+1,Data_Read+Data_in_Block
      READ(FileIn,*) PtData(i)
    END DO
  END IF
  Data_Read = Data_Read+Data_in_Block

END DO

END SUBROUTINE Read_Data

!-------------------------------------------------------------------------------
! Write the grid data in (formatted) cube format

SUBROUTINE Write_Data(Title)
CHARACTER(LEN=*), INTENT(IN) :: Title
REAL(KIND=8), DIMENSION(3) :: dx, dy, dz
INTEGER, DIMENSION(Natom) :: Atom_Number
INTEGER :: i, j
CHARACTER(LEN=2) :: Lab

! Compute delta values for each axis and get atom numbers
dx = Axis_1(:)/Net(1)
dy = Axis_2(:)/Net(2)
dz = Axis_3(:)/Net(3)
DO i=1,Natom
  Lab = Label(i)(1:2)
  CALL To_Upper(Lab)
  IF ((ICHAR(Lab(2:2)) < 97) .OR. (ICHAR(Lab(2:2)) > 122)) Lab(2:2)=' '
  DO j=LBOUND(Symbol,1),UBOUND(Symbol,1)
    IF (Lab == Symbol(j)) THEN
      Atom_Number(i) = j
      EXIT
    END IF
  END DO
END DO

! Write the header
WRITE(FileOut,1) 'File converted from MOLCAS grid format'
WRITE(FileOut,1) 'Title = '//Title
WRITE(FileOut,2) Natom, Origin(:)
WRITE(FileOut,2) NPt(1), dx(:)
WRITE(FileOut,2) NPt(2), dy(:)
WRITE(FileOut,2) NPt(3), dz(:)
DO i=1,Natom
  WRITE(FileOut,2) Atom_Number(i), 0.0, Coor(i,:)
END DO
! Write the grid data:
!  6 values per line,
!  forced new line every NPt(3) values
DO i=1,NP,NPt(3)
  WRITE(FileOut,3) PtData(i:i+NPt(3)-1)
END DO

1 FORMAT (A)
2 FORMAT (I5,4F12.6)
3 FORMAT (6ES13.5)

END SUBROUTINE Write_Data

END PROGRAM Grid2Cube
