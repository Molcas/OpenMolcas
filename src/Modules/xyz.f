************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016, Morgane Vacher                                   *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************

* Module to handle the COORD keyword in gateway/seward

#include "compiler_features.h"
#define MAXLEN 180
#define iX 1
#define iY 2
#define iZ 4
#define iXY iX+iY
#define iXZ iX+iZ
#define iYZ iY+iZ
#define iXYZ iX+iY+iZ

      Module XYZ
      Implicit None
      Private
      Type XYZAtom
        Character (Len=MAXLEN) :: Lab
        Real*8, Dimension(3) :: Coord
        Integer :: FileNum
      End Type XYZAtom
      Type(XYZAtom), Dimension(:), Allocatable :: Geom
      Character (Len=256) :: BasisAll, Symmetry
      Character (Len=256), Dimension(:,:), Allocatable :: BasisSets
      Integer, Dimension(7) :: Oper
      Integer :: FileNum=0

      Public :: Read_XYZ, Write_SewInp, Clear_XYZ, Out_Raw,
     &          Parse_Basis, Parse_Group, Symmetry

      Contains

! Read an XYZ file with molcas extensions
! Rot and Trans are transformations to apply to all coordinates in this file
      Subroutine Read_XYZ(Lu,Rot,Trans)
      Integer, Intent(In) :: Lu
      Real*8, Dimension(:,:,:), Allocatable, Intent(In) :: Rot
      Real*8, Dimension(:,:), Allocatable, Intent(In) :: Trans
      Character (Len=MAXLEN) :: Line, FName, CurrDir, Dum
      Integer :: Error, NumAt, i, Lxyz, Idx
      Logical :: Found
      Integer, External :: IsFreeUnit
      Type(XYZAtom), Dimension(:), Allocatable :: ThisGeom, TmpGeom
      Real*8, Dimension(3,5) :: Mat
      Real*8 :: Factor
#include "real.fh"
#include "constants2.fh"
#ifdef _HDF5_
#  include "mh5.fh"
#  include "Molcas.fh"
      Logical :: isH5
      Integer :: Coord_id, Attr_id, nSym, j, c
      Character (Len=LenIn), Dimension(:), Allocatable :: Labels
      Character (Len=LenIn4), Dimension(:), Allocatable :: Labels4
      Real*8, Dimension(:,:), Allocatable :: Coords
      isH5 = .False.
#endif

      Factor = One
      Read(Lu,'(A)') Line
      Line = AdjustL(Line)

      ! Try to read a number, if it fails, try to open a file
      ! Note that the slash means end-of-line in list-directed input
      NumAt = -1
      Read(Line,*,IOStat=Error) NumAt
      ! And make sure there's nothing but numbers in the first word
      ! (the above is not 100% reliable in some compilers)
      Do i=1,Len(Line)
        If (Line(i:i) .eq. ' ') Exit
        Idx = IChar(Line(i:i))
        If ((Idx .lt. IChar('0')) .or. (Idx .gt. IChar('9'))) Then
          Error = -1
          Exit
        End If
      End Do
      If (NumAt .lt. 1) Error = -1
      If (Index(Line, '/') .gt. 0) Error = -1
      Lxyz = Lu
      If (Error .ne. 0) Then
        Read(Line,'(A)') FName
        Found = .False.
        If (Index(Line, '/') .eq. 0) Then
          Call GetEnvF('CurrDir', CurrDir)
          CurrDir = Trim(CurrDir)//'/'//FName
          Call F_Inquire(CurrDir, Found)
          If (Found) FName = CurrDir
        End If
        If (.Not. Found) Call F_Inquire(FName, Found)
        If (Found) Then
#ifdef _HDF5_
          If (mh5_is_hdf5(Trim(FName))) Then
            isH5 = .True.
          Else
#endif
            Lxyz = IsFreeUnit(Lxyz)
            Call Molcas_Open(Lxyz, FName)
            Read(Lxyz,'(A)') Line
            Read(Line,*,IOStat=Error) NumAt
            If (Error .ne. 0) Then
              Write(6,*) 'Error reading file ',Trim(FName)
              Call Quit_OnUserError()
            End If
#ifdef _HDF5_
          End If
#endif
        Else
          Write(6,*) 'File ',Trim(FName),' not found!'
          Call Quit_OnUserError()
        End If
      End If
      ! Each atom has a field indicating the sequential number of the
      ! file it belongs to
      FileNum = FileNum+1

#ifdef _HDF5_
************************************************************************
* For HDF5 formatted file
************************************************************************
      If (isH5) Then
        Write(6,*) 'Reading xyz coordinates from h5 file '//Trim(FName)
        Coord_id = mh5_open_file_r(Trim(FName))
        ! check if symmetry was used
        Call mh5_fetch_attr(Coord_id,'NSYM',nSym)
        ! read numbers of atoms
        If (nSym .gt. 1) Then
          Attr_id = mh5_open_attr(Coord_id,'NATOMS_ALL')
        Else
          Attr_id = mh5_open_attr(Coord_id,'NATOMS_UNIQUE')
        End If
        Call mh5_get_attr_scalar_int(Attr_id,NumAt)
        Allocate(ThisGeom(NumAt))
        Allocate(Labels(NumAt),Coords(3,NumAt))
        ! read atom labels
        If (nSym .gt. 1) then
          Allocate(Labels4(NumAt))
          Call mh5_fetch_dset_array_str(Coord_id,'DESYM_CENTER_LABELS',
     &                                   Labels4)
          Do i=1,NumAt
            Labels(i) = Labels4(i)(1:LenIn)
          End Do
          Deallocate(Labels4)
        Else
          Call mh5_fetch_dset_array_str(Coord_id,'CENTER_LABELS',
     &                                   Labels)
        End If
        ! remove numbers from labels
        ! (otherwise there will be problems when symmetric structures
        ! are used without symmetry)
        Do i=1,NumAt
          Do j=1,LenIn
            c = IChar(Labels(i)(j:j))
            If ((c .ge. IChar('0')) .and. (c .le. IChar('9'))) Then
              Labels(i)(j:j) = ' '
            End If
          End Do
        End Do
        ! read atom coordinates
        If (nSym .gt. 1) then
          Call mh5_fetch_dset_array_real(Coord_id,
     &         'DESYM_CENTER_COORDINATES',Coords)
        Else
          Call mh5_fetch_dset_array_real(Coord_id,
     &         'CENTER_COORDINATES',Coords)
        End If
        Call mh5_close_file(Coord_id)
        ! store data
        Do i=1,NumAt
          ThisGeom(i)%Lab = Labels(i)
          ThisGeom(i)%Coord(:) = Coords(:,i)
          ThisGeom(i)%FileNum = FileNum
        End Do
        Deallocate(Labels,Coords)
      Else
#endif
************************************************************************
* For .xyz format file
************************************************************************
        Allocate(ThisGeom(NumAt))
        Read(Lxyz,'(A)',IOStat=Error) Line
        If (Error .ne. 0) Then
          Write(6,*) 'Error reading geometry'
          Call Quit_OnUserError()
        End If
        Call UpCase(Line)
        ! Units of the coordinates (default angstrom)
        If (Max(Index(Line,'BOHR'),Index(Line,'A.U.')) .le. 0) Then
          Factor = One/Angstrom
        End If
        Do i=1,NumAt
          Read(Lxyz,*,IOStat=Error) ThisGeom(i)%Lab,ThisGeom(i)%Coord(:)
          If (Error .ne. 0) Then
            Write(6,*) 'Error reading geometry'
            Call Quit_OnUserError()
          End If
          ThisGeom(i)%Coord(:) = ThisGeom(i)%Coord(:)*Factor
          ThisGeom(i)%FileNum = FileNum
        End Do
        If (Lxyz .ne. Lu) Close(Lxyz)
#ifdef _HDF5_
      End If
#endif

      ! Obtain/read transformation matrix and transform the geometry in this file
      Mat = Reshape([One, One,  One,
     &               One, Zero, Zero,
     &               Zero, One, Zero,
     &               Zero, Zero, One,
     &               Zero, Zero, Zero],
     &              Shape(Mat))
      Idx = Index(' '//Line, ' SCALE ')
      If (Idx .gt. 0) Then
        Read(Line(Idx:),*,IOStat=Error) Dum, Mat(1,1)
        Mat(2,1) = Mat(1,1)
        Mat(3,1) = Mat(1,1)
      End If
      Idx = Index(' '//Line, ' SCALEX ')
      If (Idx .gt. 0) Read(Line(Idx:),*,IOStat=Error) Dum, Mat(1,1)
      Idx = Index(' '//Line, ' SCALEY ')
      If (Idx .gt. 0) Read(Line(Idx:),*,IOStat=Error) Dum, Mat(2,1)
      Idx = Index(' '//Line, ' SCALEZ ')
      If (Idx .gt. 0) Read(Line(Idx:),*,IOStat=Error) Dum, Mat(3,1)
      Idx = Index(' '//Line, ' ROT ')
      If (Idx .gt. 0) Read(Line(Idx:),*,IOStat=Error) Dum, Mat(:,2:4)
      Idx = Index(' '//Line, ' TRANS ')
      If (Idx .gt. 0) Read(Line(Idx:),*,IOStat=Error) Dum, Mat(:,5)
      ! If Rot and Trans are given in the input, they override
      ! the inline transformations
      If (Allocated(Rot)) Then
        Mat(:,2:4) = Rot(:,:,FileNum)
      End If
      Mat(:,2:4) = Transpose(Mat(:,2:4))
      If (Allocated(Trans)) Then
        Mat(:,5) = Trans(:,FileNum)
      End If
      Mat(:,5) = Mat(:,5)*Factor
      Call TransformGeom(ThisGeom,Mat)

      ! Append the just read geometry to the general one
      If (Allocated(Geom)) Then
        Call Move_Alloc(Geom, TmpGeom)
        NumAt = Size(TmpGeom)+Size(ThisGeom)
        Allocate(Geom(NumAt))
        Geom(1:Size(TmpGeom)) = TmpGeom(:)
        Geom(Size(TmpGeom)+1:) = ThisGeom(:)
        Deallocate(TmpGeom)
        Deallocate(ThisGeom)
      Else
        Call Move_Alloc(ThisGeom, Geom)
      End If

      End Subroutine Read_XYZ

! Clear allocations
      Subroutine Clear_XYZ

      If (Allocated(Geom)) Deallocate(Geom)
      If (Allocated(BasisSets)) Deallocate(BasisSets)
      FileNum = 0

      End Subroutine Clear_XYZ

! Write an input file for seward
! Atoms belonging to GhostFiles have no charge
      Subroutine Write_SewInp(FName,GhostFiles)
      Character (Len=*), Intent(In) :: FName
      Integer, Dimension(:), Intent(In) :: GhostFiles
      Integer :: Lu, Num, i, j
      Integer, External :: IsFreeUnit
      Character (Len=MAXLEN) :: New, Old, Sym, Lab, Bas
      Logical :: Ghost

      Lu = 10
      Lu = IsFreeUnit(Lu)
      Call Molcas_Open(Lu, FName)
      ! Write symmetry
      If (Symmetry .ne. '') Write(Lu,40) 'Symmetry', Trim(Symmetry)
      ! Write the atoms, skipping symmetry-superfluous atoms
      Old = ''
      j = 0
      Do i=1,Size(Geom)
        If (Geom(i)%FileNum .eq. 0) Cycle
        j = j+1
        Call SplitLabel(Geom(i)%Lab, Sym, Num, Lab, Bas)
        ! Build an identifier for the basis
        New = Trim(Sym)//Trim(Lab)//Bas
        Ghost = Any(GhostFiles .eq. Geom(i)%FileNum)
        If (Ghost) New = Trim(New)//' *'
        ! If there is a change of basis, write the basis
        If (New .ne. Old) Then
          If (i .gt. 1) Write(Lu,10) 'End of Basis'
          Write(Lu,10) 'Basis Set'
          Write(Lu,10) FindBasis(Sym, Lab, Bas)
          If (Ghost) Write(Lu,30) 'Charge', 0.0
          Old = New
        End If
        ! Add a sequential number if the label had none
        If (Num .eq. 0) Num = j
        Write(Lu,20) Trim(Sym), Num, Geom(i)%Coord(:)
      End Do
      Write(Lu,10) 'End of Basis'
      Write(Lu,10) 'End of Coord'
      Close(Lu)
10    Format(A)
20    Format(A,I0,3(1X,ES30.20E3))
30    Format(A,/,F4.1)
40    Format(A,/,A)

      End Subroutine Write_SewInp

! Store symmetry-unique atom coordinates in an array
! Output is the number of saved coordinates
      Function Out_Raw(Array)
      Real*8, Dimension(*), Intent(InOut) :: Array
      Integer :: Out_Raw
      Integer :: i, j

      j = 0
      Do i=1,Size(Geom)
        If (Geom(i)%FileNum .eq. 0) Cycle
        Array(j+1:j+3) = Geom(i)%Coord
        j = j+3
      End Do
      Out_Raw = j

      End Function Out_Raw

! Parse the BASIS keyword to store the basis set corresponding to
! each atom + label
      Subroutine Parse_Basis(Basis)
      Character (Len=*), Intent(In) :: Basis
      Integer :: Num, Idx, IdxDot, Next, i

      ! Count number of commas
      Num = Count(Transfer(Basis, 'x', Len_Trim(Basis)) .eq. ",") + 1
      BasisAll = ''
      If (Allocated(BasisSets)) Deallocate(BasisSets)
      Allocate(BasisSets(2,Num))
      ! For each comma-separated word, split it at the first dot
      ! If the first part is empty, use it as a general basis set
      Idx = 0
      Do i=1,Num
        Next = Index(Basis(Idx+1:), ',')
        If (Next .eq. 0) Next = Len_Trim(Basis)+1
        BasisSets(2,i) = Basis(Idx+1:Next-1)
        IdxDot = Index(BasisSets(2,i), '.')
        If (IdxDot .eq. 0) Then
          BasisSets(1,i) = ''
        Else
          BasisSets(1,i) = AdjustL(BasisSets(2,i)(1:IdxDot-1))
          BasisSets(2,i)(1:IdxDot) = ''
        End If
        BasisSets(2,i) = AdjustL(BasisSets(2,i))
        If (BasisSets(1,i) .eq. '') BasisAll = BasisSets(2,i)
        Idx = Next
      End Do
      If (BasisAll .eq. '') BasisAll = 'ANO-S-MB'

      End Subroutine Parse_Basis

! Parse the GROUP keyword, detecting the symmetry,
! generating all symmetry operations and adapting
! the symmetry
      Subroutine Parse_Group(Group,Thr)
      Character (Len=*), Intent(In) :: Group
      Real*8, Intent(In) :: Thr
      Character (Len=3), Dimension(3) :: Gen
      Integer :: Error, i, j, k

      Symmetry = Group
      Call UpCase(Symmetry)
      If (Symmetry .eq. 'FULL') Then
        Call DetectSym(Thr)
      Else If ((Symmetry(1:5) .eq. 'NOSYM') .or.
     &         (Symmetry .eq. 'E') .or.
     &         (Symmetry .eq. 'C1')) Then
        Symmetry = ''
      End If
      Gen = ''
      Read(Symmetry,*,IOStat=Error) Gen
      ! Encode generators
      Oper = 0
      Do i=1,3
        If (Index(Gen(i),'X') .gt. 0) Oper(i)=Oper(i)+iX
        If (Index(Gen(i),'Y') .gt. 0) Oper(i)=Oper(i)+iY
        If (Index(Gen(i),'Z') .gt. 0) Oper(i)=Oper(i)+iZ
      End Do
      ! Create all combinations
      Oper(4) = iEOR(Oper(1),Oper(2))
      Oper(5) = iEOR(Oper(1),Oper(3))
      Oper(6) = iEOR(Oper(2),Oper(3))
      Oper(7) = iEOR(Oper(3),Oper(4))
      ! Remove duplicates
      Do i=1,7
        Do j=i+1,7
          If (Oper(j) .eq. Oper(i)) Oper(j) = 0
        End Do
      End Do
      ! Sort descending
      Do i=1,7
        Do j=i+1,7
          If (Oper(j) .ge. Oper(i)) Then
            k = Oper(j)
            Oper(j) = Oper(i)
            Oper(i) = k
          End If
        End Do
      End Do
      ! Check and adapt the symmetry
      Call AdaptSym(Thr)

      End Subroutine Parse_Group

! Private procedures follow

! Subroutine to detect the symmetry elements of the system
! The system is *not* translated or rotated
! Only elements of D2h group are tested (i.e. sign inversions of one, two or three axes)
! The result is stored in Symmetry (compatible with molcas-extra)
      Subroutine DetectSym(Thr)
      Real*8, Intent(In) :: Thr
      Logical, Dimension(7) :: Op
      Op = .False.
      Op(iX) = CheckOp(iX, Thr)
      Op(iY) = CheckOp(iY, Thr)
      Op(iZ) = CheckOp(iZ, Thr)
      Symmetry = ''
      Select Case (Count(Op))
        ! Two or three reflections: all is known
        Case (3, 2)
          Op(iXY) = Op(iX) .And. Op(iY)
          Op(iXZ) = Op(iX) .And. Op(iZ)
          Op(iYZ) = Op(iY) .And. Op(iZ)
          Op(iXYZ) = Op(iX) .And. Op(iY) .And. Op(iZ)
          If (Op(iXYZ)) Then
            Symmetry = 'x y z'
          Else
            If (Op(iXY)) Symmetry = 'xy y'
            If (Op(iXZ)) Symmetry = 'xz z'
            If (Op(iYZ)) Symmetry = 'yz z'
          End If
        ! One reflection: possibly inversion and complementary rotation
        Case (1)
          Op(iXYZ) = CheckOp(iXYZ, Thr)
          If (Op(iXYZ)) Then
            Op(iXY) = Op(iZ)
            Op(iXZ) = Op(iY)
            Op(iYZ) = Op(iX)
            If (Op(iXY)) Symmetry = 'xy xyz'
            If (Op(iXZ)) Symmetry = 'xz xyz'
            If (Op(iYZ)) Symmetry = 'yz xyz'
          Else
            If (Op(iX)) Symmetry = 'x'
            If (Op(iY)) Symmetry = 'y'
            If (Op(iZ)) Symmetry = 'z'
          End If
        ! No reflection: check rotations (only two initially)
        Case (0)
          Op(iXY) = CheckOp(iXY, Thr)
          Op(iXZ) = CheckOp(iXZ, Thr)
          If (Op(iXY)) Symmetry = 'xy'
          If (Op(iXZ)) Symmetry = Trim(Symmetry)//' xz'
          Select Case (Count(Op))
            ! Two or one rotation: the third is known
            Case (2, 1)
              Op(iYZ) = Op(iXY) .And. Op(iXZ)
            ! No rotation: check the third, and inversion if necessary
            Case (0)
              Op(iYZ) = CheckOp(iYZ, Thr)
              If (.Not. Op(iYZ)) Op(iXYZ) = CheckOp(iXYZ, Thr)
              If (Op(iXY)) Symmetry = 'yz'
              If (Op(iXYZ)) Symmetry = 'xyz'
          End Select
      End Select
      Symmetry = AdjustL(Symmetry)
      If (Symmetry .ne. '') Write(6,10) Trim(Symmetry)
      Call UpCase(Symmetry)
10    Format(6X,'Found SYMMETRY generators: ',A)
      End Subroutine DetectSym

! Function to check if a symmetry operation conserves the geometry
      Function CheckOp(Op,Thr)
      INteger, Intent(In) :: Op
      Real*8, Intent(In) :: Thr
      Logical :: CheckOp, Found
      Logical, Dimension(Size(Geom)) :: Done
      Real*8, Dimension(3) :: New
      Real*8 :: Dist
      Character (Len=MAXLEN) :: SymA, SymB, Lab, Bas
      Integer :: Num, i, j
#include "constants2.fh"
      Done = .False.
      ! For each atom, check if any of the following atoms (including itself)
      ! matches the result of the symmetry operation
      Do i=1,Size(Geom)
        If (Done(i)) Cycle
        Call SplitLabel(Geom(i)%Lab, SymA, Num, Lab, Bas)
        New = ApplySym(Op, Geom(i)%Coord)
        Found = .False.
        Do j=i,Size(Geom)
          Call SplitLabel(Geom(j)%Lab, SymB, Num, Lab, Bas)
          If (SymB .ne. SymA) Cycle
          Dist = (New(1)-Geom(j)%Coord(1))**2 +
     &           (New(2)-Geom(j)%Coord(2))**2 +
     &           (New(3)-Geom(j)%Coord(3))**2
          If (Dist*Angstrom**2 .le. Thr**2) Then
            Done(j) = .True.
            Found = .True.
            Exit
          End If
        End Do
        If (.Not. Found) Then
          CheckOp = .False.
          Return
        End If
      End Do
      CheckOp = .True.
      End Function CheckOp

! Subroutine to "fix" the symmetry, according to the threshold
! Coordinates of the atoms that are close enough are averaged
! Symmetry-superfluous atoms are marked with FileNum=0
      Subroutine AdaptSym(Thr)
      Real*8, Intent(In) :: Thr
      Real*8, Dimension(3) :: Aver, New
      Real*8 :: Dist
      Integer :: Op, nOp, Num, i, j
      Character (Len=MAXLEN) :: SymA, SymB, Lab, Bas
      Logical :: Found, Moved
#include "constants2.fh"
#include "real.fh"
      Moved = .False.
      ! Count the non-trivial operations
      nOp = Count(Oper .ne. 0)+1
      ! For each atom, find all symmetric images to average
      Do i=1,Size(Geom)
        If (Geom(i)%FileNum .eq. 0) Cycle
        Call SplitLabel(Geom(i)%Lab, SymA, Num, Lab, Bas)
        Aver = Geom(i)%Coord
        Do Op=1,7
          If (Oper(Op) .eq. 0) Exit
          Found = .False.
          Do j=i,Size(Geom)
            Call SplitLabel(Geom(j)%Lab, SymB, Num, Lab, Bas)
            If (SymB .ne. SymA) Cycle
            New = ApplySym(Oper(Op), Geom(j)%Coord)
            Dist = (New(1)-Geom(i)%Coord(1))**2 +
     &             (New(2)-Geom(i)%Coord(2))**2 +
     &             (New(3)-Geom(i)%Coord(3))**2
            If (Dist*Angstrom**2 .le. Thr**2) Then
              Found = .True.
              Aver = Aver+New
              If (j .ne. i) Geom(j)%FileNum = 0
              If (Dist .gt. Zero) Moved = .True.
              Exit
            End If
          End Do
          If (.Not. Found) Then
            Write(6,*) 'Symmetry operators do not match the geometry'
            Call Quit_OnUserError()
          End If
        End Do
        Geom(i)%Coord = Aver/Dble(nOp)
      End Do
      If (Moved)
     &  Call WarningMessage(0,
     &   'Warning! XYZ coordinates will be modified to match '//
     &   'the specified/detected symmetry. Use SYMT = 0.0 if '//
     &   'this is not desired.')
      End Subroutine AdaptSym

! Function to apply a symmetry operation to a 3D-point
! The operation is simply a possible change of sign of each axis
      Function ApplySym(Op,Coord)
      Integer, Intent(In) :: Op
      Real*8, Dimension(3), Intent(In) :: Coord
      Real*8, Dimension(3) :: ApplySym
      ApplySym = Coord
      If (iAnd(Op, iX) .gt. 0) ApplySym(1) = -ApplySym(1)
      If (iAnd(Op, iY) .gt. 0) ApplySym(2) = -ApplySym(2)
      If (iAnd(Op, iZ) .gt. 0) ApplySym(3) = -ApplySym(3)
      End Function ApplySym

! Subroutine to obtain the different components, of an atom's label
! <Sym><Num>_<Lab>.<Bas>
! Lab and Bas include the _ or .
      Subroutine SplitLabel(At,Sym,Num,Lab,Bas)
      Character (Len=*), Intent(In) :: At
      Character (Len=*), Intent(Out) :: Sym, Lab, Bas
      Integer, Intent(Out) :: Num
      Character (Len=Len(At)) :: String
      Character :: c
      Integer :: Idx, i, j
      String = At
      ! Get the Bas part
      Idx = Index(String, '.')
      If (Idx .eq. 0) Then
        Idx = Len_Trim(String)+1
        Bas = ''
      Else
        Bas = String(Idx:)
      End If
      String = String(1:Idx-1)
      ! Get the Lab part
      Idx = Index(String, '_')
      If (Idx .eq. 0) Then
        Idx = Len_Trim(String)+1
        Lab = ''
      Else
        Lab = String(Idx:)
      End If
      Sym = String(1:Idx-1)
      Num = 0
      ! Remove digits and compute Num
      Do i = Len_Trim(Sym),1,-1
        c = Sym(i:i)
        If ((IChar(c) .ge. IChar('0')) .and.
     &      (IChar(c) .le. IChar('9'))) Then
          Read(c,*) j
          Num = Num + j*10**Int(Log10(Dble(Num)))
          Sym(i:i) = ' '
        End If
      End Do
      End Subroutine SplitLabel

! Function to find the basis set that applies to a given atom
      Function FindBasis(AtSym,AtLab,AtBas)
      Character (Len=*), Intent(In) :: AtSym, AtLab, AtBas
#ifdef ALLOC_ASSIGN
      Character (Len=:), Allocatable :: FindBasis
#else
      Character (Len=MAXLEN) :: FindBasis
#endif
      Integer :: i
      ! Special case: if the Label is "MM", no basis
      If (AtLab .eq. '_MM') Then
        FindBasis = Trim(AtSym)//'...... / MM'
        Return
      End If
      ! If the atoms has a basis specified, nothing else to do
      If (AtBas .ne.  '') Then
        FindBasis = Trim(AtSym)//Trim(AtBas)
        Return
      End If
      ! Otherwise, find the basis set that matches the symbol and label
      Do i=1,Size(BasisSets,2)
        If (Trim(AtSym)//AtLab .eq. BasisSets(1,i)) Then
          FindBasis = Trim(AtSym)//'.'//Trim(BasisSets(2,i))
          Return
        End If
      End Do
      ! If none found, use the general basis
      If (BasisAll .eq. '') Then
        Write(6,*) 'No basis found for '//Trim(AtSym)//Trim(AtLab)
        Call Quit_OnUserError()
      End If
      FindBasis = Trim(AtSym)//'.'//Trim(BasisAll)
      End Function FindBasis

! Subroutine to transform a geometry according to the matrix M
      Subroutine TransformGeom(G,M)
      Type(XYZAtom), Dimension(:), Intent(InOut) :: G
      Real*8, Dimension(3,5), Intent(In) :: M
      Real*8, Dimension(3) :: Old
      Real*8, External :: DDot_
      Integer :: i, j
      ! The matrix has 5 colums: scale (1), rotate (2-4), translate (5)
      Do i=1,Size(G)
        Old(:) = M(:,1)*G(i)%Coord(:)
        Do j=1,3
          G(i)%Coord(j) = DDot_(3,Old,1,M(:,j+1),1)
        End Do
        G(i)%Coord(:) = G(i)%Coord(:)+M(:,5)
      End Do
      End Subroutine TransformGeom

      End Module XYZ
