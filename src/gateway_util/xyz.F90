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
! Copyright (C) 2016, Morgane Vacher                                   *
!               2017,2018, Ignacio Fdez. Galvan                        *
!***********************************************************************

! Module to handle the COORD keyword in gateway/seward

#include "compiler_features.h"
#define MAXLEN 180
#define iX 1
#define iY 2
#define iZ 4
#define iXY iX+iY
#define iXZ iX+iZ
#define iYZ iY+iZ
#define iXYZ iX+iY+iZ

module XYZ

implicit none
private

type XYZAtom
  character(Len=MAXLEN) :: Lab
  real*8, dimension(3) :: Coord
  integer :: FileNum
end type XYZAtom
type(XYZAtom), dimension(:), allocatable :: Geom
character(Len=256) :: BasisAll, Symmetry
character(Len=256), dimension(:,:), allocatable :: BasisSets
integer, dimension(7) :: Oper
integer :: FileNum = 0

public :: Read_XYZ, Write_SewInp, Clear_XYZ, Out_Raw, Parse_Basis, Parse_Group, Symmetry

contains

! Read an XYZ file with molcas extensions
! Rot and Trans are transformations to apply to all coordinates in this file
! If Replace=.T., the coordinates will replace existing ones, but labels will not be changed
subroutine Read_XYZ(Lu,Rot,Trans,Replace)

# ifdef _HDF5_
  use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_fetch_attr, mh5_fetch_dset, mh5_close_file
# endif

  integer, intent(In) :: Lu
  real*8, dimension(:,:,:), allocatable, intent(In) :: Rot
  real*8, dimension(:,:), allocatable, intent(In) :: Trans
  logical, optional, intent(In) :: Replace
  character(Len=MAXLEN) :: Line, FName, CurrDir
  integer :: Error, NumAt, i, Lxyz, Idx
  logical :: Found, Rep
  integer, external :: IsFreeUnit
  type(XYZAtom), dimension(:), allocatable :: ThisGeom, TmpGeom
  real*8, dimension(3,5) :: Mat
  real*8 :: Factor
# include "real.fh"
# include "constants2.fh"
# ifdef _HDF5_
# include "Molcas.fh"
  logical :: isH5
  integer :: Coord_id, nSym, j, c
  character(Len=LenIn), dimension(:), allocatable :: Labels
  character(Len=LenIn4), dimension(:), allocatable :: Labels4
  real*8, dimension(:,:), allocatable :: Coords
  isH5 = .false.
# endif

  Factor = One
  read(Lu,'(A)') Line
  Line = adjustl(Line)

  ! Try to read a number, if it fails, try to open a file
  ! Note that the slash means end-of-line in list-directed input
  NumAt = -1
  read(Line,*,IOStat=Error) NumAt
  ! And make sure there's nothing but numbers in the first word
  ! (the above is not 100% reliable in some compilers)
  do i=1,len(Line)
    if (Line(i:i) == ' ') exit
    Idx = ichar(Line(i:i))
    if ((Idx < ichar('0')) .or. (Idx > ichar('9'))) then
      Error = -1
      exit
    end if
  end do
  if (NumAt < 1) Error = -1
  if (index(Line,'/') > 0) Error = -1
  Lxyz = Lu
  if (Error /= 0) then
    read(Line,'(A)') FName
    Found = .false.
    if (index(Line,'/') == 0) then
      call GetEnvF('CurrDir',CurrDir)
      CurrDir = trim(CurrDir)//'/'//FName
      call F_Inquire(CurrDir,Found)
      if (Found) FName = CurrDir
    end if
    if (.not. Found) call F_Inquire(FName,Found)
    if (Found) then
#     ifdef _HDF5_
      if (mh5_is_hdf5(trim(FName))) then
        isH5 = .true.
      else
#     endif
        Lxyz = IsFreeUnit(Lxyz)
        call Molcas_Open(Lxyz,FName)
        read(Lxyz,'(A)') Line
        read(Line,*,IOStat=Error) NumAt
        if (Error /= 0) then
          write(6,*) 'Error reading file ',trim(FName)
          call Quit_OnUserError()
        end if
#     ifdef _HDF5_
      end if
#     endif
    else
      write(6,*) 'File ',trim(FName),' not found!'
      call Quit_OnUserError()
    end if
  end if
  ! Each atom has a field indicating the sequential number of the
  ! file it belongs to
  FileNum = FileNum+1

# ifdef _HDF5_
  if (isH5) then
    !*******************************************************************
    ! For HDF5 formatted file
    !*******************************************************************
    write(6,*) 'Reading xyz coordinates from h5 file '//trim(FName)
    Coord_id = mh5_open_file_r(trim(FName))
    ! check if symmetry was used
    call mh5_fetch_attr(Coord_id,'NSYM',nSym)
    ! read numbers of atoms
    if (nSym > 1) then
      call mh5_fetch_attr(Coord_id,'NATOMS_ALL',NumAt)
    else
      call mh5_fetch_attr(Coord_id,'NATOMS_UNIQUE',NumAt)
    end if
    allocate(ThisGeom(NumAt))
    allocate(Labels(NumAt),Coords(3,NumAt))
    ! read atom labels
    if (nSym > 1) then
      allocate(Labels4(NumAt))
      call mh5_fetch_dset(Coord_id,'DESYM_CENTER_LABELS',Labels4)
      do i=1,NumAt
        Labels(i) = Labels4(i)(1:LenIn)
      end do
      deallocate(Labels4)
    else
      call mh5_fetch_dset(Coord_id,'CENTER_LABELS',Labels)
    end if
    ! remove numbers from labels
    ! (otherwise there will be problems when symmetric structures
    ! are used without symmetry)
    do i=1,NumAt
      do j=1,LenIn
        c = ichar(Labels(i)(j:j))
        if ((c >= ichar('0')) .and. (c <= ichar('9'))) then
          Labels(i)(j:j) = ' '
        end if
      end do
    end do
    ! read atom coordinates
    if (nSym > 1) then
      call mh5_fetch_dset(Coord_id,'DESYM_CENTER_COORDINATES',Coords)
    else
      call mh5_fetch_dset(Coord_id,'CENTER_COORDINATES',Coords)
    end if
    call mh5_close_file(Coord_id)
    ! store data
    do i=1,NumAt
      ThisGeom(i)%Lab = Labels(i)
      ThisGeom(i)%Coord(:) = Coords(:,i)
      ThisGeom(i)%FileNum = FileNum
    end do
    deallocate(Labels,Coords)
  else
# endif
    !*******************************************************************
    ! For .xyz format file
    !*******************************************************************
    allocate(ThisGeom(NumAt))
    read(Lxyz,'(A)',IOStat=Error) Line
    if (Error /= 0) then
      write(6,*) 'Error reading geometry'
      call Quit_OnUserError()
    end if
    call UpCase(Line)
    ! Units of the coordinates (default angstrom)
    if (max(index(Line,'BOHR'),index(Line,'A.U.')) <= 0) then
      Factor = One/Angstrom
    end if
    do i=1,NumAt
      read(Lxyz,*,IOStat=Error) ThisGeom(i)%Lab,ThisGeom(i)%Coord(:)
      if (Error /= 0) then
        write(6,*) 'Error reading geometry'
        call Quit_OnUserError()
      end if
      ThisGeom(i)%Coord(:) = ThisGeom(i)%Coord(:)*Factor
      ThisGeom(i)%FileNum = FileNum
    end do
    if (Lxyz /= Lu) close(Lxyz)
# ifdef _HDF5_
  end if
# endif

  ! Obtain/read transformation matrix and transform the geometry in this file
  Mat = reshape([One,One,One,One,Zero,Zero,Zero,One,Zero,Zero,Zero,One,Zero,Zero,Zero],shape(Mat))
  Idx = index(' '//Line,' SCALE ')
  if (Idx > 0) then
    read(Line(Idx+5:),*,IOStat=Error) Mat(1,1)
    Mat(2,1) = Mat(1,1)
    Mat(3,1) = Mat(1,1)
  end if
  Idx = index(' '//Line,' SCALEX ')
  if (Idx > 0) read(Line(Idx+6:),*,IOStat=Error) Mat(1,1)
  Idx = index(' '//Line,' SCALEY ')
  if (Idx > 0) read(Line(Idx+6:),*,IOStat=Error) Mat(2,1)
  Idx = index(' '//Line,' SCALEZ ')
  if (Idx > 0) read(Line(Idx+6:),*,IOStat=Error) Mat(3,1)
  Idx = index(' '//Line,' ROT ')
  if (Idx > 0) read(Line(Idx+3:),*,IOStat=Error) Mat(:,2:4)
  Idx = index(' '//Line,' TRANS ')
  if (Idx > 0) read(Line(Idx+5:),*,IOStat=Error) Mat(:,5)
  ! If Rot and Trans are given in the input, they override
  ! the inline transformations
  if (allocated(Rot)) then
    Mat(:,2:4) = Rot(:,:,FileNum)
  end if
  Mat(:,2:4) = transpose(Mat(:,2:4))
  if (allocated(Trans)) then
    Mat(:,5) = Trans(:,FileNum)
  end if
  Mat(:,5) = Mat(:,5)*Factor
  call TransformGeom(ThisGeom,Mat)

  Rep = .false.
  if (present(Replace)) Rep = Replace

  if (allocated(Geom)) then
    if (Rep) then
      if (size(ThisGeom) /= size(Geom)) then
        write(6,*) 'New system size does not match previous one'
        call Quit_OnUserError()
      end if
      do i=1,size(Geom)
        Geom(i)%Coord = ThisGeom(i)%Coord
      end do
    else
      ! Append the just read geometry to the general one
      call move_alloc(Geom,TmpGeom)
      NumAt = size(TmpGeom)+size(ThisGeom)
      allocate(Geom(NumAt))
      Geom(1:size(TmpGeom)) = TmpGeom(:)
      Geom(size(TmpGeom)+1:) = ThisGeom(:)
      deallocate(TmpGeom)
      deallocate(ThisGeom)
    end if
  else
    call move_alloc(ThisGeom,Geom)
  end if

end subroutine Read_XYZ

! Clear allocations
subroutine Clear_XYZ

  if (allocated(Geom)) deallocate(Geom)
  if (allocated(BasisSets)) deallocate(BasisSets)
  FileNum = 0

end subroutine Clear_XYZ

! Write an input file for seward
! Atoms belonging to GhostFiles have no charge
subroutine Write_SewInp(FName,GhostFiles)

  character(Len=*), intent(In) :: FName
  integer, dimension(:), intent(In) :: GhostFiles
  integer :: Lu, Num, i, j
  integer, external :: IsFreeUnit
  character(Len=MAXLEN) :: New, Old, Sym, Lab, Bas
  logical :: Ghost

  Lu = 10
  Lu = IsFreeUnit(Lu)
  call Molcas_Open(Lu,FName)
  ! Write symmetry
  if (Symmetry /= '') write(Lu,40) 'Symmetry',trim(Symmetry)
  ! Write the atoms, skipping symmetry-superfluous atoms
  Old = ''
  j = 0
  do i=1,size(Geom)
    if (Geom(i)%FileNum == 0) cycle
    j = j+1
    call SplitLabel(Geom(i)%Lab,Sym,Num,Lab,Bas)
    ! Build an identifier for the basis
    New = trim(Sym)//trim(Lab)//Bas
    Ghost = any(GhostFiles == Geom(i)%FileNum)
    if (Ghost) New = trim(New)//' *'
    ! If there is a change of basis, write the basis
    if (New /= Old) then
      if (i > 1) write(Lu,10) 'End of Basis'
      write(Lu,10) 'Basis Set'
      write(Lu,10) FindBasis(Sym,Lab,Bas)
      if (Ghost) write(Lu,30) 'Charge',0.0
      Old = New
    end if
    ! Add a sequential number if the label had none
    if (Num == 0) Num = j
    write(Lu,20) trim(Sym),Num,Geom(i)%Coord(:)
  end do
  write(Lu,10) 'End of Basis'
  write(Lu,10) 'End of Coord'
  close(Lu)
10 format(A)
20 format(A,I0,3(1X,ES30.20E3))
30 format(A,/,F4.1)
40 format(A,/,A)

end subroutine Write_SewInp

! Store symmetry-unique atom coordinates in an array
! Output is the number of saved coordinates
function Out_Raw(Array)

  real*8, dimension(*), intent(InOut) :: Array
  integer :: Out_Raw
  integer :: i, j

  j = 0
  do i=1,size(Geom)
    if (Geom(i)%FileNum == 0) cycle
    Array(j+1:j+3) = Geom(i)%Coord
    j = j+3
  end do
  Out_Raw = j

end function Out_Raw

! Parse the BASIS keyword to store the basis set corresponding to
! each atom + label
subroutine Parse_Basis(Basis)

  character(Len=*), intent(In) :: Basis
  integer :: Num, Idx, IdxDot, Next, i

  ! Count number of commas
  Num = count(transfer(Basis,'x',len_trim(Basis)) == ',')+1
  BasisAll = ''
  if (allocated(BasisSets)) deallocate(BasisSets)
  allocate(BasisSets(2,Num))
  ! For each comma-separated word, split it at the first dot
  ! If the first part is empty, use it as a general basis set
  Idx = 0
  do i=1,Num
    Next = Idx+index(Basis(Idx+1:),',')
    if (Next == Idx) Next = len_trim(Basis)+1
    BasisSets(2,i) = Basis(Idx+1:Next-1)
    IdxDot = index(BasisSets(2,i),'.')
    if (IdxDot == 0) then
      BasisSets(1,i) = ''
    else
      BasisSets(1,i) = adjustl(BasisSets(2,i)(1:IdxDot-1))
      BasisSets(2,i)(1:IdxDot) = ''
    end if
    call UpCase(BasisSets(1,i))
    BasisSets(2,i) = adjustl(BasisSets(2,i))
    if (BasisSets(1,i) == '') BasisAll = BasisSets(2,i)
    Idx = Next
  end do
  if (BasisAll == '') BasisAll = 'ANO-S-MB'

end subroutine Parse_Basis

! Parse the GROUP keyword, detecting the symmetry,
! generating all symmetry operations and adapting
! the symmetry
subroutine Parse_Group(Group,Thr)

  character(Len=*), intent(In) :: Group
  real*8, intent(In) :: Thr
  character(Len=3), dimension(3) :: Gen
  integer :: Error, i, j, k

  Symmetry = Group
  call UpCase(Symmetry)
  if (Symmetry == 'FULL') then
    call DetectSym(Thr)
  else if ((Symmetry(1:5) == 'NOSYM') .or. (Symmetry == 'E') .or. (Symmetry == 'C1')) then
    Symmetry = ''
  end if
  Gen = ''
  read(Symmetry,*,IOStat=Error) Gen
  ! Encode generators
  Oper = 0
  do i=1,3
    if (index(Gen(i),'X') > 0) Oper(i) = Oper(i)+iX
    if (index(Gen(i),'Y') > 0) Oper(i) = Oper(i)+iY
    if (index(Gen(i),'Z') > 0) Oper(i) = Oper(i)+iZ
  end do
  ! Create all combinations
  Oper(4) = ieor(Oper(1),Oper(2))
  Oper(5) = ieor(Oper(1),Oper(3))
  Oper(6) = ieor(Oper(2),Oper(3))
  Oper(7) = ieor(Oper(3),Oper(4))
  ! Remove duplicates
  do i=1,7
    do j=i+1,7
      if (Oper(j) == Oper(i)) Oper(j) = 0
    end do
  end do
  ! Sort descending
  do i=1,7
    do j=i+1,7
      if (Oper(j) >= Oper(i)) then
        k = Oper(j)
        Oper(j) = Oper(i)
        Oper(i) = k
      end if
    end do
  end do
  ! Check and adapt the symmetry
  call AdaptSym(Thr)

end subroutine Parse_Group

! Private procedures follow

! Subroutine to detect the symmetry elements of the system
! The system is *not* translated or rotated
! Only elements of D2h group are tested (i.e. sign inversions of one, two or three axes)
! The result is stored in Symmetry (compatible with molcas-extra)
subroutine DetectSym(Thr)

  real*8, intent(In) :: Thr
  logical, dimension(7) :: Op

  Op = .false.
  Op(iX) = CheckOp(iX,Thr)
  Op(iY) = CheckOp(iY,Thr)
  Op(iZ) = CheckOp(iZ,Thr)
  Symmetry = ''
  select case (count(Op))
    ! Two or three reflections: all is known
    case (3,2)
      Op(iXY) = Op(iX) .and. Op(iY)
      Op(iXZ) = Op(iX) .and. Op(iZ)
      Op(iYZ) = Op(iY) .and. Op(iZ)
      Op(iXYZ) = Op(iX) .and. Op(iY) .and. Op(iZ)
      if (Op(iXYZ)) then
        Symmetry = 'x y z'
      else
        if (Op(iXY)) Symmetry = 'xy y'
        if (Op(iXZ)) Symmetry = 'xz z'
        if (Op(iYZ)) Symmetry = 'yz z'
      end if
    ! One reflection: possibly inversion and complementary rotation
    case (1)
      Op(iXYZ) = CheckOp(iXYZ,Thr)
      if (Op(iXYZ)) then
        Op(iXY) = Op(iZ)
        Op(iXZ) = Op(iY)
        Op(iYZ) = Op(iX)
        if (Op(iXY)) Symmetry = 'xy xyz'
        if (Op(iXZ)) Symmetry = 'xz xyz'
        if (Op(iYZ)) Symmetry = 'yz xyz'
      else
        if (Op(iX)) Symmetry = 'x'
        if (Op(iY)) Symmetry = 'y'
        if (Op(iZ)) Symmetry = 'z'
      end if
    ! No reflection: check rotations (only two initially)
    case (0)
      Op(iXY) = CheckOp(iXY,Thr)
      Op(iXZ) = CheckOp(iXZ,Thr)
      if (Op(iXY)) Symmetry = 'xy'
      if (Op(iXZ)) Symmetry = trim(Symmetry)//' xz'
      select case (count(Op))
        ! Two or one rotation: the third is known
        case (2,1)
          Op(iYZ) = Op(iXY) .and. Op(iXZ)
        ! No rotation: check the third, and inversion if necessary
        case (0)
          Op(iYZ) = CheckOp(iYZ,Thr)
          if (.not. Op(iYZ)) Op(iXYZ) = CheckOp(iXYZ,Thr)
          if (Op(iXY)) Symmetry = 'yz'
          if (Op(iXYZ)) Symmetry = 'xyz'
      end select
  end select
  Symmetry = adjustl(Symmetry)
  if (Symmetry /= '') write(6,10) trim(Symmetry)
  call UpCase(Symmetry)
10 format(6X,'Found SYMMETRY generators: ',A)

end subroutine DetectSym

! Function to check if a symmetry operation conserves the geometry
function CheckOp(Op,Thr)

  integer, intent(In) :: Op
  real*8, intent(In) :: Thr
  logical :: CheckOp, Found
  logical, dimension(size(Geom)) :: Done
  real*8, dimension(3) :: New
  real*8 :: Dist
  character(Len=MAXLEN) :: Sym, SymA, SymB, Lab, Bas
  integer :: Num, i, j
# include "constants2.fh"

  Done = .false.
  ! For each atom, check if any of the following atoms (including itself)
  ! matches the result of the symmetry operation
  do i=1,size(Geom)
    if (Done(i)) cycle
    call SplitLabel(Geom(i)%Lab,Sym,Num,Lab,Bas)
    SymA = trim(Sym)//trim(Lab)
    call UpCase(SymA)
    New = ApplySym(Op,Geom(i)%Coord)
    Found = .false.
    do j=i,size(Geom)
      call SplitLabel(Geom(j)%Lab,Sym,Num,Lab,Bas)
      SymB = trim(Sym)//trim(Lab)
      call UpCase(SymB)
      if (SymB /= SymA) cycle
      Dist = (New(1)-Geom(j)%Coord(1))**2+(New(2)-Geom(j)%Coord(2))**2+(New(3)-Geom(j)%Coord(3))**2
      if (Dist*Angstrom**2 <= Thr**2) then
        Done(j) = .true.
        Found = .true.
        exit
      end if
    end do
    if (.not. Found) then
      CheckOp = .false.
      return
    end if
  end do
  CheckOp = .true.

end function CheckOp

! Subroutine to "fix" the symmetry, according to the threshold
! Coordinates of the atoms that are close enough are averaged
! Symmetry-superfluous atoms are marked with FileNum=0
subroutine AdaptSym(Thr)

  real*8, intent(In) :: Thr
  real*8, dimension(3) :: Aver, New
  real*8 :: Dist
  integer :: Op, nOp, Num, i, j
  character(Len=MAXLEN) :: Sym, SymA, SymB, Lab, Bas
  logical :: Found, Moved
  logical, dimension(3) :: ZeroAxis
# include "constants2.fh"
# include "real.fh"

  Moved = .false.
  ! Count the non-trivial operations
  nOp = count(Oper /= 0)+1
  ! For each atom, find all symmetric images to average
  do i=1,size(Geom)
    if (Geom(i)%FileNum == 0) cycle
    call SplitLabel(Geom(i)%Lab,Sym,Num,Lab,Bas)
    SymA = trim(Sym)//trim(Lab)
    call UpCase(SymA)
    Aver = Geom(i)%Coord
    ZeroAxis = .false.
    do Op=1,7
      if (Oper(Op) == 0) exit
      Found = .false.
      do j=i,size(Geom)
        call SplitLabel(Geom(j)%Lab,SymB,Num,Lab,Bas)
        SymB = trim(Sym)//trim(Lab)
        call UpCase(SymB)
        if (SymB /= SymA) cycle
        New = ApplySym(Oper(Op),Geom(j)%Coord)
        Dist = (New(1)-Geom(i)%Coord(1))**2+(New(2)-Geom(i)%Coord(2))**2+(New(3)-Geom(i)%Coord(3))**2
        if (Dist*Angstrom**2 <= Thr**2) then
          Found = .true.
          Aver = Aver+New
          if (j == i) then
            if (iand(Oper(Op),iX) > 0) ZeroAxis(1) = .true.
            if (iand(Oper(Op),iY) > 0) ZeroAxis(2) = .true.
            if (iand(Oper(Op),iZ) > 0) ZeroAxis(3) = .true.
          else
            Geom(j)%FileNum = 0
          end if
          if (Dist > Zero) Moved = .true.
          exit
        end if
      end do
      if (.not. Found) then
        write(6,*) 'Symmetry operators do not match the geometry'
        call Quit_OnUserError()
      end if
    end do
    ! Make sure the axes that should be zero are exactly zero
    Geom(i)%Coord = merge([Zero,Zero,Zero],Aver/dble(nOp),ZeroAxis)
  end do
  if (Moved) call WarningMessage(0,'Warning! XYZ coordinates will be modified to match the specified/detected symmetry. '// &
                                 'Use SYMT = 0.0 if this is not desired.')

end subroutine AdaptSym

! Function to apply a symmetry operation to a 3D-point
! The operation is simply a possible change of sign of each axis
function ApplySym(Op,Coord)

  integer, intent(In) :: Op
  real*8, dimension(3), intent(In) :: Coord
  real*8, dimension(3) :: ApplySym

  ApplySym = Coord
  if (iand(Op,iX) > 0) ApplySym(1) = -ApplySym(1)
  if (iand(Op,iY) > 0) ApplySym(2) = -ApplySym(2)
  if (iand(Op,iZ) > 0) ApplySym(3) = -ApplySym(3)

end function ApplySym

! Subroutine to obtain the different components, of an atom's label
! <Sym><Num>_<Lab>.<Bas>
! Lab and Bas include the _ or .
subroutine SplitLabel(At,Sym,Num,Lab,Bas)

  character(Len=*), intent(In) :: At
  character(Len=*), intent(Out) :: Sym, Lab, Bas
  integer, intent(Out) :: Num
  character(Len=len(At)) :: String
  character :: c
  integer :: Idx, i, j, k

  String = At
  ! Get the Bas part
  Idx = index(String,'.')
  if (Idx == 0) then
    Idx = len_trim(String)+1
    Bas = ''
  else
    Bas = String(Idx:)
  end if
  String = String(1:Idx-1)
  ! Get the Lab part
  Idx = index(String,'_')
  if (Idx == 0) then
    Idx = len_trim(String)+1
    Lab = ''
  else
    Lab = String(Idx:)
  end if
  Sym = String(1:Idx-1)
  Num = 0
  k = 0
  ! Remove digits and compute Num
  do i=len_trim(Sym),1,-1
    c = Sym(i:i)
    if ((ichar(c) >= ichar('0')) .and. (ichar(c) <= ichar('9'))) then
      read(c,*) j
      Num = Num+j*10**k
      k = k+1
      Sym(i:i) = ' '
    end if
  end do

end subroutine SplitLabel

! Function to find the basis set that applies to a given atom
function FindBasis(AtSym,AtLab,AtBas)

  character(Len=*), intent(In) :: AtSym, AtLab, AtBas
# ifdef ALLOC_ASSIGN
  character(Len=:), allocatable :: FindBasis
# else
  character(Len=MAXLEN) :: FindBasis
# endif
  character(Len=len(AtSym)) :: UpSym
  character(Len=len(AtLab)) :: UpLab
  integer :: i

  ! Prepare for case-insensitive comparisons
  UpSym = AtSym
  call UpCase(UpSym)
  UpLab = AtLab
  call UpCase(UpLab)
  ! Special case: if the Label is "MM", no basis
  if (UpLab == '_MM') then
    FindBasis = trim(AtSym)//'...... / MM'
    return
  end if
  ! If the atoms has a basis specified, nothing else to do
  if (AtBas /= '') then
    FindBasis = trim(AtSym)//trim(AtBas)
    return
  end if
  ! Otherwise, find the basis set that matches the symbol and label
  do i=1,size(BasisSets,2)
    if (trim(UpSym)//UpLab == BasisSets(1,i)) then
      FindBasis = trim(AtSym)//'.'//trim(BasisSets(2,i))
      return
    end if
  end do
  ! If none found, use the general basis
  if (BasisAll == '') then
    write(6,*) 'No basis found for '//trim(AtSym)//trim(AtLab)
    call Quit_OnUserError()
  end if
  FindBasis = trim(AtSym)//'.'//trim(BasisAll)

end function FindBasis

! Subroutine to transform a geometry according to the matrix M
subroutine TransformGeom(G,M)

  type(XYZAtom), dimension(:), intent(InOut) :: G
  real*8, dimension(3,5), intent(In) :: M
  real*8, dimension(3) :: Old
  real*8, external :: DDot_
  integer :: i, j

  ! The matrix has 5 colums: scale (1), rotate (2-4), translate (5)
  do i=1,size(G)
    Old(:) = M(:,1)*G(i)%Coord(:)
    do j=1,3
      G(i)%Coord(j) = DDot_(3,Old,1,M(:,j+1),1)
    end do
    G(i)%Coord(:) = G(i)%Coord(:)+M(:,5)
  end do

end subroutine TransformGeom

end module XYZ
