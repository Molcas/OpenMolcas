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
#define iX_ 0
#define iY_ 1
#define iZ_ 2
#define iX 2**iX_
#define iY 2**iY_
#define iZ 2**iZ_
#define iXY iX+iY
#define iXZ iX+iZ
#define iYZ iY+iZ
#define iXYZ iX+iY+iZ

module XYZ

use Definitions, only: wp, iwp

implicit none
private

type XYZAtom
  character(len=MAXLEN) :: Lab
  real(kind=wp) :: Coord(3)
  integer(kind=iwp) :: FileNum
end type XYZAtom

integer(kind=iwp) :: FileNum = 0, Oper(7)
character(len=256) :: BasisAll, Symmetry
character(len=256), allocatable :: BasisSets(:,:)
type(XYZAtom), dimension(:), allocatable :: Geom

public :: Clear_XYZ, Out_Raw, Parse_Basis, Parse_Group, Read_XYZ, Symmetry, Write_SewInp

! Private extensions to mma interfaces

interface cptr2loff
  module procedure xyz_cptr2loff
end interface
interface mma_Allocate
  module procedure xyz_mma_allo_1D, xyz_mma_allo_1D_lim
end interface
interface mma_Deallocate
  module procedure xyz_mma_free_1D
end interface

contains

! Read an XYZ file with molcas extensions
! Rot and Trans are transformations to apply to all coordinates in this file
! If Replace=.T., the coordinates will replace existing ones, but labels will not be changed
subroutine Read_XYZ(Lu,Rot,Trans,Replace)

# ifdef _HDF5_
  use mh5, only: mh5_is_hdf5, mh5_open_file_r, mh5_fetch_attr, mh5_fetch_dset, mh5_close_file
# endif
  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Zero, One, Angstrom
  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: Lu
  real(kind=wp), allocatable, intent(in) :: Rot(:,:,:), Trans(:,:)
  logical(kind=iwp), optional, intent(in) :: Replace
  integer(kind=iwp) :: i, Idx, Error, Lxyz, NumAt
  real(kind=wp) :: Factor, Mat(3,5)
  logical(kind=iwp) :: Found, Rep
  character(len=MAXLEN) :: Line, FName, CurrDir
  type(XYZAtom), allocatable :: ThisGeom(:), TmpGeom(:)
  integer(kind=iwp), external :: IsFreeUnit
# ifdef _HDF5_
# include "Molcas.fh"
  integer(kind=iwp) :: c, Coord_id, j, nSym
  logical(kind=iwp) :: isH5
  real(kind=wp), allocatable :: Coords(:,:)
  character(len=LenIn4), allocatable :: Labels4(:)
  character(len=LenIn), allocatable :: Labels(:)
  isH5 = .false.
# endif

# include "macros.fh"
  unused_proc(mma_allocate(TmpGeom,[0,0]))

  Factor = One
  read(Lu,'(A)') Line
  Line = adjustl(Line)

  ! Try to read a number, if it fails, try to open a file
  ! Note that the slash means end-of-line in list-directed input
  NumAt = -1
  read(Line,*,iostat=Error) NumAt
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
        read(Line,*,iostat=Error) NumAt
        if (Error /= 0) then
          write(u6,*) 'Error reading file ',trim(FName)
          call Quit_OnUserError()
        end if
#     ifdef _HDF5_
      end if
#     endif
    else
      write(u6,*) 'File ',trim(FName),' not found!'
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
    write(u6,*) 'Reading xyz coordinates from h5 file '//trim(FName)
    Coord_id = mh5_open_file_r(trim(FName))
    ! check if symmetry was used
    call mh5_fetch_attr(Coord_id,'NSYM',nSym)
    ! read numbers of atoms
    if (nSym > 1) then
      call mh5_fetch_attr(Coord_id,'NATOMS_ALL',NumAt)
    else
      call mh5_fetch_attr(Coord_id,'NATOMS_UNIQUE',NumAt)
    end if
    call mma_allocate(ThisGeom,NumAt,label='ThisGeom')
    call mma_allocate(Labels,NumAt,label='Labels')
    call mma_allocate(Coords,3,NumAt,label='Coords')
    ! read atom labels
    if (nSym > 1) then
      call mma_allocate(Labels4,NumAt,label='Labels4')
      call mh5_fetch_dset(Coord_id,'DESYM_CENTER_LABELS',Labels4)
      do i=1,NumAt
        Labels(i) = Labels4(i)(1:LenIn)
      end do
      call mma_deallocate(Labels4)
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
    call mma_deallocate(Labels)
    call mma_deallocate(Coords)
  else
# endif
    !*******************************************************************
    ! For .xyz format file
    !*******************************************************************
    call mma_allocate(ThisGeom,NumAt,label='ThisGeom')
    read(Lxyz,'(A)',iostat=Error) Line
    if (Error /= 0) then
      write(u6,*) 'Error reading geometry'
      call Quit_OnUserError()
    end if
    call UpCase(Line)
    ! Units of the coordinates (default angstrom)
    if (max(index(Line,'BOHR'),index(Line,'A.U.')) <= 0) then
      Factor = One/Angstrom
    end if
    do i=1,NumAt
      read(Lxyz,*,iostat=Error) ThisGeom(i)%Lab,ThisGeom(i)%Coord(:)
      if (Error /= 0) then
        write(u6,*) 'Error reading geometry'
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
    read(Line(Idx+5:),*,iostat=Error) Mat(1,1)
    Mat(2,1) = Mat(1,1)
    Mat(3,1) = Mat(1,1)
  end if
  Idx = index(' '//Line,' SCALEX ')
  if (Idx > 0) read(Line(Idx+6:),*,iostat=Error) Mat(1,1)
  Idx = index(' '//Line,' SCALEY ')
  if (Idx > 0) read(Line(Idx+6:),*,iostat=Error) Mat(2,1)
  Idx = index(' '//Line,' SCALEZ ')
  if (Idx > 0) read(Line(Idx+6:),*,iostat=Error) Mat(3,1)
  Idx = index(' '//Line,' ROT ')
  if (Idx > 0) read(Line(Idx+3:),*,iostat=Error) Mat(:,2:4)
  Idx = index(' '//Line,' TRANS ')
  if (Idx > 0) read(Line(Idx+5:),*,iostat=Error) Mat(:,5)
  ! If Rot and Trans are given in the input, they override the inline transformations
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
        write(u6,*) 'New system size does not match previous one'
        call Quit_OnUserError()
      end if
      do i=1,size(Geom)
        Geom(i)%Coord = ThisGeom(i)%Coord
      end do
      call mma_deallocate(ThisGeom)
    else
      ! Append the just read geometry to the general one
      call move_alloc(Geom,TmpGeom)
      NumAt = size(TmpGeom)+size(ThisGeom)
      call mma_allocate(Geom,NumAt,label='Geom')
      Geom(1:size(TmpGeom)) = TmpGeom(:)
      Geom(size(TmpGeom)+1:) = ThisGeom(:)
      call mma_deallocate(TmpGeom)
      call mma_deallocate(ThisGeom)
    end if
  else
    call move_alloc(ThisGeom,Geom)
  end if

end subroutine Read_XYZ

! Clear allocations
subroutine Clear_XYZ()

  use stdalloc, only: mma_deallocate

  if (allocated(Geom)) call mma_deallocate(Geom)
  if (allocated(BasisSets)) call mma_deallocate(BasisSets)
  FileNum = 0

end subroutine Clear_XYZ

! Write an input file for seward
! Atoms belonging to GhostFiles have no charge
subroutine Write_SewInp(FName,GhostFiles)

  use Constants, only: Zero

  character(len=*), intent(in) :: FName
  integer(kind=iwp), intent(in) :: GhostFiles(:)
  integer(kind=iwp) :: i, j, Lu, Num
  logical(kind=iwp) :: Ghost
  character(len=MAXLEN) :: Bas, Lab, New, Old, Sym
  integer(kind=iwp), external :: IsFreeUnit

  Lu = IsFreeUnit(10)
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
      if (Ghost) write(Lu,30) 'Charge',Zero
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

  integer(kind=iwp) :: Out_Raw
  real(kind=wp), intent(inout) :: Array(*)
  integer(kind=iwp) :: i, j

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

  use fortran_strings
  use stdalloc, only: mma_allocate, mma_deallocate

  character(len=*), intent(in) :: Basis
  integer(kind=iwp) :: i, Idx, IdxDot, Next, Num

  ! Count number of commas
  Num = count(char_array(trim(Basis)) == ',')+1
  BasisAll = ''
  if (allocated(BasisSets)) call mma_deallocate(BasisSets)
  call mma_allocate(BasisSets,2,Num,label='BasisSets')
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

  character(len=*), intent(in) :: Group
  real(kind=wp), intent(in) :: Thr
  integer(kind=iwp) :: Error, i, j, k
  character(len=3) :: Gen(3)

  Symmetry = Group
  call UpCase(Symmetry)
  if (Symmetry == 'FULL') then
    call DetectSym(Thr)
  else if ((Symmetry(1:5) == 'NOSYM') .or. (Symmetry == 'E') .or. (Symmetry == 'C1')) then
    Symmetry = ''
  end if
  Gen = ''
  read(Symmetry,*,iostat=Error) Gen
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

  use Definitions, only: u6

  real(kind=wp), intent(in) :: Thr
  logical(kind=iwp) :: Op(7)

  Op = .false.
  Op(iX) = CheckOp(iX,Thr)
  Op(iY) = CheckOp(iY,Thr)
  Op(iZ) = CheckOp(iZ,Thr)
  Symmetry = ''
  select case (count(Op))
    case (3,2)
      ! Two or three reflections: all is known
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
    case (1)
      ! One reflection: possibly inversion and complementary rotation
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
    case (0)
      ! No reflection: check rotations (only two initially)
      Op(iXY) = CheckOp(iXY,Thr)
      Op(iXZ) = CheckOp(iXZ,Thr)
      if (Op(iXY)) Symmetry = 'xy'
      if (Op(iXZ)) Symmetry = trim(Symmetry)//' xz'
      select case (count(Op))
        case (2,1)
          ! Two or one rotation: the third is known
          Op(iYZ) = Op(iXY) .and. Op(iXZ)
        case (0)
          ! No rotation: check the third, and inversion if necessary
          Op(iYZ) = CheckOp(iYZ,Thr)
          if (.not. Op(iYZ)) Op(iXYZ) = CheckOp(iXYZ,Thr)
          if (Op(iXY)) Symmetry = 'yz'
          if (Op(iXYZ)) Symmetry = 'xyz'
      end select
  end select
  Symmetry = adjustl(Symmetry)
  if (Symmetry /= '') write(u6,10) trim(Symmetry)
  call UpCase(Symmetry)

  10 format(6X,'Found SYMMETRY generators: ',A)

end subroutine DetectSym

! Function to check if a symmetry operation conserves the geometry
function CheckOp(Op,Thr)

  use stdalloc, only: mma_allocate, mma_deallocate
  use Constants, only: Angstrom

  integer(kind=iwp), intent(in) :: Op
  real(kind=wp), intent(in) :: Thr
  integer(kind=iwp) :: i, j, Num
  real(kind=wp) :: Dist, New(3)
  logical(kind=iwp) :: CheckOp, Found
  character(len=MAXLEN) :: Bas, Lab, Sym, SymA, SymB
  logical(kind=iwp), allocatable :: Done(:)

  call mma_allocate(Done,size(Geom),label='Done')
  Done = .false.
  ! For each atom, check if any of the following atoms (including itself)
  ! matches the result of the symmetry operation
  Found = .false.
  outer: do i=1,size(Geom)
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
        cycle outer
      end if
    end do
    if (.not. Found) exit
  end do outer
  CheckOp = Found
  call mma_deallocate(Done)

end function CheckOp

! Subroutine to "fix" the symmetry, according to the threshold
! Coordinates of the atoms that are close enough are averaged
! Symmetry-superfluous atoms are marked with FileNum=0
subroutine AdaptSym(Thr)

  use Constants, only: Zero, Angstrom
  use Definitions, only: wp, u6

  real(kind=wp), intent(in) :: Thr
  integer(kind=iwp) :: i, j, nOp, Num, Op
  real(kind=wp) :: Aver(3), Dist, New(3)
  logical(kind=iwp) :: Found, Moved, ZeroAxis(3)
  character(len=MAXLEN) :: Bas, Lab, Sym, SymA, SymB

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
            if (btest(Oper(Op),iX_)) ZeroAxis(1) = .true.
            if (btest(Oper(Op),iY_)) ZeroAxis(2) = .true.
            if (btest(Oper(Op),iZ_)) ZeroAxis(3) = .true.
          else
            Geom(j)%FileNum = 0
          end if
          if (Dist > Zero) Moved = .true.
          exit
        end if
      end do
      if (.not. Found) then
        write(u6,*) 'Symmetry operators do not match the geometry'
        call Quit_OnUserError()
      end if
    end do
    ! Make sure the axes that should be zero are exactly zero
    Geom(i)%Coord = merge([Zero,Zero,Zero],Aver/real(nOp,kind=wp),ZeroAxis)
  end do
  if (Moved) call WarningMessage(0,'Warning! XYZ coordinates will be modified to match the specified/detected symmetry. '// &
                                 'Use SYMT = 0.0 if this is not desired.')

end subroutine AdaptSym

! Function to apply a symmetry operation to a 3D-point
! The operation is simply a possible change of sign of each axis
function ApplySym(Op,Coord)

  integer(kind=iwp), intent(in) :: Op
  real(kind=wp), intent(in) :: Coord(3)
  real(kind=wp) :: ApplySym(3)

  ApplySym = Coord
  if (btest(Op,iX_)) ApplySym(1) = -ApplySym(1)
  if (btest(Op,iY_)) ApplySym(2) = -ApplySym(2)
  if (btest(Op,iZ_)) ApplySym(3) = -ApplySym(3)

end function ApplySym

! Subroutine to obtain the different components, of an atom's label
! <Sym><Num>_<Lab>.<Bas>
! Lab and Bas include the _ or .
subroutine SplitLabel(At,Sym,Num,Lab,Bas)

  character(len=*), intent(in) :: At
  character(len=*), intent(out) :: Sym, Lab, Bas
  integer(kind=iwp), intent(out) :: Num
  integer(kind=iwp) :: i, Idx, j, k
  character(len=len(At)) :: String
  character :: c

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

  use Definitions, only: u6

  character(len=*), intent(in) :: AtSym, AtLab, AtBas
  integer(kind=iwp) :: i
  character(len=len(AtSym)) :: UpSym
  character(len=len(AtLab)) :: UpLab
# ifdef ALLOC_ASSIGN
  character(len=:), allocatable :: FindBasis
# else
  character(len=MAXLEN) :: FindBasis
# endif

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
    write(u6,*) 'No basis found for '//trim(AtSym)//trim(AtLab)
    call Quit_OnUserError()
  end if
  FindBasis = trim(AtSym)//'.'//trim(BasisAll)

end function FindBasis

! Subroutine to transform a geometry according to the matrix M
subroutine TransformGeom(G,M)

  type(XYZAtom), dimension(:), intent(inout) :: G
  real(kind=wp), intent(in) :: M(3,5)
  real(kind=wp) :: Old(3)
  integer(kind=iwp) :: i, j
  real(kind=wp), external :: DDot_

  ! The matrix has 5 colums: scale (1), rotate (2-4), translate (5)
  do i=1,size(G)
    Old(:) = M(:,1)*G(i)%Coord(:)
    do j=1,3
      G(i)%Coord(j) = DDot_(3,Old,1,M(:,j+1),1)
    end do
    G(i)%Coord(:) = G(i)%Coord(:)+M(:,5)
  end do

end subroutine TransformGeom

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define xyz_cptr2loff, xyz_mma_allo_1D, xyz_mma_allo_1D_lim, xyz_mma_free_1D
#define _TYPE_ type(XYZAtom)
#  define _FUNC_NAME_ xyz_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ xyz_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'xyz_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module XYZ
