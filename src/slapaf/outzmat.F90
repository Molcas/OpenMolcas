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
! Copyright (C) 2007, Giovanni Ghigo                                   *
!***********************************************************************

subroutine OutZMAT(nAtoms,XYZCoords,N_ZMAT)
!***********************************************************************
! Author: Giovanni Ghigo                                               *
!         Torino (Italy)  February-March 2007                          *
!                                                                      *
! This subroutine generates the nuclear coordinates in Z-matrix        *
! format using the connection indices recovered from RunFile.          *
! Remember:                                                            *
!           X dummy atoms - NAT(i)= 0 - are included in SEWARD         *
!           Z ghost atoms - NAT(i)=-1 - are NOT included in SEWARD     *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Angstrom, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, N_ZMAT
real(kind=wp), intent(inout) :: XYZCoords(3,nAtoms)
integer(kind=iwp) :: i, iAtom, iAX, iBonded, iRX, iTX, j, LuWr
real(kind=wp) :: alpha, arccos, bond, Bt(3,4), CTX(3,4), dBt(3,4,3,4), dMaxTrasl, dvec2, dXYZ(3), dXYZ2(3), prod, theta
logical(kind=iwp) :: IfTest
character(len=8) :: Label
character(len=5), allocatable :: Symbols(:)
integer(kind=iwp), allocatable :: iZmat(:,:), NAT(:)
real(kind=wp), allocatable :: ZMAT(:,:), ZMATCoords(:,:)
real(kind=wp), parameter :: ThrsTrasl = One  ! Threshold for warning

call mma_allocate(Symbols,N_ZMAT,Label='Symbols')
call mma_allocate(NAT,N_ZMAT,Label='NAT')
call mma_allocate(ZMATCoords,3,N_ZMAT,Label='ZMATCoords')
call mma_allocate(ZMAT,N_ZMAT,3,Label='ZMAT')

#ifdef _DEBUGPRINT_
IfTest = .true.
#else
IfTest = .false.
#endif

LuWr = u6
dMaxTrasl = Zero
Label = ' '
ZMAT(:,:) = Zero
ZMATCoords(:,:) = Zero
call Get_cArray('Symbol ZMAT',Symbols,N_ZMAT*5)
call mma_allocate(iZmat,3,N_ZMAT,label='iZmat')
call Get_iArray('Index ZMAT',iZmat,N_ZMAT*3)
call Get_iArray('NAT ZMAT',NAT,N_ZMAT)

if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'OutZMat - Z-Matrix Data :'
  write(LuWr,*) '          N_ZMAT=',N_ZMAT
  write(LuWr,*) '   Label  NA      i   j   k'
  do i=1,N_ZMAT
    write(LuWr,97) i,Symbols(i),NAT(i),(iZmat(j,i),j=1,3)
  end do
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'OutZMat - XYZCoords (Angstrom):'
  do i=1,nAtoms
    write(LuWr,98) i,(Angstrom*XYZCoords(j,i),j=1,3)
  end do
  write(LuWr,*)
end if

if (N_ZMAT < 3) return
if (N_ZMAT > nAtoms+3) then
  write(LuWr,'(A)') ' ZMAT cannot be defined :'
  write(LuWr,'(A)') ' Too many ghost Z atoms. '
  return
end if

! Z are ghost atoms
! A are both dummy (X) and real atoms

! Z  Z  Z
if ((NAT(1) == -1) .and. (NAT(2) == -1) .and. (NAT(3) == -1)) then
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = One/Angstrom
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-3)
  end if
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(1,3) = One/Angstrom
  else
    ZMATCoords(1,3) = XYZCoords(1,iBonded-3)
  end if
  if (iZmat(1,3) == 2) ZMATCoords(3,3) = ZMATCoords(3,2)
  ZMATCoords(:,4:N_ZMAT) = XYZCoords(:,1:N_ZMAT-3)
end if

! Z  Z  A
if ((NAT(1) == -1) .and. (NAT(2) == -1) .and. (NAT(3) >= 0)) then
  iBonded = 0
  do i=3,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = One/Angstrom
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-2)
  end if
  ZMATCoords(:,3:N_ZMAT) = XYZCoords(:,1:N_ZMAT-2)
end if

! Z  A  Z
if ((NAT(1) == -1) .and. (NAT(3) == -1) .and. (NAT(2) >= 0)) then
  dMaxTrasl = Zero
  dXYZ(:) = XYZCoords(:,1)
  do iAtom=1,N_ZMAT-2
    XYZCoords(1:2,iAtom) = XYZCoords(1:2,iAtom)-dXYZ(1:2)
  end do
  ZMATCoords(:,2) = XYZCoords(:,1)
  dMaxTrasl = maxval(abs(dXYZ(1:2)))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(2)
  if (iZmat(1,3) == 2) ZMATCoords(3,3) = ZMATCoords(3,2)
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(1,3) = One/Angstrom
  else
    ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
  end if
  ZMATCoords(:,4:N_ZMAT) = XYZCoords(:,2:N_ZMAT-2)
end if

! A  Z  Z
if ((NAT(1) >= 0) .and. (NAT(2) == -1) .and. (NAT(3) == -1)) then
  dXYZ(:) = XYZCoords(:,1)
  do iAtom=1,N_ZMAT-2
    XYZCoords(:,iAtom) = XYZCoords(:,iAtom)-dXYZ(:)
  end do
  dMaxTrasl = maxval(abs(dXYZ(:)))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(1)
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = One/Angstrom
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-2)
  end if
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(1,3) = One/Angstrom
  else
    ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
  end if
  ZMATCoords(:,4:N_ZMAT) = XYZCoords(:,2:N_ZMAT-2)
end if

! Z  A  B
if ((NAT(1) == -1) .and. (NAT(2) >= 0) .and. (NAT(3) >= 0)) then
  dXYZ(:) = XYZCoords(:,1)
  dXYZ(3) = Zero
  do iAtom=2,N_ZMAT
    ZMATCoords(:,iAtom) = XYZCoords(:,iAtom-1)-dXYZ(:)
  end do
  dMaxTrasl = maxval(abs(dXYZ(:)))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(2)
end if

! A  Z  B
if ((NAT(1) >= 0) .and. (NAT(2) == -1) .and. (NAT(3) >= 0)) then
  dXYZ = XYZCoords(:,1)
  do iAtom=1,N_ZMAT-1
    XYZCoords(:,iAtom) = XYZCoords(:,iAtom)-dXYZ(:)
  end do
  dMaxTrasl = maxval(abs(dXYZ(:)))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(1)
  iBonded = 0
  do i=3,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = One/Angstrom
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-1)
  end if
  ZMATCoords(:,3:N_ZMAT) = XYZCoords(:,2:N_ZMAT-1)
end if

! A  B  Z
if ((NAT(1) >= 0) .and. (NAT(2) >= 0) .and. (NAT(3) == -1)) then
  dMaxTrasl = Zero
  dXYZ = XYZCoords(:,1)
  do iAtom=1,N_ZMAT-1
    XYZCoords(:,iAtom) = XYZCoords(:,iAtom)-dXYZ(:)
  end do
  dMaxTrasl = maxval(abs(dXYZ(:)))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(1)
  ZMATCoords(:,2) = XYZCoords(:,2)
  if (iZmat(1,3) == 1) then
    iBonded = 0
    do i=4,N_ZMAT
      if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
    end do
    if (iBonded == 0) then
      ZMATCoords(1,3) = One/Angstrom
    else
      ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
    end if
  else
    ZMATCoords(:,3) = XYZCoords(:,2)
    ZMATCoords(1,3) = XYZCoords(1,2)+One/Angstrom
  end if
  ZMATCoords(:,4:N_ZMAT) = XYZCoords(:,3:N_ZMAT-1)
end if

! A  B  C
if ((NAT(1) >= 0) .and. (NAT(2) >= 0) .and. (NAT(3) >= 0)) ZMATCoords(:,:) = XYZCoords(:,1:N_ZMAT)

if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'OutZMat - ZMATCoords : '
  do i=1,N_ZMAT
    write(LuWr,99) i,NAT(i),(Angstrom*ZMATCoords(j,i),j=1,3)
  end do
end if

iAtom = 2
iRX = iZMAT(1,iAtom)
dXYZ(:) = ZMATCoords(:,iAtom)-ZMATCoords(:,iRX)
bond = sqrt(sum(dXYZ(:)**2))
ZMAT(iAtom,1) = bond*Angstrom

if (N_ZMAT >= 3) then
  iAtom = 3
  iRX = iZMAT(1,iAtom)
  iAX = iZMAT(2,iAtom)
  dXYZ(:) = ZMATCoords(:,iAtom)-ZMATCoords(:,iRX)
  dXYZ2(:) = ZMATCoords(:,iAX)-ZMATCoords(:,iRX)
  bond = sqrt(sum(dXYZ(:)**2))
  ZMAT(iAtom,1) = bond*Angstrom
  prod = sum(dXYZ(:)*dXYZ2(:))
  dvec2 = sqrt(sum(dXYZ2(:)**2))
  arccos = prod/(bond*dvec2)
  if (arccos > One) arccos = sign(One,arccos)
  alpha = acos(arccos)
  ZMAT(iAtom,2) = alpha/deg2rad
end if

if (N_ZMAT >= 4) then
  do iAtom=4,N_ZMAT
    if (NAT(iAtom) == -1) then
      write(LuWr,'(A)') ' ZMAT cannot be defined :'
      write(LuWr,'(A)') ' Found X atom.           '
      return
    end if
    iRX = iZMAT(1,iAtom)
    iAX = iZMAT(2,iAtom)
    iTX = iZMAT(3,iAtom)
    dXYZ(:) = ZMATCoords(:,iAtom)-ZMATCoords(:,iRX)
    dXYZ2(:) = ZMATCoords(:,iAX)-ZMATCoords(:,iRX)
    bond = sqrt(sum(dXYZ(:)**2))
    ZMAT(iAtom,1) = bond*Angstrom
    prod = sum(dXYZ(:)*dXYZ2(:))
    dvec2 = sqrt(sum(dXYZ2(:)**2))
    arccos = prod/(bond*dvec2)
    if (arccos > One) arccos = sign(One,arccos)
    alpha = acos(arccos)
    ZMAT(iAtom,2) = alpha/deg2rad
    CTX(:,1) = ZMATCoords(:,iTX)
    CTX(:,2) = ZMATCoords(:,iAX)
    CTX(:,3) = ZMATCoords(:,iRX)
    CTX(:,4) = ZMATCoords(:,iAtom)
    call Trsn(CTX,4,theta,Bt,.false.,.false.,Label,dBt,.false.)
    ZMAT(iAtom,3) = -theta/deg2rad
  end do
end if

write(LuWr,*)
write(LuWr,*) '****************************************************************'
write(LuWr,*) '* Nuclear coordinates in Z-Matrix format / Angstrom and Degree *'
write(LuWr,*) '----------------------------------------------------------------'
do i=1,N_ZMAT
  if (i == 1) write(LuWr,201) Symbols(i)
  if (i == 2) write(LuWr,202) Symbols(i),iZmat(1,i),ZMAT(i,1)
  if (i == 3) write(LuWr,203) Symbols(i),(iZmat(j,i),ZMAT(i,j),j=1,2)
  if (i >= 4) write(LuWr,204) Symbols(i),(iZmat(j,i),ZMAT(i,j),j=1,3)
end do
call mma_deallocate(iZmat)
write(LuWr,*)

call mma_deallocate(ZMAT)
call mma_deallocate(ZMATCoords)
call mma_deallocate(NAT)
call mma_deallocate(Symbols)

return

97 format(I3,1X,A,1X,I3,3X,3I4)
98 format(I3,1X,3(F12.6))
99 format(I3,1X,I3,1X,3(F12.6))
201 format(5X,A5)
202 format(5X,A5,1X,I4,F13.6)
203 format(5X,A5,1X,2(I4,F13.6))
204 format(5X,A5,1X,3(I4,F13.6))

end subroutine OutZMAT
