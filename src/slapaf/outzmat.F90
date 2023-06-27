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

implicit real*8(a-h,o-z)
real*8 XYZCoords(3,nAtoms)
#include "stdalloc.fh"
#include "SysDef.fh"
character(len=5), allocatable :: Symbols(:)
integer, allocatable :: NAT(:)
real*8, allocatable :: ZMATCoords(:,:), ZMAT(:,:)

call mma_allocate(Symbols,N_ZMAT,Label='Symbols')
call mma_allocate(NAT,N_ZMAT,Label='NAT')
call mma_allocate(ZMATCoords,3,N_ZMAT,Label='ZMATCoords')
call mma_allocate(ZMAT,N_ZMAT,3,Label='ZMAT')

call OutZMAT_Internal(nAtoms,XYZCoords,N_ZMAT,Symbols,NAT,ZMATCoords,ZMAT)

call mma_deallocate(ZMAT)
call mma_deallocate(ZMATCoords)
call mma_deallocate(NAT)
call mma_deallocate(Symbols)

return

end subroutine OutZMAT

subroutine OutZMAT_Internal(nAtoms,XYZCoords,N_ZMAT,Symbols,NAT,ZMATCoords,ZMAT)
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

implicit real*8(a-h,o-z)
implicit integer(i-n)
#include "stdalloc.fh"
#include "angstr.fh"
character(len=5) Symbols(N_ZMAT)
integer NAT(N_ZMAT)
integer, allocatable :: iZmat(:,:)
real*8 XYZCoords(3,nAtoms), ZMATCoords(3,N_ZMAT), ZMAT(N_ZMAT,3)
real*8 CTX(3,4), Bt(3,4), dBt(3,4,3,4)
logical IfTest
character(len=8) Label
real*8, parameter :: ThrsTrasl = 1.0d0  ! Threshold for warning

#ifdef _DEBUGPRINT_
IfTest = .true.
#else
IfTest = .false.
#endif

LuWr = 6
todeg = 45.0d0/atan(1.0d0)
dMaxTrasl = 0.0d0
Label = ' '
call fZero(ZMAT,N_ZMAT*3)
call fZero(ZMATCoords,N_ZMAT*3)
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
    write(LuWr,98) i,(angstr*XYZCoords(j,i),j=1,3)
  end do
  write(LuWr,*)
end if

if (N_ZMAT < 3) return
if (N_ZMAT > (nAtoms+3)) then
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
    ZMATCoords(3,2) = 1.0d0/angstr
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-3)
  end if
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(1,3) = 1.0d0/angstr
  else
    ZMATCoords(1,3) = XYZCoords(1,iBonded-3)
  end if
  if (iZmat(1,3) == 2) ZMATCoords(3,3) = ZMATCoords(3,2)
  do iAtom=4,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-3)
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-3)
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-3)
  end do
end if

! Z  Z  A
if ((NAT(1) == -1) .and. (NAT(2) == -1) .and. (NAT(3) >= 0)) then
  iBonded = 0
  do i=3,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = 1.0d0/angstr
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-2)
  end if
  do iAtom=3,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-2)
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-2)
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-2)
  end do
end if

! Z  A  Z
if ((NAT(1) == -1) .and. (NAT(3) == -1) .and. (NAT(2) >= 0)) then
  dMaxTrasl = 0.0d0
  dX = XYZCoords(1,1)
  dY = XYZCoords(2,1)
  do iAtom=1,N_ZMAT-2
    XYZCoords(1,iAtom) = XYZCoords(1,iAtom)-dX
    XYZCoords(2,iAtom) = XYZCoords(2,iAtom)-dY
  end do
  ZMATCoords(1,2) = XYZCoords(1,1)
  ZMATCoords(2,2) = XYZCoords(2,1)
  ZMATCoords(3,2) = XYZCoords(3,1)
  dMaxTrasl = max(0.0d0,abs(dX),abs(dY))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(2)
  if (iZmat(1,3) == 2) ZMATCoords(3,3) = ZMATCoords(3,2)
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(1,3) = 1.0d0/angstr
  else
    ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
  end if
  do iAtom=4,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-2)
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-2)
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-2)
  end do
end if

! A  Z  Z
if ((NAT(1) >= 0) .and. (NAT(2) == -1) .and. (NAT(3) == -1)) then
  dX = XYZCoords(1,1)
  dY = XYZCoords(2,1)
  dZ = XYZCoords(3,1)
  do iAtom=1,N_ZMAT-2
    XYZCoords(1,iAtom) = XYZCoords(1,iAtom)-dX
    XYZCoords(2,iAtom) = XYZCoords(2,iAtom)-dY
    XYZCoords(3,iAtom) = XYZCoords(3,iAtom)-dZ
  end do
  dMaxTrasl = max(0.0d0,abs(dX),abs(dY),abs(dZ))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(1)
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = 1.0d0/angstr
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-2)
  end if
  iBonded = 0
  do i=4,N_ZMAT
    if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(1,3) = 1.0d0/angstr
  else
    ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
  end if
  do iAtom=4,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-2)
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-2)
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-2)
  end do
end if

! Z  A  B
if ((NAT(1) == -1) .and. (NAT(2) >= 0) .and. (NAT(3) >= 0)) then
  dX = XYZCoords(1,1)
  dY = XYZCoords(2,1)
  do iAtom=2,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-1)-dX
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-1)-dY
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-1)
  end do
  dMaxTrasl = max(0.0d0,abs(dX),abs(dY))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(2)
end if

! A  Z  B
if ((NAT(1) >= 0) .and. (NAT(2) == -1) .and. (NAT(3) >= 0)) then
  dX = XYZCoords(1,1)
  dY = XYZCoords(2,1)
  dZ = XYZCoords(3,1)
  do iAtom=1,N_ZMAT-1
    XYZCoords(1,iAtom) = XYZCoords(1,iAtom)-dX
    XYZCoords(2,iAtom) = XYZCoords(2,iAtom)-dY
    XYZCoords(3,iAtom) = XYZCoords(3,iAtom)-dZ
  end do
  dMaxTrasl = max(0.0d0,abs(dX),abs(dY),abs(dZ))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(1)
  iBonded = 0
  do i=3,N_ZMAT
    if ((iZmat(1,i) == 2) .and. (iBonded == 0)) iBonded = i
  end do
  if (iBonded == 0) then
    ZMATCoords(3,2) = 1.0d0/angstr
  else
    ZMATCoords(3,2) = XYZCoords(3,iBonded-1)
  end if
  do iAtom=3,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-1)
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-1)
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-1)
  end do
end if

! A  B  Z
if ((NAT(1) >= 0) .and. (NAT(2) >= 0) .and. (NAT(3) == -1)) then
  dMaxTrasl = 0.0d0
  dX = XYZCoords(1,1)
  dY = XYZCoords(2,1)
  dZ = XYZCoords(3,1)
  do iAtom=1,N_ZMAT-1
    XYZCoords(1,iAtom) = XYZCoords(1,iAtom)-dX
    XYZCoords(2,iAtom) = XYZCoords(2,iAtom)-dY
    XYZCoords(3,iAtom) = XYZCoords(3,iAtom)-dZ
  end do
  dMaxTrasl = max(0.0d0,abs(dX),abs(dY),abs(dZ))
  if (dMaxTrasl >= ThrsTrasl) write(LuWr,'(A)') '  Warning: MaxTrasl >= ThrsTrasl for atom ',Symbols(1)
  ZMATCoords(1,2) = XYZCoords(1,2)
  ZMATCoords(2,2) = XYZCoords(2,2)
  ZMATCoords(3,2) = XYZCoords(3,2)
  if (iZmat(1,3) == 1) then
    iBonded = 0
    do i=4,N_ZMAT
      if ((iZmat(1,i) == 3) .and. (iBonded == 0)) iBonded = i
    end do
    if (iBonded == 0) then
      ZMATCoords(1,3) = 1.0d0/angstr
    else
      ZMATCoords(1,3) = XYZCoords(1,iBonded-2)
    end if
  else
    ZMATCoords(1,3) = XYZCoords(1,2)+1.0d0/angstr
    ZMATCoords(2,3) = XYZCoords(2,2)
    ZMATCoords(3,3) = XYZCoords(3,2)
  end if
  do iAtom=4,N_ZMAT
    ZMATCoords(1,iAtom) = XYZCoords(1,iAtom-1)
    ZMATCoords(2,iAtom) = XYZCoords(2,iAtom-1)
    ZMATCoords(3,iAtom) = XYZCoords(3,iAtom-1)
  end do
end if

! A  B  C
if ((NAT(1) >= 0) .and. (NAT(2) >= 0) .and. (NAT(3) >= 0)) then
  call dCopy_(N_ZMAT*3,XYZCoords,1,ZMATCoords,1)
end if

if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'OutZMat - ZMATCoords : '
  do i=1,N_ZMAT
    write(LuWr,99) i,NAT(i),(angstr*ZMATCoords(j,i),j=1,3)
  end do
end if

iAtom = 2
iRX = iZMAT(1,iAtom)
dX = ZMATCoords(1,iAtom)-ZMATCoords(1,iRX)
dY = ZMATCoords(2,iAtom)-ZMATCoords(2,iRX)
dZ = ZMATCoords(3,iAtom)-ZMATCoords(3,iRX)
bond = sqrt(dX*dX+dY*dY+dZ*dZ)
ZMAT(iAtom,1) = bond*angstr

if (N_ZMAT >= 3) then
  iAtom = 3
  iRX = iZMAT(1,iAtom)
  iAX = iZMAT(2,iAtom)
  dX = ZMATCoords(1,iAtom)-ZMATCoords(1,iRX)
  dY = ZMATCoords(2,iAtom)-ZMATCoords(2,iRX)
  dZ = ZMATCoords(3,iAtom)-ZMATCoords(3,iRX)
  dX2 = ZMATCoords(1,iAX)-ZMATCoords(1,iRX)
  dY2 = ZMATCoords(2,iAX)-ZMATCoords(2,iRX)
  dZ2 = ZMATCoords(3,iAX)-ZMATCoords(3,iRX)
  bond = sqrt(dX*dX+dY*dY+dZ*dZ)
  ZMAT(iAtom,1) = bond*angstr
  prod = dX*dX2+dY*dY2+dZ*dZ2
  dvec2 = sqrt(dX2*dX2+dY2*dY2+dZ2*dZ2)
  arccos = prod/(bond*dvec2)
  if (arccos > 1.0d0) arccos = sign(1.0d0,arccos)
  alpha = acos(arccos)
  ZMAT(iAtom,2) = alpha*todeg
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
    dX = ZMATCoords(1,iAtom)-ZMATCoords(1,iRX)
    dY = ZMATCoords(2,iAtom)-ZMATCoords(2,iRX)
    dZ = ZMATCoords(3,iAtom)-ZMATCoords(3,iRX)
    dX2 = ZMATCoords(1,iAX)-ZMATCoords(1,iRX)
    dY2 = ZMATCoords(2,iAX)-ZMATCoords(2,iRX)
    dZ2 = ZMATCoords(3,iAX)-ZMATCoords(3,iRX)
    bond = sqrt(dX*dX+dY*dY+dZ*dZ)
    ZMAT(iAtom,1) = bond*angstr
    prod = dX*dX2+dY*dY2+dZ*dZ2
    dvec2 = sqrt(dX2*dX2+dY2*dY2+dZ2*dZ2)
    arccos = prod/(bond*dvec2)
    if (arccos > 1.0d0) arccos = sign(1.0d0,arccos)
    alpha = acos(arccos)
    ZMAT(iAtom,2) = alpha*todeg
    do ii=1,3
      CTX(ii,1) = ZMATCoords(ii,iTX)
      CTX(ii,2) = ZMATCoords(ii,iAX)
      CTX(ii,3) = ZMATCoords(ii,iRX)
      CTX(ii,4) = ZMATCoords(ii,iAtom)
    end do
    call Trsn(CTX,4,theta,Bt,.false.,.false.,Label,dBt,.false.)
    ZMAT(iAtom,3) = -theta*todeg
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

201 format(5X,A5)
202 format(5X,A5,1X,I4,F13.6)
203 format(5X,A5,1X,2(I4,F13.6))
204 format(5X,A5,1X,3(I4,F13.6))

97 format(I3,1X,A,1X,I3,3X,3I4)
98 format(I3,1X,3(F12.6))
99 format(I3,1X,I3,1X,3(F12.6))

end subroutine OutZMAT_Internal
