!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine CoSys(Cent,R,xyz)

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 Cent(3,3), r(3), xyz(3,2)
integer iComp(3)
logical Linear

iRout = 221
iPrint = nPrint(iRout)
Lu = 6
if (iPrint >= 99) then
  call RecPrt('CoSys: Cent',' ',Cent,3,3)
end if

! Check if linear

Co = Zero
Crap = Zero
RR1 = Zero
RR2 = Zero
ThrAcos = 1.0D-6
do i=1,3
  Co = Co+(Cent(i,1)-Cent(i,2))*(Cent(i,3)-Cent(i,2))
  RR1 = RR1+(Cent(i,1)-Cent(i,2))**2
  RR2 = RR2+(Cent(i,3)-Cent(i,2))**2
end do
RR1 = sqrt(RR1)
RR2 = sqrt(RR2)
Co = Co/(RR1*RR2)
do i=1,3
  Crap = Crap+((Cent(i,3)-Cent(i,2))/RR2-sign(One,Co)*(Cent(i,1)-Cent(i,2))/RR1)**2
end do
Crap = sqrt(Crap)
if (iPrint >= 99) then
  write(6,*) 'Co=',Co
  write(6,*) 'Crap=',Crap
end if
if (Crap < 1.0D-6) then
  Si = Crap
  if (Co < Zero) then
    Fi = Pi-ArSin(Si)
  else
    Fi = ArSin(Si)
  end if
else
  if ((Co > One) .and. (Co < One+ThrAcos)) Co = One
  if ((Co < -One) .and. (Co > -One-ThrAcos)) Co = -One
  if ((Co > One+ThrAcos) .or. (Co < -One-ThrAcos)) then
    call WarningMessage(2,'Error in CoSys')
    write(6,*) 'Error in cosys: arcos(',Co,')'
    call Abend()
  end if
  Fi = ArCos(Co)
  Si = sqrt(One-Co**2)
end if
if (iPrint >= 99) then
  write(6,*) 'Fi,Pi=',Fi,Pi
  write(6,*) 'Pi-Fi=',Pi-Fi
end if

Linear = abs(Si) < 1.0D-13

! Form reference axis

RR = Zero
do i=1,3
  R(i) = Cent(i,3)-Cent(i,1)
  RR = RR+R(i)**2
end do
RR = sqrt(RR)
if ((RR1 >= RR2) .and. (RR1 >= RR)) then
  do i=1,3
    R(i) = Cent(i,1)-Cent(i,2)
  end do
  RR = RR1
else if (RR2 >= RR) then
  do i=1,3
    R(i) = Cent(i,3)-Cent(i,2)
  end do
  RR = RR2
end if
call DScal_(3,One/RR,R,1)

if (iPrint >= 99) then
  write(6,*) 'Linear=',Linear
  write(6,*) 'RR=',RR
  call RecPrt('R',' ',R,3,1)
end if

99 continue
if (Linear) then

  nComp = 0
  do i=1,3
    if (R(i) == Zero) then
      nComp = nComp+1
      iComp(nComp) = i
    end if
  end do

  ! Compute the WDC B-matrix

  call dcopy_(6,[Zero],0,xyz,1)
  if (nComp == 0) then
    !write(6,*) ' Case nComp == 0'

    xyz(1,1) = R(1)
    xyz(2,1) = R(2)
    xyz(3,1) = -R(3)
    r12 = R(1)**2+R(2)**2-R(3)**2
    r11 = R(1)**2+R(2)**2+R(3)**2
    xyz(1,1) = xyz(1,1)-(r12/r11)*R(1)
    xyz(2,1) = xyz(2,1)-(r12/r11)*R(2)
    xyz(3,1) = xyz(3,1)-(r12/r11)*R(3)
    r2 = sqrt(xyz(1,1)**2+xyz(2,1)**2+xyz(3,1)**2)
    xyz(1,1) = xyz(1,1)/r2
    xyz(2,1) = xyz(2,1)/r2
    xyz(3,1) = xyz(3,1)/r2
    xyz(1,2) = R(2)*xyz(3,1)-R(3)*xyz(2,1)
    xyz(2,2) = R(3)*xyz(1,1)-R(1)*xyz(3,1)
    xyz(3,2) = R(1)*xyz(2,1)-R(2)*xyz(1,1)
    r2 = sqrt(xyz(1,2)**2+xyz(2,2)**2+xyz(3,2)**2)
    xyz(1,2) = xyz(1,2)/r2
    xyz(2,2) = xyz(2,2)/r2
    xyz(3,2) = xyz(3,2)/r2

  else if (nComp == 1) then
    !write(6,*) ' Case nComp == 1'

    i = iComp(1)
    xyz(i,1) = One
    xyz(1,2) = R(2)*xyz(3,1)-R(3)*xyz(2,1)
    xyz(2,2) = R(3)*xyz(1,1)-R(1)*xyz(3,1)
    xyz(3,2) = R(1)*xyz(2,1)-R(2)*xyz(1,1)
    r2 = sqrt(xyz(1,2)**2+xyz(2,2)**2+xyz(3,2)**2)
    xyz(1,2) = xyz(1,2)/r2
    xyz(2,2) = xyz(2,2)/r2
    xyz(3,2) = xyz(3,2)/r2

  else if (nComp == 2) then
    !write(6,*) ' Case nComp == 2'

    i = iComp(1)
    xyz(i,1) = One
    i = iComp(2)
    xyz(i,2) = One

  else
    call WarningMessage(2,'Error in CoSys')
    write(Lu,*) ' CoSys: nComp == 3'
    call Abend()
  end if

else     ! Non-linear

  ! Form the cross product R12xR32
  RR = Zero
  do i=1,3
    j = mod(i,3)+1
    k = mod(i+1,3)+1
    R21j = Cent(j,1)-Cent(j,2)
    R21k = Cent(k,1)-Cent(k,2)
    R23j = Cent(j,3)-Cent(j,2)
    R23k = Cent(k,3)-Cent(k,2)
    xyz(i,2) = R21j*R23k-R21k*R23j
    RR = RR+xyz(i,2)**2
  end do
  if (RR == Zero) then
    Linear = .true.
    if (iPrint >= 99) write(6,*) 'Linear=',Linear
    Go To 99
  end if
  call DScal_(3,One/sqrt(RR),xyz(1,2),1)
  if (iPrint >= 99) write(6,*) 'RR=',RR

  RR = Zero
  do i=1,3
    j = mod(i,3)+1
    k = mod(i+1,3)+1
    xyz(i,1) = xyz(j,2)*R(k)-xyz(k,2)*R(j)
    RR = RR+xyz(i,1)**2
  end do
  call DScal_(3,One/sqrt(RR),xyz(1,1),1)
  if (iPrint >= 99) write(6,*) 'RR=',RR
  if (iPrint >= 99) then
    call RecPrt('xyz',' ',xyz,3,2)
  end if

end if
if (iPrint >= 99) then
  call RecPrt(' Reference Axis',' ',R,3,1)
  call RecPrt(' Perpendicular Axes',' ',xyz,3,2)
end if

return

end subroutine CoSys
