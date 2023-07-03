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

#include "compiler_features.h"
#ifdef _Test_Numerical_
subroutine OutofP1(xyz,nCent,Teta,Bt,ldB)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 Bt(3,nCent), xyz(3,nCent), R42(3), R43(3), C14X(3,3), BR14X(3,3), dBR14X(3,3,3,3)
logical ldB
character*8 Label

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! First some diagnostics

! 4-->1 (Bond)
RX1 = xyz(1,1)-xyz(1,4)
RY1 = xyz(2,1)-xyz(2,4)
RZ1 = xyz(3,1)-xyz(3,4)
R41KV = RX1*RX1+RY1*RY1+RZ1*RZ1
Q41 = sqrt(R41KV)
e41x = RX1/Q41
e41y = RY1/Q41
e41z = RZ1/Q41
! 4-->2 (Bond in plane)
RX2 = xyz(1,2)-xyz(1,4)
RY2 = xyz(2,2)-xyz(2,4)
RZ2 = xyz(3,2)-xyz(3,4)
R42KV = RX2*RX2+RY2*RY2+RZ2*RZ2
Q42 = sqrt(R42KV)
e42x = RX2/Q42
e42y = RY2/Q42
e42z = RZ2/Q42
! 4-->3 (Bond in plane)
RX3 = xyz(1,3)-xyz(1,4)
RY3 = xyz(2,3)-xyz(2,4)
RZ3 = xyz(3,3)-xyz(3,4)
R43KV = RX3*RX3+RY3*RY3+RZ3*RZ3
Q43 = sqrt(R43KV)
e43x = RX3/Q43
e43y = RY3/Q43
e43z = RZ3/Q43

! Get the angle between e43 and e42

CosFi1 = e43x*e42x+e43y*e42y+e43z*e42z

Fi1 = ArCos(CosFi1)
if (abs(CosFi1) > One) call RecPrt('xyz(1)',' ',xyz,3,4)
dFi1 = 180.d0*Fi1/Pi
if ((dFi1 > 177.5d0) .or. (dFi1 < 2.5d0)) then
  write(6,*) 'Warning: auxiliary Angle close to end of range'
end if

! Dirty exit! This happens when an earlier structure is ill defined.

if (abs(Fi1-Pi) < 1.0D-13) then
  Teta = 0.0d0
  call FZero(Bt,3*nCent)
  return
end if

! Get the angle between e41 and e43

CosFi2 = e41x*e43x+e41y*e43y+e41z*e43z

Fi2 = ArCos(CosFi2)
if (abs(CosFi2) > One) call RecPrt('xyz(2)',' ',xyz,3,4)
dFi2 = 180.d0*Fi2/Pi
if ((dFi2 > 177.5d0) .or. (dFi2 < 2.5d0)) then
  write(6,*) 'Warning: auxiliary Angle close to end of range'
end if

! Get the angle between e41 and e42

CosFi3 = e41x*e42x+e41y*e42y+e41z*e42z

Fi3 = ArCos(CosFi3)
if (abs(CosFi3) > One) call RecPrt('xyz(3)',' ',xyz,3,4)
dFi3 = 180.d0*Fi3/Pi
if ((dFi3 > 177.5d0) .or. (dFi3 < 2.5d0)) then
  write(6,*) 'Warning: auxiliary Angle close to end of range'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('xyz',' ',xyz,3,nCent)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! The first two centers are trivially

call dcopy_(3,xyz(1,1),1,C14X(1,1),1)
call dcopy_(3,xyz(1,4),1,C14X(1,2),1)

! The 3rd is

R42(1) = xyz(1,2)-xyz(1,4)
R42(2) = xyz(2,2)-xyz(2,4)
R42(3) = xyz(3,2)-xyz(3,4)
R43(1) = xyz(1,3)-xyz(1,4)
R43(2) = xyz(2,3)-xyz(2,4)
R43(3) = xyz(3,3)-xyz(3,4)
C14X(1,3) = R42(2)*R43(3)-R42(3)*R43(2)
C14X(2,3) = R42(3)*R43(1)-R42(1)*R43(3)
C14X(3,3) = R42(1)*R43(2)-R42(2)*R43(1)

! Exit if 2-3-4 are collinear
! (equivalent to the above check, but this is more concrete)

if ((C14X(1,3)**2+C14X(2,3)**2+C14X(3,3)**2) < 1.0D-10) then
  Teta = 0.0d0
  call FZero(Bt,3*nCent)
  return
end if
C14X(1,3) = C14X(1,3)+xyz(1,4)
C14X(2,3) = C14X(2,3)+xyz(2,4)
C14X(3,3) = C14X(3,3)+xyz(3,4)

mCent = 3
call Bend(C14X,mCent,Teta,BR14X,.false.,.false.,Label,dBR14X,ldB)

Teta = Teta-Pi/Two
dTeta = 180.d0*Teta/Pi
if ((dTeta > 87.5d0) .or. (dTeta < -87.5d0)) then
  write(6,*) 'Warning: Out of plane angle close to end of range'
end if

! Compute the WDC matrix

do ix=1,3
  iy = mod(ix+1,4)+(ix+1)/4
  iz = mod(iy+1,4)+(iy+1)/4

  Bt(ix,1) = -BR14X(ix,1)
  Bt(ix,2) = R43(iz)*BR14X(iy,3)-R43(iy)*BR14X(iz,3)
  Bt(ix,3) = -R42(iz)*BR14X(iy,3)+R42(iy)*BR14X(iz,3)

  Bt(ix,4) = -(Bt(ix,1)+Bt(ix,2)+Bt(ix,3))

end do
#ifdef _DEBUGPRINT_
call RecPrt('Outofp: R43',' ',R43,1,3)
call RecPrt('Outofp: R43',' ',R42,1,3)
call RecPrt('Outofp: BR14X',' ',BR14X,3,3)
call RecPrt('Outofp: B matrix',' ',Bt,3,nCent)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine OutofP1
#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(OutofP1)

#endif
