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
subroutine OutofP0(xyz,nCent,Teta,ldB)

use Constants, only: deg2rad
use Constants, only: Zero, One, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: nCent
real(kind=wp) :: xyz(3,nCent), Teta
logical(kind=iwp) :: ldB
integer(kind=iwp) :: mCent
real(kind=wp) :: BR14X(3,3), C14X(3,3), CosFi1, CosFi2, CosFi3, dBR14X(3,3,3,3), dFi1, dFi2, dFi3, e41x, e41y, e41z, e42x, e42y, &
                 e42z, e43x, e43y, e43z, Fi1, Fi2, Fi3, Q41, Q42, Q43, R41KV, R42(3), R42KV, R43(3), R43KV, RX1, RX2, RX3, RY1, &
                 RY2, RY3, RZ1, RZ2, RZ3
character(len=8) :: Label
real(kind=wp), external :: ArCos

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
dFi1 = Fi1/deg2rad
if ((dFi1 > 177.5_wp) .or. (dFi1 < 2.5_wp)) write(u6,*) 'Warning: auxiliary Angle close to end of range'

! Dirty exit! This happens when an earlier structure is ill defined.

if (abs(Fi1-Pi) < 1.0e-13_wp) then
  Teta = Zero
  return
end if

! Get the angle between e41 and e43

CosFi2 = e41x*e43x+e41y*e43y+e41z*e43z

Fi2 = ArCos(CosFi2)
if (abs(CosFi2) > One) call RecPrt('xyz(2)',' ',xyz,3,4)
dFi2 = Fi2/deg2rad
if ((dFi2 > 177.5_wp) .or. (dFi2 < 2.5_wp)) write(u6,*) 'Warning: auxiliary Angle close to end of range'

! Get the angle between e41 and e42

CosFi3 = e41x*e42x+e41y*e42y+e41z*e42z

Fi3 = ArCos(CosFi3)
if (abs(CosFi3) > One) call RecPrt('xyz(3)',' ',xyz,3,4)
dFi3 = Fi3/deg2rad
if ((dFi3 > 177.5_wp) .or. (dFi3 < 2.5_wp)) write(u6,*) 'Warning: auxiliary Angle close to end of range'
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

C14X(:,1) = xyz(:,1)
C14X(:,2) = xyz(:,4)

! The 3rd is

R42(:) = xyz(:,2)-xyz(:,4)
R43(:) = xyz(:,3)-xyz(:,4)
C14X(1,3) = R42(2)*R43(3)-R42(3)*R43(2)
C14X(2,3) = R42(3)*R43(1)-R42(1)*R43(3)
C14X(3,3) = R42(1)*R43(2)-R42(2)*R43(1)

! Exit if 2-3-4 are collinear
! (equivalent to the above check, but this is more concrete)

if ((C14X(1,3)**2+C14X(2,3)**2+C14X(3,3)**2) < 1.0e-10_wp) then
  Teta = Zero
  return
end if
C14X(:,3) = C14X(:,3)+xyz(:,4)

mCent = 3
call Bend(C14X,mCent,Teta,BR14X,.false.,.false.,Label,dBR14X,ldB)

Teta = Teta-Half*Pi
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine OutofP0
#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(OutofP0)

#endif
