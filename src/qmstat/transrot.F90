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

!----------------------------------------------------------------------*
! Compute rotation matrix for this specific solvent.                   *
!----------------------------------------------------------------------*
subroutine TransRot(Cordst,i,Rot,xt,yt,zt,Ax,Ay,Az)

use qmstat_global, only: nCent
use Constants, only: Zero, One, Ten
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: Cordst(3,3), Rot(3,3), xt, yt, zt, Ax, Ay, Az
integer(kind=iwp) :: i
integer(kind=iwp) :: IFLAG
real(kind=wp) :: A, ANORM, DELR, DELX, DELY, DELZ, TAL, XA, XH, XO, YA, YH, YO, ZA, ZH, ZO
#include "warnings.h"

XO = CORDST(1,1)-AX
YO = CORDST(2,1)-AY
ZO = CORDST(3,1)-AZ
XH = CORDST(1,2)-AX
YH = CORDST(2,2)-AY
ZH = CORDST(3,2)-AZ
XA = CORDST(1,3)-AX
YA = CORDST(2,3)-AY
ZA = CORDST(3,3)-AZ
DELX = (XH+XA)/2.-XO
DELY = (YH+YA)/2.-YO
DELZ = (ZH+ZA)/2.-ZO
DELR = DELX*DELX+DELY*DELY+DELZ*DELZ
DELR = DELR-1.225449_wp
!This is a check of the water geometry.
!If we enter here, something is wrong.
if (abs(DELR) > 1.0e-4_wp) then
  write(u6,*) 'Molecule',((i-1)/nCent)+1
  write(u6,*) ' WARNING IN TRANSROT ','delr',delr
  write(u6,*) ' O',XO,YO,ZO
  write(u6,*) ' H',XH,YH,ZH
  write(u6,*) ' A',XA,YA,ZA
  call Quit(_RC_GENERAL_ERROR_)
end if
XT = XO+0.3_wp/1.107_wp*DELX
YT = YO+0.3_wp/1.107_wp*DELY
ZT = ZO+0.3_wp/1.107_wp*DELZ
ROT(1,3) = (XO-XT)/0.3_wp
ROT(2,3) = (YO-YT)/0.3_wp
ROT(3,3) = (ZO-ZT)/0.3_wp
ROT(1,2) = (XH-XA)/2.86_wp
ROT(2,2) = (YH-YA)/2.86_wp
ROT(3,2) = (ZH-ZA)/2.86_wp
ANORM = 1./sqrt(ROT(1,3)**2+ROT(2,3)**2+ROT(3,3)**2)
ROT(1,3) = ROT(1,3)*ANORM
ROT(2,3) = ROT(2,3)*ANORM
ROT(3,3) = ROT(3,3)*ANORM
ANORM = 1./sqrt(ROT(1,2)**2+ROT(2,2)**2+ROT(3,2)**2)
ROT(1,2) = ROT(1,2)*ANORM
ROT(2,2) = ROT(2,2)*ANORM
ROT(3,2) = ROT(3,2)*ANORM
ROT(1,1) = 1.-ROT(1,3)**2-ROT(1,2)**2
if (ROT(1,1) < 0.) ROT(1,1) = Zero
ROT(1,1) = sqrt(ROT(1,1))
ROT(2,1) = 1.-ROT(2,2)**2-ROT(2,3)**2
if (ROT(2,1) < 0.) ROT(2,1) = Zero
ROT(2,1) = sqrt(ROT(2,1))
ROT(3,1) = 1.-ROT(3,3)**2-ROT(3,2)**2
if (ROT(3,1) < 0.) ROT(3,1) = Zero
ROT(3,1) = sqrt(ROT(3,1))

IFLAG = 0
do
  TAL = ROT(1,1)*ROT(1,2)+ROT(2,1)*ROT(2,2)+ROT(3,1)*ROT(3,2)
  ROT(1,1) = ROT(1,1)-TAL*ROT(1,2)
  ROT(2,1) = ROT(2,1)-TAL*ROT(2,2)
  ROT(3,1) = ROT(3,1)-TAL*ROT(3,2)
  TAL = ROT(1,1)*ROT(1,3)+ROT(2,1)*ROT(2,3)+ROT(3,1)*ROT(3,3)
  ROT(1,1) = ROT(1,1)-TAL*ROT(1,3)
  ROT(2,1) = ROT(2,1)-TAL*ROT(2,3)
  ROT(3,1) = ROT(3,1)-TAL*ROT(3,3)
  A = One/sqrt(ROT(1,1)**2+ROT(2,1)**2+ROT(3,1)**2)
  ROT(1,1) = ROT(1,1)*A
  ROT(2,1) = ROT(2,1)*A
  ROT(3,1) = ROT(3,1)*A
  IFLAG = IFLAG+1
  if (IFLAG > 3) then
    write(u6,*) ' STOP IN TRANSROT'
    call Quit(_RC_GENERAL_ERROR_)
  end if
  if (A <= Ten) exit
end do

return

end subroutine TransRot
