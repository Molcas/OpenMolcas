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

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "warnings.h"
dimension Rot(3,3), Cordst(MxPut*MxCen,3)

XO = CORDST(I,1)-AX
YO = CORDST(I,2)-AY
ZO = CORDST(I,3)-AZ
XH = CORDST(I+1,1)-AX
YH = CORDST(I+1,2)-AY
ZH = CORDST(I+1,3)-AZ
XA = CORDST(I+2,1)-AX
YA = CORDST(I+2,2)-AY
ZA = CORDST(I+2,3)-AZ
DELX = (XH+XA)/2.-XO
DELY = (YH+YA)/2.-YO
DELZ = (ZH+ZA)/2.-ZO
DELR = DELX*DELX+DELY*DELY+DELZ*DELZ
DELR = DELR-1.225449d0
!This is a check of the water geometry.
!If we enter here, something is wrong.
if (abs(DELR) > .0001) then
  write(6,*) 'Molecule',((i-1)/5)+1
  write(6,*) ' WARNING IN TRANSROT ','delr',delr
  write(6,*) ' O',XO,YO,ZO
  write(6,*) ' H',XH,YH,ZH
  write(6,*) ' A',XA,YA,ZA
  call Quit(_RC_GENERAL_ERROR_)
end if
XT = XO+.3d0/1.107d0*DELX
YT = YO+.3d0/1.107d0*DELY
ZT = ZO+.3d0/1.107d0*DELZ
ROT(1,3) = (XO-XT)/.3d0
ROT(2,3) = (YO-YT)/.3d0
ROT(3,3) = (ZO-ZT)/.3d0
ROT(1,2) = (XH-XA)/2.86d0
ROT(2,2) = (YH-YA)/2.86d0
ROT(3,2) = (ZH-ZA)/2.86d0
ANORM = 1./sqrt(ROT(1,3)**2+ROT(2,3)**2+ROT(3,3)**2)
ROT(1,3) = ROT(1,3)*ANORM
ROT(2,3) = ROT(2,3)*ANORM
ROT(3,3) = ROT(3,3)*ANORM
ANORM = 1./sqrt(ROT(1,2)**2+ROT(2,2)**2+ROT(3,2)**2)
ROT(1,2) = ROT(1,2)*ANORM
ROT(2,2) = ROT(2,2)*ANORM
ROT(3,2) = ROT(3,2)*ANORM
ROT(1,1) = 1.-ROT(1,3)**2-ROT(1,2)**2
if (ROT(1,1) < 0.) ROT(1,1) = 0d0
ROT(1,1) = sqrt(ROT(1,1))
ROT(2,1) = 1.-ROT(2,2)**2-ROT(2,3)**2
if (ROT(2,1) < 0.) ROT(2,1) = 0d0
ROT(2,1) = sqrt(ROT(2,1))
ROT(3,1) = 1.-ROT(3,3)**2-ROT(3,2)**2
if (ROT(3,1) < 0.) ROT(3,1) = 0d0
ROT(3,1) = sqrt(ROT(3,1))

IFLAG = 0
732 continue
TAL = ROT(1,1)*ROT(1,2)+ROT(2,1)*ROT(2,2)+ROT(3,1)*ROT(3,2)
ROT(1,1) = ROT(1,1)-TAL*ROT(1,2)
ROT(2,1) = ROT(2,1)-TAL*ROT(2,2)
ROT(3,1) = ROT(3,1)-TAL*ROT(3,2)
TAL = ROT(1,1)*ROT(1,3)+ROT(2,1)*ROT(2,3)+ROT(3,1)*ROT(3,3)
ROT(1,1) = ROT(1,1)-TAL*ROT(1,3)
ROT(2,1) = ROT(2,1)-TAL*ROT(2,3)
ROT(3,1) = ROT(3,1)-TAL*ROT(3,3)
A = 1d0/sqrt(ROT(1,1)**2+ROT(2,1)**2+ROT(3,1)**2)
ROT(1,1) = ROT(1,1)*A
ROT(2,1) = ROT(2,1)*A
ROT(3,1) = ROT(3,1)*A
IFLAG = IFLAG+1
if (IFLAG > 3) then
  write(6,*) ' STOP IN TRANSROT'
  call Quit(_RC_GENERAL_ERROR_)
end if
if (A > 10.) GO TO 732

return

end subroutine TransRot
