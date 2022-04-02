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
use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Cordst(3,3), Ax, Ay, Az
integer(kind=iwp), intent(in) :: i
real(kind=wp), intent(out) :: Rot(3,3), xt, yt, zt
integer(kind=iwp) :: IFLAG
real(kind=wp) :: A, ANORM, DELR, DEL(3), TAL, XA(3), XH(3), XO(3)
#include "warnings.h"

XO(1) = CORDST(1,1)-AX
XO(2) = CORDST(2,1)-AY
XO(3) = CORDST(3,1)-AZ
XH(1) = CORDST(1,2)-AX
XH(2) = CORDST(2,2)-AY
XH(3) = CORDST(3,2)-AZ
XA(1) = CORDST(1,3)-AX
XA(2) = CORDST(2,3)-AY
XA(3) = CORDST(3,3)-AZ
DEL = (XH+XA)*Half-XO
DELR = DEL(1)**2+DEL(2)**2+DEL(3)**2
DELR = DELR-1.225449_wp
!This is a check of the water geometry.
!If we enter here, something is wrong.
if (abs(DELR) > 1.0e-4_wp) then
  write(u6,*) 'Molecule',((i-1)/nCent)+1
  write(u6,*) ' WARNING IN TRANSROT ','delr',delr
  write(u6,*) ' O',XO
  write(u6,*) ' H',XH
  write(u6,*) ' A',XA
  call Quit(_RC_GENERAL_ERROR_)
end if
XT = XO(1)+0.3_wp/1.107_wp*DEL(1)
YT = XO(2)+0.3_wp/1.107_wp*DEL(2)
ZT = XO(3)+0.3_wp/1.107_wp*DEL(3)
ROT(1,3) = (XO(1)-XT)/0.3_wp
ROT(2,3) = (XO(2)-YT)/0.3_wp
ROT(3,3) = (XO(3)-ZT)/0.3_wp
ROT(:,2) = (XH-XA)/2.86_wp
ANORM = One/sqrt(ROT(1,3)**2+ROT(2,3)**2+ROT(3,3)**2)
ROT(:,3) = ROT(:,3)*ANORM
ANORM = One/sqrt(ROT(1,2)**2+ROT(2,2)**2+ROT(3,2)**2)
ROT(:,2) = ROT(:,2)*ANORM
ROT(1,1) = One-ROT(1,3)**2-ROT(1,2)**2
if (ROT(1,1) < Zero) ROT(1,1) = Zero
ROT(1,1) = sqrt(ROT(1,1))
ROT(2,1) = One-ROT(2,2)**2-ROT(2,3)**2
if (ROT(2,1) < Zero) ROT(2,1) = Zero
ROT(2,1) = sqrt(ROT(2,1))
ROT(3,1) = One-ROT(3,3)**2-ROT(3,2)**2
if (ROT(3,1) < Zero) ROT(3,1) = Zero
ROT(3,1) = sqrt(ROT(3,1))

IFLAG = 0
do
  TAL = ROT(1,1)*ROT(1,2)+ROT(2,1)*ROT(2,2)+ROT(3,1)*ROT(3,2)
  ROT(:,1) = ROT(:,1)-TAL*ROT(:,2)
  TAL = ROT(1,1)*ROT(1,3)+ROT(2,1)*ROT(2,3)+ROT(3,1)*ROT(3,3)
  ROT(:,1) = ROT(:,1)-TAL*ROT(:,3)
  A = One/sqrt(ROT(1,1)**2+ROT(2,1)**2+ROT(3,1)**2)
  ROT(:,1) = ROT(:,1)*A
  IFLAG = IFLAG+1
  if (IFLAG > 3) then
    write(u6,*) ' STOP IN TRANSROT'
    call Quit(_RC_GENERAL_ERROR_)
  end if
  if (A <= Ten) exit
end do

return

end subroutine TransRot
