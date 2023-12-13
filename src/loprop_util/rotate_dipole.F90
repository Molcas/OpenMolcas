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

subroutine Rotate_Dipole(rMP,EC,nij,nElem,nPert,ij,ii,jj,Dipole_Rot_A,Dipole_Rot_B,Dipole_Rot_AB,R_A,R_B)

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nij, nElem, nPert, ij, ii, jj
real(kind=wp), intent(in) :: rMP(nij,0:nElem-1,0:nPert-1), EC(3,nij)
real(kind=wp), intent(out) :: Dipole_Rot_A, Dipole_Rot_B, Dipole_Rot_AB, R_A, R_B
integer(kind=iwp) :: nL
real(kind=wp) :: E_Bond(3), QQ, R_Q(3,3), R_temp, R_test(3), rMu_A(3), rMu_AB(3), rMu_B(3), T(3,3), T_R, Tmp(3,3), x_R, y_R, z_R

! The three dipole components for the bond
!write(u6,*) 'Test 0'
!call xFlush(u6)
!?? x_my = rMP(ij,1,0)
!?? y_my = rMP(ij,2,0)
!?? z_my = rMP(ij,3,0)
! Vector from ii to jj.
x_R = EC(1,ii)-EC(1,jj)
y_R = EC(2,ii)-EC(2,jj)
z_R = EC(3,ii)-EC(3,jj)
T_R = sqrt(x_R**2+y_R**2+z_R**2)
!?? write(u6,*) 'My=',x_my,y_my,z_my
!?? write(u6,*) 'R=',x_R,y_R,z_R
! Unit vector in the direction from ii to jj
E_Bond(1) = x_r/T_R
E_Bond(2) = y_r/T_R
E_Bond(3) = z_r/T_R
! Get transformation matrix (T) that rotates the coordinate system to align with the bond
nL = 1
call GS(E_bond,nL,T,3,.true.,.false.)
call RecPrt('T-matrix',' ',T,3,3)

call RecPrt('EC(*,ij) origional',' ',EC(1,ij),1,3)
call RecPrt('EC(*,ii) origional',' ',EC(1,ii),1,3)
call RecPrt('EC(*,jj) origional',' ',EC(1,jj),1,3)
! Rotate EC(1,ij) to the new coordinate system
call DGEMM_('T','N',3,1,3,One,T,3,EC(1,ij),3,Zero,R_test,3)
call RecPrt('EC(Bond system)',' ',R_Test,1,3)

! Rotate the dipole moment for the bond and ii and jj to the new coordinate system
tmp(1,1) = rMP(ij,1,0)
tmp(2,1) = rMP(ij,2,0)
tmp(3,1) = rMP(ij,3,0)
call DGEMM_('T','N',3,1,3,One,T,3,tmp,nij,Zero,rMu_AB,3)
call RecPrt('rMu_AB',' ',rMu_AB,1,3)
tmp(1,1) = rMP(ii,1,0)
tmp(2,1) = rMP(ii,2,0)
tmp(3,1) = rMP(ii,3,0)
call DGEMM_('T','N',3,1,3,One,T,3,tmp,nij,Zero,rMu_A,3)
call RecPrt('rMu_A',' ',rMu_A,1,3)
tmp(1,1) = rMP(jj,1,0)
tmp(2,1) = rMP(jj,2,0)
tmp(3,1) = rMP(jj,3,0)
call DGEMM_('T','N',3,1,3,One,T,3,tmp,nij,Zero,rMu_B,3)
call RecPrt('rMu_B',' ',rMu_B,1,3)

R_Q(1,1) = rMP(ij,4,0)
R_Q(1,2) = Half*rMP(ij,5,0)
R_Q(1,3) = Half*rMP(ij,6,0)
R_Q(2,2) = rMP(ij,7,0)
R_Q(2,3) = Half*rMP(ij,8,0)
R_Q(3,3) = rMP(ij,9,0)
R_Q(2,1) = R_Q(1,2)
R_Q(3,1) = R_Q(1,3)
R_Q(3,2) = R_Q(2,3)
call RecPrt('R_Q',' ',R_Q,3,3)
call DGEMM_('N','T',3,3,3,One,R_Q,3,T,3,Zero,tmp,3)
call DGEMM_('N','N',3,3,3,One,T,3,tmp,3,Zero,R_Q,3)
call RecPrt('R_Q',' ',R_Q,3,3)
QQ = R_Q(3,3)-Half*(R_Q(1,1)+R_Q(2,2))
R_Test(3) = R_Q(3,3)/(Two*QQ)+R_Test(3)
call RecPrt('EC(Bond system) New',' ',R_Test,1,3)

call DGEMM_('N','N',3,3,3,One,T,3,R_Test,3,Zero,tmp,3)
call RecPrt('EC New',' ',tmp,1,3)
write(u6,*)
! Store the relevant data
Dipole_Rot_A = rMu_A(3)
Dipole_Rot_B = rMu_B(3)
Dipole_Rot_AB = rMu_AB(3)
!R_A = Zero
!R_B = T_R
! Rotate EC(1,ii) to the new coordinate system
call DGEMM_('T','N',3,1,3,One,T,3,EC(1,ii),3,Zero,R_test,3)
call RecPrt('EC(ii)',' ',R_Test,1,3)
R_temp = R_Test(3)
! Rotate EC(1,jj) to the new coordinate system
call DGEMM_('T','N',3,1,3,One,T,3,EC(1,jj),3,Zero,R_test,3)
call RecPrt('EC(jj)',' ',R_Test,1,3)
! Put A in origo, and R_B the appropriate distance away.
R_B = R_Test(3)-R_temp
R_A = Zero
write(u6,*) 'Dipoles = ',Dipole_Rot_A,Dipole_Rot_B,Dipole_Rot_AB

return

end subroutine Rotate_Dipole
