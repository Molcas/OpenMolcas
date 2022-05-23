************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Rotate_Dipole(rMP, EC, nij, nElem, nPert, ij, ii, jj,
     &                         Dipole_Rot_A,Dipole_Rot_B,Dipole_Rot_AB,
     &                         R_A,R_B)
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 rMP(nij,0:nElem-1,0:nPert-1), EC(3,nij), E_Bond(3)
      Real*8 T(3,3), R_test(3), R_Q(3,3), Tmp(3,3),
     &       rMu_AB(3), rMu_A(3), rMu_B(3)
*
* The three dipole components for the bond
c      Print *,'Test 0'
c      Call xFlush(6)
C??      x_my=rMP(ij,1,0)
C??      y_my=rMP(ij,2,0)
C??      z_my=rMP(ij,3,0)
* Vector from ii to jj.
      x_R =   EC(1,ii) - EC(1,jj)
      y_R =   EC(2,ii) - EC(2,jj)
      z_R =   EC(3,ii) - EC(3,jj)
      T_R =Sqrt(x_R**2+y_R**2+z_R**2)
C??      Write (6,*) 'My=',x_my,y_my,z_my
C??      Write (6,*) 'R=',x_R,y_R,z_R
* Unit vector in the direction from ii to jj
      E_Bond(1)=x_r/T_R
      E_Bond(2)=y_r/T_R
      E_Bond(3)=z_r/T_R
* Get transformation matrix (T) that rotates the coordinate system to align with the bond
      Call GS(E_bond,1,T,3,.True.,.False.)
      Call RecPrt('T-matrix',' ',T,3,3)
*
      Call RecPrt('EC(*,ij) origional',' ',EC(1,ij),1,3)
      Call RecPrt('EC(*,ii) origional',' ',EC(1,ii),1,3)
      Call RecPrt('EC(*,jj) origional',' ',EC(1,jj),1,3)
* Rotate EC(1,ij) to the new coordinate system
      Call DGEMM_('T','N',
     &            3,1,3,
     &            1.0d0,T,3,
     &            EC(1,ij),3,
     &            0.0d0,R_test,3)
      Call RecPrt('EC(Bond system)',' ',R_Test,1,3)
*
* Rotate the dipole moment for the bond and ii and jj to the new coordinate system
      tmp(1,1)=rMP(ij,1,0)
      tmp(2,1)=rMP(ij,2,0)
      tmp(3,1)=rMP(ij,3,0)
      Call DGEMM_('T','N',
     &            3,1,3,
     &            1.0d0,T,3,
     &            tmp,nij,
     &            0.0d0,rMu_AB,3)
      Call RecPrt('rMu_AB',' ',rMu_AB,1,3)
      tmp(1,1)=rMP(ii,1,0)
      tmp(2,1)=rMP(ii,2,0)
      tmp(3,1)=rMP(ii,3,0)
      Call DGEMM_('T','N',
     &            3,1,3,
     &            1.0d0,T,3,
     &            tmp,nij,
     &            0.0d0,rMu_A,3)
      Call RecPrt('rMu_A',' ',rMu_A,1,3)
      tmp(1,1)=rMP(jj,1,0)
      tmp(2,1)=rMP(jj,2,0)
      tmp(3,1)=rMP(jj,3,0)
      Call DGEMM_('T','N',
     &            3,1,3,
     &            1.0d0,T,3,
     &            tmp,nij,
     &            0.0d0,rMu_B,3)
      Call RecPrt('rMu_B',' ',rMu_B,1,3)
*
      R_Q(1,1)=rMP(ij,4,0)
      R_Q(1,2)=half*rMP(ij,5,0)
      R_Q(1,3)=half*rMP(ij,6,0)
      R_Q(2,2)=rMP(ij,7,0)
      R_Q(2,3)=half*rMP(ij,8,0)
      R_Q(3,3)=rMP(ij,9,0)
      R_Q(2,1)=R_Q(1,2)
      R_Q(3,1)=R_Q(1,3)
      R_Q(3,2)=R_Q(2,3)
      Call RecPrt('R_Q',' ',R_Q,3,3)
      Call DGEMM_('N','T',
     &            3,3,3,
     &            1.0d0,R_Q,3,
     &            T,3,
     &            0.0d0,tmp,3)
      Call DGEMM_('N','N',
     &            3,3,3,
     &            1.0d0,T,3,
     &            tmp,3,
     &            0.0d0,R_Q,3)
      Call RecPrt('R_Q',' ',R_Q,3,3)
      QQ=R_Q(3,3)-Half*(R_Q(1,1)+R_Q(2,2))
      R_Test(3)=R_Q(3,3)/(Two*QQ) + R_Test(3)
      Call RecPrt('EC(Bond system) New',' ',R_Test,1,3)
*
      Call DGEMM_('N','N',
     &            3,3,3,
     &            1.0d0,T,3,
     &            R_Test,3,
     &            0.0d0,tmp,3)
      Call RecPrt('EC New',' ',tmp,1,3)
      Write (6,*)
* Store the relavant data
      Dipole_Rot_A  = rMu_A(3)
      Dipole_Rot_B  = rMu_B(3)
      Dipole_Rot_AB = rMu_AB(3)
C      R_A           = 0.0D0
C      R_B           = T_R
* Rotate EC(1,ii) to the new coordinate system
      Call DGEMM_('T','N',
     &            3,1,3,
     &            1.0d0,T,3,
     &            EC(1,ii),3,
     &            0.0d0,R_test,3)
      Call RecPrt('EC(ii)',' ',R_Test,1,3)
      R_temp        = R_Test(3)
* Rotate EC(1,jj) to the new coordinate system
      Call DGEMM_('T','N',
     &            3,1,3,
     &            1.0d0,T,3,
     &            EC(1,jj),3,
     &            0.0d0,R_test,3)
      Call RecPrt('EC(jj)',' ',R_Test,1,3)
* Put A in origo, and R_B the appropriate distance away.
      R_B           = R_Test(3) - R_temp
      R_A           = 0.0D0
      write(6, *)'Dipoles = ',Dipole_Rot_A,Dipole_Rot_B,Dipole_Rot_AB
*
      Return
      End
