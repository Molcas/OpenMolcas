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
      Subroutine Dynamic_Properties(Temp,nAtoms,rMP,nij,nPert,
     &                              nElem,Delta,EC,Polar,
     &                              iANr,Bond_Threshold,ChPol,
     &                              ChPolBB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Temp(nij), A(3), B(3), EC(3,nij), Polar(6,nij),
     &       rMP(nij,0:nElem-1,0:nPert-1), ChPol(6, nij),
     &       ChPolBB(6, nij)
      Integer iAnr(nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
C     Call RecPrt('rMP',' ',rMP,nij*nElem,nPert)
      Write (6,*)
      Write (6,*)  ' D y n a m i c  P r o p e r t i e s'
      Write (6,*)
      Write (6,*)  ' Properties computed with FFPT'
      Write (6,*)

      Do iPol = 1, 6
         Do iAtom = 1, nAtoms
            Do jAtom = 1, iAtom
               ij=iAtom*(iAtom-1)/2+jAtom
               ChPol(iPol, ij) = 0.0
               ChPolBB(iPol, ij) = 0.0
            End Do
         End Do
      End Do


*
*     iPol: index vector for polarizability
*           (1,2,3,4,5,6)=(xx,yx,yy,zx,zy,zz)
*
      Do iPol = 1, 6
         Call FZero(Temp,nij)
C        Write (6,*)
         Do iAtom = 1, nAtoms
            ii = iAtom*(iAtom+1)/2
            call dcopy_(3,EC(1,ii),1,A,1)
            Do jAtom = 1, iAtom
               jj = jAtom*(jAtom+1)/2
               call dcopy_(3,EC(1,jj),1,B,1)
*
               ij=iAtom*(iAtom-1)/2+jAtom
*                                                                      *
************************************************************************
*                                                                      *
*              Polarizabilities: alpha(iAtom,jAtom,iCar,jCar)          *
*                                                                      *
************************************************************************
*                                                                      *
*              iCar, jCar: index of cartesian for each perturbation
*                        1=x, 2=y, 3=z
               iCar = Int((One+sqrt(Eight*DBLE(iPol)-Three))/Two)
               jCar = iPol   -iCar*(iCar-1)/2
*
*              iPert, jPert: index vector to actuall perturbation
*                    (1,2,3,4,5,6)=(+dx,-dx,+dy,-dy,+dz,-dz)
               iPert=(jCar-1)*2+1
               jPert=iPert+1
*
*              Contribution due to change of localized dipole moment
*
*------------- mu(iAtom,jAtom,iCar,+dF(jCar))
*------------- mu(iAtom,jAtom,iCar,-dF(jCar))
*
               Pol1a=(rMP(ij,iCar,iPert)-rMP(ij,iCar,jPert))/(Two*Delta)
*
               iPert_=(iCar-1)*2+1
               jPert_=iPert_+1
               Pol1b=(rMP(ij,jCar,iPert_)-rMP(ij,jCar,jPert_))
     &              /(Two*Delta)
               Pol1 = Half * (Pol1a+Pol1b)
C              Write (*,*) 'Pol1',ij,iCar,iPert,jPert
C              Write (*,*) rMP(ij,iCar,iPert),rMP(ij,iCar,jPert)
*
*              Contribution due to change of localized charges
*
               If (iAtom .ne. jAtom) Then
                  Rij_iCar=B(iCar)-A(iCar)
C                 Write (*,*) rMP(ij,0,iPert),rMP(ij,0,jPert)
                  Pol2= (rMP(ij,0,iPert)-rMP(ij,0,jPert))
     &                * Rij_iCar/(Two*Delta)
               Else
                  Pol2=Zero
               End If
*
C              Write (*,*) Pol1, Pol2

               Temp(ij) = Temp(ij) + Pol1 + Pol2
               Polar(iPol,ij)=Temp(ij)
*                                                                      *
************************************************************************
*              Compute Charge contrib                                  *
************************************************************************
*                                                                      *
               ChPol(iPol,ij)=ChPol(iPol,ij) + Pol2
               ChPolBB(iPol,ij)=ChPolBB(iPol,ij) + Pol2
************************************************************************


*                                                                      *
************************************************************************
*                                                                      *
            End Do   ! jAtom
         End Do      ! iAtom
*                                                                      *
************************************************************************
*                                                                      *
      End Do   ! iPol
*                                                                      *
************************************************************************
*                                                                      *
*            Move the polarizabilities to the atoms if needed          *
*                                                                      *
************************************************************************
*                                                                      *
      Call Move_Polar(Polar,EC,nAtoms,nij,iANr,Bond_Threshold)
      Call Move_Polar(ChPol,EC,nAtoms,nij,iANr,Bond_Threshold)

      Return

      End
