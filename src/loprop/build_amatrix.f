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
      Subroutine Build_AMatrix(nAtoms,iANr,AMatrix,AInvMatrix,EC,
     &                         nij,Alpha)
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "constants.fh"
      Real*8 A(3), B(3), AMatrix(nAtoms,nAtoms),
     &       AInvMatrix(nAtoms,nAtoms), EC(3,nij)
      Integer iANr(nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
      Do iAtom = 1, nAtoms
         ii = iAtom*(iAtom+1)/2
         call dcopy_(3,EC(1,ii),1,A,1)
         R_BS_i=Bragg_Slater(iANr(iAtom))
         Do jAtom = 1, iAtom-1
            jj = jAtom*(jAtom+1)/2
            call dcopy_(3,EC(1,jj),1,B,1)
            R_BS_j=Bragg_Slater(iANr(jAtom))
            rij2=(A(1)-B(1))**2
     &          +(A(2)-B(2))**2
     &          +(A(3)-B(3))**2
            rij02=((R_BS_i+R_BS_j))**2
            AMatrix(iAtom,jAtom)= Exp(-Alpha*(rij2/rij02))/Two
            AMatrix(jAtom,iAtom)= Exp(-Alpha*(rij2/rij02))/Two
         End Do
*
         temp=zero
         Do jAtom = 1, nAtoms
            jj = jAtom*(jAtom+1)/2
            call dcopy_(3,EC(1,jj),1,B,1)
            R_BS_j=Bragg_Slater(iANr(jAtom))
            rij2=(A(1)-B(1))**2
     &          +(A(2)-B(2))**2
     &          +(A(3)-B(3))**2
            rij02=((R_BS_i+R_BS_j))**2
            If (jAtom.ne.iAtom) Then
               temp =temp- Exp(-Alpha*(rij2/rij02))/Two
            End If
         End Do
         AMatrix(iAtom,iAtom)=Temp
*
      End Do
c      Call RecPrt('A-matrix','(5G12.5)',AMatrix,nAtoms,nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
      Shift=Zero
      Do i = 1, nAtoms
         Do j = 1, i
            Shift=Max(Abs(AMatrix(i,j)),Shift)
         End Do
      End Do
      Shift=Two*Shift
C     Write (*,*) 'Shift=',Shift
      Call DaXpY_(nAtoms**2,Shift,One,0,AMatrix,1)
c      Call RecPrt('A-matrix(Shifted)','(5G12.5)',AMatrix,nAtoms, nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
      Call MINV(AMatrix,AInvMatrix,ISING,DET,nAtoms)
c      Call RecPrt('AInv-matrix',' ',AIntMatrix,nAtoms,nAtoms)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
