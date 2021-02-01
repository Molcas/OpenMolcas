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
      Subroutine Strtch(xyz,nCent,Avst,B,lWrite,Label,dB,ldB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 B(3,nCent), xyz(3,nCent), dB(3,nCent,3,nCent), R(3)
      Logical lWrite, ldB
      Character*8 Label
#include "angstr.fh"
*
      R(1)=xyz(1,2)-xyz(1,1)
      R(2)=xyz(2,2)-xyz(2,1)
      R(3)=xyz(3,2)-xyz(3,1)
      R2=R(1)**2+R(2)**2+R(3)**2
      RR=Sqrt(R2)
      Avst=RR
*
      aRR=RR*angstr
      If (lWrite) Write (6,'(1X,A,A,2(F10.6,A))') Label,
     &      ' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
      If (aRR.LT.1d-6) then
         Call WarningMessage(2,'Abend in Strtch')
         Write (6,*) '***************** ERROR **********************'
         Write (6,*) ' Short (or negative) distance for coordinate: '
         Write (6,'(1X,A,A,2(F10.6,A))') Label,
     &      ' : Bond Length=',aRR,' / Angstrom',RR,' / bohr'
         Write (6,*) '**********************************************'
         Write (6,*)
         Call Quit_OnUserError()
      EndIf
*
*---- Compute the WDC B-matrix.
*
      B(1,1)=-R(1)/RR
      B(2,1)=-R(2)/RR
      B(3,1)=-R(3)/RR
*
*---- Renormalize
*
      xRR=Sqrt(B(1,1)**2+B(2,1)**2+B(3,1)**2)
      B(1,1)=B(1,1)/xRR
      B(2,1)=B(2,1)/xRR
      B(3,1)=B(3,1)/xRR
*
*.... Utilize translational invariance.
      B(1,2)=-B(1,1)
      B(2,2)=-B(2,1)
      B(3,2)=-B(3,1)
*
*---- Compute the cartesian derivative of the B-matrix.
*
      If (ldB) Then
*
         Do i=1,3
            Do j=1,i
               If (i.eq.j) Then
                  dB(i,1,j,1)= (One-B(j,1)*B(i,1))/RR
               Else
                  dB(i,1,j,1)= (-B(j,1)*B(i,1))/RR
               End If
               dB(j,1,i,1)= dB(i,1,j,1)
*
               dB(i,2,j,1)=-dB(i,1,j,1)
               dB(j,1,i,2)= dB(i,2,j,1)
*
               dB(i,1,j,2)=-dB(i,1,j,1)
               dB(j,2,i,1)= dB(i,1,j,2)
*
               dB(i,2,j,2)=-dB(i,2,j,1)
               dB(j,2,i,2)=dB(i,2,j,2)
            End Do
         End Do
*
      End If
      Return
      End
