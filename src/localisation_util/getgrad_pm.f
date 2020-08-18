************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine GetGrad_PM(nAtoms,nOrb2Loc,PA,GradNorm,Rmat,Debug)
C
C     Thomas Bondo Pedersen, December 2005.
C
C     Purpose: compute the gradient of the Pipek-Mezey functional.
C
      Implicit Real*8 (a-h,o-z)
      Real*8 Rmat(nOrb2Loc,nOrb2Loc)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Logical Debug
#include "WrkSpc.fh"

      RMat(:,:)=0.0D0
      Do iAtom = 1,nAtoms
         Do j = 1,nOrb2Loc
            Rjj = PA(j,j,iAtom)
            Do i = 1,nOrb2Loc
               Rmat(i,j) = Rmat(i,j) +
     &                     PA(i,j,iAtom)*Rjj
            End Do
         End Do
      End Do

      GradNorm = 0.0d0
      Do i = 1,nOrb2Loc-1
         Do j = i+1,nOrb2Loc
            GradNorm = GradNorm
     &               + (Rmat(i,j)-Rmat(j,i))**2
         End Do
      End Do
      GradNorm = 4.0d0*sqrt(GradNorm)

      If (Debug) Then
         Fun = 0.0d0
         Do i = 1,nOrb2Loc
            Fun = Fun + Rmat(i,i)
         End Do
         Write(6,*) 'GetGrad_PM: functional = Tr(R) = ',Fun
      End If

      End
