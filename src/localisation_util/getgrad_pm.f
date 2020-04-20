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
#define _TEST2_
#ifdef _TEST2_
      SubRoutine GetGrad_PM(nAtoms,nOrb2Loc,iTab_Ptr,PA,GradNorm,Rmat,
     &                      Debug)
#else
      SubRoutine GetGrad_PM(nAtoms,nOrb2Loc,iTab_Ptr,GradNorm,Rmat,
     &                      Debug)
#endif
C
C     Thomas Bondo Pedersen, December 2005.
C
C     Purpose: compute the gradient of the Pipek-Mezey functional.
C
      Implicit Real*8 (a-h,o-z)
      Integer iTab_Ptr(nAtoms)
      Real*8 Rmat(nOrb2Loc,nOrb2Loc)
#ifdef _TEST2_
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
#endif
      Logical Debug
#include "WrkSpc.fh"

*     Call FZero(Rmat,nOrb2Loc**2)
#ifdef _TEST2_
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
*     Call RecPrt('RMAT',' ',RMAT,nOrb2Loc,nOrb2Loc)
#else
      RMat(:,:)=0.0D0
      Do iAtom = 1,nAtoms
         ip0 = iTab_Ptr(iAtom) - 1
         Do j = 1,nOrb2Loc
            Rjj = Work(ip0+nOrb2Loc*(j-1)+j)
            Do i = 1,nOrb2Loc
               Rmat(i,j) = Rmat(i,j) +
     &                     Work(ip0+nOrb2Loc*(j-1)+i)*Rjj
            End Do
         End Do
      End Do
*     Call RecPrt('RMAT',' ',RMAT,nOrb2Loc,nOrb2Loc)
#endif

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
