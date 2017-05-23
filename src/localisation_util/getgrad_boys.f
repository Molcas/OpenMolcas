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
      SubRoutine GetGrad_Boys(nOrb2Loc,ipLbl,nComp,Rmat,GradNorm,Debug)
C
C     Thomas Bondo Pedersen, December 2005.
C
C     Purpose: compute R-matrix and gradient norm for Boys functional.
C
      Implicit Real*8 (a-h,o-z)
      Integer ipLbl(nComp)
      Real*8  Rmat(nOrb2Loc,nOrb2Loc)
      Logical Debug
#include "WrkSpc.fh"

      Call FZero(Rmat,nOrb2Loc**2)
      Do iComp = 1,nComp
         ip0 = ipLbl(iComp) - 1
         Do j = 1,nOrb2Loc
            Rjj = Work(ip0+nOrb2Loc*(j-1)+j)
            Do i = 1,nOrb2Loc
               Rmat(i,j) = Rmat(i,j) +
     &                     Work(ip0+nOrb2Loc*(j-1)+i)*Rjj
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
         Write(6,*) 'GetGrad_Boys: functional = Tr(R) = ',Fun
      End If

      End
