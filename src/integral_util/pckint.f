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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine PckInt(abab,mZeta,nab,ab,rKappa,Mode,Zeta,nZeta,
     &                  qKappa)
************************************************************************
*                                                                      *
* Object: to keep the diagonal angular indices of a integral batch.    *
*         The integrals are also stripped of the prefactor due to      *
*         the product of gaussians. In case of numerical different-    *
*         iation the prefactor will be due to the undifferentiated     *
*         charge densities.                                            *
*                                                                      *
* Called from: k2Loop                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             April '92                                                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "real.fh"
      Real*8 abab(mZeta,nab,nab), ab(nZeta,nab), rKappa(mZeta),
     &       Zeta(mZeta), qKappa(mZeta)
      Logical Mode
*
      iRout=399
      iPrint=nPrint(399)
*
      If (Mode) Then
*--------Integrals
         Do iab = 1, nab
            Do iZeta = 1, mZeta
               ab(iZeta,iab) =
     &             Sqrt(Sqrt(Two*Zeta(iZeta))*Abs(abab(iZeta,iab,iab)))
     &                       / rKappa(iZeta)
            End Do
         End Do
      Else
*--------Integrals for numerical estimation of the gradient.
         Do iab = 1, nab
            Do iZeta = 1, mZeta
               ab(iZeta,iab) = Sqrt(Two*Zeta(iZeta))*
     &                       abab(iZeta,iab,iab)
     &                       / (rKappa(iZeta)*qKappa(iZeta))
            End Do
         End Do
      End If
      If (iPrint.ge.99) Then
         Write (6,*) 'nZeta,mZeta=',nZeta,mZeta
         Call RecPrt(' abab','(5G20.10)',abab,mZeta,nab**2)
         Call RecPrt(' rKappa','(5G20.10)',rKappa,mZeta,1)
         Call RecPrt(' Zeta  ','(5G20.10)',Zeta  ,mZeta,1)
         Do iab = 1, nab
            Call RecPrt(' ab ','(5G20.10)',ab(1,iab),mZeta,1)
         End Do
      End If
*
      Return
      End
