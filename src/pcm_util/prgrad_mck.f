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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine PrGrad_mck(Label,Grad,nGrad,Names,iPrint)
************************************************************************
*                                                                      *
* Object: to print set gradient with respect to the symmetrical dis-   *
*         placements.                                                  *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             October '91                                              *
************************************************************************
      use Symmetry_Info, only: lIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
      Real*8 Grad(nGrad)
      Real*8 CGrad(3,MxAtom)
      Character CNames(MxAtom)*(LENIN5)
      Character Label*(*), Names(nGrad)*(LENIN6)
      Character Namei*(LENIN5)
*
*     Call qEnter('PrGrad')
*
      Write (6,*)
      Call Banner(Label,1,Len(Label)+30)
      Write (6,*)
      If (iPrint.EQ.4) then
         Call TrGrd_Alaska_(CGrad,CNames,Grad,nGrad,iCen)
         Write (6,'(1x,A,A)') ' Irreducible representation: ',lIrrep(0)
         Write (6,'(1x,A)')
     &              '--------------------------------------------------'
         Write (6,'(1x,A)')
     &              '                    X           Y           Z     '
         Write (6,'(1x,A)')
     &              '--------------------------------------------------'
         Do iGrad = 1, iCen
            TempX = CGrad(1,iGrad)
            TempY = CGrad(2,iGrad)
            TempZ = CGrad(3,iGrad)
            Namei = CNames(iGrad)
            Write (6,'(2X,A,3X,3F12.6)') Namei, TempX, TempY, TempZ
         End Do
         Write (6,'(1x,A)')
     &              '--------------------------------------------------'
      else
*
*        Modified by Luca De Vico november 2005 Teokem
*        I need to print the full gradient vector
*        to use it for constrained optimizations
*        with TRANSVERSE option in SLAPAF MODULE
*
*         mGrad=Min(21,nGrad)
*
         mGrad=nGrad
         Write (6,'(15x,A,A)') ' Irreducible representation: ',lIrrep(0)
         Write (6,*)
         Do iGrad = 1, mGrad
            Temp = Grad(iGrad)
            If (Abs(Temp).lt.1.0D-15) Temp = Zero
            Write (6,'(16X,A,15X,E15.7)') Names(iGrad), Temp
         End Do
*
*         If (nGrad.gt.21) Then
*            Write (6,*)
*            Write (6,*) '   ... list is truncated ... '
*            Write (6,*)
*         End If
*
*        End of modfications
*
      EndIf
      Write (6,*)
*
*     Call qExit('PrGrad')
      Return
      End
