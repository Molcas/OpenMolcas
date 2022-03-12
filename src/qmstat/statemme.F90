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
!
!-- This is just an interface for the state transformation. Either
!   we use usual AO-basis route, or we take the reduced MO-basis
!   route.
!
      Subroutine StateMME(MoOrNot,nAObas,nMObas,nState,nTyp,iCi,iBigT   &
     &                   ,iMME,iCent,ipAvRed,Cha,Dip,Qua)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6),iCent(MxBas**2)
      Dimension Cha(MxStOT,MxQCen),Dip(MxStOT,3,MxQCen)
      Dimension Qua(MxStOT,6,MxQCen)
      Logical MoOrNot

      If(.not.MoOrNot) then
        Call StateMMEao(nAObas,nState,nTyp,iBigT,iMME,iCent,Cha,Dip,Qua)
      Else
        Call StateMMEmo(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME,iCent  &
     &                 ,ipAvRed,Cha,Dip,Qua)
      Endif

      Return
      End
