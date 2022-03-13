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

! This is just an interface for the state transformation. Either
! we use usual AO-basis route, or we take the reduced MO-basis route.
subroutine StateMME(MoOrNot,nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME,iCent,ipAvRed,Cha,Dip,Qua)

use Index_Functions, only: nTri3_Elem
use Definitions, only: wp, iwp

implicit none
#include "maxi.fh"
logical(kind=iwp) :: MoOrNot
integer(kind=iwp) :: nAObas, nMObas, nState, nTyp, iCi, iBigT, iMME(nTri3_Elem(MxMltp)), iCent(MxBas**2), ipAvRed
real(kind=wp) :: Cha(MxStOT,MxQCen), Dip(MxStOT,3,MxQCen), Qua(MxStOT,6,MxQCen)

if (.not. MoOrNot) then
  call StateMMEao(nAObas,nState,nTyp,iBigT,iMME,iCent,Cha,Dip,Qua)
else
  call StateMMEmo(nAObas,nMObas,nState,nTyp,iCi,iBigT,iMME,iCent,ipAvRed,Cha,Dip,Qua)
end if

return

end subroutine StateMME
