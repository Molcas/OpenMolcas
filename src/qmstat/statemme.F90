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
subroutine StateMME(MoOrNot,nAObas,nMObas,nState,nTyp,MME,iCent,Cha,Dip,Qua)

use qmstat_global, only: MxMltp
use Index_Functions, only: nTri3_Elem, nTri_Elem
use Data_Structures, only: Alloc1DArray_Type
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: MoOrNot
integer(kind=iwp), intent(in) :: nAObas, nMObas, nState, nTyp, iCent(nTri_Elem(nAObas))
type(Alloc1DArray_Type), intent(in) :: MME(nTri3_Elem(MxMltp))
real(kind=wp), intent(inout) :: Cha(nTri_Elem(nState),*), Dip(nTri_Elem(nState),3,*), Qua(nTri_Elem(nState),6,*)

if (.not. MoOrNot) then
  call StateMMEao(nAObas,nState,nTyp,MME,iCent,Cha,Dip,Qua)
else
  call StateMMEmo(nAObas,nMObas,nState,nTyp,MME,iCent,Cha,Dip,Qua)
end if

return

end subroutine StateMME
