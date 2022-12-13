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

subroutine Get_Cmo(CMO,nCMO)

use RunFile_data, only: lw
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nCMO
real(kind=wp), intent(out) :: CMO(nCMO)
integer(kind=iwp) :: mCMO
logical(kind=iwp) :: Found
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: nBas(0:7) = -1, nSym = -1
#endif
character(len=lw) :: Label

Label = 'Last orbitals'
call qpg_dArray(Label,Found,mCmo)
if (.not. Found) Label = 'Guessorb'
call Get_dArray_chk(Label,CMO,nCMO)

return

end subroutine Get_Cmo
