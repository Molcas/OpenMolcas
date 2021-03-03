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

subroutine Print_Input(Title,nSym,nBas,wSet,nSet)

use Definitions, only: wp, iwp, u6

implicit none
character(len=72), intent(in) :: Title
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nSet
real(kind=wp), intent(in) :: wSet(nSet)
integer(kind=iwp) :: iS, iSym
integer(kind=iwp), parameter :: lPaper = 132

call Banner(Title,1,lPaper-7)
write(u6,*)
write(u6,*)
write(u6,*) 'Number of symmetries:',nSym
write(u6,*) 'Basis functions:',(nBas(iSym),iSym=1,nSym)
write(u6,*) 'Normalized weights:',(wSet(iS),iS=1,nSet)

return

end subroutine Print_Input
