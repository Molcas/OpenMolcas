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

subroutine PrDiOp(Text,nSym,nBas,XInt)
!***********************************************************************
!                                                                      *
!     Object: Print a diagonal block matrix                            *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Text
integer(kind=iwp), intent(in) :: nSym, nBas(nSym)
real(kind=wp), intent(in) :: XInt(*)
integer(kind=iwp) :: iOff, iSym, lText, nBs

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Loop over pairs of symmetry labels. Skip all offdiagonal         *
!     blocks.                                                          *
!                                                                      *
!----------------------------------------------------------------------*

lText = min(120,len(Text))
write(u6,'(6X,A)') Text(1:lText)
iOff = 1
do iSym=1,nSym
  nBs = nBas(iSym)
  if (nBs /= 0) then
    write(u6,'(6X,A,I2)') 'Symmetry species',iSym
    call TriPrt(' ',' ',XInt(iOff),nBs)
  end if
  iOff = iOff+nBs*(nBs+1)/2
end do

!----------------------------------------------------------------------*
!     Terminate procedure                                              *
!----------------------------------------------------------------------*

return

end subroutine PrDiOp
