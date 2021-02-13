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

implicit real*8(A-H,O-Z)
character*(*) Text
dimension XInt(*)
integer nBas(*)

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Loop over pairs of symmetry labels. Skip all offdiagonal         *
!     blocks.                                                          *
!                                                                      *
!----------------------------------------------------------------------*

lText = min(120,len(Text))
write(6,'(6X,A)') Text(1:lText)
iOff = 0
do iSym=1,nSym
  nBs = nBas(iSym)
  if (nBs /= 0) then
    write(6,'(6X,A,I2)') 'Symmetry species',iSym
    call TriPrt(' ',' ',XInt(iOff+1),nBs)
  end if
end do

!----------------------------------------------------------------------*
!     Terminate procedure                                              *
!----------------------------------------------------------------------*

return

end
