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

subroutine DELF(FNAM,INUM1,INUM2)

use Definitions, only: iwp

implicit none
character(len=6), intent(in) :: FNAM
integer(kind=iwp), intent(in) :: INUM1, INUM2
integer(kind=iwp) :: I, LU
character(len=8) :: FN
integer(kind=iwp), external :: IsFreeUnit

FN(1:6) = FNAM
do I=inum1,inum2
  write(fn(7:8),'(I2.2)') I
  !write(u6,*) 'File ',FN,' to be deleted'
  LU = IsFreeUnit(8)
  call Molcas_Open(LU,fn)
  close(LU,status='DELETE')
end do

return

end subroutine DELF
