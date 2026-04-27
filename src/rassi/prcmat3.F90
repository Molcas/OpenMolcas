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

subroutine PRCMAT3(NSS,SMATR,SMATI,DIR)
! Write out spin matrix elements in parsable format

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSS, DIR
real(kind=wp), intent(in) :: SMATR(NSS,NSS), SMATI(NSS,NSS)
integer(kind=iwp) :: f_iostat, ISS, JSTA, LU
logical(kind=iwp) :: is_error
character(len=200) :: FILENAME
character :: DIRECTION
integer(kind=iwp), external :: IsFreeUnit

write(DIRECTION,'(I1)') DIR
FILENAME = 'spin-'//DIRECTION//'.txt'
Lu = IsFreeUnit(88)
call molcas_open_ext2(Lu,FILENAME,'SEQUENTIAL','FORMATTED',f_iostat,.false.,1,'REPLACE',is_error)
write(Lu,*) '#NROW NCOL REAL IMAG'
do JSTA=1,NSS
  do ISS=1,NSS
    write(Lu,'(I6,1X,I6,A1,ES25.16,A1,ES25.16)') ISS,JSTA,' ',SMATR(ISS,JSTA),' ',SMATI(ISS,JSTA)
  end do
end do
close(Lu)

end subroutine PRCMAT3
