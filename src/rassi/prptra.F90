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

subroutine PRPTRA(ND,NCPL,TRA)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ND, NCPL
real(kind=wp), intent(in) :: TRA(ND,NCPL)
integer(kind=iwp) :: ICPL, ID, IEND, ISTA

if ((ND < 0) .or. (NCPL < 0)) then
  call WarningMessage(2,'Program bug: Erroneous call to PRPTRA.')
  write(u6,*) 'PRPTRA error: Wrong arguments.'
  write(u6,*) 'PRPTRA: ND,NCPPL=',ND,NCPL
  call ABEND()
end if
if ((ND == 0) .or. (NCPL == 0)) then
  call WarningMessage(1,'Program bug? Strange call to PRPCSF.')
  write(u6,*) 'PRPTRA warning: Strange arguments.'
  write(u6,*) 'PRPTRA: ND,NCPPL=',ND,NCPL
else
  do ISTA=1,NCPL,5
    IEND = min(NCPL,ISTA+4)
    write(u6,*)
    write(u6,'(8x,5(I8,8X))') (ICPL,ICPL=ISTA,IEND)
    do ID=1,ND
      write(u6,'(1x,5F16.8)') (TRA(ID,ICPL),ICPL=ISTA,IEND)
    end do
  end do
end if

end subroutine PRPTRA
