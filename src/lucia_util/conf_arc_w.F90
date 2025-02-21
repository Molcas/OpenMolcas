!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine CONF_ARC_W(IOCC_MIN,IOCC_MAX,NORB,NEL,IVERTEXW,IARCW)
! Obtain arcweights for single and double occupied arcs
! from vertex weights
!
! Jeppe Olsen, October 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NORB, IOCC_MIN(NORB), IOCC_MAX(NORB), NEL, IVERTEXW(NORB+1,NEL+1)
integer(kind=iwp), intent(out) :: IARCW(NORB,NEL,2)
integer(kind=iwp) :: I, J, K

IARCW(:,:,:) = 0
! IARCW(I,J,K) is weight of arc with occupation K ending at (I,J)
! IARCW(I,J,K) = Sum(J-K < L <= J)   IVERTEXW(I-1,L)
do I=1,NORB
  do J=1,NEL
    if ((IOCC_MIN(I) <= J) .and. (J <= IOCC_MAX(I))) then
      do K=1,NEL
        if (K == 1) IARCW(I,J,K) = IVERTEXW(I-1+1,J+1)
        if ((K == 2) .and. (J >= 2)) IARCW(I,J,K) = IVERTEXW(I-1+1,J+1)+IVERTEXW(I-1+1,J-1+1)
      end do
    end if
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Arc weights for single occupied arcs'
call IWRTMA(IARCW(1,1,1),NORB,NEL,NORB,NEL)
write(u6,*) ' Arc weights for double occupied arcs'
call IWRTMA(IARCW(1,1,2),NORB,NEL,NORB,NEL)
#endif

end subroutine CONF_ARC_W
