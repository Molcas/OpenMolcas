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
subroutine CONF_VERTEX_W(IOCC_MIN,IOCC_MAX,NORB,NEL,IVERTEXW)
! Obtain vertex weights for configuration graph
!
! Jeppe Olsen, October 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NORB, IOCC_MIN(NORB), IOCC_MAX(NORB), NEL
integer(kind=iwp), intent(out) :: IVERTEXW(NORB+1,NEL+1)
integer(kind=iwp) :: IEL, IORB

!write(u6,*) ' CONF_VERTEX : NORB, NEL = ',NORB,NEL
IVERTEXW(:,:) = 0
IVERTEXW(1,1) = 1
do IORB=1,NORB
  !Jesper DO IEL=0,NEL
  do IEL=IOCC_MIN(IORB),IOCC_MAX(IORB)
    !Jesper if ((IOCC_MIN(IORB) <= IEL) .and. (IEL <= IOCC_MAX(IORB))) then

    if (IEL == 0) IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)

    if (IEL == 1) IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)+IVERTEXW(IORB-1+1,IEL+1-1)

    if (IEL >= 2) IVERTEXW(IORB+1,IEL+1) = IVERTEXW(IORB-1+1,IEL+1)+IVERTEXW(IORB-1+1,IEL+1-1)+IVERTEXW(IORB-1+1,IEL+1-2)

    !Jesper end if
  end do
end do
! Check whether a configuration has enough singly occupied orbitals
!Jesper ??? call REDUCE_VERTEX_WEIGHTS(IVERTEXW(NORB+1,NEL+1),IVERTEXW,NEL,NORB,IOCC_MIN,IOCC_MAX)

#ifdef _DEBUGPRINT_
write(u6,*) ' Vertex weights as an (NORB+1)*(NEL+1) matrix'
call IWRTMA(IVERTEXW,NORB+1,NEL+1,NORB+1,NEL+1)
#endif

end subroutine CONF_VERTEX_W
