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
subroutine CONF_GRAPH(IOCC_MIN,IOCC_MAX,NORB,NEL,IARCW,NCONF,ISCR)
! A group of configurations is described by the
! accumulated min and max, IOCC_MIN and IOCC_MAX.
!
! Find arcweights of corresponding graph and total number
! of configurations (all symmetries)
!
! IARCW(I,J,K) : gives weight of arc ending at vertex (I,J)
!                with occupation K (=1,2)
! ISCR : Length should be (NORB+1)*(NEL+1)
!
! Jeppe Olsen, Oct. 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NORB, IOCC_MIN(NORB), IOCC_MAX(NORB), NEL
integer(kind=iwp), intent(out) :: IARCW(NORB,NEL,2), NCONF, ISCR(NORB+1,NEL+1)

! Set up vertex weights
call CONF_VERTEX_W(IOCC_MIN,IOCC_MAX,NORB,NEL,ISCR)
NCONF = ISCR(NORB+1,NEL+1)
! Obtain arcweights from vertex weights
!write(u6,*) ' CONF_GRAPH, NORB, NEL = ',NORB,NEL
call CONF_ARC_W(IOCC_MIN,IOCC_MAX,NORB,NEL,ISCR,IARCW)

#ifdef _DEBUGPRINT_
write(u6,*) ' IOCMIN and IOCMAX'
call IWRTMA(IOCC_MIN,1,NORB,1,NORB)
call IWRTMA(IOCC_MAX,1,NORB,1,NORB)
write(u6,*) ' Arcweights for single occupied arcs'
call IWRTMA(IARCW(:,:,1),NORB,NEL,NORB,NEL)
write(u6,*) ' Arcweights for double occupied arcs'
call IWRTMA(IARCW(:,:,2),NORB,NEL,NORB,NEL)
write(u6,*) ' Total number of configurations ',NCONF
#endif

end subroutine CONF_GRAPH
