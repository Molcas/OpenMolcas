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

!#define _DEBUGPRINT_
subroutine SIGVST(ISGVST,NSMST)
! Obtain ISGVST(ISM) : Symmetry of sigma v on string of symmetry ism

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NSMST
integer(kind=iwp), intent(out) :: ISGVST(NSMST)
integer(kind=iwp) :: ISM

ISGVST(:) = [(-ISM+2,ISM=1,NSMST)]

#ifdef _DEBUGPRINT_
write(u6,*) ' ISGVST array'
write(u6,*) ' ============'
call IWRTMA(ISGVST,1,NSMST,1,NSMST)
#endif

end subroutine SIGVST
