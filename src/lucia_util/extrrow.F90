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
! Copyright (C) 1996, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine EXTRROW(INMAT,IROW,NROW,NCOL,IOUTVEC)
! Extract row IROW from integer matrix INMAT
!
! Jeppe Olsen, Winter 1996

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NROW, NCOL, INMAT(NROW,NCOL), IROW
integer(kind=iwp), intent(out) :: IOUTVEC(NCOL)

IOUTVEC(:) = INMAT(IROW,:)

#ifdef _DEBUGPRINT_
write(u6,*) ' Output vector from EXTRROW'
write(u6,*) ' Extracted ROW ',IROW
call IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
#endif

end subroutine EXTRROW
