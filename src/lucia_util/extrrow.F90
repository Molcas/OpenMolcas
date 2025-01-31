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

subroutine EXTRROW(INMAT,IROW,NROW,NCOL,IOUTVEC)
! Extract row IROW from integer matrix INMAT
!
! Jeppe Olsen, Winter 1996

implicit real*8(A-H,O-Z)
dimension INMAT(NROW,NCOL)
dimension IOUTVEC(NCOL)

do ICOL=1,NCOL
  IOUTVEC(ICOL) = INMAT(IROW,ICOL)
end do

NTEST = 0
if (NTEST >= 100) then
  write(6,*) ' Output vector from EXTRROW'
  write(6,*) ' Extracted ROW ',IROW
  call IWRTMA(IOUTVEC,1,NCOL,1,NCOL)
end if

end subroutine EXTRROW
