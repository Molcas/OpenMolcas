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
! Copyright (C) 1989, Per Ake Malmqvist                                *
!***********************************************************************
!****************************************************************
!  PROGRAM RASSI        PER-AAKE MALMQVIST
!  SUBROUTINE WRMAT     IBM-3090 RELEASE 89 01 31
!  PRINT OUT A SYMMETRY-BLOCKED MATRIX. THE COMBINED SYMMETRY
!  OF ROWS AND COLUMNS IS ISY12. NDIM1 GIVES THE NUMBER OF
!  OF ROWS WITHIN EACH SYMMETRY TYPE, SIMILAR NDIM2, COLUMNS.
!****************************************************************

subroutine WRMAT(TEXT,ISY12,NDIM1,NDIM2,NMAT,XMAT)

use Symmetry_Info, only: MUL, nIrrep
use Definitions, only: wp, iwp, u6

implicit none
character(len=*) :: TEXT
integer(kind=iwp) :: ISY12, NDIM1(nIrrep), NDIM2(nIrrep), NMAT
real(kind=wp) :: XMAT(NMAT)
integer(kind=iwp) :: ISTA, ISY1, ISY2, NN

ISTA = 1
write(u6,'(/,1X,A,/)') TEXT
do ISY1=1,nIrrep
  ISY2 = MUL(ISY1,ISY12)
  NN = NDIM1(ISY1)*NDIM2(ISY2)
  if (NN /= 0) then
    write(u6,*)
    write(u6,'(A,2I2)') ' SYMMETRY LABELS OF ROWS/COLS:',ISY1,ISY2
    call WRMAT1(NDIM1(ISY1),NDIM2(ISY2),XMAT(ISTA))
  end if
  ISTA = ISTA+NN
end do
write(u6,*)
write(u6,*) repeat('*',80)

end subroutine WRMAT
