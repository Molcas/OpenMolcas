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

subroutine CHO_CHKINT(XINT,DIAG,ISYM,NERR,TOL,REPORT)
!
! Purpose: check diagonals in qualified integral columns.

use Cholesky, only: iiBstR, IndRed, iQuAB, LuPri, nnBstR, nQual
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: XINT(*), DIAG(*), TOL
integer(kind=iwp), intent(in) :: ISYM
integer(kind=iwp), intent(out) :: NERR
logical(kind=iwp), intent(in) :: REPORT
integer(kind=iwp) :: I, II, IK, JJ, KK
real(kind=wp) :: DF
character(len=*), parameter :: SECNAM = 'CHO_CHKINT'

NERR = 0
do I=1,NQUAL(ISYM)
  II = IQUAB(I,ISYM)
  JJ = INDRED(II,2)
  IK = II-IIBSTR(ISYM,2)
  KK = NNBSTR(ISYM,2)*(I-1)+IK
  DF = DIAG(JJ)-XINT(KK)
  if (abs(DF) > TOL) then
    NERR = NERR+1
    if (REPORT) then
      write(LUPRI,*) SECNAM,': diag error: ',DIAG(JJ),XINT(KK)
      write(LUPRI,*) '            diagonal elm    : ',JJ,' (rs1) ',II,' (rs2)'
      write(LUPRI,*) '            integral row,col: ',IK,I
    end if
  end if
end do

end subroutine CHO_CHKINT
