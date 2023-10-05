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

subroutine CHO_MCA_CALCINT_4(XINT,LINT,ISHLCD,ISHLAB)
!
! Purpose: calculate qualified integral columns from
!          shell quadruple (ISHLC ISHLD|ISHLA ISHLB).
!
! Version 4: avoid storage of full shell quadruple in interface to
!            seward; get qualified directly as in Version 2 and 3!
!            Changes from Version 3:
!            - only one shell quadruple is computed (not an entire
!              column).

use Cholesky, only: INF_IN2, IPRINT, iSP2F, LuPri, NCOLAB, nSym, TINTEG
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LINT, ISHLCD, ISHLAB
real(kind=wp), intent(inout) :: XINT(LINT)
integer(kind=iwp) :: i, ILOC, IRC, ISHLA, ISHLB, ISHLC, ISHLD, NAB(8)
real(kind=wp) :: C1, C2, W1, W2
logical(kind=iwp), parameter :: LOCDBG = .false.
character(len=*), parameter :: SECNAM = 'CHO_MCA_CALCINT_4'

! Set mapping from shell pair AB to qualified columns.
! ----------------------------------------------------

IRC = 0
ILOC = 2
call CHO_SETSHP2Q_2(IRC,ILOC,ISHLAB,NAB)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_SETSHP2Q_2 returned ',IRC
  call CHO_QUIT('Error termination in '//SECNAM,IRC)
end if

! Print.
! ------

if (IPRINT >= INF_IN2) then
  call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
  call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
  NCOLAB = sum(NAB(1:NSYM))
  write(LUPRI,'(/,A,I5,1X,I5,A,I5,1X,I5,A,I9,A)') 'Calculating shell quadruple (',ISHLC,ISHLD,'|',ISHLA,ISHLB,'):',NCOLAB, &
                                                  ' columns have been qualified'
  write(LUPRI,'(89A)') ('=',i=1,89)
end if

! Set mapping from shell pair CD to reduced set.
! ----------------------------------------------

IRC = 0
ILOC = 2
call CHO_SETSHP2RS_2(IRC,ILOC,ISHLCD,NAB)
if (IRC /= 0) then
  write(LUPRI,*) SECNAM,': CHO_SETSHP2RS_2 returned ',IRC
  call CHO_QUIT('Error termination in '//SECNAM,IRC)
end if

! Calculate integrals.
! --------------------

call CWTIME(C1,W1)
call CHO_MCA_INT_1(ISHLCD,ISHLAB,XINT,LINT,LOCDBG .or. (IPRINT >= 100))
call CWTIME(C2,W2)
TINTEG(1,1) = TINTEG(1,1)+C2-C1
TINTEG(2,1) = TINTEG(2,1)+W2-W1

end subroutine CHO_MCA_CALCINT_4
