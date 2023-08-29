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

subroutine CHO_MCA_DIAGINT(ISHLA,ISHLB,SCR,LSCR)
!
! Purpose: call Seward to calculate diagonal shell (AB|AB).

#ifdef _DEBUGPRINT_
use Cholesky, only: CutInt, ThrInt
#endif
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ISHLA, ISHLB, LSCR
real(kind=wp), intent(out) :: SCR(LSCR)
#include "itmax.fh"
#ifdef _DEBUGPRINT_
character(len=*), parameter :: SECNAM = 'CHO_MCA_DIAGINT'
#endif
external :: Integral_WrOut_Cho_diag

SCR(:) = Zero

#ifdef _DEBUGPRINT_
CUTINT1 = CutInt
THRINT1 = ThrInt
#endif

call EVAL_IJKL(ISHLA,ISHLB,ISHLA,ISHLB,SCR,LSCR,Integral_WrOut_Cho_diag)

#ifdef _DEBUGPRINT_
CUTINT2 = CutInt
THRINT2 = ThrInt
if ((CUTINT2 /= CUTINT1) .or. (THRINT2 /= THRINT1)) then
  write(LUPRI,*) SECNAM,': CutInt before Eval_Ints_: ',CUTINT1
  write(LUPRI,*) SECNAM,': CutInt after  Eval_Ints_: ',CUTINT2
  write(LUPRI,*) SECNAM,': ThrInt before Eval_Ints_: ',THRINT1
  write(LUPRI,*) SECNAM,': ThrInt after  Eval_Ints_: ',THRINT2
  call CHO_QUIT('Integral prescreening error detected in '//SECNAM,102)
end if
#endif

end subroutine CHO_MCA_DIAGINT
