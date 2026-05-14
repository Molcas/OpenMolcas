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

subroutine CHO_INTCHK_REG(LABEL,ISHLCD,ISHLAB)
!
! Purpose: register a shell quadruple (CD|AB) for minimal integral
!          check using LABEL to keep track of its origin.

use Cholesky, only: iChkq, iSP2F, nChkq, nnShl
use Definitions, only: iwp

implicit none
character(len=8), intent(in) :: LABEL
integer(kind=iwp), intent(in) :: ISHLCD, ISHLAB
integer(kind=iwp) :: ID, ISHLA, ISHLB, ISHLC, ISHLD
character(len=8) :: LAB
character(len=*), parameter :: SECNAM = 'CHO_INTCHK_REG'

! Check shell pair index.
! -----------------------

if ((ISHLCD < 1) .or. (ISHLCD > NNSHL)) call CHO_QUIT('Shell index error 1 in '//SECNAM,103)
if ((ISHLAB < 1) .or. (ISHLAB > NNSHL)) call CHO_QUIT('Shell index error 2 in '//SECNAM,103)

! Registration.
! -------------

call CHO_INVPCK(ISP2F(ISHLCD),ISHLC,ISHLD,.true.)
call CHO_INVPCK(ISP2F(ISHLAB),ISHLA,ISHLB,.true.)
LAB = LABEL
call CHO_INTCHK_ID_OF(LAB,ID,1)
if ((ID < 1) .or. (ID > NCHKQ)) ID = NCHKQ+1 ! junk yard
ICHKQ(1,ID) = ISHLC
ICHKQ(2,ID) = ISHLD
ICHKQ(3,ID) = ISHLA
ICHKQ(4,ID) = ISHLB

end subroutine CHO_INTCHK_REG
