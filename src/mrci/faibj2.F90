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

subroutine faibj2(IFTA,IFTB,ICOUP1,ICOUP,INDA,INDB,MYSYM,INTSYM,NYSYM,NSIJ,MYL,NYL,FACS,IPOA,IPOB,INMY,INNY,INDX,iTYP)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IFTA, IFTB, ICOUP1, ICOUP, INDA, INDB, MYSYM, INTSYM(*), NYSYM, NSIJ, MYL, NYL, IPOA(9), IPOB(9), INMY, INNY, &
                     INDX(*), iTYP
real(kind=wp) :: FACS
#include "mrci.fh"
integer(kind=iwp), external :: JSUNP

IFTA = 0
IFTB = 0
!goto (109,110,111,112,113),ITYP
if (ITYP == 1) then
  INDA = IRC(2)+ICOUP1
  INDB = IRC(2)+ICOUP
  IFTA = 1
  IFTB = 1
end if
if (ITYP == 2) then
  INDA = IRC(3)+ICOUP1
  INDB = IRC(3)+ICOUP
end if
if (ITYP == 3) then
  INDA = IRC(2)+ICOUP1
  INDB = IRC(3)+ICOUP
  IFTA = 1
end if
if (ITYP == 4) then
  INDA = IRC(3)+ICOUP1
  INDB = IRC(2)+ICOUP
  IFTB = 1
end if
if (ITYP == 5) then
  INDA = IRC(1)+ICOUP1
  INDB = IRC(1)+ICOUP
end if
!vv : unroll inline function to make GCC compiler works proper..
MYSYM = JSUNP(INTSYM,INDA)
NYSYM = MUL(MYSYM,NSIJ)
MYL = MUL(MYSYM,LSYM)
NYL = MUL(NYSYM,LSYM)
FACS = One
call IPO(IPOA,NVIR,MUL,NSYM,MYL,IFTA)
call IPO(IPOB,NVIR,MUL,NSYM,NYL,IFTB)
INMY = INDX(INDA)+1
INNY = INDX(INDB)+1

return

end subroutine faibj2
