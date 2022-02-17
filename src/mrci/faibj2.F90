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

use mrci_global, only: IRC, LSYM, NSYM, NVIR
use Symmetry_Info, only: Mul
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: IFTA, IFTB, INDA, INDB, MYSYM, NYSYM, MYL, NYL, IPOA(9), IPOB(9), INMY, INNY
integer(kind=iwp), intent(in) :: ICOUP1, ICOUP, INTSYM(*), NSIJ, INDX(*), iTYP
real(kind=wp) :: FACS
integer(kind=iwp), external :: JSUNP

IFTA = 0
IFTB = 0
select case (ITYP)
  case (1)
    INDA = IRC(2)+ICOUP1
    INDB = IRC(2)+ICOUP
    IFTA = 1
    IFTB = 1
  case (2)
    INDA = IRC(3)+ICOUP1
    INDB = IRC(3)+ICOUP
  case (3)
    INDA = IRC(2)+ICOUP1
    INDB = IRC(3)+ICOUP
    IFTA = 1
  case (4)
    INDA = IRC(3)+ICOUP1
    INDB = IRC(2)+ICOUP
    IFTB = 1
  case (5)
    INDA = IRC(1)+ICOUP1
    INDB = IRC(1)+ICOUP
end select
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
