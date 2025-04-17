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
! Copyright (C) 1989, Jeppe Olsen                                      *
!***********************************************************************

subroutine GETCNF_MCLR(KCNF,KTYP,K,ICONF,IREFSM,NEL)
! Obtain configuration number K.
! Occupation in KCNF
! Type in KTYP
!
! Jeppe Olsen, summer of 89

use MCLR_Data, only: MINOP, NCNATS, NTYP
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: KCNF(*)
integer(kind=iwp), intent(out) :: KTYP
integer(kind=iwp), intent(in) :: K, ICONF(*), IREFSM, NEL
integer(kind=iwp) :: ICNFB1, ICNFB2, JCL, JOCC, JOP, JTYP, KADD, KREL, NJCNF

ICNFB1 = 1
ICNFB2 = 1
KTYP = 0
do JTYP=1,NTYP
  JOP = JTYP-1+MINOP
  JCL = (NEL-JOP)/2
  JOCC = JOP+JCL

  NJCNF = NCNATS(JTYP,IREFSM)
  if ((K >= ICNFB1) .and. (K <= ICNFB1+NJCNF-1)) then
    KREL = K-ICNFB1+1
    KADD = (KREL-1)*JOCC
    KTYP = JTYP
    KCNF(1:JOCC) = ICONF(ICNFB2+KADD:ICNFB2+KADD+JOCC-1)
  end if

  ICNFB1 = ICNFB1+NJCNF
  ICNFB2 = ICNFB2+NJCNF*JOCC
end do

end subroutine GETCNF_MCLR
