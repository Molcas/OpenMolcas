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

!#define _DEBUGPRINT_
subroutine WEIGHT_LUCIA(Z,NEL,NORB1,NORB2,NORB3,MNRS1,MXRS1,MNRS3,MXRS3,ISCR)
! construct vertex weights
!
! Reverse lexical ordering is used for restricted space

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(out) :: Z(*), ISCR(*)
integer(kind=iwp), intent(in) :: NEL, NORB1, NORB2, NORB3, MNRS1, MXRS1, MNRS3, MXRS3
integer(kind=iwp) :: KLFREE, KLMAX, KLMIN, KW, NORB

NORB = NORB1+NORB2+NORB3

#ifdef _DEBUGPRINT_
write(u6,*) ' >>>> WEIGHT <<<<<'
write(u6,*) ' NORB1 NORB2 NORB3 ',NORB1,NORB2,NORB3
write(u6,*) ' NEL MNRS1 MXRS1 MNRS3 MXRS3'
write(u6,*) NEL,MNRS1,MXRS1,MNRS3,MXRS3
#endif

KLFREE = 1
KLMAX = KLFREE
KLFREE = KLFREE+NORB

KLMIN = KLFREE
KLFREE = KLFREE+NORB

KW = KLFREE
!KLFREE = KW+(NEL+1)*(NORB+1)
! Max and min arrays for strings
call RSMXMN(ISCR(KLMAX),ISCR(KLMIN),NORB1,NORB2,NORB3,NEL,MNRS1,MXRS1,MNRS3,MXRS3)
! Arc weights
call GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL)

end subroutine WEIGHT_LUCIA
