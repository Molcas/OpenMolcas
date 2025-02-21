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
subroutine WEIGHT_SPGP(Z,NORBTP,NELFTP,NORBFTP,ISCR)
! construct vertex weights for given supergroup
!
! Reverse lexical ordering is used

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(_OUT_) :: Z(*), ISCR(*)
integer(kind=iwp), intent(in) :: NORBTP, NELFTP(NORBTP), NORBFTP(NORBTP)
integer(kind=iwp) :: KLFREE, KLMAX, KLMIN, KW, NEL, NORB

NORB = sum(NORBFTP)
NEL = sum(NELFTP)

#ifdef _DEBUGPRINT_
write(u6,*) ' Subroutine WEIGHT_SPGP in action'
write(u6,*) ' ================================'
write(u6,*) 'NELFTP'
call IWRTMA(NELFTP,1,NORBTP,1,NORBTP)
#endif

KLFREE = 1
KLMAX = KLFREE
KLFREE = KLFREE+NORB

KLMIN = KLFREE
KLFREE = KLFREE+NORB

KW = KLFREE
KLFREE = KW+(NEL+1)*(NORB+1)
! Max and min arrays for strings
call MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),NORBTP,NORBFTP,NELFTP)
! Arc weights
call GRAPW(ISCR(KW),Z,ISCR(KLMIN),ISCR(KLMAX),NORB,NEL)

end subroutine WEIGHT_SPGP
