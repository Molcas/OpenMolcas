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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine REO_PTDET(NOPEN,NALPHA,IZ_PTDET,IREO_PTDET,ILIST_PTDET,NLIST_PTDET,ISCR)
! A list(ILIST_PTDET) of prototype determinants with NOPEN unpaired electrons and
! NALPHA alpha electrons is given.
!
! Obtain 1 : Z matrix for this set of prototype dets
!        2 : Reorder array going from lexical order to
!            the order specified by the prototype dets
!            given in ILIST_PTDET. The reordering goes
!            from lexical order to actual order.
!            Prototype determinants not included in
!            ILIST_PTDET are given zero address
!
! ISCR : Min length : 2*NOPEN + (NALPHA+1)*(NOPEN+1)
!
! Jeppe Olsen, December 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NOPEN, NALPHA, NLIST_PTDET, ILIST_PTDET(NOPEN,NLIST_PTDET)
integer(kind=iwp), intent(out) :: IZ_PTDET(NOPEN,NALPHA)
integer(kind=iwp), intent(_OUT_) :: IREO_PTDET(*), ISCR(*)
integer(kind=iwp) :: IDUM(1), ILEX, JPTDT, KLMAX, KLMIN, KLW, NTOT_PTDET
integer(kind=iwp), external :: IBINOM, IZNUM_PTDT

! 1 : Set up lexical order array for prototype determinants
!     (alpha considered as occupied electron)
KLMIN = 1
KLMAX = KLMIN+NOPEN
KLW = KLMAX+NOPEN
call MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),1,[NOPEN],[NALPHA])

! Arc weights

call GRAPW(ISCR(KLW),IZ_PTDET,ISCR(KLMIN),ISCR(KLMAX),NOPEN,NALPHA)

! Reorder array

! Total number of prototype determinants
NTOT_PTDET = 0
if ((NALPHA >= 0) .and. (NALPHA <= NOPEN)) then
  NTOT_PTDET = IBINOM(NOPEN,NALPHA)
else
  NTOT_PTDET = 0
end if
IREO_PTDET(1:NTOT_PTDET) = 0

do JPTDT=1,NLIST_PTDET
  !write(u6,*) ' JPTDT = ',JPTDT
  ! Lexical address of prototype determiant JPTDT
  if (NALPHA == 0) then
    ILEX = 1
  else
    ILEX = IZNUM_PTDT(ILIST_PTDET(:,JPTDT),NOPEN,NALPHA,IZ_PTDET,IDUM,0)
  end if
  IREO_PTDET(ILEX) = JPTDT
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Reorder array for prototype determinants'
call IWRTMA(IREO_PTDET,1,NTOT_PTDET,1,NTOT_PTDET)
#endif

end subroutine REO_PTDET
