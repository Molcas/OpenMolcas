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
! Jeppe Olsen, December 2001

use Definitions, only: u6

implicit real*8(A-H,O-Z)
! Input
integer ILIST_PTDET(NOPEN,*)
! Output
integer IZ_PTDET(NOPEN,NALPHA), IREO_PTDET(*)
! Local scratch : Min length : 2*NOPEN + (NALPHA+1)*(NOPEN+1)
integer ISCR(*)
integer NOPEN_ARR(1), NALPHA_ARR(1), IDUM(1)
NTEST = 0

! 1 : Set up lexical order array for prototype determinants
!     (alpha considered as occupied electron)
KLMIN = 1
KLMAX = KLMIN+NOPEN
KLW = KLMAX+NOPEN
NOPEN_ARR(1) = NOPEN
NALPHA_ARR(1) = NALPHA
call MXMNOC_SPGP(ISCR(KLMIN),ISCR(KLMAX),1,NOPEN_ARR,NALPHA_ARR,NTEST)
NOPEN = NOPEN_ARR(1)
NALPHA = NALPHA_ARR(1)

! Arc weights

call GRAPW(ISCR(KLW),IZ_PTDET,ISCR(KLMIN),ISCR(KLMAX),NOPEN,NALPHA,NTEST)
NOPEN_ARR(1) = NOPEN
NALPHA_ARR(1) = NALPHA

! Reorder array

! Total number of prototype determinants
NTOT_PTDET = 0
if ((NALPHA >= 0) .and. (NALPHA <= NOPEN)) then
  NTOT_PTDET = IBION_LUCIA(NOPEN,NALPHA)
else
  NTOT_PTDET = 0
end if
IZERO = 0
call ISETVC(IREO_PTDET,IZERO,NTOT_PTDET)

do JPTDT=1,NLIST_PTDET
  !write(u6,*) ' JPTDT = ',JPTDT
  ! Lexical address of prototype determiant JPTDT
  if (NALPHA == 0) then
    ILEX = 1
  else
    ILEX = IZNUM_PTDT(ILIST_PTDET(1,JPTDT),NOPEN,NALPHA,IZ_PTDET,IDUM,0)
  end if
  IREO_PTDET(ILEX) = JPTDT
end do

if (NTEST >= 100) then
  write(u6,*) ' Reorder array for prototype determinants'
  call IWRTMA(IREO_PTDET,1,NTOT_PTDET,1,NTOT_PTDET)
end if

end subroutine REO_PTDET
