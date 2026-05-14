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

subroutine SETSXCI()

use sxci, only: IDXCI, IDXSX
use gas_data, only: NGAS, NGSSH
use general_data, only: NSYM
use Molcas, only: MxGAS
use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: I, IGAS, IGSSH, IOFF_GSSH(mxgas), ISTOT, ISYM, NGSSHT

!---------------------------------------------------------
!--  SET INDEX VECTORS FOR CI/SX INTEGRAL ORDERING
!---------------------------------------------------------

NGSSHT = 0
do IGAS=1,NGAS
  IOFF_GSSH(IGAS) = NGSSHT
  NGSSHT = NGSSHT+sum(NGSSH(IGAS,1:NSYM))
end do
ISTOT = 0
do ISYM=1,NSYM
  do IGAS=1,NGAS
    do IGSSH=1,NGSSH(IGAS,ISYM)
      IOFF_GSSH(IGAS) = IOFF_GSSH(IGAS)+1
      ISTOT = ISTOT+1
      IDXCI(ISTOT) = IOFF_GSSH(IGAS)
    end do
  end do
end do

do I=1,ISTOT
  IDXSX(IDXCI(I)) = I
end do

#ifdef _DEBUGPRINT_
write(u6,'(1X,A)') 'REORDERING VECTOR FOR CI'
write(u6,'(1X,12I5)') (IDXCI(I),I=1,ISTOT)
write(u6,'(1X,A)') 'REORDERING VECTOR FOR SX'
write(u6,'(1X,12I5)') (IDXSX(I),I=1,ISTOT)
#endif

end subroutine SETSXCI
