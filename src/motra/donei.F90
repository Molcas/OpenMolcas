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

subroutine DONEI(DLT,DSQ,CMO)

#include "intent.fh"

use motra_global, only: Debug, iPrint, nBas, nFro, nSym
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(_OUT_) :: DLT(*), DSQ(*)
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp) :: IB, IJ, ISTLT, ISTSQ, ISYM, JB, NB, NF

ISTSQ = 0
ISTLT = 0
do ISYM=1,NSYM
  NF = NFRO(ISYM)
  NB = NBAS(ISYM)
  if (NB*NF > 0) call DGEMM_('N','T',NB,NB,NF,One,CMO(ISTSQ+1),NB,CMO(ISTSQ+1),NB,Zero,DSQ(ISTSQ+1),NB)
  call DSCAL_(NB*NB,Two,DSQ(ISTSQ+1),1)
  IJ = ISTLT
  do IB=1,NB
    do JB=1,IB
      IJ = IJ+1
      DLT(IJ) = Two*DSQ(ISTSQ+JB+(IB-1)*NB)
    end do
    DLT(IJ) = Half*DLT(IJ)
  end do
  ISTSQ = ISTSQ+NB*NB
  ISTLT = ISTLT+NB*(NB+1)/2
end do

if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'Frozen one-body density matrix in AO basis'
  ISTLT = 1
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(u6,'(6X,A,I2)') 'symmetry species:',ISYM
      call TRIPRT(' ',' ',DLT(ISTLT),NB)
      ISTLT = ISTLT+NB*(NB+1)/2
    end if
  end do
end if

return

end subroutine DONEI
