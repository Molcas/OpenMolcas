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

implicit real*8(A-H,O-Z)
#include "motra_global.fh"
real*8 CMO(*)
dimension DSQ(*), DLT(*)

ISTSQ = 0
ISTLT = 0
do ISYM=1,NSYM
  NF = NFRO(ISYM)
  NB = NBAS(ISYM)
  if (NB*NF > 0) call DGEMM_('N','T',NB,NB,NF,1.0d0,CMO(ISTSQ+1),NB,CMO(ISTSQ+1),NB,0.0d0,DSQ(ISTSQ+1),NB)
  call DSCAL_(NB*NB,2.0d0,DSQ(ISTSQ+1),1)
  IJ = ISTLT
  do IB=1,NB
    do JB=1,IB
      IJ = IJ+1
      DLT(IJ) = 2.0d0*DSQ(ISTSQ+JB+(IB-1)*NB)
    end do
    DLT(IJ) = 0.5d0*DLT(IJ)
  end do
  ISTSQ = ISTSQ+NB*NB
  ISTLT = ISTLT+NB*(NB+1)/2
end do

if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(6,'(6X,A)') 'Frozen one-body density matrix in AO basis'
  ISTLT = 1
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(6,'(6X,A,I2)') 'symmetry species:',ISYM
      call TRIPRT(' ',' ',DLT(ISTLT),NB)
      ISTLT = ISTLT+NB*(NB+1)/2
    end if
  end do
end if

return

end subroutine DONEI
