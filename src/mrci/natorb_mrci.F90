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

subroutine NATORB_MRCI(CMO,DMO,CNO,OCC,SCR)

implicit real*8(A-H,O-Z)
dimension CMO(NCMO), DMO(NBTRI), CNO(NCMO), OCC(NBAST)
dimension SCR((NBMAX*(NBMAX+1))/2)
#include "SysDef.fh"
#include "mrci.fh"

call DCOPY_(NBAST,[0.0d00],0,OCC,1)
call DCOPY_(NCMO,CMO,1,CNO,1)
! LOOP OVER SYMMETRY LABELS
! Present index of end of processed CMO block:
IECMO = 0
! Present end of index to MO-s, to be translated by ICH array:
IEO = 0
! Present end of index to basis functions:
IEB = 0
do ISYM=1,NSYM
  NO = NORB(ISYM)
  NIAV = NISH(ISYM)+NASH(ISYM)+NVIR(ISYM)
  NB = NBAS(ISYM)
  if (NB == 0) Go To 100
  ISB = IEB+1
  IEB = IEB+NB
  ! ORBITALS PRE-FROZEN IN MOTRA, OR FROZEN IN MRCI:
  NF = NFMO(ISYM)+NFRO(ISYM)
  NBF = NB*NF
  ISCMO = IECMO+1
  IECMO = IECMO+NBF
  ! (DO NOTHING WITH THE FROZEN ORBITALS)
  if (NF > 0) call DCOPY_(NF,[2.0d00],0,OCC(ISB),1)
  IEO = IEO+NFRO(ISYM)
  ! ORBITALS EXPLICITLY USED IN CI:
  NBO = NB*NIAV
  ISCMO = IECMO+1
  IECMO = IECMO+NBO
  ISO = IEO+1
  IEO = IEO+NIAV
  ! TRANSFER SYMMETRY BLOCK OF DMO TO TRIANGULAR SCRATCH MATRIX:
  I12 = 0
  do I=ISO,IEO
    IO1 = ICH(I)
    do J=ISO,I
      IO2 = ICH(J)
      IO12 = (IO1*(IO1-1))/2+IO2
      if (IO1 < IO2) IO12 = (IO2*(IO2-1))/2+IO1
      I12 = I12+1
      SCR(I12) = DMO(IO12)
    end do
  end do
  ! DIAGONALIZE AND TRANSFORM ORBITALS:
  call JACOB(SCR,CNO(ISCMO),NIAV,NB)
  ! PICK OCCUP NR FROM DIAGONAL:
  II = 0
  do I=1,NIAV
    II = II+I
    OCC(ISB+NF-1+I) = SCR(II)
  end do
  ! ORDER BY DECREASING NATURAL OCCUPANCY:
  NN = NO-NFRO(ISYM)
  do I=1,NN-1
    OMAX = OCC(ISB+NF-1+I)
    IMAX = I
    do J=I+1,NN
      OC = OCC(ISB+NF-1+J)
      if (OMAX >= OC) goto 30
      IMAX = J
      OMAX = OC
30    continue
    end do
    if (IMAX == I) goto 40
    OCC(ISB+NF-1+IMAX) = OCC(ISB+NF-1+I)
    OCC(ISB+NF-1+I) = OMAX
    ISTA1 = ISCMO+NB*(I-1)
    ISTA2 = ISCMO+NB*(IMAX-1)
    call DCOPY_(NB,CNO(ISTA1),1,SCR,1)
    call DCOPY_(NB,CNO(ISTA2),1,CNO(ISTA1),1)
    call DCOPY_(NB,SCR,1,CNO(ISTA2),1)
40  continue
  end do
  ! ORBITALS PRE-DELETED IN MOTRA OR DELETED IN MRCI:
  ND = NDMO(ISYM)+NDEL(ISYM)
  NBD = NB*ND
  ISCMO = IECMO+1
  IECMO = IECMO+NBD
  IEO = IEO+NDEL(ISYM)
  ! (DO NOTHING WITH THE DELETED ORBITALS)
100 continue
end do

return

end subroutine NATORB_MRCI
