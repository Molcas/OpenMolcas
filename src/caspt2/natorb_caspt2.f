************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1994, Per Ake Malmqvist                                *
************************************************************************
*--------------------------------------------*
* 1994  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE NATORB_CASPT2(DMAT,CMO,OCC,CNAT)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
      DIMENSION DMAT(*),CMO(*),OCC(*),CNAT(*)
C Given DMAT, symmetry-blocked array of triangular
C density matrices in MO basis, and symmetry-blocked
C array CMO of MO coefficients, return array of
C natural occupation numbers and MO coefficients of
C natural orbitals.

      CALL QENTER('NATORB')

      IDMAT=0
      IOCC=0
      ICMO=0
      DO ISYM=1,NSYM
        NF=NFRO(ISYM)
        NO=NORB(ISYM)
        ND=NDEL(ISYM)
        NB=NBAS(ISYM)
C  Frozen orbitals:
        IF(NF.GT.0) THEN
          CALL DCOPY_(NF,[2.0D00],0,OCC(IOCC+1),1)
          IOCC=IOCC+NF
          CALL DCOPY_(NB*NF,CMO(ICMO+1),1,CNAT(ICMO+1),1)
          ICMO=ICMO+NB*NF
        END IF
C Inactive, active, and secondary orbitals:
        IF(NO.GT.0) THEN
          NTMP=(NO*(NO+1))/2
          CALL GETMEM('TMP','ALLO','REAL',LTMP,NTMP)
          CALL DCOPY_(NB*NO,CMO(ICMO+1),1,CNAT(ICMO+1),1)
C For correct order, change sign.
          CALL DYAX(NTMP,-1.0D0,DMAT(IDMAT+1),1,WORK(LTMP),1)
          CALL NIDiag(WORK(LTMP),CNAT(ICMO+1),NO,NB,0)
          CALL JACORD(WORK(LTMP),CNAT(ICMO+1),NO,NB)
          CALL VEIG(NO,WORK(LTMP),OCC(IOCC+1))
C Change back to positive sign.
          CALL DSCAL_(NO,-1.0D0,OCC(IOCC+1),1)
          IDMAT=IDMAT+NTMP
          IOCC=IOCC+NO
          ICMO=ICMO+NB*NO
          CALL GETMEM('TMP','FREE','REAL',LTMP,NTMP)
        END IF
C Deleted orbitals:
        IF(ND.GT.0) THEN
          CALL DCOPY_(ND,[0.0D0],0,OCC(IOCC+1),1)
          IOCC=IOCC+ND
          CALL DCOPY_(NB*ND,CMO(ICMO+1),1,CNAT(ICMO+1),1)
          ICMO=ICMO+NB*ND
        END IF
      END DO

      CALL QEXIT('NATORB')

      RETURN
      END
