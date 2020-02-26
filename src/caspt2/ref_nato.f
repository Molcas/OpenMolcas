************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE REF_NATO(DREF,CMO,OCC,CNAT)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
      DIMENSION DREF(*),CMO(*),OCC(*),CNAT(*)
* Purpose: compute natural orbitals and natural occupation numbers
* for the reference wave function.
* Given DREF, a triangular density matrix
* in active/active MO basis, and symmetry-blocked array CMO of MO
* coefficients, return array of  natural occupation numbers and MO
* coefficients of  natural orbitals. Frozen, inactive and virtual
C orbitals are copied unchanged.

      CALL QENTER('REF_NATO')

      IDREF=0
      IOCC=0
      ICMO=0
      DO ISYM=1,NSYM
        NF=NFRO(ISYM)
        NI=NISH(ISYM)
        NA=NASH(ISYM)
        NB=NBAS(ISYM)
C Frozen and inactive orbitals:
        NFI=NF+NI
        IF(NFI.GT.0) THEN
          CALL DCOPY_(NFI,[2.0D00],0,OCC(IOCC+1),1)
          IOCC=IOCC+NFI
          CALL DCOPY_(NB*NFI,CMO(ICMO+1),1,CNAT(ICMO+1),1)
          ICMO=ICMO+NB*NFI
        END IF
C Active orbitals:
        IF(NA.GT.0) THEN
          NTMP=(NA*(NA+1))/2
          CALL GETMEM('TMP','ALLO','REAL',LTMP,NTMP)
          CALL DCOPY_(NB*NA,CMO(ICMO+1),1,CNAT(ICMO+1),1)
C For correct ordering, change sign.
          LIJ=LTMP
          DO I=1,NA
           II=I+IDREF
           DO J=1,I
            JJ=J+IDREF
            WORK(LIJ)=-DREF((II*(II-1))/2+JJ)
            LIJ=LIJ+1
           END DO
          END DO
          CALL NIDiag(WORK(LTMP),CNAT(ICMO+1),NA,NB,0)
          CALL JACORD(WORK(LTMP),CNAT(ICMO+1),NA,NB)
          CALL VEIG(NA,WORK(LTMP),OCC(IOCC+1))
          CALL GETMEM('TMP','FREE','REAL',LTMP,NTMP)
C Change back to positive sign.
          CALL DSCAL_(NA,-1.0D0,OCC(IOCC+1),1)
* Certain CAS or RAS wave functions can legitimately have
* occupation numbers that are exactly 0 or 2. These may become
* inappropriate by rounding. Fix that as well.
          DO I=1,NA
           OC=OCC(IOCC+I)
           IF(OC.LT.0.0D0) OC=0.0D0
           IF(OC.GT.2.0D0) OC=2.0D0
           OCC(IOCC+I)=OC
          END DO
          IDREF=IDREF+NA
          IOCC=IOCC+NA
          ICMO=ICMO+NB*NA
        END IF
C Secondary and deleted orbitals:
        NSD=NB-(NFI+NA)
        IF(NSD.GT.0) THEN
          CALL DCOPY_(NSD,[0.0D0],0,OCC(IOCC+1),1)
          IOCC=IOCC+NSD
          CALL DCOPY_(NB*NSD,CMO(ICMO+1),1,CNAT(ICMO+1),1)
          ICMO=ICMO+NB*NSD
        END IF
      END DO

      CALL QEXIT('REF_NATO')

      RETURN
      END
