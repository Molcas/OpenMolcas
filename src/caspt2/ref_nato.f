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
      use definitions, only: iwp, wp
      use constants, only: Zero, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM,NFRO,NISH,NASH,NBAS

      IMPLICIT None

      real(kind=wp), intent(In):: DREF(*),CMO(*)
      real(kind=wp), intent(out):: OCC(*),CNAT(*)

      real(kind=wp), ALLOCATABLE:: TMP(:)
      integer(kind=iwp) IDREF, IOCC, ICMO, ISYM, NF, NI, NA, NB, I, II,
     &                  J, JJ, LIJ, NFI, NSD, NTMP
      real(kind=wp) OC

* Purpose: compute natural orbitals and natural occupation numbers
* for the reference wave function.
* Given DREF, a triangular density matrix
* in active/active MO basis, and symmetry-blocked array CMO of MO
* coefficients, return array of  natural occupation numbers and MO
* coefficients of  natural orbitals. Frozen, inactive and virtual
C orbitals are copied unchanged.


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
        IF(NFI>0) THEN
          CALL DCOPY_(NFI,[Two],0,OCC(IOCC+1),1)
          IOCC=IOCC+NFI
          CALL DCOPY_(NB*NFI,CMO(ICMO+1),1,CNAT(ICMO+1),1)
          ICMO=ICMO+NB*NFI
        END IF
C Active orbitals:
        IF(NA>0) THEN
          NTMP=(NA*(NA+1))/2
          CALL mma_allocate(TMP,NTMP,LABEL='TMP')
          CALL DCOPY_(NB*NA,CMO(ICMO+1),1,CNAT(ICMO+1),1)
C For correct ordering, change sign.
          LIJ=1
          DO I=1,NA
           II=I+IDREF
           DO J=1,I
            JJ=J+IDREF
            TMP(LIJ)=-DREF((II*(II-1))/2+JJ)
            LIJ=LIJ+1
           END DO
          END DO
          CALL NIDiag(TMP,CNAT(ICMO+1),NA,NB)
          CALL JACORD(TMP,CNAT(ICMO+1),NA,NB)
          CALL VEIG(NA,TMP,OCC(IOCC+1))
          CALL mma_deallocate(TMP)
C Change back to positive sign.
          CALL DSCAL_(NA,-One,OCC(IOCC+1),1)
* Certain CAS or RAS wave functions can legitimately have
* occupation numbers that are exactly 0 or 2. These may become
* inappropriate by rounding. Fix that as well.
          DO I=1,NA
           OC=OCC(IOCC+I)
           IF(OC<Zero) OC=Zero
           IF(OC>Two) OC=Two
           OCC(IOCC+I)=OC
          END DO
          IDREF=IDREF+NA
          IOCC=IOCC+NA
          ICMO=ICMO+NB*NA
        END IF
C Secondary and deleted orbitals:
        NSD=NB-(NFI+NA)
        IF(NSD>0) THEN
          CALL DCOPY_(NSD,[Zero],0,OCC(IOCC+1),1)
          IOCC=IOCC+NSD
          CALL DCOPY_(NB*NSD,CMO(ICMO+1),1,CNAT(ICMO+1),1)
          ICMO=ICMO+NB*NSD
        END IF
      END DO

      END SUBROUTINE REF_NATO
