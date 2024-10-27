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
      SUBROUTINE HEFVAL(IST,JST,DVALUE)
      use caspt2_global, only: LUCIEX, IDTCEX
      use EQSOLV
      use stdalloc, only: mma_allocate, mma_deallocate
      IMPLICIT NONE
C Apart from input call parameters, we need two vectors stored on
C LUSOLV. Vector nr IVECC (presently=2) contains the contravariant
C elements of the solution to the CASPT2 equations.
C IVECW is the number (presently=6) of the vector on LUSOLV
C where a contravariant representation of the CASPT2 Right-Hand Side
C vector is stored. This depends on the MOs used, but is actually
C the same for all the root states.

#include "caspt2.fh"
#include "pt2_guga.fh"

      INTEGER IST,JST
      REAL*8 DVALUE

      INTEGER I
      INTEGER NTG1,NTG2,NTG3
      INTEGER IDCI
      REAL*8 OVL,DUMMY(1)
      REAL*8, ALLOCATABLE:: TG1(:), TG2(:), TG3(:), CI1(:), CI2(:)

C We evaluate the effective Hamiltonian matrix element in two steps.

      NTG1=NASHT**2
      NTG2=NASHT**4
      NTG3=(NTG1*(NTG1+1)*(NTG1+2))/6
C Note: Need proper allocation even if unused, sinced allocated
C arrays are in subroutine parameter lists of MKTG3, HCOUP.
      NTG1=MAX(1,NTG1)
      NTG2=MAX(1,NTG2)
      NTG3=MAX(1,NTG3)
      CALL mma_allocate(TG1,NTG1,Label='TG1')
      CALL mma_allocate(TG2,NTG2,Label='TG2')
      CALL mma_allocate(TG3,NTG3,Label='TG3')
      TG1(1)=0.0D0
      TG2(1)=0.0D0
      TG3(1)=0.0D0

      CALL mma_allocate(CI1,MXCI,Label='CI1')
      CALL mma_allocate(CI2,MXCI,Label='CI2')
      IF(ISCF.EQ.0) THEN
C Read root vectors nr. IST and JST from LUCI.
        IDCI=IDTCEX
        DO I=1,NSTATE
          IF(I.EQ.IST) THEN
            CALL DDAFILE(LUCIEX,2,CI1,NCONF,IDCI)
            IF(I.EQ.JST) THEN
              CALL DCOPY_(NCONF,CI1,1,CI2,1)
            END IF
          ELSE IF(I.EQ.JST) THEN
            CALL DDAFILE(LUCIEX,2,CI2,NCONF,IDCI)
          ELSE
            CALL DDAFILE(LUCIEX,0,DUMMY,NCONF,IDCI)
          END IF
        END DO
      END IF

      CALL MKTG3(STSYM,STSYM,CI1,CI2,OVL,TG1,TG2,NTG3,TG3)
      CALL mma_deallocate(CI1)
      CALL mma_deallocate(CI2)

      CALL HCOUP(IVECW,IVECC,OVL,TG1,TG2,TG3,DVALUE)

      CALL mma_deallocate(TG1)
      CALL mma_deallocate(TG2)
      CALL mma_deallocate(TG3)

      END SUBROUTINE HEFVAL
