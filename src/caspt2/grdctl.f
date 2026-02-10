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
* Copyright (C) 1998,2006, Per Ake Malmqvist                           *
*               2019, Stefano Battaglia                                *
************************************************************************
*--------------------------------------------*
* 1998, 2006, Per Ake Malmqvist              *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
* 2006 update: Use RAS1..RAS3
*--------------------------------------------*
      SUBROUTINE GRDCTL(HEFF)
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      use caspt2_global, only: TAT, TORB
      use caspt2_global, only: LUCIEX, IDTCEX
      use EQSOLV, only: iVecC, iVecW
      use stdalloc, only: mma_allocate,mma_deallocate
      use caspt2_module, only: NSTATE,E2CORR,E2TOT,IFMSCOUP,ISCF,JSTATE,
     &                         NCONF,NSYM,ORBIN,STSYM,NISH,NRAS1,NRAS2,
     &                         NRAS3,NSSH,NAES
      IMPLICIT None
      real(kind=wp), intent(inout):: HEFF(NSTATE,NSTATE)

      real(kind=wp), Allocatable :: CI(:),SGM(:),BRACI(:)
      integer(kind=iwp) ID,I,IOFF1,IOFF2,ISYM,NI,NR1,NR2,NR3,NS,J,IJ,JI,
     &                  ISTART,ISTATE,ITO,ITOEND,ITOSTA,NSG
      real(kind=wp), external:: DDot_

C Purpose: Compute three sets of quantities, needed by MCLR, used to
C compute forces and derivatives.


C--------------------------------------------------------------------
C (A):
C The active/active density matrix contribution <Psi_1| Etu |Psi_1>
C and all the rest was already computed in the PRPCTL section.

C--------------------------------------------------------------------
C (B):
C Compute also (Proj_CAS)*(Ham)*(Wave op) | Psi_0 >
C and (Proj_CAS)*((Wave op)**(dagger))*(Ham) | Psi_0 >.
C Read the Psi_0 wave function:
      call mma_allocate(CI,NCONF,Label='GRDCI')
      IF(ISCF.EQ.0) THEN
C This is an ordinary CASSCF or RASSCF calculation.
       ID=IDTCEX
       DO I=1,JSTATE-1
         CALL DDAFILE(LUCIEX,0,CI,NCONF,ID)
       END DO
       CALL DDAFILE(LUCIEX,2,CI,NCONF,ID)
      ELSE
       CI(1)=One
      END IF

      IF(ORBIN.EQ.'TRANSFOR') THEN
C Read, and transpose, the active orbital transformation matrices
        TAT(:)=Zero
        IOFF1=0
        IOFF2=0
        DO ISYM=1,NSYM
          NI=NISH(ISYM)
          NR1=NRAS1(ISYM)
          NR2=NRAS2(ISYM)
          NR3=NRAS3(ISYM)
          NS=NSSH(ISYM)
* Skip inactive transformation matrix:
          IOFF1=IOFF1+NI**2
* Copy RAS1 transformation matrix transposed to TAT:
          DO I=1,NR1
            DO J=1,NR1
              IJ=I+NR1*(J-1)
              JI=J+NR1*(I-1)
              TAT(IOFF2+JI)=TORB(IOFF1+IJ)
            END DO
          END DO
          IOFF1=IOFF1+NR1**2
          IOFF2=IOFF2+NR1**2
* Copy RAS2 transformation matrix transposed to TAT:
          DO I=1,NR2
            DO J=1,NR2
              IJ=I+NR2*(J-1)
              JI=J+NR2*(I-1)
              TAT(IOFF2+JI)=TORB(IOFF1+IJ)
            END DO
          END DO
          IOFF1=IOFF1+NR2**2
          IOFF2=IOFF2+NR2**2
* Copy RAS2 transformation matrix transposed to TAT:
          DO I=1,NR3
            DO J=1,NR3
              IJ=I+NR3*(J-1)
              JI=J+NR3*(I-1)
              TAT(IOFF2+JI)=TORB(IOFF1+IJ)
            END DO
          END DO
          IOFF1=IOFF1+NR3**2
          IOFF2=IOFF2+NR3**2
* Skip virtual transformation matrix:
          IOFF1=IOFF1+NS**2
        END DO
      END IF

      NSG=NCONF
      call mma_allocate(SGM,NSG,Label='GRDSGM')
C Compute (Proj_CAS)*(Ham)*(Wave op) acting on Psi_0:
      CALL DCOPY_(NCONF,[Zero],0,SGM,1)
      IF(ISCF.EQ.0) THEN
       CALL W1TW2(IVECW,IVECC,CI,SGM)
      ELSE
       SGM(1)=E2TOT*CI(1)
      END IF
      IF(IFMSCOUP) THEN
C Multi-State HEFF:
C Loop over all the other root states, and compute off-diagonal
C effective Hamiltonian matrix elements.
        call mma_allocate(BRACI,NCONF,Label='BRACI')
        ID=IDTCEX
        DO ISTATE=1,NSTATE
          IF(ISTATE.NE.JSTATE) THEN
            CALL DDAFILE(LUCIEX,2,BRACI,NCONF,ID)
          ELSE
            CALL DDAFILE(LUCIEX,0,BRACI,NCONF,ID)
          END IF
        END DO
        call mma_deallocate(BRACI)
      END IF
C-End of Multi-State insert -----------------------------------------

      IF(ORBIN.EQ.'TRANSFOR') THEN
C Transform SGM to use original MO:
        ITOEND=0
        DO ISYM=1,NSYM
         NI=NISH(ISYM)
         NR1=NRAS1(ISYM)
         NR2=NRAS2(ISYM)
         NR3=NRAS3(ISYM)
         NS=NSSH(ISYM)
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NR1**2+NR2**2+NR3**2
*         ITO=ITOSTA+NI**2
         ITO=ITOSTA
         IF(NR1.GT.0) THEN
           ISTART=NAES(ISYM)+1
           CALL TRACI_RPT2(ISTART,NR1,TAT(ITO),STSYM,
     &                                           NSG,SGM)
         END IF
         ITO=ITO+NR1**2
         IF(NR2.GT.0) THEN
           ISTART=NAES(ISYM)+NR1+1
           CALL TRACI_RPT2(ISTART,NR2,TAT(ITO),STSYM,
     &                                          NSG,SGM)
         END IF
         ITO=ITO+NR2**2
         IF(NR3.GT.0) THEN
           ISTART=NAES(ISYM)+NR1+NR2+1
           CALL TRACI_RPT2(ISTART,NR1,TAT(ITO),STSYM,
     &                                          NSG,SGM)
         END IF
        END DO
      END IF

C-Multi-State insert ------------------------------------------------
C In Multi-State calculations, it is very convenient to use the opportunity
C here, when the SGM wave function is available, to compute the coupling
C elements in Heff simply by contraction with the bra CAS functions.
      IF(IFMSCOUP) THEN
C Multi-State HEFF:
C Loop over all the other root states, and compute off-diagonal
C effective Hamiltonian matrix elements. Given that the SGM
C wave function is available, this is a very efficient way of
C computing the multi-state coupling elements.
        call mma_allocate(BRACI,NCONF,Label='BRACI')
        ID=IDTCEX
        DO ISTATE=1,NSTATE
          IF(ISTATE.NE.JSTATE) THEN
            CALL DDAFILE(LUCIEX,2,BRACI,NCONF,ID)
            HEFF(ISTATE,JSTATE)=HEFF(ISTATE,JSTATE) +
     &      DDOT_(NCONF,BRACI,1,SGM,1)
          ELSE
            CALL DDAFILE(LUCIEX,0,BRACI,NCONF,ID)
            HEFF(JSTATE,JSTATE)=HEFF(JSTATE,JSTATE)+E2CORR
          END IF
        END DO
        call mma_deallocate(BRACI)
      END IF
C-End of Multi-State insert -----------------------------------------

C Similar, but (Proj_CAS)*((Wave op)**(dagger))*(Ham) | Psi_0 >.
      CALL DCOPY_(NCONF,[Zero],0,SGM,1)
      IF(ISCF.EQ.0) THEN
       CALL W1TW2(IVECC,IVECW,CI,SGM)
      ELSE
       SGM(1)=E2TOT*CI(1)
      END IF

      IF(ORBIN.EQ.'TRANSFOR') THEN
C Transform SGM to use original MO:
        ITOEND=0
        DO ISYM=1,NSYM
         NI=NISH(ISYM)
         NR1=NRAS1(ISYM)
         NR2=NRAS2(ISYM)
         NR3=NRAS3(ISYM)
         NS=NSSH(ISYM)
         ITOSTA=ITOEND+1
         ITOEND=ITOEND+NR1**2+NR2**2+NR3**2
*         ITO=ITOSTA+NI**2
         ITO=ITOSTA
         IF(NR1.GT.0) THEN
           ISTART=NAES(ISYM)+1
           CALL TRACI_RPT2(ISTART,NR1,TAT(ITO),STSYM,
     &                                          NSG,SGM)
         END IF
         ITO=ITO+NR1**2
         IF(NR2.GT.0) THEN
           ISTART=NAES(ISYM)+NR1+1
           CALL TRACI_RPT2(ISTART,NR2,TAT(ITO),STSYM,
     &                                          NSG,SGM)
         END IF
         ITO=ITO+NR2**2
         IF(NR3.GT.0) THEN
           ISTART=NAES(ISYM)+NR1+NR2+1
           CALL TRACI_RPT2(ISTART,NR1,TAT(ITO),STSYM,
     &                                          NSG,SGM)
         END IF
        END DO
      END IF


C--------------------------------------------------------------------
C (C):
C Compute also < Psi_0 | Epq | Psi_1 > and < Psi_0 | Epqrs | Psi_1 >.
C Actually, the first was already computed in PRPCTL.

C ... Unfinished code.
C--------------------------------------------------------------------

      call mma_deallocate(SGM)
      call mma_deallocate(CI)

      END SUBROUTINE GRDCTL
