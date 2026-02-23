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
* Copyright (C) 1998, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE STINI(JSTATE)
#ifdef _DMRG_
      use qcmaquis_interface, only:qcmaquis_interface_set_state
      use iso_c_binding, only: c_int
      use caspt2_module, only: DMRG
#endif
      use caspt2_global, only:iPrGlb
      use caspt2_global, only: DREF, PREF
      use PrintLevel, only: DEBUG, USUAL
      use caspt2_module, only: CPUFG3, ERef, nAshT, EASUM,
     &                         TIOFG3, EPSA, mState, RefEne,
     &                         CPUSIN, TIOSIN
      use pt2_guga, only: iAdr10, CLab10
      use constants, only: Zero
      use definitions, only: iwp, wp, u6
      IMPLICIT NONE
      integer(kind=iwp), intent(in):: JSTATE

      CHARACTER(LEN=50)  STLNE2
C     timers
      REAL(kind=wp) CPU0,CPU1,CPU,  CPTF0, CPTF11, CPE,
     &       TIO0,TIO1,TIO, TIOTF0,TIOTF11,TIOE
C     indices
      INTEGER(kind=iwp) :: I,J,IFTEST=0
************************************************************************
      CALL TIMING(CPTF0,CPE,TIOTF0,TIOE)
************************************************************************

      Write(STLNE2,'(A,I0)')
     &                'Compute H0 matrices for state ',MSTATE(JSTATE)
      Call StatusLine('CASPT2: ',STLNE2)
      IF(IPRGLB.GE.USUAL) THEN
        WRITE(u6,'(20A4)')('****',I=1,20)
        WRITE(u6,'(A,I4)')
     &   ' Compute H0 matrices for state ',MSTATE(JSTATE)
        WRITE(u6,'(20A4)')('----',I=1,20)
      END IF

* Reinitialize labels for saving density matrices on disk.
* The fields IADR10 and CLAB10 are kept in the module pt2_guga.F90
      IADR10(:,1)=-1
      IADR10(:,2)=0
      CLAB10(:)='   EMPTY'
      IADR10(1,1)=0

#ifdef _DMRG_
      if (DMRG) then
        ! set state number here because in poly1 we have no reference
        ! to which state we are computing
        if (iPrGlb >= DEBUG) then
          write (u6,*) 'STINI setting DMRG state number to ',
     &                mstate(jstate)-1
        endif
        ! Convert to the root number despite having
        ! set only the checkpoint file paths for the desired state(s)
        call qcmaquis_interface_set_state(int(mstate(jstate)-1,c_int))
      end if
#endif
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(u6,*)' STINI calling POLY3...'
      END IF
      CALL TIMING(CPU0,CPU,TIO0,TIO)

      CALL POLY3(1)

      CALL TIMING(CPU1,CPU,TIO1,TIO)
      CPUFG3=CPU1-CPU0
      TIOFG3=TIO1-TIO0
      IF (IPRGLB.GE.DEBUG) THEN
        WRITE(u6,*)' STINI back from POLY3.'
      END IF

* GETDPREF: Restructure GAMMA1 and GAMMA2, as DREF and PREF arrays.
      CALL GETDPREF(DREF,SIZE(DREF),PREF,SIZE(PREF))

      IF ( IFTEST.NE.0 ) THEN
        WRITE(u6,*)' DREF for state nr. ',MSTATE(JSTATE)
        DO I=1,NASHT
          WRITE(u6,'(1x,14f10.6)')(DREF((I*(I-1))/2+J),J=1,I)
        END DO
        WRITE(u6,*)
      END IF

      EREF=REFENE(JSTATE)
* With new DREF, recompute EASUM:
      EASUM=Zero
      DO I=1,NASHT
        EASUM=EASUM+EPSA(I)*DREF((I*(I+1))/2)
      END DO

      IF(IPRGLB.GE.USUAL) THEN
       WRITE(u6,'(20A4)')('----',I=1,20)
       WRITE(u6,'(A)')' H0 matrices have been computed.'
       WRITE(u6,*)
      ENDIF
************************************************************************
      CALL TIMING(CPTF11,CPE,TIOTF11,TIOE)
      CPUSIN=CPTF11-CPTF0
      TIOSIN=TIOTF11-TIOTF0
************************************************************************
      END SUBROUTINE STINI
