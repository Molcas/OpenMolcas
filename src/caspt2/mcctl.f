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
*               2019, Stefano Battaglia                                *
************************************************************************
*--------------------------------------------*
* 1998  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
*--------------------------------------------*
      SUBROUTINE MCCTL(HEFF,NSTATE,JSTATE)
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: E2Corr, NLYRoot, mState
      use constants, only: Zero
      use definitions, only: iwp,wp,u6
      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: NSTATE, JSTATE
      REAL(kind=wp), intent(inout):: HEFF(NSTATE,NSTATE)

      INTEGER(kind=iwp) ISTATE
      REAL(kind=wp) DVALUE
      REAL(kind=wp) TOTCPU1, TOTWALL1, TOTCPU2, TOTWALL2
      Character(len=160) string
      real(kind=wp), allocatable :: cpu_timing(:), wall_timing(:)


      Call mma_allocate( cpu_timing,nstate,'timing in mcctl')
      Call mma_allocate(wall_timing,nstate,'timing in mcctl')
      cpu_timing(:)=Zero
      wall_timing(:)=Zero

C The ket state is JSTATE.
C Loop over the bra states

      DO ISTATE=1,NSTATE

        Write(string,'(A,I0,A,I0,A,I0)')
     &     'Multistate coupling between state ',ISTATE,' and ',JSTATE,
     &     ' out of ',NSTATE
        Call StatusLine('CASPT2: MCCTL: ',string)
        TOTCPU1=Zero; TOTWALL1=Zero; TOTCPU2=Zero; TOTWALL2=Zero;
        Call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        IF(ISTATE==JSTATE) THEN
          HEFF(ISTATE,JSTATE)=HEFF(ISTATE,JSTATE)+E2CORR
        ELSE
C Compute the effective Hamiltonian:
          CALL HEFVAL(ISTATE,JSTATE,DVALUE)
          HEFF(ISTATE,JSTATE)=HEFF(ISTATE,JSTATE)+DVALUE
        END IF

        Call CWTIME(TOTCPU2,TOTWALL2)
        cpu_timing(istate)=TOTCPU2-TOTCPU1
        wall_timing(istate)=TOTWALL2-TOTWALL1

      END DO

      IF ((IPRGLB.GE.VERBOSE).OR.(NLYROOT.NE.0)) THEN
       WRITE(u6,*)
       WRITE(u6,*) 'Hamiltonian Effective Couplings'
       WRITE(u6,*) '-------------------------------'
       WRITE(u6,*)
       WRITE(u6,'(10X,6X,A3,I4,A3)') ' | ', MSTATE(JSTATE), ' > '
       DO ISTATE=1,NSTATE
        WRITE(u6,'(A3,I4,A3,ES22.14)')
     &   ' < ',MSTATE(ISTATE),' | ', HEFF(ISTATE,JSTATE)
       ENDDO
      ENDIF

      IF ((IPRGLB.GE.VERBOSE).OR.(NLYROOT.NE.0)) THEN
       WRITE(u6,*)
       WRITE(u6,'(A,I4,A)')
     &           'Time spent for multi-state couplings for root ',
     &            MSTATE(JSTATE),':'
       WRITE(u6,*) '----------------- CPU TIME  -------- WALL TIME'
       DO ISTATE=1,NSTATE
        WRITE(u6,'(A3,I4,A,F18.3,2x,F18.3)') ' < ',MSTATE(ISTATE),' |',
     &    cpu_timing(istate),  wall_timing(istate)
       ENDDO
      ENDIF
      Call mma_deallocate(cpu_timing)
      Call mma_deallocate(wall_timing)

      END SUBROUTINE MCCTL
