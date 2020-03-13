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
      SUBROUTINE MCCTL(HEFF)
      IMPLICIT NONE
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
#include "SysDef.fh"
#include "stdalloc.fh"
      REAL*8 HEFF(NSTATE,NSTATE)

      INTEGER ISTATE
      REAL*8 DVALUE
      REAL*8 TOTCPU1, TOTWALL1, TOTCPU2, TOTWALL2
      Character(len=160) string
      real*8, allocatable :: cpu_timing(:), wall_timing(:)

      CALL QENTER('MCCTL')

      Call mma_allocate( cpu_timing,nstate,'timing in mcctl')
      Call mma_allocate(wall_timing,nstate,'timing in mcctl')
      Call dcopy_(nstate,[0.0d0],0, cpu_timing,1)
      Call dcopy_(nstate,[0.0d0],0,wall_timing,1)
C The ket state is JSTATE.
C Loop over the bra states
      DO ISTATE=1,NSTATE
        Write(string,'(A,I4,A,I4,A,I4)')
     &     'Multistate coupling between state',ISTATE,' and',JSTATE,
     &     ' out of ',NSTATE
        Call StatusLine('CASPT2: MCCTL: ',trim(string))
        TOTCPU1=0.0d0; TOTWALL1=0.0d0; TOTCPU2=0.0d0; TOTWALL2=0.0d0;
        Call CWTIME(TOTCPU1,TOTWALL1) !start clock for total time

        IF(ISTATE.EQ.JSTATE) THEN
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
       WRITE(6,*)
       WRITE(6,*) 'Hamiltonian Effective Couplings'
       WRITE(6,*) '-------------------------------'
       WRITE(6,*)
       WRITE(6,'(10X,6X,A3,I4,A3)') ' | ', MSTATE(JSTATE), ' > '
       DO ISTATE=1,NSTATE
        WRITE(6,'(A3,I4,A3,ES22.14)')
     &   ' < ',MSTATE(ISTATE),' | ', HEFF(ISTATE,JSTATE)
       ENDDO
      ENDIF

      IF ((IPRGLB.GE.VERBOSE).OR.(NLYROOT.NE.0)) THEN
       WRITE(6,*)
       WRITE(6,'(A,I4,A)')
     &           'Time spent for multi-state couplings for root ',
     &            MSTATE(JSTATE),':'
       WRITE(6,*) '----------------- CPU TIME  -------- WALL TIME'
       DO ISTATE=1,NSTATE
        WRITE(6,'(A3,I4,A,F18.3,2x,F18.3)') ' < ',MSTATE(ISTATE),' |',
     &    cpu_timing(istate),  wall_timing(istate)
       ENDDO
      ENDIF
      Call mma_deallocate(cpu_timing)
      Call mma_deallocate(wall_timing)

      CALL QEXIT('MCCTL')
      RETURN
      END
