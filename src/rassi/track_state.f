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
      SUBROUTINE Track_State(OVLP)
      IMPLICIT NONE
#include "Molcas.fh"
#include "cntrl.fh"
#include "real.fh"
#include "prgm.fh"
      INTEGER iState,initState,newState
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      INTEGER j
#endif
      REAL*8 MaxOv,ThisOv
      REAL*8 ovlp(nstate,nstate)
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='Track_State')



*     Check that there are 2 JOB files, with the same number of states
      IF (nJob.ne.2) THEN
        Call SysAbendMsg('Track_State',
     &  'The number of JOB files should be 2.','')
      END IF
      IF (nStat(2).ne.nStat(1)) THEN
        Call SysAbendMsg('Track_State',
     &  'The number of states in the JOB files should be the same.','')
      END IF

*     Get the root to track
      CALL Get_iScalar('Relax CASSCF root',initState)

*     Find the root in the first JOB (current states) that has maximum
*     overlap with the tracked root in the second JOB (previous states)
#ifdef _DEBUGPRINT_
      WRITE(6,*)
      WRITE(6,*)'OVERLAP MATRIX FOR THE ORIGINAL STATES:'
      WRITE(6,*)
      DO iState=1,nState
        WRITE(6,'(5(1X,F15.8))')(Ovlp(iState,j),j=1,iState)
      END DO
#endif
      IF (IPGLOB.ge.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*) 'Initial root: ',initState
        WRITE(6,*) 'Overlaps with current states:'
      END IF
      MaxOv=Zero
      newState=0
      DO iState=1,nStat(1)
        ThisOv=Ovlp(iState,initState+nStat(1))
        IF (IPGLOB.ge.USUAL) THEN
          WRITE(6,'(I5,1X,F15.8)') iState,ThisOv
        END IF
        IF (ABS(ThisOv).gt.MaxOv) THEN
          MaxOv=ABS(ThisOv)
          newState=iState
        END IF
      END DO
      IF (IPGLOB.ge.USUAL) THEN
        WRITE(6,*) 'New root: ',newState
      END IF

*     If no state is found, something wrong has happened
      IF (newState.eq.0) THEN
        Call SysAbendMsg('Track_State',
     &  'No overlaps!','')
      END IF

*     Store the new state for geometry optimization
      IF (newState.ne.initState) THEN
        CALL Put_iScalar('Relax CASSCF root',newState)
        CALL Put_iScalar('Relax Original root',newState)
        CALL Put_iScalar('NumGradRoot',newState)
      END IF

      END
