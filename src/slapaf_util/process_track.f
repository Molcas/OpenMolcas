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
* Copyright (C) 2016, Ignacio Fdez. Galvan                             *
************************************************************************
      SUBROUTINE Process_Track()
      use Slapaf_Info, only: RootMap
      use Slapaf_parameters, only: Request_RASSI
      IMPLICIT NONE
#include "print.fh"
#include "nadc.fh"
#include "real.fh"
#include "stdalloc.fh"
      INTEGER :: nOv,nRoots,i
      INTEGER, DIMENSION(:), ALLOCATABLE :: OldMap,RootIdx
      INTEGER, DIMENSION(2) :: MaxId
      REAL*8, DIMENSION(:), ALLOCATABLE :: Ovlp
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: Overlaps
      LOGICAL :: Found,Done
      CHARACTER(LEN=8) :: Method
*
      CALL Get_cArray('Relax Method',Method,8)
      IF ((Method .NE. 'CASSCF')  .AND.
     &    (Method .NE. 'RASSCF')  .AND.
     &    (Method .NE. 'CASSCFSA').AND.
     &    (Method .NE. 'RASSCFSA').AND.
     &    (Method .NE. 'CASPT2')  .AND.
     &    (Method .NE. 'RASPT2')  ) THEN
      CALL WarningMessage(2,'Error in Process_Track')
        WRITE(6,*) '***************** ERROR ********************'
        WRITE(6,*) ' The TRACK keyword can only be used with'
        WRITE(6,*) ' states computed by the RASSCF or CASPT2'
        WRITE(6,*) ' programs.'
        WRITE(6,*) '********************************************'
        CALL Quit_OnUserError()
      END IF
*
* Find the number of roots, and whether state overlaps are available
*
      nRoots=1
      CALL Qpg_iScalar('Number of roots',Found)
      IF (Found) CALL Get_iScalar('Number of roots',nRoots)
      CALL qpg_dArray('State Overlaps',Found,nOv)
      IF (Found.AND.(nOv.EQ.(2*nRoots)**2)) THEN
        Request_RASSI = .FALSE.
      ELSE
        Request_RASSI = .TRUE.
      END IF
*
* Make sure that the root mapping is only done once per iteration
*
      CALL Qpg_iScalar('Track Done',Done)
      CALL Get_lScalar('Track Done',Done)
      IF (Request_RASSI .OR. Done) RETURN
*
* Modify the root map according to the overlaps:
*
*   RootMap(i) = N
*
* where i is the original root number (at iter=1), and N is the
* root number at the current iteration.
*
      CALL mma_allocate(Ovlp,nOv)
      CALL mma_allocate(Overlaps,nRoots,nRoots)
      CALL get_dArray('State Overlaps',Ovlp,nOv)
      DO i=1,nRoots
        Overlaps(:,i)=Ovlp((2*i-1)*nRoots+1:2*i*nRoots)
      END DO
      CALL mma_deallocate(Ovlp)
      IF (nPrint(1).GE.5) THEN
        CALL RecPrt('Overlaps with previous states','',
     &              Overlaps,nRoots,nRoots)
        CALL mma_allocate(OldMap,nRoots)
        CALL iCopy(nRoots,RootMap,1,OldMap,1)
      END IF
      CALL mma_allocate(RootIdx,nRoots)
      DO i=1,nRoots
        MaxId=MAXLOC(ABS(Overlaps))
        RootIdx(MaxId(1))=MaxId(2)
        Overlaps(MaxId(1),:)=Zero
        Overlaps(:,MaxId(2))=Zero
      END DO
      DO i=1,nRoots
        RootMap(i)=RootIdx(RootMap(i))
      END DO
      CALL Put_iArray('Root Mapping',RootMap,nRoots)
      IF (nPrint(1).GE.5) THEN
        WRITE(6,*)
        WRITE(6,100) 'Root map'
        WRITE(6,*)
        WRITE(6,100) 'Original  Prev.  This'
        WRITE(6,100) '  root    iter.  iter.'
        WRITE(6,100) '----------------------'
        DO i=1,nRoots
          WRITE(6,101) i,OldMap(i),RootMap(i)
        END DO
        WRITE(6,100) '----------------------'
        WRITE(6,*)
        CALL mma_deallocate(OldMap)
      END IF
100   FORMAT(3X,A)
101   FORMAT(3X,I6,1X,I6,1X,I6)

      CALL mma_deallocate(RootIdx)
      CALL mma_deallocate(Overlaps)
*
* Update the RunFile for automatic gradients
*
      CALL Qpg_iScalar('Relax CASSCF root',Found)
      IF (Found) THEN
        Call Get_iScalar('Relax CASSCF root',i)
        IF (RootMap(i).NE.i) THEN
          CALL Put_iScalar('Relax CASSCF root',RootMap(i))
          CALL Put_iScalar('Relax Original root',RootMap(i))
        END IF
      END IF
      CALL Qpg_iScalar('NumGradRoot',Found)
      IF (Found) THEN
        Call Get_iScalar('NumGradRoot',i)
        IF (RootMap(i).NE.i)
     &    CALL Put_iScalar('NumGradRoot',RootMap(i))
      END IF
      CALL Put_lScalar('Track Done',.TRUE.)
*
      RETURN
*
      END SUBROUTINE Process_Track
