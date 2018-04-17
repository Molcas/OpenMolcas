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
      SUBROUTINE Write_QMMM(Coord,nAtIn,Iter)

      IMPLICIT NONE

#include "constants2.fh"
#include "Molcas.fh"
#include "periodic_table.fh"
#include "qmmm.fh"
#include "stdalloc.fh"

      INTEGER, INTENT(IN) :: Iter,nAtIn
      REAL*8, DIMENSION(3*nAtIn,Iter), INTENT(IN) :: Coord

      INTEGER :: i,iAtIn,iAtNum,iAtOut,iAtTot,iFirst,k,Lu_XYZ,nAtOut
      INTEGER :: nAtTot
      INTEGER, DIMENSION(:), ALLOCATABLE :: AT
      REAL*8, DIMENSION(:), ALLOCATABLE :: Charge
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CoordMMO
      CHARACTER(LEN=16) :: FileName
      CHARACTER(LEN=LENIN) :: Symbol
      CHARACTER(LEN=LENIN), DIMENSION(:), ALLOCATABLE :: LabMMO
      LOGICAL :: Found,isMMI,isMMO,isQM

      CALL Qpg_dArray('MMO Coords',Found,nAtOut)
      IF (Found) THEN
         nAtOut = nAtOut/3
         nAtTot = nAtIn+nAtOut
         CALL mma_allocate(Charge,nAtIn)
         CALL mma_allocate(CoordMMO,3,nAtOut)
         CALL mma_allocate(LabMMO,nAtOut)
         CALL mma_allocate(AT,nAtTot)
         CALL Get_dArray('Nuclear charge',Charge,nAtIn)
         CALL Get_dArray('MMO Coords',CoordMMO,3*nAtOut)
         CALL Get_cArray('MMO Labels',LabMMO,LENIN*nAtOut)
         CALL Get_iArray('Atom Types',AT,nAtTot)

! Write one file per iteration and one for the current geometry
         DO k = 1,2
            IF (k==1) THEN
               WRITE(FileName,'(I16.4)') Iter
               FileName = 'QMMMITXYZ.'//TRIM(ADJUSTL(FileName))
            ELSE
               FileName='QMMMENDXYZ'
            END IF
            Lu_XYZ = 1
            CALL molcas_open(Lu_XYZ,FileName)
            WRITE(Lu_XYZ,'(I6)') nAtTot
            WRITE(Lu_XYZ,*)
            iAtIn = 1
            iAtOut = 1
100         FORMAT (A6,3(1X,F12.6))
            DO iAtTot = 1,nAtTot
               isQM = AT(iAtTot)==QM
               isMMI = AT(iAtTot)==MMI
               isMMO = AT(iAtTot)==MMO
               IF (isQM.OR.isMMI) THEN
                  iAtNum = NINT(Charge(iAtIn))
                  iFirst = INDEX(PTab(iAtNum),' ')+1
                  Symbol = PTab(iAtNum)(iFirst:2)
                  WRITE(Lu_XYZ,100) Symbol,
     &                 (Coord(3*(iAtIn-1)+i,Iter)*Angstrom,i=1,3)
                  iAtIn = iAtIn+1
               ELSE IF (isMMO) THEN
                  Symbol = LabMMO(iAtOut)
                  IF (INDEX(Symbol,'_')>1) THEN
                     Symbol=Symbol(1:INDEX(Symbol,'_')-1)
                  END IF
                  WRITE(Lu_XYZ,100) Symbol,
     &                 (CoordMMO(i,iAtOut)*Angstrom,i=1,3)
                  iAtOut = iAtOut+1
               END IF
            END DO
            CLOSE(Lu_XYZ)
         END DO
         CALL mma_deallocate(Charge)
         CALL mma_deallocate(CoordMMO)
         CALL mma_deallocate(LabMMO)
         CALL mma_deallocate(AT)
      END IF

      RETURN
      END
