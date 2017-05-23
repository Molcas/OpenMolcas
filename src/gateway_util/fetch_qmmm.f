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
#ifdef _GROMACS_
      SUBROUTINE Fetch_QMMM(CastMM,nCastMM)

      IMPLICIT NONE

#include "Molcas.fh"
#include "periodic_table.fh"
#include "qmmm.fh"
#include "real.fh"
#include "stdalloc.fh"

      INTEGER, INTENT(IN) :: nCastMM
      INTEGER, DIMENSION(nCastMM), INTENT(IN) :: CastMM

      INTEGER :: iAtGMX,iAtNmbGMX,iAtOut,iCastMM,iFirst,iGrpGMX,iLast
      INTEGER :: iOk,ipCR,ipGMS,iXYZ,LuWr,LuXYZ,nAtGMX,nAtIn,nAtOut
      INTEGER, DIMENSION(:), ALLOCATABLE :: AT
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CoordGMX,CoordMMO
      CHARACTER(LEN=256) :: LogFileName, Message, TPRFileName
      CHARACTER(LEN=LENIN) :: Symbol
      CHARACTER(LEN=LENIN), DIMENSION(:), ALLOCATABLE :: LabMMO
      LOGICAL :: Exist

      INTEGER, EXTERNAL :: init_commrec,isFreeUnit,mmslave_copyx
      INTEGER, EXTERNAL :: mmslave_get_atomnumber,mmslave_get_group_id
      INTEGER, EXTERNAL :: mmslave_init,mmslave_natoms,mmslave_read_tpr

      LuWr = 6

! Initialize Gromacs mmslave
      ipCR = init_commrec()
      CALL prgmtranslate('GMX.LOG',LogFileName,iLast)
      LogFileName(iLast+1:iLast+1) = CHAR(0)
      ipGMS = mmslave_init(%val(ipCR),LogFileName)

! Tell Gromacs to read tpr file
      TPRFileName = TPRDefName
      iLast = LEN_TRIM(TPRFileName)
      CALL f_inquire(TPRFileName,Exist)
      IF (.NOT.Exist) THEN
         Message = 'File '//TPRFileName(1:iLast)//' not found'
         CALL WarningMessage(2,Message)
         CALL Quit_OnUserError()
      END IF
      TPRFileName(iLast+1:iLast+1) = CHAR(0)
      iOk = mmslave_read_tpr(TPRFileName,%val(ipGMS))
      IF (iOk.NE.1) THEN
         Message = 'Error reading tpr file'
         CALL WarningMessage(2,Message)
         CALL Quit_OnUserError()
      END IF

! Fetch coordinates from Gromacs
      nAtGMX = mmslave_natoms(%val(ipGMS))
      CALL mma_allocate(CoordGMX,3,nAtGMX)
      CALL dcopy_(3*nAtGMX,Zero,0,CoordGMX,1)
      iOk = mmslave_copyX(%val(ipGMS),%val(nAtGMX),CoordGMX)
      IF (iOk.NE.1) THEN
         Message = 'Fetch_QMMM: mmslave_copyx is not ok'
         CALL WarningMessage(2,Message)
         CALL Abend()
      END IF

! Find out and count atom types
      nAtIn = 0
      nAtOut = 0
      CALL mma_allocate(AT,nAtGMX)
      DO iAtGMX = 1,nAtGMX
         iGrpGMX = mmslave_get_group_id(%val(ipGMS),%val(iAtGMX-1))
         IF (iGrpGMX==QMGMX) THEN
            AT(iAtGMX) = QM
         ELSE IF (iGrpGMX==MMGMX) THEN
            AT(iAtGMX) = MMO
         ELSE
            Message = 'Fetch_QMMM: unknown index group'
            CALL WarningMessage(2,Message)
            CALL Abend()
         END IF
      END DO
      DO iCastMM = 1,nCastMM
         IF (AT(CastMM(iCastMM))==QM) THEN
            Message = 'Attempting to cast atom of QM type'
            CALL WarningMessage(2,Message)
            CALL Quit_OnUserError()
         ELSE IF (AT(CastMM(iCastMM))==MMI) THEN
            Message = 'Attempting to cast atom of MMI type'
            CALL WarningMessage(2,Message)
            CALL Quit_OnUserError()
         END IF
         AT(CastMM(iCastMM)) = MMI
      END DO
      DO iAtGMX = 1,nAtGMX
         IF (AT(iAtGMX)==QM.OR.AT(iAtGMX)==MMI) THEN
            nAtIn = nAtIn+1
         ELSE
            nAtOut = nAtOut+1
         END IF
      END DO

! Put QM and inner MM atoms in xyz file and outer MM atoms on runfile
      LuXYZ = 1
      LuXYZ = isFreeUnit(LuXYZ)
      CALL molcas_open(LuXYZ,'GMX.XYZ')
      WRITE(LuXYZ,'(i5,/)') nAtIn
      CALL mma_allocate(CoordMMO,3,nAtOut)
      CALL mma_allocate(LabMMO,nAtOut)
      iAtOut = 1
      DO iAtGMX = 1,nAtGMX
         iAtNmbGMX = mmslave_get_atomnumber(%val(ipGMS),%val(iAtGMX-1))
         iFirst = INDEX(PTab(iAtNmbGMX),' ')+1
         IF (AT(iAtGMX)==QM) THEN
            Symbol = PTab(iAtNmbGMX)(iFirst:2)
         ELSE
            Symbol = PTab(iAtNmbGMX)(iFirst:2)//'_MM'
         END IF
         IF (AT(iAtGMX)==QM.OR.AT(iAtGMX)==MMI) THEN
            WRITE(LuXYZ,'(a8,3E24.15E3)') Symbol(1:5),
     &           (CoordGMX(iXYZ,iAtGMX)*NmToAng,iXYZ=1,3)
         ELSE
            DO iXYZ = 1,3
               CoordMMO(iXYZ,iAtOut) = CoordGMX(iXYZ,iAtGMX)
            END DO
            LabMMO(iAtOut) = Symbol
            iAtOut = iAtOut+1
         END IF
      END DO
      CLOSE(LuXYZ)
      CALL Put_iArray('Atom Types',AT,nAtGMX)
      CALL dscal_(3*nAtOut,One/AuToNm,CoordMMO,1)
      CALL Put_dArray('MMO Coords',CoordMMO,3*nAtOut)
      CALL Put_cArray('MMO Labels',LabMMO,LENIN*nAtOut)

! Clean up
      CALL mma_deallocate(CoordGMX)
      CALL mma_deallocate(AT)
      CALL mma_deallocate(CoordMMO)
      CALL mma_deallocate(LabMMO)
      CALL mmslave_done(%val(ipGMS))

      RETURN
      END
#elif defined (NAGFOR)
      SUBROUTINE empty_Fetch_QMMM()
      END
#endif
