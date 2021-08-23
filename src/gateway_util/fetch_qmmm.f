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

      USE, INTRINSIC :: iso_c_binding, only: c_int, c_loc, c_ptr
      IMPLICIT NONE

#include "Molcas.fh"
#include "periodic_table.fh"
#include "qmmm.fh"
#include "real.fh"
#include "stdalloc.fh"

      INTEGER, INTENT(IN) :: nCastMM
      INTEGER, DIMENSION(nCastMM), INTENT(IN) :: CastMM

      INTEGER :: iAtGMX,iAtNmbGMX,iAtOut,iCastMM,iFirst,iGrpGMX,iLast
      INTEGER :: iOk,iXYZ,LuWr,LuXYZ,nAtGMX,nAtIn,nAtOut
      TYPE(c_ptr) :: ipCR, ipGMS
      INTEGER, DIMENSION(:), ALLOCATABLE :: AT
      REAL*8, DIMENSION(:,:), ALLOCATABLE :: CoordGMX,CoordMMO
      CHARACTER(LEN=256) :: LogFileName, Message, TPRFileName
      CHARACTER(LEN=LENIN) :: Symbol
      CHARACTER(LEN=LENIN), DIMENSION(:), ALLOCATABLE :: LabMMO
      LOGICAL :: Exist

      INTEGER, EXTERNAL :: isFreeUnit
      INTERFACE
        SUBROUTINE mmslave_done(gms) BIND(C,NAME='mmslave_done_')
          USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr
          TYPE(c_ptr), VALUE :: gms
        END SUBROUTINE mmslave_done
        FUNCTION mmslave_init(cr,log_) BIND(C,NAME='mmslave_init_')
          USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_ptr
          TYPE(c_ptr) :: mmslave_init
          TYPE(c_ptr), VALUE :: cr
          CHARACTER(kind=c_char) :: log_(*)
        END FUNCTION mmslave_init
        FUNCTION mmslave_get_atomnumber(gms,id)
     &           BIND(C,NAME='mmslave_get_atomnumber_')
          USE, INTRINSIC :: iso_c_binding, only: c_int, c_ptr
          INTEGER(kind=c_int) :: mmslave_get_atomnumber
          TYPE(c_ptr), VALUE :: gms
          INTEGER(kind=c_int), VAlUE :: id
        END FUNCTION mmslave_get_atomnumber
        FUNCTION mmslave_get_group_id(gms,id)
     &           BIND(C,NAME='mmslave_get_group_id_')
          USE, INTRINSIC :: iso_c_binding, only: c_int, c_ptr
          INTEGER(kind=c_int) :: mmslave_get_group_id
          TYPE(c_ptr), VALUE :: gms
          INTEGER(kind=c_int), VAlUE :: id
        END FUNCTION mmslave_get_group_id
        FUNCTION mmslave_natoms(gms) BIND(C,NAME='mmslave_natoms_')
          USE, INTRINSIC :: iso_c_binding, only: c_int, c_ptr
          INTEGER(kind=c_int) :: mmslave_natoms
          TYPE(c_ptr), VALUE :: gms
        END FUNCTION mmslave_natoms
        FUNCTION mmslave_read_tpr(tpr,gms)
     &           BIND(C,NAME='mmslave_read_tpr_')
          USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_int, c_ptr
          INTEGER(kind=c_int) :: mmslave_read_tpr
          CHARACTER(kind=c_char) :: tpr(*)
          TYPE(c_ptr), VALUE :: gms
        END FUNCTION mmslave_read_tpr
        FUNCTION init_commrec() BIND(C,NAME='init_commrec_')
          USE, INTRINSIC :: iso_c_binding, ONLY: c_ptr
          TYPE(c_ptr) :: init_commrec
        END FUNCTION init_commrec
      END INTERFACE

      LuWr = 6

! Initialize Gromacs mmslave
      ipCR = init_commrec()
      CALL prgmtranslate('GMX.LOG',LogFileName,iLast)
      LogFileName(iLast+1:iLast+1) = CHAR(0)
      ipGMS = mmslave_init(ipCR,LogFileName)

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
      iOk = mmslave_read_tpr(TPRFileName,ipGMS)
      IF (iOk.NE.1) THEN
         Message = 'Error reading tpr file'
         CALL WarningMessage(2,Message)
         CALL Quit_OnUserError()
      END IF

! Fetch coordinates from Gromacs
      nAtGMX = mmslave_natoms(ipGMS)
      CALL mma_allocate(CoordGMX,3,nAtGMX)
      CALL dcopy_(3*nAtGMX,Zero,0,CoordGMX,1)
      iOk = mmslave_copyX_wrapper(ipGMS,INT(nAtGMX,kind=c_int),CoordGMX)
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
         iGrpGMX = mmslave_get_group_id(ipGMS,INT(iAtGMX-1,kind=c_int))
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
         iAtNmbGMX = mmslave_get_atomnumber(ipGMS,
     &                                      INT(iAtGMX-1,kind=c_int))
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
      CALL mmslave_done(ipGMS)

      RETURN

      CONTAINS

      FUNCTION mmslave_copyx_wrapper(gms,natoms,x)
      INTEGER :: mmslave_copyx_wrapper
      TYPE(c_ptr) :: gms
      INTEGER(kind=c_int) :: natoms
      REAL*8, TARGET :: x(*)
      INTERFACE
        FUNCTION mmslave_copyx(gms,natoms,x)
     &           BIND(C,NAME='mmslave_copyx_')
          USE, INTRINSIC :: iso_c_binding, only: c_int, c_ptr
          INTEGER(kind=c_int) :: mmslave_copyx
          TYPE(c_ptr), VALUE :: gms, x
          INTEGER(kind=c_int), VALUE :: natoms
        END FUNCTION mmslave_copyx
      END INTERFACE
      mmslave_copyx_wrapper = mmslave_copyx(gms,natoms,c_loc(x(1)))
      END FUNCTION mmslave_copyx_wrapper

      END
#elif defined (NAGFOR)
      SUBROUTINE empty_Fetch_QMMM()
      END
#endif
