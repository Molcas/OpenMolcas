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
      SUBROUTINE CHO_MAXABSDIAG(DIAG,IRED,DGMAX)
C
C     Purpose: set max. abs. DIAG (reduced set IRED) in each symmetry, and
C              return the global max. abs. in DGMAX.
C
      use ChoSwp, only: IndRed
#include "implicit.fh"
      DIMENSION DIAG(*)
#include "cholesky.fh"
#include "choptr.fh"

      CHARACTER*14 SECNAM
      PARAMETER (SECNAM = 'CHO_MAXABSDIAG')

      INTEGER AB, AB1, AB2

      LOGICAL LOCDBG
#if defined (_DEBUGPRINT_)
      PARAMETER (LOCDBG = .TRUE.)
#else
      PARAMETER (LOCDBG = .FALSE.)
#endif

      IF (CHO_1CENTER) THEN ! specialization for 1-center approximation
         CALL CHO_MAXABSDIAG_1C(DIAG,IRED,DGMAX)
         RETURN
      END IF

      IF (IRED .EQ. 1) THEN
         DO ISYM = 1,NSYM
            IF (NNBSTR(ISYM,IRED) .LT. 1) THEN
               DIAMAX(ISYM) = 0.0D0
            ELSE
               DIAMAX(ISYM) = ABS(DIAG(IIBSTR(ISYM,IRED)+1))
               AB1 = IIBSTR(ISYM,IRED) + 2
               AB2 = IIBSTR(ISYM,IRED) + NNBSTR(ISYM,IRED)
               DO AB = AB1,AB2
                  DIAMAX(ISYM) = MAX(DIAMAX(ISYM),ABS(DIAG(AB)))
               END DO
            END IF
            DIAMAXT(ISYM)=DIAMAX(ISYM)
         END DO
      ELSE IF ((IRED.EQ.2) .OR. (IRED.EQ.3)) THEN
         DO ISYM = 1,NSYM
            IF (NNBSTR(ISYM,IRED) .LT. 1) THEN
               DIAMAX(ISYM) = 0.0D0
            ELSE
               AB = INDRED(IIBSTR(ISYM,IRED)+1,IRED)
               DIAMAX(ISYM) = ABS(DIAG(AB))
               AB1 = IIBSTR(ISYM,IRED) + 2
               AB2 = IIBSTR(ISYM,IRED) + NNBSTR(ISYM,IRED)
               DO IAB = AB1,AB2
                  AB = INDRED(IAB,IRED)
                  DIAMAX(ISYM) = MAX(DIAMAX(ISYM),ABS(DIAG(AB)))
               END DO
            END IF
            IF (NNBSTR(ISYM,1) .LT. 1) THEN
               DIAMAXT(ISYM) = 0.0D0
            ELSE
               DIAMAXT(ISYM) = ABS(DIAG(IIBSTR(ISYM,1)+1))
               AB1 = IIBSTR(ISYM,1) + 2
               AB2 = IIBSTR(ISYM,1) + NNBSTR(ISYM,1)
               DO AB = AB1,AB2
                  DIAMAXT(ISYM) = MAX(DIAMAXT(ISYM),ABS(DIAG(AB)))
               END DO
            END IF
         END DO
      ELSE
         WRITE(LUPRI,*) SECNAM,': unknown reduced set, IRED = ',IRED
         CALL CHO_QUIT('Unknown reduced set in '//SECNAM,104)
      END IF

      DGMAX = DIAMAX(1)
      DO ISYM = 2,NSYM
         DGMAX = MAX(DGMAX,DIAMAX(ISYM))
      END DO

      IF (LOCDBG) THEN
         WRITE(LUPRI,*) SECNAM,': in reduced set ',IRED,':'
         WRITE(LUPRI,*) 'DIAMAX  = ',(DIAMAX(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,*) 'DIAMAXT = ',(DIAMAXT(ISYM),ISYM=1,NSYM)
         WRITE(LUPRI,*) 'DGMAX   = ',DGMAX
      END IF

      END
      SubRoutine Cho_MaxAbsDiag_1C(Diag,iLoc,DGMax)
C
C     Specialization for 1-Center approximation: only find max for
C     1-center diagonals.
C
      use ChoArr, only: iSP2F, iAtomShl
      use ChoSwp, only: nnBstRSh, iiBstRSh, IndRed
#include "implicit.fh"
      Real*8 Diag(*)
#include "cholesky.fh"
#include "choptr.fh"

      Character*17 SecNam
      Parameter (SecNam = 'Cho_MaxAbsDiag_1C')

      Logical LocDbg
#if defined (_DEBUGPRINT_)
      Parameter (LocDbg = .true.)
#else
      Parameter (LocDbg = .false.)
#endif

      If (iLoc .eq. 1) Then
         Do iSym = 1,nSym
            DiaMax(iSym) = 0.0d0
            Do iShlAB = 1,nnShl
               Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
               If (iAtomShl(iShlA) .eq. iAtomShl(iShlB)) Then
                  i1 = iiBstR(iSym,1) + iiBstRSh(iSym,iShlAB,1) + 1
                  i2 = i1 + nnBstRSh(iSym,iShlAB,1) - 1
                  Do i = i1,i2
                     DiaMax(iSym)=max(DiaMax(iSym),Diag(i))
                  End Do
               End If
            End Do
            DiaMaxT(iSym)=DiaMax(iSym)
         End Do
      Else If (iLoc.eq.2 .or. iLoc.eq.3) Then
         Do iSym = 1,nSym
            DiaMax(iSym) = 0.0d0
            Do iShlAB = 1,nnShl
               Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
               If (iAtomShl(iShlA) .eq. iAtomShl(iShlB)) Then
                  i1 = iiBstR(iSym,iLoc) + iiBstRSh(iSym,iShlAB,iLoc)
     &               + 1
                  i2 = i1 + nnBstRSh(iSym,iShlAB,iLoc) - 1
                  Do i = i1,i2
                     DiaMax(iSym)=max(DiaMax(iSym),Diag(IndRed(i,iLoc)))
                  End Do
               End If
            End Do
            DiaMaxT(iSym) = 0.0d0
            Do iShlAB = 1,nnShl
               Call Cho_InvPck(iSP2F(iShlAB),iShlA,iShlB,.true.)
               If (iAtomShl(iShlA) .eq. iAtomShl(iShlB)) Then
                  i1 = iiBstR(iSym,1) + iiBstRSh(iSym,iShlAB,1) + 1
                  i2 = i1 + nnBstRSh(iSym,iShlAB,1) - 1
                  Do i = i1,i2
                     DiaMaxT(iSym)=max(DiaMaxT(iSym),Diag(i))
                  End Do
               End If
            End Do
         End Do
      Else
         Write(LuPri,*) SecNam,': unknown reduced set, iLoc = ',iLoc
         Call Cho_Quit('Unknown reduced set in '//SecNam,104)
      End If

      DGMax = DiaMax(1)
      Do iSym = 2,nSym
         DGMax = max(DGMax,DiaMax(iSym))
      End Do

      If (LocDbg) Then
         Write(LuPri,*) SecNam,': in reduced set ',iLoc,':'
         Write(LuPri,*) 'DiaMax  = ',(DiaMax(iSym),iSym=1,nSym)
         Write(LuPri,*) 'DiaMaxT = ',(DiaMaxT(iSym),iSym=1,nSym)
         Write(LuPri,*) 'DGMax   = ',DGMax
      End If

      End
