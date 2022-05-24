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
      SUBROUTINE CHO_MCA_DRV()
C
C     Purpose: MOLCAS interface to Cholesky decomposition driver.
C
#include "implicit.fh"
#include "cholesky.fh"
#include "choptr2.fh"
#include "WrkSpc.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_MCA_DRV')

      LOGICAL INDEXATION, DOFOCK, DOGRAD
      LOGICAL VERBOSE, FREEK2

      CALL STATUSLINE('Seward: ','Cholesky decomposition of ERIs')

C     Initialize integral program (this does some memory
C     allocations; thus, DO NOT move this.
C     --------------------------------------------------

#if defined (_DEBUGPRINT_)
      CALL CHO_PRESCR(CUTINT1,THRINT1)
#endif

      Call Set_Basis_Mode('Valence')
      Call Setup_iSD()
      NSHELL     = -1
      INDEXATION = .TRUE.
      THRAO      = 0.0D0
      DOFOCK     = .FALSE.
      DOGRAD     = .FALSE.
      CALL SETUP_INTS(NSHELL,INDEXATION,THRAO,DOFOCK,DOGRAD)

#if defined (_DEBUGPRINT_)
      CALL CHO_PRESCR(CUTINT2,THRINT2)
      WRITE(LUPRI,*) SECNAM,': CutInt before Setup_Ints: ',CUTINT1
      WRITE(LUPRI,*) SECNAM,': CutInt after  Setup_Ints: ',CUTINT2
      WRITE(LUPRI,*) SECNAM,': ThrInt before Setup_Ints: ',THRINT1
      WRITE(LUPRI,*) SECNAM,': ThrInt after  Setup_Ints: ',THRINT2
      IF (CUTINT2.NE.CUTINT1 .OR. THRINT2.NE.THRINT1) THEN
         CALL CHO_QUIT('Initialization error in '//SECNAM,102)
      END IF
#endif

C     Start the Cholesky decomposition program.
C     -----------------------------------------

      ICODE = 0
      CALL CHO_DRV(ICODE)
      IF (ICODE .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': decomposition driver returned code ',
     &                  ICODE
         CALL CHO_QUIT('Decomposition failed!',104)
      END IF

C     Finalize integral program.
C     --------------------------

      VERBOSE = .FALSE.
      FREEK2  = .TRUE.
      CALL TERM_INTS(VERBOSE,FREEK2)

C     Halt execution if requested.
C     ----------------------------

      IF (HALTIT) THEN
         WRITE(LUPRI,*) SECNAM,': halting execution after ',
     &                  'decomposition as requested...'
         CALL GASYNC
         CALL CHO_QUIT('End of Test (in '//SECNAM//')',100)
      END IF

      CALL GASYNC
      Call Free_iSD()

      If (l_mySP .gt. 0) Then
         Call GetMem('mySP','Free','Inte',ip_mySP,l_mySP)
         l_mySP = 0
      End If
      Call Cho_X_dealloc(irc)

      END
