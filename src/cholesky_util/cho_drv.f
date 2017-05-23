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
* Copyright (C) 2003-2010, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine Cho_Drv(iReturn)
C
C     Thomas Bondo Pedersen, April 2010.
C
C     Purpose: driver for the Cholesky decomposition of two-electron
C              integrals. On entry, the integral program must have been
C              initialized and the relevant index info (#irreps, basis
C              functions, shells, etc.) must have been set up.
C
      Implicit None
      Integer iReturn
#include "cholesky.fh"

      iReturn=0
      If (Cho_DecAlg .eq. 5) Then ! parallel two-step algorithm
         Call Cho_Drv_ParTwoStep(iReturn)
      Else
         Call Cho_Drv_(iReturn)
      End If

      End
C
C=======================================================================
C
      SUBROUTINE CHO_DRV_(IRETURN)
C
C     Thomas Bondo Pedersen, 2003-2010.
C
C     Purpose: driver for the Cholesky decomposition of two-electron
C              integrals. On entry, the integral program must have been
C              initialized and the relevant index info (#irreps, basis
C              functions, shells, etc.) must have been set up.
C
C     NOTE: this is the "old" version prior to the parallel two-step
C           algorithm. It is still used for the "old" algorithms.
C
C     Return codes, IRETURN:
C
C        0 -- successful execution
C        1 -- decomposition failed
C        2 -- memory has been out of bounds
C
#include "implicit.fh"
#include "cho_para_info.fh"
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "chosubscr.fh"
#include "WrkSpc.fh"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
      LOGICAL SKIP_PRESCREEN, ALLOC_BKM
      PARAMETER (SKIP_PRESCREEN=.FALSE., ALLOC_BKM=.TRUE.)

      LOGICAL LCONV

      CHARACTER*8 SECNAM
      PARAMETER (SECNAM = 'CHO_DRV_')

      PARAMETER (DUMTST = 0.123456789D0, DUMTOL = 1.0D-15)

      CALL QENTER('_DRV_')
#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [ENTER]')
#endif

C     Start overall timing.
C     ---------------------

      IF (IPRINT .GE. INF_TIMING) CALL CHO_TIMER(TCPU0,TWALL0)

C     Set return code.
C     ----------------

      IRETURN = 0

C     Make a dummy allocation.
C     ------------------------

      l_START = 1
      CALL CHO_MEM('DRVDUM','ALLO','REAL',ip_START,l_START)
      WORK(ip_START) = DUMTST

C     INITIALIZATION.
C     ===============

      ISEC = 1
      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
      END IF
      CALL CHO_INIT(SKIP_PRESCREEN,ALLOC_BKM)
      CALL CHO_GASYNC()
      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
         CALL CHO_PRTTIM('Cholesky initialization',
     &                   TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                   TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                   1)
      END IF
#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_INIT]')
      IRC = 0
      CALL CHO_DUMP(IRC,LUPRI)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_DUMP returned ',IRC
         CALL CHO_QUIT('[1] Error detected in CHO_DUMP',103)
      END IF
#endif

C     GET DIAGONAL.
C     =============

      ISEC = 2
      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
         WRITE(LUPRI,'(/,A)')
     &   '***** Starting Cholesky diagonal setup *****'
         CALL CHO_FLUSH(LUPRI)
      END IF
      CALL CHO_GETDIAG(KDIAG,LCONV)
      CALL CHO_GASYNC()
      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
         CALL CHO_PRTTIM('Cholesky diagonal setup',
     &                   TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                   TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                   1)
      END IF
#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_GETDIAG]')
      IRC = 0
      CALL CHO_DUMP(IRC,LUPRI)
      IF (IRC .NE. 0) THEN
         WRITE(LUPRI,*) SECNAM,': CHO_DUMP returned ',IRC
         CALL CHO_QUIT('[2] Error detected in CHO_DUMP',103)
      END IF
      CALL CHO_PRINTLB() ! print vector dimension on each node
#endif

C     DECOMPOSITION.
C     ==============

      ISEC = 3
      IF (LCONV) THEN
         IF (RSTCHO) THEN
            WRITE(LUPRI,'(//,10X,A,A,A,//)')
     &      '***** ',SECNAM,': restarted calculation converged. *****'
            CALL CHO_DZERO(TIMSEC(1,ISEC),4)
         ELSE
            WRITE(LUPRI,'(A,A)')
     &      SECNAM,': logical error: converged but not restart?!?!'
            CALL CHO_QUIT('Error in '//SECNAM,103)
         END IF
      ELSE
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
            WRITE(LUPRI,'(/,A)')
     &      '***** Starting Cholesky decomposition *****'
            CALL CHO_FLUSH(LUPRI)
         END IF
         CALL CHO_P_SETADDR()
         IF (CHO_SSCREEN) THEN
            CALL CHO_SUBSCR_INIT()
         END IF
         CALL CHO_DECDRV(WORK(KDIAG))
         CALL CHO_GASYNC()
         IF (CHO_DECALG.EQ.2) THEN
            ! generate vectors from map
            CALL CHO_P_OPENVR(2) ! close files
            CALL CHO_P_OPENVR(1) ! re-open (problem on dec-alpha)
            IF (IPRINT .GE. INF_TIMING) THEN
               CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
               CALL CHO_PRTTIM('Cholesky map generation',
     &                         TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                         TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                         2)
            END IF
            IRC = 0
            CALL CHO_X_GENVEC(IRC,WORK(KDIAG))
            CALL CHO_GASYNC()
            IF (IRC .NE. 0) THEN
               WRITE(LUPRI,'(A,A)')
     &         SECNAM,': decomposition failed!'
               WRITE(LUPRI,'(A,A,I9)')
     &         SECNAM,': CHO_X_GENVEC returned ',IRC
               IRETURN = 1
               CALL CHO_QUIT('Error',104)
            END IF
            IF (IPRINT .GE. INF_TIMING) THEN
               CALL CHO_TIMER(TC,TW)
               CALL CHO_PRTTIM('Cholesky vector generation',
     &                         TC,TIMSEC(2,ISEC),
     &                         TW,TIMSEC(4,ISEC),
     &                         2)
            END IF
         END IF
         IF (CHO_SSCREEN) THEN
            CALL CHO_SUBSCR_FINAL()
         END IF
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
            CALL CHO_PRTTIM('Cholesky decomposition',
     &                      TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                      TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                      1)
         END IF
      END IF
#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER DECOMPOSITION]')
#endif

C     CHECK DIAGONAL.
C     ===============

      ISEC = 4
      IF (LCONV) THEN
         CALL CHO_DZERO(TIMSEC(1,ISEC),4)
      ELSE
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
            WRITE(LUPRI,'(/,A)')
     &      '***** Starting Cholesky diagonal check *****'
            CALL CHO_FLUSH(LUPRI)
         END IF
         CALL CHO_MEM('Cho.Rst','MAX ','REAL',KWRK,LWRK)
         CALL CHO_RESTART(WORK(KDIAG),WORK(KWRK),LWRK,.TRUE.,LCONV)
         CALL CHO_GASYNC()
         CALL CHO_MEM('Cho.Rst','FREE','REAL',KWRK,LWRK)
         IF (.NOT. LCONV) THEN
            WRITE(LUPRI,'(A,A)')
     &      SECNAM,': Decomposition failed!'
            IRETURN = 1
            CALL CHO_QUIT('Decomposition failed!',104)
         END IF
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
            CALL CHO_PRTTIM('Cholesky diagonal check',
     &                      TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                      TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                      1)
         END IF
#if defined (_DEBUG_)
         CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER DIAGONAL CHECK]')
#endif
      END IF

C     PARALLEL RUNS: WRITE GLOBAL DIAGONAL TO DISK.
C     =============================================

      CALL CHO_P_WRDIAG()

C     CHECK INTEGRALS.
C     ================

      ISEC = 5
      IF (CHO_INTCHK) THEN
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
            WRITE(LUPRI,'(/,A)')
     &      '***** Starting Cholesky integral check *****'
            CALL CHO_FLUSH(LUPRI)
         END IF
         CALL CHO_DBGINT()
         CALL CHO_GASYNC()
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
            CALL CHO_PRTTIM('Cholesky integral check',
     &                      TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                      TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                      1)
         END IF
#if defined (_DEBUG_)
         CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER INTEGRAL CHECK]')
#endif
      ELSE
         CALL CHO_DZERO(TIMSEC(1,ISEC),4)
      END IF

C     REORDER VECTORS.
C     ================

      ISEC = 6
      IF (CHO_REORD) THEN
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
            WRITE(LUPRI,'(/,A)')
     &      '***** Starting vector reordering *****'
            CALL CHO_FLUSH(LUPRI)
         END IF
         LIRS1F = NNBSTRT(1)*3
         CALL CHO_MEM('rs1tof','ALLO','INTE',KIRS1F,LIRS1F)
         CALL CHO_MEM('mx.reo','MAX ','REAL',KWRK,LWRK)
         CALL CHO_REOVEC(IWORK(KIRS1F),3,NNBSTRT(1),WORK(KWRK),LWRK)
         CALL CHO_GASYNC()
         CALL CHO_MEM('mx.reo','FREE','REAL',KWRK,LWRK)
         CALL CHO_MEM('rs1tof','FREE','INTE',KIRS1F,LIRS1F)
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
            CALL CHO_PRTTIM('Vector reordering',
     &                      TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                      TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                      1)
         END IF
#if defined (_DEBUG_)
         CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER VECTOR REORDERING]')
#endif
      ELSE
         CALL CHO_DZERO(TIMSEC(1,ISEC),4)
      END IF

C     FAKE PARALLEL: DISTRIBUTE VECTORS.
C     Note: after this section, InfVec(*,3,*) is changed
C     => vector disk addresses are screwed up!!
C     ==================================================

      ISEC = 7
      IF (CHO_FAKE_PAR .AND. NPROCS.GT.1 .AND. Is_Real_Par()) THEN
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
            WRITE(LUPRI,'(/,A)')
     &      '***** Starting vector distribution *****'
            CALL CHO_FLUSH(LUPRI)
         END IF
         CALL CHO_PFAKE_VDIST()
         CALL CHO_P_WRRSTC(XNPASS)
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
            CALL CHO_PRTTIM('Vector distribution',
     &                      TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                      TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                      1)
         END IF
#if defined (_DEBUG_)
         CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_PFAKE_VDIST]')
#endif
      ELSE
         CALL CHO_DZERO(TIMSEC(1,ISEC),4)
      END IF

C     FINALIZATIONS.
C     ==============

      ISEC = 8
      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
         WRITE(LUPRI,'(/,A)')
     &   '***** Starting Cholesky finalization *****'
         CALL CHO_FLUSH(LUPRI)
      END IF
      CALL CHO_TRCIDL_FINAL()
      CALL CHO_FINAL(.True.)
      CALL CHO_GASYNC()
      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
         CALL CHO_PRTTIM('Cholesky finalization',
     &                   TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                   TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                   1)
      END IF
#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_FINAL]')
#endif

C     STATISTICS.
C     ===========

      IF (IPRINT .GE. 1) THEN
         ISEC = 9
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(1,ISEC),TIMSEC(3,ISEC))
            WRITE(LUPRI,'(/,A)')
     &      '***** Starting Cholesky statistics *****'
            CALL CHO_FLUSH(LUPRI)
         END IF
         CALL CHO_P_STAT()
         CALL CHO_GASYNC()
         IF (IPRINT .GE. INF_TIMING) THEN
            CALL CHO_TIMER(TIMSEC(2,ISEC),TIMSEC(4,ISEC))
            CALL CHO_PRTTIM('Cholesky statistics',
     &                      TIMSEC(2,ISEC),TIMSEC(1,ISEC),
     &                      TIMSEC(4,ISEC),TIMSEC(3,ISEC),
     &                      1)
         END IF
      END IF
#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [AFTER CHO_STAT]')
#endif

C     Close vector and reduced storage files as well as restart files.
C     Deallocate all memory (using flush) and test bound.
C     Print total timing.
C     ----------------------------------------------------------------

      CALL CHO_P_OPENVR(2)

      TST = DUMTST - WORK(ip_START)
      IF (ABS(TST) .GT. DUMTOL) THEN
         WRITE(LUPRI,*) SECNAM,': memory has been out of bounds!!!'
         CALL CHO_FLUSH(LUPRI)
         IRETURN = 2
      END IF
      CALL CHO_MEM('DRVDUM','FLUS','REAL',ip_START,l_START)

      IF (IPRINT .GE. INF_TIMING) THEN
         CALL CHO_TIMER(TCPU1,TWALL1)
         CALL CHO_PRTTIM('Cholesky procedure',TCPU1,TCPU0,
     &                   TWALL1,TWALL0,1)
      END IF

#if defined (_DEBUG_)
      CALL CHO_PRTMAXMEM('CHO_DRV_ [EXIT]')
#endif
      CALL QEXIT('_DRV_')

      END
