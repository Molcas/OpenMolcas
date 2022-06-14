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
      SUBROUTINE CHO_GETDIAG(LCONV)
C
C     Purpose: get diagonal in first reduced set. On exit, DIAG
C              points to the diagonal and flag LCONV tells
C              if the diagonal is converged.
C
      use ChoArr, only: iSP2F, MySP, n_MySP, iSimRI
      use ChoSwp, only: IndRSh, IndRSh_Hidden
      use ChoSwp, only: IndRed, IndRed_Hidden
      use ChoSwp, only: Diag, Diag_Hidden
#include "implicit.fh"
      LOGICAL LCONV
#include "cholesky.fh"
#include "choprint.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      CHARACTER*11 SECNAM
      PARAMETER (SECNAM = 'CHO_GETDIAG')

      LOGICAL LOCDBG, DODUMMY, SYNC
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYLST(8)

      Integer, Allocatable:: KIBUF(:)
      Real*8, Allocatable:: KBUF(:), KSCR(:), KWRK(:)

      IF (RSTDIA) THEN

C        Set mySP list.
C        Always the trivial list for serial runs (restart not possible
C        in parallel runs).
C        -------------------------------------------------------------

         N_MYSP = NNSHL
         l_MySP=0
         If (Allocated(MYSP)) l_MySP=SIZE(MySP)
         IF (l_MYSP .EQ. NNSHL) THEN
            DO ISP = 1,N_MYSP
               MYSP(ISP) = ISP
            END DO
         ELSE
            CALL CHO_QUIT('MYSP allocation error in '//SECNAM,101)
         END IF

C        Read index array NNBSTRSH and set IIBSTRSH etc.
C        -----------------------------------------------

         CALL CHO_RSTD_GETIND1()

C        Allocate mapping arrays between reduced sets.
C        ---------------------------------------------

         MMBSTRT  = NNBSTRT(1)

         Call mma_allocate(IndRed_Hidden,NNBSTRT(1),3,
     &                     Label='IndRed_Hidden')
         IndRed => IndRed_Hidden

         Call mma_allocate(IndRSh_Hidden,NNBSTRT(1),
     &                     Label='IndRSh_Hidden')
         IndRSh => IndRSh_Hidden

C        Read mapping arrays.
C        --------------------

         CALL CHO_RSTD_GETIND2()

C        Check reduced to full shell pair mapping with the one on disk.
C        --------------------------------------------------------------

         NERR = -1
         CALL CHO_RSTD_CHKSP2F(iSP2F,SIZE(iSP2F),NERR)
         IF (NERR .NE. 0) THEN
            WRITE(LUPRI,*) SECNAM,': ',NERR,' errors detected in ',
     &                     'reduced-to-full shell pair mapping!'
            CALL CHO_QUIT('SP2F error in '//SECNAM,102)
         END IF

C        Allocation: diagonal.
C        ---------------------

         NEEDR = 1
         NEEDI = 4*NEEDR

         CALL mma_allocate(Diag_Hidden,NNBSTRT(1),Label='Diag_Hidden')
         Call mma_allocate(KBUF,NEEDR,Label='KBUF')
         Call mma_allocate(KIBUF,NEEDI,Label='KIBUF')

C        Read diagonal.
C        --------------

         CALL CHO_GETDIAG1(Diag_Hidden,KBUF,KIBUF,NEEDR,NDUMP)

         Call mma_deallocate(KIBUF)
         Call mma_deallocate(KBUF)

      ELSE

C        Calculate diagonal and get 1st reduced set.
C        -------------------------------------------

         Call mma_maxDBLE(LMAX)
         LMAX = LMAX/2 - MX2SH
         IF (LMAX .LT. 5*LBUF) THEN
            LBUF = MAX(LMAX/5,1)
         END IF

         LSCR  = MX2SH
         NEEDR = LBUF
         NEEDI = 4*LBUF

         Call mma_allocate(KBUF,NEEDR,Label='KBUF')
         Call mma_allocate(KSCR,LSCR ,Label='KSCR')
         Call mma_allocate(KIBUF,NEEDI,Label='KIBUF')

         NDUMP = 0

         CALL CHO_CALCDIAG(KBUF,KIBUF,LBUF,KSCR,LSCR,NDUMP)

         Call mma_deallocate(KIBUF)
         Call mma_deallocate(KBUF)
         Call mma_deallocate(KSCR)

C        Allocate diagonal and mapping array between reduced sets.
C        Reallocate buffer.
C        ---------------------------------------------------------

         MMBSTRT  = NNBSTRT(1)
         Call mma_allocate(IndRed_Hidden,NNBSTRT(1),3,
     &                     Label='IndRed_Hidden')
         IndRed => IndRed_Hidden
         Call mma_allocate(IndRSh_Hidden,NNBSTRT(1),
     &                     Label='IndRSh_Hidden')
         IndRSh => IndRSh_Hidden
         CALL mma_allocate(Diag_Hidden,NNBSTRT(1),Label='Diag_Hidden')

         NEEDR = LBUF
         NEEDI = 4*LBUF
         Call mma_allocate(KBUF,NEEDR,Label='KBUF')
         Call mma_allocate(KIBUF,NEEDI,Label='KIBUF')

C        Get diagonal in first reduced set.
C        ----------------------------------

         CALL CHO_GETDIAG1(Diag_Hidden,KBUF,KIBUF,LBUF,NDUMP)

C        Deallocate back to and including buffer.
C        ----------------------------------------

         Call mma_deallocate(KIBUF)
         Call mma_deallocate(KBUF)

      END IF

C     Set local and global info. On exit, Diag points to the local
C     diagonal.
C     -------------------------------------------------------------

      CALL CHO_P_SETGL()

C     Write local diagonal to disk.
C     -----------------------------

      IOPT = 1
      CALL CHO_IODIAG(Diag,IOPT)

C     Allocate memory for iscratch array for reading vectors.
C     -------------------------------------------------------

      DODUMMY = .NOT.(CHO_IOVEC.EQ.1 .OR. CHO_IOVEC.EQ.2 .OR.
     &                CHO_IOVEC.EQ.3 .OR. CHO_IOVEC.EQ.4 .OR.
     &                (FRAC_CHVBUF.GT.0.0D0 .AND. FRAC_CHVBUF.LT.1.0D0))
      CALL CHO_ALLO_ISCR(DODUMMY)

C     Initialize reduced set dimension(s) used for reading vectors.
C     -------------------------------------------------------------

      CALL CHO_INIRSDIM()

C     For RI simulation, zero 1-center diagonals smaller than
C     THR_SIMRI. Indices of zeroed diagonals are stored in ip_ISIMRI.
C     ---------------------------------------------------------------

      IF (CHO_SIMRI) THEN
         Call mma_allocate(iSimRI,NNBSTRT(1),Label='iSimRI')
         CALL CHO_SIMRI_Z1CDIA(Diag,THR_SIMRI,ISIMRI)
      END IF

C     Update diagonal if restart, else just do analysis.
C     --------------------------------------------------

      LCONV = .FALSE.
      IF (RSTCHO) THEN
         Call mma_maxDBLE(LWRK)
         Call mma_allocate(KWRK,LWRK,Label='KWRK')
         IF (LOCDBG) THEN
            WRITE(LUPRI,*) SECNAM,': restart diagonal:'
            DO ISYM = 1,NSYM
               ISYLST(ISYM) = ISYM
            END DO
            CALL CHO_PRTDIA(Diag,ISYLST,NSYM,1)
         END IF
         CALL CHO_RESTART(Diag,KWRK,LWRK,.FALSE.,LCONV)
         Call mma_deallocate(KWRK)
         IPRTRED = 2  ! print flag for cho_prtred
      ELSE
         IF (IPRINT .GE. INF_PASS) THEN
            BIN1 = 1.0D2
            STEP = 1.0D-1
            NBIN = 18
            SYNC = .FALSE.
            CALL CHO_P_ANADIA(Diag,SYNC,BIN1,STEP,NBIN,.TRUE.)
         END IF
         IPRTRED = 1  ! print flag for cho_prtred
      END IF
      IF (IPRINT .GE. INF_PASS) THEN
         CALL CHO_P_PRTRED(IPRTRED)
      END IF

      END
