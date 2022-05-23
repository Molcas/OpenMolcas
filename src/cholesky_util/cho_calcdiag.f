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
      SUBROUTINE CHO_CALCDIAG(BUF,IBUF,LENBUF,SCR,LENSCR,
     &                        IIBSTRSH,NNBSTRSH,MSYM,MMSHL,NDUMP)
C
C     Purpose: shell-driven calculation of the integral diagonal and
C              setup of the first reduced set.
C
#include "implicit.fh"
      DIMENSION BUF(LENBUF), SCR(LENSCR)
      INTEGER   IBUF(4,LENBUF)
      INTEGER   IIBSTRSH(MSYM,MMSHL,3), NNBSTRSH(MSYM,MMSHL,3)
#include "cholesky.fh"
#include "choprint.fh"
#include "choorb.fh"
#include "choptr.fh"
#include "choptr2.fh"
#include "WrkSpc.fh"

      CHARACTER*12 SECNAM
      PARAMETER (SECNAM = 'CHO_CALCDIAG')

      PARAMETER (INFO_DEBUG = 4, INFO_INSANE = 10)

      DIMENSION SCRMAX(8), XNCD(1)

      INTEGER  CHO_ISAOSH
      EXTERNAL CHO_ISAOSH

      MULD2H(I,J)=IEOR(I-1,J-1)+1
      ITRI(I,J)=MAX(I,J)*(MAX(I,J)-3)/2+I+J
      IBASSH(I,J)=IWORK(ip_IBASSH-1+NSYM*(J-1)+I)
      NBASSH(I,J)=IWORK(ip_NBASSH-1+NSYM*(J-1)+I)
      NBSTSH(I)=IWORK(ip_NBSTSH-1+I)
      IATOMSHL(I)=IWORK(ip_IATOMSHL-1+I)
      ISP2F(I)=IWORK(ip_iSP2F-1+I)

C     Check dimensions.
C     -----------------

      IF (LENBUF .LT. LBUF) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,': LENBUF >= LBUF required!'
         WRITE(LUPRI,'(1X,A,I10)')    'LENBUF = ',LENBUF
         WRITE(LUPRI,'(1X,A,I10,/)')  'LBUF   = ',LBUF
         CALL CHO_QUIT('Buffer error in '//SECNAM,102)
      END IF
      LSCR = MX2SH
      IF (LENSCR .LT. LSCR) THEN
         WRITE(LUPRI,'(//,1X,A,A)') SECNAM,
     &   ': LENSCR >= MX2SH required!'
         WRITE(LUPRI,'(1X,A,I10)')    'LENSCR = ',LENSCR
         WRITE(LUPRI,'(1X,A,I10,/)')  'MX2SH  = ',LSCR
         CALL CHO_QUIT('Scratch space error in '//SECNAM,102)
      END IF

C     Open scratch files.
C     -------------------

      IUNIT = -1
      CALL CHO_OPEN(IUNIT,'_CHO_DIASCR2')
      JUNIT = -1
      CALL CHO_OPEN(JUNIT,'_CHO_DIASCR1')
      REWIND(IUNIT)
      REWIND(JUNIT)

C     Make JUNIT available outside this routine.
C     ------------------------------------------

      LUSCR = JUNIT

C     Initialize abs. max. diag. array.
C     ---------------------------------

      CALL CHO_DZERO(DIAMAX,NSYM)

C     Allocate array for storing 10 most negative diagonals
C     (there should be none, of course, but they do show up)
C     ------------------------------------------------------

      l_NegCalcDiag=10
      Call GetMem('NCD','Allo','Real',ip_NegCalcDiag,l_NegCalcDiag)
      Call Cho_dZero(Work(ip_NegCalcDiag),l_NegCalcDiag)
      n_NegCalcDiag=0

C     Calculate diagonal in loop over shell-pairs.
C     CHO_NO2CENTER on: skip all 2-center diagonals.
C     ----------------------------------------------

      IF (CHO_NO2CENTER .AND. l_IATOMSHL.LT.NSHELL) THEN
         CALL CHO_QUIT(SECNAM//': iAtomShl not allocated correctly!',
     &                 103)
      END IF

      XLDIAG = 0.0D0
      ICOUNT = 0
      NDUMP  = 0
      N_MYSP = 0
      IOPT = 2
      CALL CHO_P_DISTRIB_SP(IOPT,IWORK(ip_MYSP),N_MYSP)
      CALL GETMEM('MAXIDA','MAX','REAL',KINTD,LINTD)
      CALL XSETMEM_INTS(LINTD) ! set memory for seward
      DO I_MYSP = 0,N_MYSP-1

         ISAB = IWORK(ip_MYSP+I_MYSP)

         ISHLAB = ISP2F(ISAB)
         CALL CHO_INVPCK(ISHLAB,ISHLA,ISHLB,.TRUE.)

         IF (CHO_NO2CENTER) THEN
            IF (IATOMSHL(ISHLA) .NE. IATOMSHL(ISHLB)) THEN
               GO TO 1 ! cycle loop
            END IF
         END IF

         NUMA = NBSTSH(ISHLA)
         NUMB = NBSTSH(ISHLB)
         IF (ISHLA .EQ. ISHLB) THEN
            NUMAB = NUMA*(NUMA + 1)/2
         ELSE
            NUMAB = NUMA*NUMB
         END IF

         SHA = ISHLA
         SHB = ISHLB

         CALL CHO_MCA_DIAGINT(ISHLA,ISHLB,SCR,NUMAB)

         IF (IPRINT .GE. INFO_INSANE) THEN
            IF ((ISHLA.EQ.1) .AND. (ISHLB.EQ.1)) THEN
               IF (CHO_PRESCREEN) THEN
                  CALL CHO_HEAD(SECNAM//': Prescreened Diagonal',
     &                          '=',80,LUPRI)
               ELSE
                  CALL CHO_HEAD(SECNAM//': Unscreened Diagonal',
     &                          '=',80,LUPRI)
               END IF
            END IF
            WRITE(LUPRI,'(/,2X,A,I10,1X,I10,1X,I10)')
     &      'Diagonal shell block A,B,AB = ',ISHLA,ISHLB,
     &                                       ITRI(ISHLA,ISHLB)
            IF (ISHLA .EQ. ISHLB) THEN
               CALL CHO_OUTPAK(SCR,NUMA,1,LUPRI)
            ELSE
               CALL CHO_OUTPUT(SCR,1,NUMA,1,NUMB,NUMA,NUMB,1,LUPRI)
            END IF
         END IF

         IF (ISHLA .EQ. ISHLB) THEN
            DO IA = 1,NUMA
               ISYMA  = CHO_ISAOSH(IA,ISHLA)
               DO IB = 1,IA
                  ISYMB  = CHO_ISAOSH(IB,ISHLB)
                  ISYMAB = MULD2H(ISYMB,ISYMA)
                  IAB    = ITRI(IA,IB)
                  DIAAB  = SCR(IAB)
                  IF (DIAAB .LT. 0.0D0) THEN
                     n_NegCalcDiag=n_NegCalcDiag+1
                     Call UpdateMostNegative(l_NegCalcDiag,
     *                                       Work(ip_NegCalcDiag),DIAAB)
                  END IF
                  IF (DIAAB .GT. THRDIAG) THEN
                     DIAMAX(ISYMAB) = MAX(DIAMAX(ISYMAB),DIAAB)
                     ICOUNT = ICOUNT + 1
                     BUF(ICOUNT)    = SCR(IAB)
                     IBUF(1,ICOUNT) = ISAB
                     IBUF(2,ICOUNT) = IAB
                     IBUF(3,ICOUNT) = ISYMAB
                     IBUF(4,ICOUNT) = IAB
                     IF (ICOUNT .EQ. LBUF) THEN
                        CALL CHO_WRBUF(LBUF,BUF,IBUF,LBUF,IUNIT)
                        XLDIAG = XLDIAG + DBLE(LBUF)
                        ICOUNT = 0
                        NDUMP  = NDUMP + 1
                     END IF
                  END IF
               END DO
            END DO
         ELSE
            DO ISYMB = 1,NSYM
               DO IBB = 1,NBASSH(ISYMB,ISHLB)
                  IB = IBASSH(ISYMB,ISHLB) + IBB
                  DO ISYMA = 1,NSYM
                     DO IAA = 1,NBASSH(ISYMA,ISHLA)
                        IA = IBASSH(ISYMA,ISHLA) + IAA
                        ISYMAB = MULD2H(ISYMA,ISYMB)
                        IAB    = NUMA*(IB - 1) + IA
                        DIAAB  = SCR(IAB)
                        IF (DIAAB .LT. 0.0D0) THEN
                           n_NegCalcDiag=n_NegCalcDiag+1
                           Call UpdateMostNegative(l_NegCalcDiag,
     *                                       Work(ip_NegCalcDiag),DIAAB)
                        END IF
                        IF (DIAAB .GT. THRDIAG) THEN
                           DIAMAX(ISYMAB) = MAX(DIAMAX(ISYMAB),DIAAB)
                           ICOUNT = ICOUNT + 1
                           BUF(ICOUNT)    = SCR(IAB)
                           IBUF(1,ICOUNT) = ISAB
                           IBUF(2,ICOUNT) = IAB
                           IBUF(3,ICOUNT) = ISYMAB
                           IBUF(4,ICOUNT) = IAB
                           IF (ICOUNT .EQ. LBUF) THEN
                              CALL CHO_WRBUF(LBUF,BUF,IBUF,LBUF,
     &                                       IUNIT)
                              XLDIAG = XLDIAG + DBLE(LBUF)
                              ICOUNT = 0
                              NDUMP  = NDUMP + 1
                           END IF
                        END IF
                     END DO
                  END DO
               END DO
            END DO
         END IF

    1    CONTINUE ! for cycling loop

      END DO
      IF (ICOUNT .GT. 0) THEN ! flush buffer
         IF (ICOUNT .GT. LBUF) THEN
            CALL CHO_QUIT('Logical error in '//SECNAM,103)
         END IF
         CALL CHO_WRBUF(ICOUNT,BUF,IBUF,LBUF,IUNIT)
         XLDIAG = XLDIAG + DBLE(ICOUNT)
         ICOUNT = 0
         NDUMP  = NDUMP + 1
      END IF
      CALL XRLSMEM_INTS ! release memory (seward)
      CALL CHO_GADGOP(DIAMAX,NSYM,'max') ! sync abs. max. diag.
      n_NegCalcDiag_local=n_NegCalcDiag
      XNCD(1)=dble(n_NegCalcDiag)
      CALL CHO_GADGOP(XNCD,1,'+')
      n_NegCalcDiag=int(XNCD(1))
      If (n_NegCalcDiag.gt.0) Then
         Call WarningMessage(1,
     *          'WARNING: negative integral diagonal elements computed')
         Write(LuPri,'(3X,A,A)')
     *   'All negative integral diagonal elements have been',
     *   ' removed (zeroed) - they are considered irrelevant!'
         Write(LuPri,'(3X,A,I10)')
     *   'Number of negative elements computed:   ',n_NegCalcDiag
         Write(LuPri,'(3X,A,I10)')
     *   'Number of negative elements (this node):',n_NegCalcDiag_local
         If (n_NegCalcDiag_local.gt.0) Then
            ll=min(n_NegCalcDiag_local,l_NegCalcDiag)
            Write(LuPri,'(I5,A)')
     *      ll,' most negative elements (this node):'
            Write(LuPri,'(1P,10D12.4)')
     *      (Work(ip_NegCalcDiag+i),i=0,ll-1)
         End If
         Call CHO_GADGOP(Work(ip_NegCalcDiag),1,'min')
         Write(LuPri,'(3X,A,1P,D12.4)')
     *   'Most negative element overall: ',Work(ip_NegCalcDiag)
      End If
      Call GetMem('NCD','Free','Real',ip_NegCalcDiag,l_NegCalcDiag)
      ip_NegCalcDiag=-99999999
      l_NegCalcDiag=0

      IF (IPRINT .GE. INFO_DEBUG) THEN
         CALL CHO_HEAD(SECNAM//': Diagonal Info','=',80,LUPRI)
         XXX   = DBLE(NBAST)
         XMDIA = XXX*(XXX + 1.0D0)/2.0D0
         XLDIA = XLDIAG
         SAVD  = 1.0D2*(XMDIA - XLDIA)/XMDIA
         WRITE(LUPRI,'(/,2X,A,1P,D15.6)')
     &   'Screening threshold for initial diagonal: ',THRDIAG
         WRITE(LUPRI,'(2X,A,F15.1,/,2X,A,F15.1)')
     &   'Dimension of unscreened initial diagonal: ',XMDIA,
     &   'Dimension of   screened initial diagonal: ',XLDIA
         WRITE(LUPRI,'(2X,A,7X,F8.3,A)')
     &   'Saving from screening                   : ',SAVD,'%'
         DO ISYM = 1,NSYM
            WRITE(LUPRI,'(2X,A,I2,12X,A,1P,D15.6)')
     &      'Maximum diagonal, symmetry',ISYM,': ',DIAMAX(ISYM)
         END DO
         WRITE(LUPRI,'(2X,A,5X,I10)')
     &   'Number of negative diagonals computed   : ',n_NegCalcDiag
      END IF

C     Read through the file to get first reduced set.
C     -----------------------------------------------

      CALL CHO_IZERO(NNBSTRSH(1,1,1),NSYM*NNSHL)

      REWIND(IUNIT)
      REWIND(JUNIT)
      IF (SCDIAG) THEN ! screen diagonal
         DEL1 = THRCOM*THRCOM/DAMP(1)
         DO ISYM = 1,NSYM
            IF (ABS(DIAMAX(ISYM)) .GT. 0.0D0) THEN
               SCRMAX(ISYM) = DEL1/DIAMAX(ISYM)
            ELSE
               SCRMAX(ISYM) = 1.0D15
            END IF
         END DO
         IF (CHO_USEABS) THEN
            DO IDUMP = 1,NDUMP
               CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
               IF (IDUMP .EQ. NDUMP) THEN
                  CALL CHO_CLOSE(IUNIT,'DELETE')
               END IF
               DO L = 1,LENGTH
                  DIAGAB = BUF(L)
                  ISYMAB = IBUF(3,L)
                  IF (ABS(DIAGAB) .LT. SCRMAX(ISYMAB)) THEN
                     BUF(L)    = 0.0D0
                     IBUF(2,L) = -1
                  ELSE
                     ISHLAB = IBUF(1,L)
                     NNBSTRSH(ISYMAB,ISHLAB,1) =
     &                                     NNBSTRSH(ISYMAB,ISHLAB,1) + 1
                     IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
                  END IF
               END DO
               CALL CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
            END DO
         ELSE
            DO IDUMP = 1,NDUMP
               CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
               IF (IDUMP .EQ. NDUMP) THEN
                  CALL CHO_CLOSE(IUNIT,'DELETE')
               END IF
               DO L = 1,LENGTH
                  DIAGAB = BUF(L)
                  ISYMAB = IBUF(3,L)
                  IF (DIAGAB .LT. SCRMAX(ISYMAB)) THEN
                     BUF(L)    = 0.0D0
                     IBUF(2,L) = -1
                  ELSE
                     ISHLAB = IBUF(1,L)
                     NNBSTRSH(ISYMAB,ISHLAB,1) =
     &                                     NNBSTRSH(ISYMAB,ISHLAB,1) + 1
                     IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
                  END IF
               END DO
               CALL CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
            END DO
         END IF
      ELSE ! no screening at all
         DO IDUMP = 1,NDUMP
            CALL CHO_RDBUF(LENGTH,BUF,IBUF,LBUF,IUNIT)
            IF (IDUMP .EQ. NDUMP) THEN
               CALL CHO_CLOSE(IUNIT,'DELETE')
            END IF
            DO L = 1,LENGTH
               ISHLAB = IBUF(1,L)
               ISYMAB = IBUF(3,L)
               NNBSTRSH(ISYMAB,ISHLAB,1) = NNBSTRSH(ISYMAB,ISHLAB,1) + 1
               IBUF(2,L) = NNBSTRSH(ISYMAB,ISHLAB,1)
            END DO
            CALL CHO_WRBUF(LENGTH,BUF,IBUF,LBUF,JUNIT)
         END DO
      END IF

      CALL CHO_GAIGOP(NNBSTRSH(1,1,1),NSYM*NNSHL,'+') ! sync
      CALL CHO_SETREDIND(IIBSTRSH,NNBSTRSH,NSYM,NNSHL,1)

      END
      Subroutine UpdateMostNegative(n,X,Val)
      Implicit None
      Integer n
      Real*8  X(n)
      Real*8  Val

      Integer i, j

      If (Val.ge.X(n)) Return
      i=0
      Do While (i.lt.n)
         i=i+1
         If (Val.lt.X(i)) Then
            Do j=n,i+1,-1
               X(j)=X(j-1)
            End Do
            X(i)=Val
            i=n+1 ! break while loop
         End If
      End Do

      End
