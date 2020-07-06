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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
************************************************************************
      SUBROUTINE TAB2(NREF,IOCR,nIOCR,L0,L1,L2,L3,INTNUM,LV,LSYM,
     *ICIALL,IFCORE,ICOR,NONE,JONE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IOCR(nIOCR),L0(*),L1(*),L2(*),L3(*),ICOR(*),
     *JONE(*)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      DIMENSION IORB(MXVERT)
      DIMENSION K00(MXVERT),K11(MXVERT),K22(MXVERT),K33(MXVERT)
      DIMENSION L00(MXVERT),L11(MXVERT),L22(MXVERT),L33(MXVERT)
*
      CALL QENTER('TAB2')
      IEL=2
      IF(IFIRST.NE.0)IEL=1
*
C     NUMBER OF ACTIVE ELECTRONS
      NAC=N-2*NIORB
*
C     UPPER LIMIT FOR NUMBER OF ELECTRONS IN ACTIVE SPACE
      NACU=NAC+IEL
*
      IUT=0
      IB(1)=INT(2*S)
      IA(1)=INT(N-2*S)/2
      IJ(LN+1)=0
      IJ(LN)=1
      NIJ=1
      IJR=1
      IJS=2
      IJRL=IJR
      IORB(1)=0
      DO 10 II=1,LN
         IIM=LN-II+1-LV
         IAC=IIM-NIORB-1
         IIM2=(IIM-1)*2
*
C        S=0
*
16       INUM=N-2*IA(IJR)-IB(IJR)
         IF(INTNUM.EQ.0.OR.IIM.LE.0)GO TO 109
         IF(IAC.GE.0)GO TO 409
         IF(IB(IJR).EQ.0)GO TO 109
         IF(INUM+IIM2.EQ.N-2)GO TO 11
         GO TO 109
409      IF(IB(IJR).GT.IAC+3)GO TO 11
109      NIJ=NIJ+1
         IA(NIJ)=IA(IJR)
         IB(NIJ)=IB(IJR)
         IORB(NIJ)=IORB(IJR)+2
         IF(IIM.LE.0)GO TO 11
         CALL CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
         IF(ISTOP.EQ.1)NIJ=NIJ-1
*
C        S=1
*
11       IF(IB(IJR).EQ.0)GO TO 12
         IF(INTNUM.EQ.0.OR.IIM.LE.0)GO TO 112
         IF(IAC.LT.0)GO TO 112
         IF(INUM+1.GT.NACU)GO TO 12
         IF(INUM+1.NE.NACU)GO TO 112
         IBS=IB(IJR)-1
         IF(IBS.NE.0.AND.IBS.NE.2)GO TO 12
112      NIJ=NIJ+1
         IA(NIJ)=IA(IJR)
         IB(NIJ)=IB(IJR)-1
         IORB(NIJ)=IORB(IJR)+1
         IF(IIM.LE.0)GO TO 12
         CALL CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
         IF(ISTOP.EQ.1)NIJ=NIJ-1
*
C        S=2
*
12       IF(IA(IJR).EQ.0)GO TO 13
         IF(INTNUM.EQ.0.OR.IIM.LE.0)GO TO 113
         IF(IAC.LT.0)GO TO 413
         IF(IB(IJR)+1.GT.IAC+3)GO TO 13
         IF(INUM+1.GT.NACU)GO TO 13
         IF(INUM+1.NE.NACU)GO TO 113
         IBS=IB(IJR)+1
         IF(IBS.NE.0.AND.IBS.NE.2)GO TO 13
         GO TO 113
413      IF(IB(IJR).GE.2)GO TO 13
113      NIJ=NIJ+1
         IA(NIJ)=IA(IJR)-1
         IB(NIJ)=IB(IJR)+1
         IORB(NIJ)=IORB(IJR)+1
         IF(IIM.LE.0)GO TO 13
         CALL CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
         IF(ISTOP.EQ.1)NIJ=NIJ-1
*
C        S=3
*
13       IF(IA(IJR).EQ.0)GO TO 14
         IF(INTNUM.EQ.0)GO TO 114
         IF(IAC.LT.0)GO TO 114
         IF(IB(IJR).GT.IAC+3)GO TO 14
         IF(INUM+2.GT.NACU)GO TO 14
         IF(INUM+2.NE.NACU)GO TO 114
         IBS=IB(IJR)
         IF(IBS.NE.0.AND.IBS.NE.2)GO TO 14
114      NIJ=NIJ+1
         IA(NIJ)=IA(IJR)-1
         IB(NIJ)=IB(IJR)
         IORB(NIJ)=IORB(IJR)
         IF(IIM.LE.0)GO TO 14
         CALL CHEL(IA(NIJ),IB(NIJ),IIM,IEL,ISTOP)
         IF(ISTOP.EQ.1)NIJ=NIJ-1
*
14       IF(IJR.EQ.IJRL)GO TO 15
         IJR=IJR+1
         GO TO 16
C        DELETE VERTICES
15       NIJ1=NIJ-1
         IN=IJS
         IUT=IJS
         IF(NIJ1.LT.IJS)GO TO 21
         DO 20 IJD=IJS,NIJ1
            JJ1=NIJ-IJD+IJS-1
            J=JJ1+1
            DO 25 K=IJS,JJ1
               IF(IA(J).NE.IA(K))GO TO 25
               IF(IB(J).NE.IB(K))GO TO 25
               GO TO 26
25          CONTINUE
            GO TO 20
26          IA(J)=-1
            IB(J)=-1
20       CONTINUE
C        PACK VERTICES
         IJS1=IJS+1
         DO 30 J=IJS1,NIJ
            IF(IA(J).NE.-1)GO TO 31
            IF(IB(J).NE.-1)GO TO 31
            IN=IN+1
            GO TO 30
31          IN=IN+1
            IUT=IUT+1
            IA(IUT)=IA(IN)
            IB(IUT)=IB(IN)
            IORB(IUT)=IORB(IN)
30       CONTINUE
C        ORDER VERTICES
         IUT1=IUT-1
         IF(IUT1.LT.IJS)GO TO 21
         DO 41 J=IJS,IUT1
            J11=J+1
            DO 42 K=J11,IUT
               IF (IA(J)-IA(K).LT.0) THEN
                  GO TO 43
               ELSE IF (IA(J)-IA(K).EQ.0) THEN
                  GO TO 44
               ELSE
                  GO TO 42
               END IF
44             IF(IB(J).GT.IB(K))GO TO 42
43             IAT=IA(J)
               IBT=IB(J)
               IA(J)=IA(K)
               IB(J)=IB(K)
               IA(K)=IAT
               IB(K)=IBT
42          CONTINUE
41       CONTINUE
21       IF(II.NE.LN)IJ(LN-II)=IUT
         IJR=IJS
         IJS=IUT+1
         IJRL=IUT
         NIJ=IUT
10    CONTINUE
      IF(N.NE.2)GO TO 35
      IUT=IUT+1
      IA(IUT)=0
      IB(IUT)=0
      IA(IUT-1)=0
      IB(IUT-1)=1
      IA(IUT-2)=0
      IB(IUT-2)=2
35    JJ2=0
      DO 40 II=1,LN
         I=LN-II+1-LV
         IIM2=(I-1)*2
         I0=I+LV
         JJ1=IJ(I0+1)+1
         JJ2=IJ(I0)
         J3=JJ2+1
         IF(I0.NE.1)J4=IJ(I0-1)
         IF(I0.EQ.1)J4=IUT
C        DETERMINE CASE DOWN
         DO 50 J=JJ1,JJ2
            IA1=IA(J)
            IB1=IB(J)
            INUM=2*IA1+IB1
            K00(J)=0
            K11(J)=0
            K22(J)=0
            K33(J)=0
            DO 60 JJ=J3,J4
               IF(IA1.EQ.IA(JJ))GO TO 61
               IF((IA1-IA(JJ)).NE.1)GO TO 60
               IF(IB1.EQ.IB(JJ))GO TO 62
               IF((IB(JJ)-IB1).NE.1)GO TO 60
               IF(I.GT.NIORB.OR.INTNUM.EQ.0)GO TO 59
               IF(I.LE.0)GO TO 59
               IF(IB1.GE.2)GO TO 60
59             K22(J)=JJ
               GO TO 60
62             K33(J)=JJ
               GO TO 60
61             IF(IB1.EQ.IB(JJ))GO TO 63
               IF((IB1-IB(JJ)).NE.1)GO TO 60
               K11(J)=JJ
               GO TO 60
63             IF(I.GT.NIORB.OR.INTNUM.EQ.0)GO TO 64
               IF(I.LE.0)GO TO 64
               IF(IB1.EQ.0)GO TO 64
               IF(INUM-IIM2.EQ.2)GO TO 60
               IF(IA1.LT.I-1)GO TO 60
64             K00(J)=JJ
60          CONTINUE
50       CONTINUE
C        DETERMINE CASE UP
         DO 70 J=J3,J4
            IA1=IA(J)
            IB1=IB(J)
            INUM=2*IA1+IB1
            L00(J)=0
            L11(J)=0
            L22(J)=0
            L33(J)=0
            DO 80 JJ=JJ1,JJ2
               IF(IA(JJ).EQ.IA1)GO TO 81
               IF((IA(JJ)-IA1).NE.1)GO TO 80
               IF(IB(JJ).EQ.IB1)GO TO 82
               IF((IB1-IB(JJ)).NE.1)GO TO 80
               IF(I.GT.NIORB.OR.INTNUM.EQ.0)GO TO 79
               IF(I.LE.0)GO TO 79
               IF(IB1.GE.3)GO TO 80
79             L22(J)=JJ
               GO TO 80
82             L33(J)=JJ
               GO TO 80
81             IF(IB(JJ).EQ.IB1)GO TO 83
               IF((IB(JJ)-IB1).NE.1)GO TO 80
               L11(J)=JJ
               GO TO 80
83             IF(I.GT.NIORB.OR.INTNUM.EQ.0)GO TO 84
               IF(I.LE.0)GO TO 84
               IF(IB1.EQ.0)GO TO 84
               IF(INUM-IIM2.EQ.2)GO TO 80
               IF(IA1.LT.I-1)GO TO 80
84             L00(J)=JJ
80          CONTINUE
70       CONTINUE
40    CONTINUE
      IV0=IUT
      IV1=IUT-1
      IV2=IUT-2
      IV3=IUT-3
      K00(IUT)=0
      K11(IUT)=0
      K22(IUT)=0
      K33(IUT)=0
      K00(IUT+1)=0
      K11(IUT+1)=0
      K22(IUT+1)=0
      K33(IUT+1)=0
      IF(ICIALL.NE.0)CALL CIALL(LSYM,NREF,IOCR,nIOCR,L00,L11,L22,L33,LV)
      CALL DELTAB(NREF,IOCR,L0,L1,L2,L3,INTNUM,LV,IFCORE,ICOR,
     *NONE,JONE,K00,K11,K22,K33,L00,L11,L22,L33)
      DO 103 I=1,ILIM
         ISTA=(I-1)*MXVERT
         IF(IPRINT.GE.5)WRITE(IW,101)
101      FORMAT(///,6X,'TAB2',//,8X,'J',8X,'A',3X,'B',7X,
     *   'K0',2X,'K1',2X,'K2',2X,'K3',2X,'L0',2X,'L1',2X,'L2',2X,'L3',/)
         IF(IPRINT.GE.5)WRITE(IW,100)(J,IA(J),IB(J),K0(ISTA+J),
     *   K1(ISTA+J),K2(ISTA+J),K3(ISTA+J),L0(ISTA+J),L1(ISTA+J),
     *   L2(ISTA+J),L3(ISTA+J),J=1,IUT)
100   FORMAT(6X,I3,5X,2I4,5X,8I4)
103   CONTINUE
      IBMAX=0
      DO 450 J=1,IUT
         IF(IB(J).GT.IBMAX)IBMAX=IB(J)
450   CONTINUE
      IUT1=IUT-1
      DO 511 IL=1,ILIM
         ISTA=(IL-1)*MXVERT
         IX(ISTA+1)=1
         DO 410 II=1,LN
            IF(II.EQ.LN)GO TO 420
            I=LN-II
            IJL=IJ(I)
            IJS=IJ(I+1)+1
            GO TO 430
420         IJL=IUT
            IJS=IUT-3
            IF(IFIRST.NE.0)IJS=IUT-1
430         DO 440 J=IJS,IJL
            ISUM=0
            IF(L0(ISTA+J).EQ.0)GO TO 441
            ISUM=ISUM+IX(ISTA+L0(ISTA+J))
441         IF(L1(ISTA+J).EQ.0)GO TO 442
            IY(ISTA+L1(ISTA+J),1)=ISUM
            ISUM=ISUM+IX(ISTA+L1(ISTA+J))
442         IF(L2(ISTA+J).EQ.0)GO TO 443
            IY(ISTA+L2(ISTA+J),2)=ISUM
            ISUM=ISUM+IX(ISTA+L2(ISTA+J))
443         IF(L3(ISTA+J).EQ.0)GO TO 444
            IY(ISTA+L3(ISTA+J),3)=ISUM
            ISUM=ISUM+IX(ISTA+L3(ISTA+J))
444         IX(ISTA+J)=ISUM
440         CONTINUE
410      CONTINUE
511   CONTINUE
      DO 104 I=1,ILIM
         ISTA=(I-1)*MXVERT
         IF(IPRINT.GE.5)WRITE(IW,102)
102      FORMAT(///,6X,'INDEX TABLE',//,8X,'J',8X,'Y1',3X,'Y2',
     *   3X,'Y3',9X,'X',/)
         IF(IPRINT.GE.5)WRITE(IW,200)(J,IY(ISTA+J,1),IY(ISTA+J,2),
     *   IY(ISTA+J,3),IX(ISTA+J),J=1,IUT)
200      FORMAT(6X,I3,5X,3I5,5X,I5)
104   CONTINUE
      IF(IPRINT.GE.2)WRITE(IW,210)IUT
210   FORMAT(/,6X,'NUMBER OF VERTICES',I10)
      WRITE(IW,214)
214   FORMAT(///,6X,'INTERNAL CONFIGURATIONS (FORMAL)')

      IF (IFIRST.EQ.0) THEN
* This is a normal calculation, both singles and doubles included.
      WRITE(IW,215)(IX(IUT+1-ITTT+MXVERT*(ITTT-1)),ITTT=1,4)
215   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,
     */,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,
     */,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
      ELSE
* ''FIRST'' keyword has been given. Then this is just a so-called
* first-order CI, i.e., singles only.
      WRITE(IW,215)(IX(IUT+1-ITTT+MXVERT*(ITTT-1)),ITTT=1,2)
c216   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
c     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
      END IF

      IRC(1)=IX(IUT)
      DO 8 I=2,ILIM
         ISTA=(I-1)*MXVERT
         IRC(I)=IX(ISTA+IUT+1-I)+IRC(I-1)
8     CONTINUE
      ISUM=IRC(ILIM)
      If (ISUM.GT.LIX) Then
         Write(6,*) 'Tab2: ISUM.GT.LIX'
         Write(6,*) 'ISUM,LIX=',ISUM,LIX
         Call QTrace
         Call Abend
      End If
      DO 5 I=1,10
         FBB=I-1
         FB=FBB/(FBB+1)
         BS1(I)=SQRT(FB)
         FB=(FBB+2)/(FBB+1)
         BS2(I)=SQRT(FB)
         IF(I.GT.1)BS3(I)=D1/BS1(I)
         BS4(I)=D1/BS2(I)
         FB=FBB*FBB-1
         IF(I.GT.1)BL1(I)=SQRT(FB)/FBB
         FB=(FBB+2)**2-1
         BL2(I)=SQRT(FB)/(FBB+2)
5     CONTINUE
C     PUT ZEROS IN VECTORS
      DO 6 I=1,LN
         COUP(I)=D0
         COUP1(I)=D0
6     CONTINUE
      IN=0
      DO 305 I=1,LN
         II=LN-I+1
         IJS=IJ(II+1)+1
         IJL=IJ(II)
         IJFS=IJF(II+1)+1
         IJFL=IJF(II)
         DO 320 IIJ=IJS,IJL
            DO 330 IIJF=IJFS,IJFL
               IF(IA(IIJ).NE.IAF(IIJF))GO TO 330
               IF(IB(IIJ).NE.IBF(IIJF))GO TO 330
               IPO(IIJ)=IIJF
               GO TO 320
330         CONTINUE
320      CONTINUE
305   CONTINUE
      IPO(IUT)=IJF(1)+1
      IF(IPRINT.GE.10)WRITE(IW,350)(IPO(J),J=1,IUT)
350   FORMAT(/,6X,5I5)
      JMAX=0
      LN1=LN+1
      DO 400 I=2,LN1
         I1=I-1
         JL=IJ(I1)-IJ(I)
         IF(JL.LT.JMAX)GO TO 400
         JMAX=JL
400   CONTINUE
*
      IF (IPRINT.GE.2) THEN
         WRITE(IW,411)JMAX
         WRITE(IW,412)IBMAX, MAXB
411      FORMAT(/,6X,'NUMBER OF VERTICES IN ONE ROW',I6,
     &     6X,'(PRESENT LIMIT 31)')
         IF(JMAX.GT.31) Call Abend
412      FORMAT(6X,'MAXIMUM B VALUE',I20,6X,'(PRESENT LIMIT ',I5,')')
      END IF
      IF(IBMAX.GT.MAXB) Call Abend
*
      CALL QEXIT('TAB2')
      RETURN
      END
