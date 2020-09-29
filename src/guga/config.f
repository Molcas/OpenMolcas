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
      SUBROUTINE CONFIG(NREF,IOCR,nIOCR,L0,L1,L2,L3,JSYM,JSY,INTNUM,
     &                  LSYM,JJS,ISO,LV,IFCORE,ICOR,NONE,JONE,JREFX,
     &                  NFREF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IOCR(nIOCR),L0(*),L1(*),L2(*),L3(*),JSYM(*),JSY(*),
     &          JJS(*),ISO(*),ICOR(*),JONE(*),JREFX(*)
#include "real_guga.fh"
#include "integ.fh"
      COMMON/CNSTS/D0,D1,D2
      COMMON/D/JNDX(500 000)
      DIMENSION IOC(55),ISP(55)
*
      JSYL=30000
      JSYLL=3000
      JRX=9000
      IBS=0
      IEL=2
      LSYM=1
      NFREF=0
      LNS=NIORB+LV+1
* Compute wave function symmetry
      IRR=0
      DO 5 I=LNS,LN
        IRR=IRR+1
        IF(IOCR(IRR).EQ.1)LSYM=MUL(LSYM,NSM(I))
5     CONTINUE
      WRITE(IW,'(6X,A,I3)') 'WAVE-FUNCTION SYMMETRY LABEL:', LSYM

* Initialize arrays JJS, JNDX, ICASE, JSY
      DO 3 I=1,18
        JJS(I)=0
3     CONTINUE
C CONSTRUCT JNDX
* ILIM=2 or ILIM=4 was set by input: Normally 4, but 2 if keyword FIRST
* has been given. ILIM is in integ.fh
      INTOT=IRC(ILIM)
      DO 1 I=1,INTOT
        JNDX(I)=0
1     CONTINUE
      DO 2 J=1,MXCASE
        ICASE(J)=0
2     CONTINUE
      DO 4 I=1,JSYLL
        JSY(I)=0
4     CONTINUE
      JND=0
      LMN=0

      DO 10 IIJ=1,ILIM
* Special cases:
      JRC(IIJ)=LMN
      IF(N.eq.0 .and. IIJ.gt.1) goto 10
      IF(N.eq.1 .and. IIJ.gt.2) goto 10
      IF(N.eq.2 .and. ISPIN.eq.1 .and. IIJ.eq.3) goto 10
      IF(N.eq.2 .and. ISPIN.eq.3 .and. IIJ.eq.4) goto 10
* ISTA is actually an offset, not a start point.
      ISTA=(IIJ-1)*MXVERT
* Element in row IJJ=IV0+1-IIJ of the tables is the top vertex, for
* IIJ=1 (Valence), 2 (Singles), 3 (Doubles T), 4 (Doubles S).
* It is found as element ISTA+IJJ of arrays L0, L1 etc.
      IJJ=IV0+1-IIJ
      KM=1
      J2(KM)=IJJ
11    KM=KM+1
      IWAY(KM)=0
12    KM1=KM-1
* At KM, trying for a way down.
      IF(L0(ISTA+J2(KM1)).EQ.0.OR.IWAY(KM).GE.1)GO TO 14
* IWAY is 0, and the next vertex L0(ISTA+J2(KM1)) is actually there:
      J2(KM)=L0(ISTA+J2(KM1))
      IWAY(KM)=1
      IOC(KM1)=0
      ISP(KM1)=0
      GO TO 20

14    IF(L1(ISTA+J2(KM1)).EQ.0.OR.IWAY(KM).GE.2)GO TO 15
* IWAY is 1, and the next vertex is actually there:
      J2(KM)=L1(ISTA+J2(KM1))
      IWAY(KM)=2
      IOC(KM1)=1
      ISP(KM1)=1
      GO TO 20

15    IF(L2(ISTA+J2(KM1)).EQ.0.OR.IWAY(KM).GE.3)GO TO 16
* IWAY is 2, and the next vertex is actually there:
      J2(KM)=L2(ISTA+J2(KM1))
      IWAY(KM)=3
      IOC(KM1)=1
      ISP(KM1)=2
      GO TO 20

16    IF(L3(ISTA+J2(KM1)).EQ.0.OR.IWAY(KM).GE.4)GO TO 17
* IWAY is 3, and the next vertex is actually there:
      J2(KM)=L3(ISTA+J2(KM1))
      IWAY(KM)=4
      IOC(KM1)=2
      ISP(KM1)=3
      GO TO 20

17    KM=KM-1
* No more ways to try from this vertex. Up to higher vertex:
      IF(KM.EQ.1)GO TO 10
* If KM was 1, then finish this IIJ value and do next one.
* Else, goto 12, to take new route from the higher vertex.
      GO TO 12

20    CONTINUE
* While trying to find a new wav through the graph, we found
* a feasible edge from a vertex J2(KM) to a vertex at level KM1.
      IF(KM1.EQ.NIORB+LV)IBS=IB(J2(KM))
      IF(KM.NE.LN+1)GO TO 11
* KM has reached top of the graph.

      JND=JND+1
* A formally legal way was found. Should it be accepted?
      NSJ=1
      INHOLE=0
      DO 110 I=1,LN
       IF(IOC(I).EQ.1)NSJ=MUL(NSJ,NSM(I))
       IF(I.LE.NIORB+LV.AND.I.GT.LV)INHOLE=INHOLE+2-IOC(I)
110   CONTINUE
C     STRIKE OUT INTERNAL CONFIGURATIONS
      IPART=0
      IF(JND.GT.IRC(1))IPART=IPART+1
      IF(JND.GT.IRC(2))IPART=IPART+1
      IF(IPART.EQ.0.AND.NSJ.NE.LSYM)GO TO 12
C     TEST IF TO HIGHLY EXCITED
C     TEST ALSO IF REFERENCE STATE
      IFEXC=0
      IFREF=0
      JJ1=0
      DO 111 IREF=1,NREF
      JHOLE=0
      JPART=IPART
      DO 112 I=1,LN
      IF(I.GT.LV)GO TO 250
      IDIF=IOC(I)
      GO TO 251
250   IF(I.GT.NIORB+LV)GO TO 252
      IDIF=IOC(I)-2
      GO TO 251
252   JJ1=JJ1+1
      IF(IOC(I).EQ.IOCR(JJ1))GO TO 112
      IDIF=IOC(I)-IOCR(JJ1)
251   IF(IDIF.GT.0)GO TO 114
      JHOLE=JHOLE-IDIF
      GO TO 112
114   JPART=JPART+IDIF
112   CONTINUE
      If (JPART.NE.JHOLE) Then
         Write (6,*) 'Config: JPART.NE.JHOLE'
         Write (6,*) 'JPART,JHOLE=',JPART,JHOLE
         Call Abend
      End If
      IF(JPART.LE.IEL)IFEXC=1
      IF(JPART.EQ.0)IFREF=1
111   CONTINUE
      IF(IFEXC.EQ.0)GO TO 12
      IF(IPART.NE.2.OR.INTNUM.EQ.0)GO TO 115
C     INTERACTING SPACE
      IF(INHOLE.EQ.2.AND.IBS.NE.0)GO TO 12
C     NO CORE-CORE CORRELATION
115   IF(IFCORE.EQ.0)GO TO 116
      NCORR=0
      DO 117 I=1,LN
      IF(ICOR(I).EQ.0)GO TO 117
      NCORR=NCORR+2-IOC(I)
117   CONTINUE
      IF(NCORR.GT.1)GO TO 12
C     SINGLY OCCUPIED ORBITALS
116   IF(NONE.EQ.0)GO TO 118
      DO 119 I=1,NONE
      IF(IOC(JONE(I)).NE.1)GO TO 12
119   CONTINUE
118   LMN=LMN+1
      L=JND
      IND=LMN
      JNDX(L)=IND
      IF(IIJ.EQ.1) THEN
C CONSTRUCT INDEX LIST FOR REFERENCE STATES
        If (LMN.GT.JRX) Then
           Write (6,*) 'Config: LMN.GT.JRX'
           Write (6,*) 'LMN,JRX=',LMN,JRX
           Write (6,*) 'This error is almost certainly that the problem'
           Write (6,*) ' at hand requires larger arrays than GUGA is'
           Write (6,*) ' compiled for. Please check your input against'
           Write (6,*) ' the manual. If you are certain that this'
           Write (6,*) ' calculation should be possible, please report'
           Write (6,*) ' this as a bug to the Molcas developers. Use'
           Write (6,*) ' the link ''http://www.teokem.lu.se/molcas'' '
           Write (6,*) ' and then ''Entry for Users'' and'
           Write (6,*) ' ''Report a new bug''.'
           Call Abend
        End If
        JREFX(LMN)=0
        IF(IFREF.EQ.0)GO TO 210
        NFREF=NFREF+1
        JREFX(LMN)=NFREF
      END IF

 210  CONTINUE
      JRC(IIJ)=LMN
      IF(LMN.GT.JSYL)GO TO 985
      JSYM(LMN)=NSJ
      IF(IIJ.LE.2)GO TO 985
      NSJ1=NSJ+1
      IF(IIJ.EQ.3)JJS(NSJ1)=JJS(NSJ1)+1
      IF(IIJ.EQ.4)JJS(NSJ1+9)=JJS(NSJ1+9)+1

 985  CONTINUE
      M=(LMN-1)*LN
      If (M+LN.GT.ISPA) Then
         Write (6,*) 'Config: M+LN.GT.ISPA'
         Write (6,*) 'M,LN,ISPA=',M,LN,ISPA
         Write (6,*) 'This error may be due to a bug.'
         Write (6,*) ' Please check your input against the manual.'
         Write (6,*) ' If there is no input errors, please report'
         Write (6,*) ' this as a bug to the Molcas developers. Use'
         Write (6,*) ' the link ''http://www.teokem.lu.se/molcas'' '
         Write (6,*) ' and then ''Entry for Users'' and'
         Write (6,*) ' ''Report a new bug''.'
         Call Abend
      End If
      DO 130 K=1,LN
      ISO(M+K)=ISP(K)
130   CONTINUE
      GO TO 12
 10   CONTINUE

      JRC(ILIM)=LMN
      WRITE(IW,214)
214   FORMAT(//,6X,'INTERNAL CONFIGURATIONS (REAL)')
      IX1=JRC(1)
      IX2=JRC(2)-JRC(1)
      IF(IFIRST.NE.0)GO TO 205

*      IF(N.EQ.2) JRC(3)=JRC(2)

      IX3=JRC(3)-JRC(2)
      IX4=JRC(4)-JRC(3)
      WRITE(IW,215)IX1,IX2,IX3,IX4
215   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7,
     */,6X,'NUMBER OF TRIPLET COUPLED DOUBLES',I7,
     */,6X,'NUMBER OF SINGLET COUPLED DOUBLES',I7)
      If (LMN.GT.JSYL) Then
         Write (6,*) 'Config: LMN.GT.JSYL'
         Write (6,*) 'LMN,JSYL=',LMN,JSYL
         Call Abend()
      End If
      If (IX1.GE.8192.OR.IX2.GE.8192.OR.
     &    IX3.GE.8192.OR.IX4.GE.8192) Then
         Write (6,*) 'Config: IX?.GE.8192'
         Write (6,*) 'IX1,IX2,IX3,IX4=',IX1,IX2,IX3,IX4
         Write (6,*) 'This error is almost certainly that the problem'
         Write (6,*) ' at hand requires larger arrays than GUGA is'
         Write (6,*) ' compiled for. Please check your input against'
         Write (6,*) ' the manual. If you are certain that this'
         Write (6,*) ' calculation should be possible, please report'
         Write (6,*) ' this as a bug to the Molcas developers. Use'
         Write (6,*) ' the link ''http://www.teokem.lu.se/molcas'' '
         Write (6,*) ' and then ''Entry for Users'' and'
         Write (6,*) ' ''Report a new bug''.'
         Call Abend()
      End If
C     SORT BY SYMMETRY
      IF(NSYM.EQ.1)GO TO 410
      ITU=2
400   IRC1=IRC(ITU)+1
      IRC2=IRC(ITU+1)
      JRC1=JRC(ITU)+1
      JRC2=JRC(ITU+1)
      JRC21=JRC2-1
      IF(JRC21.LT.JRC1)GO TO 401
      DO 402 I=JRC1,JRC21
      I1=I+1
      DO 403 J=I1,JRC2
      IF(JSYM(J).GE.JSYM(I))GO TO 403
      ITEMP=JSYM(I)
      JSYM(I)=JSYM(J)
      JSYM(J)=ITEMP
      M1=(I-1)*LN
      M2=(J-1)*LN
      DO 404 K=1,LN
      ITEMP=ISO(M1+K)
      ISO(M1+K)=ISO(M2+K)
      ISO(M2+K)=ITEMP
404   CONTINUE
      DO 420 K=IRC1,IRC2
      IF(JNDX(K).NE.I)GO TO 421
      JNDX(K)=J
      GO TO 420
421   IF(JNDX(K).NE.J)GO TO 420
      JNDX(K)=I
420   CONTINUE
403   CONTINUE
402   CONTINUE
401   IF(ITU.EQ.3)GO TO 406
      ITU=3
      GO TO 400

406   CONTINUE
* JJS(2..9): JJS(I+1)=Nr of internal triplet states per symmetry I.
* JJS(11..18): JJS(I+10)=Nr of internal singlet states per symmetry I.
      WRITE(IW,407)(JJS(I+1),I=1,NSYM)
407   FORMAT(/6X,'INTERNAL TRIPLET STATES PER SYMMETRY:',6X,8I5)
      WRITE(IW,408)(JJS(I+10),I=1,NSYM)
408   FORMAT( 6X,'INTERNAL SINGLET STATES PER SYMMETRY:',6X,8I5)
      DO 405 I=2,NSYM
       I1=I+1
       JJS(I1)=JJS(I)+JJS(I1)
       JJS(I1+9)=JJS(I+9)+JJS(I1+9)
405   CONTINUE
* Now, JJS is changed to contain instead the corresponding cumulative sum
* summed over the symmetries.
      GO TO 410

205   WRITE(IW,216)IX1,IX2
216   FORMAT(/,6X,'NUMBER OF VALENCE STATES',I16,
     */,6X,'NUMBER OF DOUBLET COUPLED SINGLES',I7)
      If (IX1.GE.8192.OR.IX2.GE.8192) Then
         Write (6,*) 'Config: IX?.GE.8192'
         Write (6,*) 'IX1,IX2=',IX1,IX2
         Call Abend()
      End If

C     PACK OCCUPATION AND SYMMETRY VECTORS FOR CI
410   LMN0=JRC(ILIM)*LN
CPAM97      M1=(LMN0+29)/30
      M1=(LMN0+14)/15
      If (M1.GT.MXCASE) Then
         Write (6,*) 'Config: M1.GT.MXCASE'
         Write (6,*) 'M1,MXCASE=',M1,MXCASE
         Call Abend()
      End If
      M=0
      DO 411 L=1,LMN
      DO 412 K=1,LN
      M=M+1
      MND=ISO(M)
C     IOCC((M+29)/30)=OR(IOCC((M+29)/30),
C    1SHIFT(MND,2*((M+29)/30*30-M)))
CPAM97      QOCC((M+29)/30)=PACK(QOCC((M+29)/30), MND, 2*M-(2*M-1)/60*60, 2)
      CALL ICPCK(ICASE,M,MND)
412   CONTINUE
CPAM97      NSJ=JSYM(L)-1
CPAM97      JSY((L+9)/10)=IOR(JSY((L+9)/10),ISHFT(NSJ,29-3*MOD(L-1,10)))
      CALL JSPCK(JSY,L,JSYM(L))
411   CONTINUE
*
      Return
      End
