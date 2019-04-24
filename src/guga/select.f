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
      SUBROUTINE CI_SELECT(INDOUT,ICAD,IBUFL,L0,L1,L2,L3,
     &                     KBUF,NTPB,NBINS)
      IMPLICIT REAL*8 (A-H,O-Z)
      External Int8, Int2, int1, int4
#include "SysDef.fh"
      DIMENSION INDOUT(*),ICAD(*),IBUFL(*),L0(*),L1(*),L2(*),L3(*)
#include "real_guga.fh"
#include "integ.fh"
#include "files_guga.fh"
      COMMON/CNSTS/D0,D1,D2
#include "addr_guga.fh"
*
      CALL QENTER('SELECT')
      CALL JTIME(IST)
      CALL CI_SELECT_INTERNAL(INDOUT,IBUFL)
      CALL QEXIT('SELECT')
      RETURN
*
*     This is to allow type punning without an explicit interface
      CONTAINS
      SUBROUTINE CI_SELECT_INTERNAL(INDOUT,IBUFL)
      USE ISO_C_BINDING
      INTEGER, TARGET :: INDOUT(*),IBUFL(*)
      REAL*8, POINTER :: dINDOUT(:),dIBUFL(:)
      CALL C_F_POINTER(C_LOC(INDOUT),dINDOUT,[1])
      CALL C_F_POINTER(C_LOC(IBUFL),dIBUFL,[1])
*
      KBUF2=(RTOI+1)*KBUF+2
      ID=0
      DO 5 IREC=1,NBINS
         IBUFL(IREC)=0
         ICAD(IREC)=ID
         INDOUT(ID+KBUF2)=-1
         ID=ID+KBUF2
5     CONTINUE
      II=0
      IID=0
      JJ=0
      JJD=0
      JTYP=0
C     DIAGONAL ELEMENTS
      IADD11=0
      IAD10(3)=IADD10
      IOUT=0
      NMAT=0
      IDIAG=1
      DO 205 M3=1,LN
      DO 200 M4=1,M3
      DO 300 M1=M3,LN
      M2MIN=1
      IF(M1.EQ.M3)M2MIN=M4
      DO 310 M2=M2MIN,M1
      IF(M1.NE.M3.OR.M2.NE.M4)GO TO 310
      IF(M1.EQ.M2)GO TO 310
      NA=ICH(M1)
      NB=ICH(M2)
      IF(NA.GE.NB)GO TO 131
      NSAVE=NA
      NA=NB
      NB=NSAVE
131   NC=ICH(M3)
      ND=ICH(M4)
      IF(NC.GE.ND)GO TO 130
      NSAVE=NC
      NC=ND
      ND=NSAVE
130   IF(NA.GT.NC)GO TO 132
      IF(NA.EQ.NC)GO TO 133
      NSAV1=NA
      NSAV2=NB
      NA=NC
      NB=ND
      NC=NSAV1
      ND=NSAV2
      GO TO 132
133   IF(NB.GT.ND)GO TO 132
      NSAVE=NB
      NB=ND
      ND=NSAVE
132   CALL INT7(ND,NB,NA,IDIAG,dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB)
310   CONTINUE
300   CONTINUE
200   CONTINUE
205   CONTINUE
      CALL AIAI(dINDOUT,INDOUT,ICAD,IBUFL,KBUF,NTPB,NBINS)
      CALL EMPTY(dIBUFL,IBUFL(RtoI*kBuf+1),ICAD,dINDOUT,KBUF,NTPB)
      ICOP1(nCOP+1)=IOUT
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+IOUT
      ICOP1(nCOP+1)=-1
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      WRITE(IW,601)NMAT
601   FORMAT(/6X,'COEFFICIENTS FOR DIAG',I9)
      CALL JTIME(IFIN)
      ITIM=IFIN-IST
      WRITE(IW,701)ITIM
701   FORMAT(6X,'TIME FOR DIAG',I17)
      IAD10(4)=IADD10
      IST=IFIN
      JTYP=1
CFUE start modification
C     IF(IFIRST.EQ.0)CALL AI(JTYP,INDOUT,L0,L1,L2,L3)
      IF(IFIRST.EQ.0) THEN
        CALL AI(JTYP,INDOUT,L0,L1,L2,L3)
      ELSE
        ICOP1(nCOP+1)=0
        CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
        CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
        ICOP1(nCOP+1)=-1
        CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
        CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      END IF
CFUE end modification
      CALL JTIME(IFIN)
      ITIM=IFIN-IST
      WRITE(IW,702)ITIM
702   FORMAT(6X,'TIME FOR ABCI',I17)
      IAD10(5)=IADD10
      IST=IFIN
      IOUT=0
      NMAT=0
      IDIAG=0
      DO 405 M3=1,LN
      NSC=NSM(ICH(M3))
      DO 400 M4=1,M3
      NSD=NSM(ICH(M4))
      NSCD=MUL(NSC,NSD)
      DO 500 M1=M3,LN
      NSA=NSM(ICH(M1))
      NSABC=MUL(NSA,NSCD)
      M2MIN=1
      IF(M1.EQ.M3)M2MIN=M4
      DO 510 M2=M2MIN,M1
      NA=ICH(M1)
      NB=ICH(M2)
      NSB=NSM(NB)
      IF(NSB.NE.NSABC)GO TO 510
      IF(NA.GE.NB)GO TO 231
      NSAVE=NA
      NA=NB
      NB=NSAVE
231   NC=ICH(M3)
      ND=ICH(M4)
      IF(NC.GE.ND)GO TO 230
      NSAVE=NC
      NC=ND
      ND=NSAVE
230   IF(NA.GT.NC)GO TO 232
      IF(NA.EQ.NC)GO TO 233
      NSAV1=NA
      NSAV2=NB
      NA=NC
      NB=ND
      NC=NSAV1
      ND=NSAV2
      GO TO 232
233   IF(NB.GT.ND)GO TO 232
      NSAVE=NB
      NB=ND
      ND=NSAVE
232   IOUT=IOUT+1
      ICOP1(IOUT)=0
      IF(IOUT.LT.NBUF)GO TO 460
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
460   IOUT=IOUT+1
*      IND1=NA+2**8*NB
*      IND2=IND1+2**16*NC
*      ICOP1(IOUT)=IND2+2**24*ND
      IND1=IOR(NA,ISHFT(NB,8))
      IND2=IOR(IND1,ISHFT(NC,16))
      ICOP1(IOUT)=IOR(IND2,ISHFT(ND,24))
      IF(IOUT.LT.NBUF)GO TO 511
      ICOP1(nCOP+1)=NBUF
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+NBUF
      IOUT=0
511   IF(NA.NE.NC)GO TO 520
      IF(NB.EQ.ND)GO TO 525
      IF(NB.EQ.NC)GO TO 524
      DO 101 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT61(ND,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
      CALL INT62(ND,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
101   CONTINUE
      GO TO 510
524   DO 108 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT9(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
108   CONTINUE
      GO TO 510
525   IF(NA.EQ.NB)GO TO 510
      CALL INT7(ND,NB,NA,IDIAG,dINDOUT,INDOUT,ICAD,IBUFL,
     *KBUF,NTPB)
      GO TO 510
520   IF(NB.NE.ND)GO TO 530
      IF(NC.EQ.ND)GO TO 529
      CALL INT5(ND,NC,NA)
      GO TO 510
529   DO 109 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT9(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
109   CONTINUE
      GO TO 510
530   IF(NB.NE.NC)GO TO 535
      DO 102 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT4(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
102   CONTINUE
      GO TO 510
535   IF(NA.EQ.NB)GO TO 546
      IF(NC.NE.ND)GO TO 547
      DO 103 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT8(NB,NA,NC,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
103   CONTINUE
      GO TO 510
546   IF(NC.EQ.ND)GO TO 510
      DO 104 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT8(ND,NC,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
104   CONTINUE
      GO TO 510
547   IF(NB.GT.ND)GO TO 540
      DO 105 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT3(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
105   CONTINUE
      GO TO 510
540   IF(NB.GT.NC)GO TO 545
      DO 106 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT2(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
106   CONTINUE
      GO TO 510
545   DO 107 ITT=1,ILIM
      IT1=(ITT-1)*MXVERT
      IT2=IT1
      CALL INT1(ND,NC,NB,NA,IT1,IT2,II,IID,JJ,JJD,JTYP,INDOUT,
     *L0,L1,L2,L3)
107   CONTINUE
510   CONTINUE
500   CONTINUE
400   CONTINUE
405   CONTINUE
      ICOP1(nCOP+1)=IOUT
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT=NMAT+IOUT
      ICOP1(nCOP+1)=-1
      CALL dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      CALL iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      WRITE(IW,600)NMAT
600   FORMAT(/6X,'COEFFICIENTS FOR IJKL',I9)
      CALL JTIME(IFIN)
      ITIM=IFIN-IST
      WRITE(IW,703)ITIM
703   FORMAT(6X,'TIME FOR IJKL',I17)
      IAD10(6)=IADD10
      IST=IFIN
      NMAT=0
      CALL AIBJ(L0,L1,L2,L3,INDOUT)
      CALL JTIME(IFIN)
      ITIM=IFIN-IST
      WRITE(IW,704)ITIM
704   FORMAT(6X,'TIME FOR AIBJ',I17)
      IAD10(7)=IADD10
      IST=IFIN
      CALL AIJK(INDOUT,L0,L1,L2,L3)
      CALL JTIME(IFIN)
      ITIM=IFIN-IST
      WRITE(IW,705)ITIM
705   FORMAT(6X,'TIME FOR AIJK',I17)
      IAD10(8)=IADD10
      IST=IFIN
      CALL ONEEL_GUGA()
      IAD10(9)=IADD10
      JTYP=0
      CALL AI(JTYP,INDOUT,L0,L1,L2,L3)
      CALL JTIME(IFIN)
      ITIM=IFIN-IST
      WRITE(IW,706)ITIM
706   FORMAT(6X,'TIME FOR ONEEL',I16)
      NULLIFY(dINDOUT,dIBUFL)
      RETURN
      END SUBROUTINE CI_SELECT_INTERNAL
*
      END
