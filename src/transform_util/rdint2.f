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
* Copyright (C) 1998, Jun-ya Hasegawa                                  *
*               2004, Giovanni Ghigo                                   *
************************************************************************
      SUBROUTINE RDINT2(IPRX,DoTCVA)
C
C     SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. TEST SECTION
C
C     THIS SUBROUTINE READS AND CHECKS THE RESULT OF THE SECOND ORDER
C     TWO-ELECTRON TRANSFORMATION PROGRAM TRA2. IT CAN BE CALLED BY
C     TR2CTL IMMEDIATELY AFTER THE CALL TO TRA2
C      IPRX=0 Do not print Coulomb nor Exchange Integrals
C      IPRX=1 Print only Coulomb Integrals
C      IPRX=2 Print all Integrals
C      IPRX=3 Print only Exchange Integrals
C
*     Modified for caspt2 by Jun-ya Hasegawa(98/08/18)
*     Modified again by G. Ghigo - September-December 2004
C
C     A,B are MO indices, counting only non-frozen and non-deleted.
C     T,U are occupied MO indices, only non-frozen and non-deleted.
C
C     <AB/TU> ARE ALWAYS GENERATED
C     EXCHANGE INTEGRALS <AT/BU> ARE GENERATED AS FOLLOWS:
C     <AT/BU> IF ISP.GE.ISR
C     <AT/UB> IF ISP.GT.ISS AND ISP.NE.ISQ
C     <TA/BU> IF ISQ.GT.ISR AND ISP.NE.ISQ
C     <TA/UB> IF ISQ.GE.ISS AND ISP.NE.ISQ
C
C     IAD2M CONTAINS START ADRESS FOR EACH TYPE OF INTEGRALS:
C     IAD2M(1,iSymIJAB)   COULOMB INTEGRALS <AB|TU>
C     IAD2M(2,iSymIJAB)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T > SYM U
C     IAD2M(3,iSymIJAB)   EXCHANGE INTEGRALS <AB|TU> FOR SYM T < SYM U
C     THE LAST ADRESS IS ZERO IF SYM T = SYM U

      Implicit Integer (i-n)
      Implicit Real*8 (A-H,O-Z)

      Parameter (MxFro = 50)
      Parameter (MxDel = 50)

#include "rasdim.fh"
#include "caspt2.fh"
#include "trafo.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"

      Dimension IAD2M(3*36*36)
      Logical DoTCVA
      Logical Found

C GG-Dec04  The following informations must be passed to the Cholesky
C transformation section through RunFile. COMMON blocks could not be
C used due to several conflicts.
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iArray('nFroPT',nFro,nSym)
      Call Get_iArray('nDelPT',nDel,nSym)
      Call Get_iArray('nIsh',nIsh,nSym)
* PAM 2007      Call Get_iArray('nAsh',nAsh,nSym)
* Replaced by:
      Do i=1,nSym
        nAsh(i)=0
      End Do
      Call qpg_iArray('nAsh',Found,nData)
      If(Found .and. nData.eq.nSym) Then
        Call Get_iArray('nAsh',nAsh,nSym)
      End If
* End of replacement
      Do i=1,nSym
        nOrb(i) = nBas(i) - nFro(i) - nDel(i)
        nIsh(i) = nIsh(i) - nFro(i)
        nOsh(i) = nIsh(i) + nAsh(i)
        nSsh(i) = nOrb(i) - nOsh(i)
      EndDo

      Write(6,*)
      Write(6,*) 'SECOND ORDER TWO-ELECTRON TRANFORMATION PROGRAM. ',
     &' TEST SECTION:'
      Write(6,*) ' A,B are MO-symmetry indices, counting only',
     &'  non-frozen and non-deleted.'
      Write(6,*) ' I,J are occupied MO-symmetry indices, only',
     &'  non-frozen and non-deleted.'
      Write(6,*) ' i,j are occupied MO indices'
c      Write(6,*) ' <AB/TU> COULOMB INTEGRALS ARE ALWAYS GENERATED.'
c      Write(6,*) ' <AB|TU> EXCHANGE INTEGRALS ARE GENERATED AS FOLLOWS:'
c      Write(6,*) '   <AT/BU> IF ISP.GE.ISR'
c      Write(6,*) '   <AT/UB> IF ISP.GT.ISS AND ISP.NE.ISQ'
c      Write(6,*) '   <TA/BU> IF ISQ.GT.ISR AND ISP.NE.ISQ'
c      Write(6,*) '   <TA/UB> IF ISQ.GE.ISS AND ISP.NE.ISQ'
c      Write(6,*) ' IAD2M(1,iSymIJAB)  COULOMB INTEGRALS <AB|TU>'
c      Write(6,*) ' IAD2M(2,iSymIJAB)  EXCHANGE INTEGRALS <AB|TU> FOR SYM
c     &  T > SYM U'
c      Write(6,*) ' IAD2M(3,iSymIJAB)  EXCHANGE INTEGRALS <AB|TU> FOR SYM
c     &  T < SYM U'
c      Write(6,*) ' THE LAST ADRESS IS ZERO IF SYM T = SYM U'
      Write(6,*)
      Write(6,'(A,8I3)')'        Symmetries :',(i,i=1,nSym)
      Write(6,*)
      Write(6,'(A,8I3)')'           Frozen  :',(nFro(i),i=1,nSym)
      Write(6,'(A,8I3)')'      Inactive (I) :',(nIsh(i),i=1,nSym)
      Write(6,'(A,8I3)')'        Active (A) :',(nAsh(i),i=1,nSym)
      Write(6,'(A,8I3)')'     Secondary (S) :',(nSsh(i),i=1,nSym)
      Write(6,'(A,8I3)')'          Deleted  :',(nDel(i),i=1,nSym)
      Write(6,*)
      Write(6,'(A,8I3)')'  Total correlated :',(nOrb(i),i=1,nSym)
      Call XFlush(6)

C
C     READ ADDRESS RECORD ON UNIT LUINTM
C
      IAD13=0
      LIADUT=3888
      CALL iDAFILE(LUINTM,2,IAD2M,LIADUT,IAD13)
C
C     LOOP OVER QUADRUPLES OF SYMMETRIES (iSymI,iSymJ,iSymA,iSymB)
C
      iSymIJAB=0
      LTotInt=0
      LTotEx1=0
      LTotEx2=0
      DO 104 iSymI=1,NSYM
       nOrbI=nOrb(iSymI)
       nOccI=nOsh(iSymI)
       DO 103 iSymJ=1,iSymI
        nOrbJ=nOrb(iSymJ)
        nOccJ=nOsh(iSymJ)
        iSymIJ=MUL(iSymI,iSymJ)
        DO 102 iSymA=1,NSYM
         nOrbA=nOrb(iSymA)
         nOccA=nOsh(iSymA)
         iSymIJA=MUL(iSymIJ,iSymA)
         ISR=iSymA
         DO 101 iSymB=1,iSymA
          nOrbB=nOrb(iSymB)
          nOccB=nOsh(iSymB)
          iSymIJAB=iSymIJAB+1

          IF(iSymIJA.NE.iSymB) GO TO 101
          IF((nOccI*nOccJ).EQ.0) GO TO 101
          IF(nOrbI*nOrbJ*nOrbA*nOrbB.EQ.0) GO TO 101
C
C         FIND ADDRESSES FOR THIS SYMMETRY BLOCK
C
          IADC =IAD2M(3*iSymIJAB-2)
          IADX1=IAD2M(3*iSymIJAB-1)
          IADX2=IAD2M(3*iSymIJAB)
          WRITE(6,1000) iSymA,iSymB,iSymI,iSymJ
1000      FORMAT(/1X,'SYMMETRY BLOCK < A B | I J >',4I4)

          IF(IADC.EQ.0) THEN
C           WRITE(6,*)
           WRITE(6,1100)
1100       FORMAT(1X,'NO COULOMB INTEGRALS FOR THIS SYMMETRY BLOCK')
          ELSE
           WRITE(6,1200) IADC
1200       FORMAT(1X,'ADDRESS FOR COULOMB INTEGRALS',I8)
           IAD13C=IADC
          ENDIF
          IF(IADX1.EQ.0) THEN
C           WRITE(6,*)
           WRITE(6,1110)
1110       FORMAT(1X,'NO EXCHAN1 INTEGRALS FOR THIS SYMMETRY BLOCK')
          ELSE
           WRITE(6,1210) IADX1
1210       FORMAT(1X,'ADDRESS FOR EXCHAN1 INTEGRALS',I8)
           IAD131=IADX1
          ENDIF
          IF(IADX2.EQ.0) THEN
C           WRITE(6,*)
           WRITE(6,1120)
1120       FORMAT(1X,'NO EXCHAN2 INTEGRALS FOR THIS SYMMETRY BLOCK')
          ELSE
           WRITE(6,1220) IADX2
1220       FORMAT(1X,'ADDRESS FOR EXCHAN2 INTEGRALS',I8)
           IAD132=IADX2
          ENDIF
          LREC=nOrbA*nOrbB
          IF(iSymA.EQ.iSymB) LREC=(nOrbA**2+nOrbA)/2
          LRECX=nOrbA*nOrbB
          IF(.NOT.DoTCVA) LRECX=(nOrbA-nOccA)*(nOrbB-nOccB)
          LInt=0
          LEx1=0
          LEx2=0
          DO 10 NI=1,nOccI
           NUM=nOccJ
           IF(iSymI.EQ.iSymJ) NUM=NI
            DO 11 NJ=1,NUM

C     THE LOOP ABOVE OVER T AND U RECOVERS ONE BLOCK OF INTEGRALS <AB|TU>
C     FOR EACH PAIR T,U. TO PROCESS ONE BLOCK FOR ALL A AND B THE
C     FOLLOWING LOOP STRUCTURE IS USED:

            IF (IADC.NE.0) THEN
             LInt=LInt+LREC
             If (IPRX.GT.0. and. IPRX.LT.3) then
              Call GetMem('Tmp','ALLO','REAL',iTmp,LREC)
              CALL dDAFILE(LUINTM,2,WORK(iTmp),LREC,IAD13C)
          WRITE(6,1300) NI,NJ,IAD13C-LREC,(WORK(I),I=iTmp,iTmp+LREC-1)
              Call GetMem('Tmp','FREE','REAL',iTmp,LREC)
1300          FORMAT(/1X,'<AB|IJ> COULOMB INTEGR.S FOR |ij> PAIR',2I3,
     &      '  DiskAdd=',I8  /(8F10.6))
             EndIf
            ENDIF

C      THE EXCHANGE INTEGRALS OF TYPE 1 ,(AT|BU) ARE PROCESSED AS THE
C      COULOMB INTEGRALS. IF NST.NE.NSU THERE ARE ALSO EXCHANGE
C      INTEGRALS OF TYPE 2, (AU|BT). THE ORDERING IS STILL T,U AND A,B
C      BUT T IS NOW THE FOURTH INDEX AND U THE SECOND
C      EXCHANGE INTEGRALS ARE ALWAYS QUADRATIC IN A,B

            IF (IADX1.NE.0) THEN
             LEx1=LEx1+LRECX
             If (IPRX.GT.1) then
              Call GetMem('Tmp','ALLO','REAL',iTmp,LRECX)
              CALL dDAFILE(LUINTM,2,WORK(iTmp),LRECX,IAD131)
          WRITE(6,1310) NI,NJ,IAD131-LRECX,(WORK(I),I=iTmp,iTmp+LRECX-1)
              Call GetMem('Tmp','FREE','REAL',iTmp,LRECX)
1310          FORMAT(/1X,'EXCHAN1 INTEGRALS FOR |ij> PAIR',2I3,
     &      '  DiskAdd=',I8  /(8F10.6))
             EndIf
            ENDIF

            IF (IADX2.NE.0) THEN
             LEx2=LEx2+LRECX
             If (IPRX.GT.1) then
              Call GetMem('Tmp','ALLO','REAL',iTmp,LRECX)
              CALL dDAFILE(LUINTM,2,WORK(iTmp),LRECX,IAD132)
          WRITE(6,1320) NI,NJ,IAD132-LRECX,(WORK(I),I=iTmp,iTmp+LRECX-1)
              Call GetMem('Tmp','FREE','REAL',iTmp,LRECX)
1320          FORMAT(/1X,'EXCHAN2 INTEGRALS FOR |ij> PAIR',2I3,
     &      '  DiskAdd=',I8  /(8F10.6))
             EndIf
            ENDIF
11          CONTINUE
10        CONTINUE

          Write(6,*)
          Write(6,1010) LInt,LEx1,LEx2
1010      Format(3X,'LCou=',I8,' , LEx1=',I8,' , LEx2=',I8)
          LTotInt=LTotInt+LInt
          LTotEx1=LTotEx1+LEx1
          LTotEx2=LTotEx2+LEx2
C
C     ALL INTEGRALS FOR SYMMETRY BLOCK iSymI,iSymJ,iSymA,iSymB ARE READ
C
101      CONTINUE
102     CONTINUE
103    CONTINUE
104   CONTINUE
      Write(6,*)
      Write(6,2000) LTotInt,LTotEx1,LTotEx2
2000  Format(3X,'LTotCou=',I8,' , LTotEx1=',I8,' , LTotEx2=',I8)
      Write(6,*) '   LTotTot=',LTotInt+LTotEx1+LTotEx2
      Write(6,*)
C
      RETURN
      END
