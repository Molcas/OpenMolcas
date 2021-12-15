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
      SUBROUTINE BJAI(IAD,EPSI,EPSE,E2BJAI,VECL2)

      IMPLICIT REAL*8 (A-H,O-Z)

#include "mxdim.fh"
#include "corbinf.fh"
#include "files_mbpt2.fh"
#include "WrkSpc.fh"
      DIMENSION EPSI(*),EPSE(*),IAD(3888)
      Logical DoCholesky, Debug
      Data Debug/.False./
c      Data Debug/.True./   ! CGG
      Dimension dVectorA(255) ! CGG
      Dimension dVectorB(255) ! CGG

#include "SysDef.fh"

      MUL(I,J)=1+IEOR(I-1,J-1)


      SKAL2=-9999999.9
      IAD13=0
      LIADUT=3888

      Call DecideOnCholesky(DoCholesky)
      If (Debug) then
      Write(6,*)
      Write(6,'(A,8I3)')'      nOcc:',(nOcc(i),i=1,nSym)
      Write(6,'(A,8I3)')'      nExt:',(nExt(i),i=1,nSym)
      Write(6,'(A,8I3)')'      nOrb:',(nOrb(i),i=1,nSym)
      EndIf

      CALL iDAFILE(LUINTM,2,IAD,LIADUT,IAD13)

      VECL2=1.0D+00
      E2BJAI=0.0D+00
      ISPQRS=0
      NRI=0
      DO 10 iSymI=1,nSym
       nOccI=nOcc(iSymI)
       nOrbI=nOrb(iSymI)
       NRJ=0
       DO 20 iSymJ=1,iSymI
        nOccJ=nOcc(iSymJ)
        nOrbJ=nOrb(iSymJ)
        NRA=0
        DO 30 iSymA=1,nSym
         nExtA=nExt(iSymA)
         nOrbA=nOrb(iSymA)
         NRB=0
         DO 40 iSymB=1,iSymA
          nExtB=nExt(iSymB)
          nOrbB=nOrb(iSymB)
          ISPQRS=ISPQRS+1
          IF (nOccI*nOccJ*nExtA*nExtB.EQ.0) GOTO 41
          IF (MUL(iSymI,iSymJ).EQ.MUL(iSymA,iSymB)) THEN
           IAD1=IAD(3*ISPQRS-1)
           IAD2=IAD(3*ISPQRS)

      If(Debug) then
      Write(6,*)
      Write(6,*)' [BJAI] Integrals <A,B|I,J> : ',iSymA,iSymB,iSymI,iSymJ
      EndIf
*    ---   iSymI.NE.iSymJ   ---
           IF (iSymI.NE.iSymJ) THEN
            LAB = nOrbA * nOrbB
            If (DoCholesky) LAB=nExtA * nExtB
            LAB1= nExtA * nExtB
            LB  = nExtB
            Call GetMem('INT1','ALLO','REAL',LINT1,LAB)
            Call GetMem('INT2','ALLO','REAL',LINT2,LAB)
            Call GetMem('AIBJ','ALLO','REAL',LAIBJ,LAB1) ! AB|IJ
            Call GetMem('AJBI','ALLO','REAL',LAJBI,LAB1) ! AB|JI
            DO 50 iI=1,nOccI
             DO 60 iJ=1,nOccJ
      If(Debug) then
      Write(6,*)
      Write(6,*) ' *  i,j = ',iI,iJ
      EndIf
              Call dDAFile(LUINTM,2,Work(LINT1),LAB,IAD1)
              Call dDAFile(LUINTM,2,Work(LINT2),LAB,IAD2)
      If (Debug .and. DoCholesky ) then
      Call RecPrt('Int1:','(8F10.6)',Work(lInt1),nExtB,nExtA)
      Call RecPrt('Int2:','(8F10.6)',Work(lInt2),nExtB,nExtA)
      EndIf
      If (Debug .and. .NOT.DoCholesky ) then
      Call RecPrt('Int1:','(8F10.6)',Work(lInt1),nOrbB,nOrbA)
      Call RecPrt('Int2:','(8F10.6)',Work(lInt2),nOrbB,nOrbA)
      EndIf
              EDIFIJ=EPSI(NRI+iI)+EPSI(NRJ+iJ)
              IADA=nOcc(iSymA) * nOrbB + nOcc(iSymB)
              If (DoCholesky) IADA=0  ! CGG
              IADB=LINT2+IADA
              IADA=LINT1+IADA
              I=0
       If(Debug) then
       Write(6,*)
       EndIf
              DO 70 iA=0,nExtA-1
               DO 80 iB=0,nExtB-1
                A=Work(IADA+iB)
                B=Work(IADB+iB)
       If(Debug) dVectorA(iB+1)=A
       If(Debug) dVectorB(iB+1)=B
                Work(LAIBJ+I)=A+B
                Work(LAJBI+I)=A-B
                I=I+1
80             CONTINUE
       If(Debug) then
       Write(6,'(A,I3,6F10.6)') 'A:',iA+1,(dVectorA(j),j=1,nExtB)
       Write(6,'(A,I3,6F10.6)') 'B:',iA+1,(dVectorB(j),j=1,nExtB)
       Write(6,*) ' -------'
       EndIf
               If (DoCholesky) then
                IADA=IADA + nExtB
                IADB=IADB + nExtB
               else
                IADA=IADA + nOrbB
                IADB=IADB + nOrbB
               EndIf
70            CONTINUE
              I=0
              DO 71 iA=1,nExtA
               DO 72 iB=1,nExtB
                EDIF=EPSE(NRA+iA)+EPSE(NRB+iB)-EDIFIJ
                Work(LINT1+I)=Work(LAIBJ+I)/EDIF
                Work(LINT2+I)=Work(LAJBI+I)/EDIF
                I=I+1
72             CONTINUE
71            CONTINUE
              SKAL1=DDOT_(LAB1,Work(LINT1),1,Work(LAIBJ),1)
              SKAL2=DDOT_(LAB1,Work(LINT2),1,Work(LAJBI),1)
              E2BJAI=E2BJAI-SKAL1-3.0D+00*SKAL2
              CKAL1=DDOT_(LAB1,Work(LINT1),1,Work(LINT1),1)
              CKAL2=DDOT_(LAB1,Work(LINT2),1,Work(LINT2),1)
              VECL2=VECL2+CKAL1+3.0D+00*CKAL2
60           CONTINUE
50          CONTINUE
            Call GetMem('INT1','FREE','REAL',LINT1,LAB)
            Call GetMem('INT2','FREE','REAL',LINT2,LAB)
            Call GetMem('AIBJ','FREE','REAL',LAIBJ,LAB1)
            Call GetMem('AJBI','FREE','REAL',LAJBI,LAB1)
           ENDIF

*    ---   iSymI.EQ.iSymJ   ---
           IF (iSymI.EQ.iSymJ) THEN
            LA  = nExtA
            LAB = nOrbA**2
            If (DoCholesky) LAB=nExtA**2
            LAA=(LA*(LA+1))/2
            Call GetMem('INT','ALLO','REAL',LINT,2*LAB)
            Call GetMem('AIBJ','ALLO','REAL',LAIBJ,2*LAA)
            Call GetMem('AJBI','ALLO','REAL',LAJBI,2*LAA)
            IADAB=LINT + nOcc(iSymA) * (nOrbA+1)
            If (DoCholesky) IADAB=LINT
            DO 90 iI=1,nOccI
             DO 100 iJ=1,iI
      If(Debug) then
      Write(6,*)
      Write(6,*) ' *  i,j = ',iI,iJ
      EndIf
              Call dDAFile(LUINTM,2,Work(LINT),LAB,IAD1)
      If (Debug .and. DoCholesky )
     &Call RecPrt('Int:','(8F10.6)',Work(lInt),nExtA,nExtA)
      If (Debug .and. .NOT.DoCholesky )
     &Call RecPrt('Int:','(8F10.6)',Work(lInt),nOrbA,nOrbA)
              EDIFIJ=EPSI(NRI+iI)+EPSI(NRJ+iJ)
              I=0
              IADA=IADAB
              DO 110 iA=0,nExtA-1
               IADB=IADAB
       If(Debug) then
       Write(6,*) ' iA=',iA+1,'  IADA:'
       Write(6,'(8F10.6)') (WORK(IADA+j),j=0,iA)
       EndIf
               DO 120 iB=0,iA
                A=Work(IADA+iB)
                B=Work(IADB+iA)
                Work(LAIBJ+I)=A+B
                Work(LAJBI+I)=A-B
                I=I+1
                If (DoCholesky) then
                  IADB=IADB+nExtA
                else
                  IADB=IADB+nOrbA
                EndIf
120            CONTINUE
               If (DoCholesky) then
                 IADA=IADA+nExtA
               else
                 IADA=IADA+nOrbA
               EndIf
110           CONTINUE
              I=0
              DO 130 iA=0,LA-1
               DO 140 iB=0,iA
                EDIF=EPSE(NRA+iA+1)+EPSE(NRA+iB+1)-EDIFIJ
                IF (iA.EQ.iB) EDIF=2.0D+00*EDIF
                Work(LINT+I)=Work(LAIBJ+I)/EDIF
                Work(LINT+LAA+I)=Work(LAJBI+I)/EDIF
                IF ((iA.EQ.iB).AND.(iI.NE.iJ)) THEN
                 VECL2=VECL2+2.0D+00*Work(LINT+I)**2
                ELSEIF ((iA.NE.iB).AND.(iI.EQ.iJ)) THEN
                 VECL2=VECL2+0.5D+00*Work(LINT+I)**2
                ELSE
                 VECL2=VECL2+Work(LINT+I)**2
                ENDIF
                I=I+1
140            CONTINUE
130           CONTINUE
              SKAL1=DDOT_(LAA,Work(LINT),1,Work(LAIBJ),1)
              IF (iI.NE.iJ) THEN
              SKAL2=3.0D+00*DDOT_(LAA,Work(LINT+LAA),1,Work(LAJBI),1)
              ENDIF
              IF (iI.EQ.iJ) E2BJAI=E2BJAI-0.5*SKAL1
              IF (iI.NE.iJ) E2BJAI=E2BJAI-SKAL1-SKAL2
              IF (iI.NE.iJ) THEN
               CKAL2=DDOT_(LAA,Work(LINT+LAA),1,Work(LINT+LAA),1)
               VECL2=VECL2+3.0D+00*CKAL2
              ENDIF
100          CONTINUE
90          CONTINUE
            Call GetMem('INT','FREE','REAL',LINT,LAB+LA)
            Call GetMem('AIBJ','FREE','REAL',LAIBJ,LAA)
            Call GetMem('AJBI','FREE','REAL',LAJBI,LAA)
           ENDIF
          ENDIF
41        NRB=NRB+nExtB
40       CONTINUE
         NRA=NRA+nExtA
30      CONTINUE
        NRJ=NRJ+nOccJ
20     CONTINUE
       NRI=NRI+nOccI
10    CONTINUE
      VECL2=sqrt(1.0D+00/VECL2)
      RETURN
      END
