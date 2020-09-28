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
      SUBROUTINE SXHAM(D,P,PA,FP,SXN,F1,F2,DIA,G,H,HDIAG,DF,DDIAG)
C
C RASSCF VERSION IBM-3090: SX-SECTION
C
C Objective: Compute the auxiliary matrices G,H,F1,F2,and SXN
C defined as:
C
C G(pq)=sum(r,s)(2*P(pqrs)-D(pq)*D(rs))*FP(rs)
C H(pq)=G(pq)+DF(pq)+DF(qp)  where
C DF(pq)=sum(r)D(pr)*FP(rq)
C F1 is the occupied part of the SX Fock matrix FP
C F2 is the active-external part of the same matrix
C DIA is the occupied part of the density matrix(squared)
C SXN contains normalization factors for the SX states
C HDIAG is the diagonal of the SX Hamiltonian
C All matrices are blocked according to symmetry
C Each block having the dimension
C G, DIA, AND F1:     (NIO+NAO)**2
C F2:                 (NAO+NEO)**2
C H:                  NAO*NAE
C SXN:                (NIO+NAO)*(NAO+NEO)
C HDIAG:              NROOT+(NIO+NAO)*(NAO+NEO)
C
C These matrices are used to simplify the calculation
C of the sigma vector (see SIGVEC for further details).
C
C Called from SXCTL
C
C Subroutine calls: none
C
C ********** IBM-3090 RELEASE 89 01 23 **********
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "warnings.fh"
#include "rasscf.fh"
#include "general.fh"
#include "gas.fh"
#include "output_ras.fh"
      Parameter (ROUTINE='SXHAM   ')
      DIMENSION D(*),P(*),PA(*),FP(*),SXN(*),F1(*),F2(*),DIA(*),
     *          G(*),H(*),HDIAG(*),DF(*),DDIAG(*)
C -- THRA: THRESHOLD FOR WARNING, ACTIVE OCC NO CLOSE TO 0 OR 2.
      DATA THRA/1.D-06/
      DIMENSION P2Act(1)
C Local print level (if any)
      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
C
C Loop over all symmetry blocks
C
      ISTBM=0
      ISTD=0
      ISTFP=0
      ISTIA=0
      ISTAE=0
      ISTH=0
      IASHI=0
      ISTZ=0
      IX1=0
      DO 100 ISYM=1,NSYM
       IX=IX1+NFRO(ISYM)
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
c      NR1O=NRS1(ISYM)
c      NR2O=NRS2(ISYM)
c      NR3O=NRS3(ISYM)
c      NAO=NR1O+NR2O+NR3O
       NIA=NIO+NAO
       NEO=NSSH(ISYM)
       NAE=NAO+NEO
       NO=NORB(ISYM)
       IF(NIA.EQ.0.OR.NAE.EQ.0) GO TO 98
C
C Form the diagonal of the density matrix for all orbitals
C
       justone=0
       DO 10 NP=1,NO
        IF(NP.LE.NIO) DDIAG(NP)=2.0D0
        IF(NP.GT.NIA) DDIAG(NP)=0.0D0
        IF(NP.GT.NIO.AND.NP.LE.NIA) THEN
         DDIAG(NP)=D(ISTD+ITRI(NP-NIO+1))
         If ( DDIAG(NP).LT.THRA .AND. ISCF.EQ.0 ) Then
          If (JustOne.eq.0) Then
           Call WarningMessage(1,'Problems in orbital optimization.')
           JustOne=JustOne+1
          End If
          Write(LF,'(1x,a,i2,a,i4,a,d14.6)')' Warning: In symmetry ',
     &     ISYM,', orbital p=',NP,' has diagonal density matrix ' //
     &     'element D(p,p) close to zero. D(p,p)=',DDIAG(NP)
          End If
         If ( 2.0D0-DDIAG(NP).LT.THRA .AND. ISCF.EQ.0 ) Then
          If (JustOne.eq.0) Then
           Call WarningMessage(1,'Problems in orbital optimization.')
           JustOne=JustOne+1
          End If
          Write(LF,'(1x,a,i2,a,i4,a,d14.6)')' Warning: In symmetry ',
     &     ISYM,', orbital p=',NP,' has diagonal density matrix ' //
     &     'element D(p,p) close to two. (2 - D(p,p))=',2.0D0-DDIAG(NP)
         End If
        ENDIF
10     CONTINUE
C
C Compute the normalization constants:
C SXN(pq)=4*(P(prrp)-P(prpr))+D(rr)+D(pp)
C
       NPR=ISTBM
C NP loops over active and secondary indices.
       DO 15 NP=NIO+1,NO
C NR loops over inactive and active indices.
       DO 16 NR=1,NIA
        NPR=NPR+1
        SXN(NPR)=1.0D0
        IF(NR.GE.NP) GO TO 95
        DRR=DDIAG(NR)
        DPP=DDIAG(NP)
        PRPR=-2.0D0*DPP
        IF(NP.GT.NIA) PRPR=0.0D0
        IF(NR.GT.NIO.AND.NP.LE.NIA) THEN
         NT=NP-NIO+IASHI
         NU=NR-NIO+IASHI
         NTU=ITRI(NT)+NU
         NTUTU=ITRI(NTU+1)
         PRPR=-4.0D0*PA(NTUTU)
        ENDIF
        SXNRM2=PRPR+DRR+DPP
        IF(SXNRM2.LT.-1.0D-12) THEN
          Call WarningMessage(1,'Negative norm occured in SXHAM.')
          Write(LF,*)'SXHAM Error: Negative SXNRM2.'
          Write(LF,*)' Symmetry block ISYM:',ISYM
          Write(LF,*)'      Orbitals NP,NR:',NP,NR
          Write(LF,*)'                PRPR:',PRPR
          Write(LF,*)'                 DRR:',DRR
          Write(LF,*)'                 DPP:',DPP
          Write(LF,*)'      Norm**2 SXNRM2:',SXNRM2
         Write(LF,*)
         Write(LF,*)' The squared norm of a trial vector has been'
         Write(LF,*)' computed to be negative.'
         Write(LF,*)' This is possible only for some severe malfunction'
         Write(LF,*)' of the rasscf program. Please issue a bug report.'
         Write(LF,*)
          Call Quit(_RC_GENERAL_ERROR_)
        ENDIF
C
        IF(SXNRM2.LT.1.0D-12) THEN
          SXN(NPR)=0.0D0
        ELSE
          SXN(NPR)=1.0D0/sqrt(SXNRM2)
        END IF
95      CONTINUE
16     CONTINUE
15     CONTINUE
C
C Form the occupied (F1) and active-external (F2) part of FP
C and the occupied part DIA of the density matrix D
C
       IPQ=ISTFP
       CALL FZERO(DIA(ISTIA+1),NIA**2)
       DO 20 NP=1,NO
       DO 21 NQ=1,NP
        IPQ=IPQ+1
        IF(NP.LE.NIA.AND.NQ.LE.NIA) THEN
         F1(ISTIA+NIA*(NP-1)+NQ)=FP(IPQ)
         F1(ISTIA+NIA*(NQ-1)+NP)=FP(IPQ)
        ENDIF
        IF(NP.GT.NIO.AND.NQ.GT.NIO)THEN
         F2(NAE*(NP-NIO-1)+NQ-NIO+ISTAE)=FP(IPQ)
         F2(NAE*(NQ-NIO-1)+NP-NIO+ISTAE)=FP(IPQ)
        ENDIF
        IF(NP.LE.NIO.AND.NP.EQ.NQ) DIA(ISTIA+NIA*(NP-1)+NP)=2.0D0
        IF((NP.GT.NIO.AND.NP.LE.NIA).AND.NQ.GT.NIO)THEN
         NT=NP-NIO
         NU=NQ-NIO
         NTU=ITRI(NT)+NU+ISTD
         DIA(ISTIA+NIA*(NP-1)+NQ)=D(NTU)
         DIA(ISTIA+NIA*(NQ-1)+NP)=D(NTU)
        ENDIF
21     CONTINUE
20     CONTINUE
       IF(IPRLEV.ge.DEBUG) THEN
         write(6,*) 'DIA in SXHAM:'
         write(6,*) (DIA(ISTIA+i),i=1,NIA**2)
       END IF
C
C Form the matrix DF (row index active, column index all orbitals)
C
       IF(NAO.NE.0) THEN
        CALL DGEMM_('N','N',
     &              NAO,NAE,NAO,
     &              1.0d0,DIA(ISTIA+NIO*NIA+NIO+1),NIA,
     &              F2(ISTAE+1),NAE,
     &              0.0d0,DF(NAO*NIO+1),NAO)
        IF(NIO.NE.0) THEN
         CALL DGEMM_('N','N',
     &               NAO,NIO,NAO,
     &               1.0d0,DIA(ISTIA+NIO*NIA+NIO+1),NIA,
     &               F1(ISTIA+NIO+1),NIA,
     &               0.0d0,DF,NAO)
        ENDIF
       ENDIF
C
C compute the G matrix:
C G(ij)=-2*FP(ij)
C G(ti)=G(it)=-DF(ti)
C G(tu)=sum(vx)(2P(tuvx)-D(tu)D(vx))F(vx)
C
       If (ExFac.ne.1.0D0) Then
          Call Get_Temp('nP2Act  ',P2Act,1)
          nP2Act=Int(P2Act(1))
          Call Get_Temp('P2_RAW  ',P,nP2Act)
       End If
       DO 30 NP=1,NIA
       DO 31 NQ=1,NP
        IPQ=ISTIA+NIA*(NP-1)+NQ
        IQP=ISTIA+NIA*(NQ-1)+NP
        IF(NP.LE.NIO) THEN
         G(IPQ)=-2.0D0*F1(IPQ)
         G(IQP)=G(IPQ)
        ELSE IF(NP.GT.NIO.AND.NQ.LE.NIO) THEN
         G(IPQ)=-DF(NAO*(NQ-1)+NP-NIO)
         G(IQP)=G(IPQ)
        ELSE IF(NQ.GT.NIO) THEN
         DTU=DIA(ISTIA+NIA*(NP-1)+NQ)
         GTU=0.0D0
         NTT=NP-NIO+IASHI
         NUT=NQ-NIO+IASHI
         NTUT=ITRI(NTT)+NUT
C
         IASHJ=0
         ISTFPJ=0
         NVX=0
         DO 24 JSYM=1,NSYM
          NAOJ=NASH(JSYM)
          NIOJ=NISH(JSYM)
          IF(NAOJ.EQ.0) GO TO 23
          DO 22 NV=1,NAOJ
           NVT=NV+IASHJ
          DO 25 NX=1,NV
           NXT=NX+IASHJ
           NVX=NVX+1
           NVXT=ITRI(NVT)+NXT
           NTUVX=ITRI(MAX(NTUT,NVXT))+MIN(NTUT,NVXT)
           FAC=2.0D0
           FACD=2.0D0
           IF(NV.EQ.NX) FACD=1.0D0
           IF(NTUT.LT.NVXT) THEN
            IF(NVT.NE.NXT.AND.NTT.EQ.NUT) FAC=4.0D0
            IF(NVT.EQ.NXT.AND.NTT.NE.NUT) FAC=1.0D0
           ENDIF
           NVXF=ISTFPJ+ITRI(NV+NIOJ)+NX+NIOJ
           GTU=GTU+FP(NVXF)*(FAC*P(NTUVX)-FACD*DTU*D(NVX))
25        CONTINUE
22        CONTINUE
          IASHJ=IASHJ+NAOJ
23        ISTFPJ=ISTFPJ+ITRI(NORB(JSYM)+1)
24       CONTINUE
         G(IPQ)=GTU
         G(IQP)=GTU
        ENDIF
31     CONTINUE
30     CONTINUE
       If (ExFac.ne.1.0D0) Call Get_Temp('P2_KS   ',P,nP2Act)
C
C FORM THE H MATRIX (only the nae*nao block is needed)
C H(tu)=G(tu)+DF(tu)+DF(ut)
C H(at)=DF(ta)
C H(ab) and H(ta) are not included.
C
       IF(NAO.NE.0) THEN
        IPQ=ISTH
        DO 40 NP=1,NAO
        DO 41 NQ=1,NAE
         IPQ=IPQ+1
         IF(NQ.LE.NAO) THEN
          H(IPQ)=G(ISTIA+NIA*(NP+NIO-1)+NQ+NIO)
     *          +DF(NAO*(NP+NIO-1)+NQ)+DF(NAO*(NQ+NIO-1)+NP)
         ELSE
          H(IPQ)=DF(NAO*(NQ+NIO-1)+NP)
         ENDIF
41      CONTINUE
40      CONTINUE
       ENDIF
C
C Diagonal elements of the SX Hamiltonian
C
       NPR=ISTBM
       DO 50 NP=1,NAE
        FPP=F2(ISTAE+NAE*(NP-1)+NP)
        HPP=0.0D0
        DPP=0.0D0
        IF(NP.LE.NAO) THEN
         HPP=H(ISTH+NAE*(NP-1)+NP)
         DPP=DIA(ISTIA+NIA*(NP+NIO-1)+NP+NIO)
        ENDIF
        NRR=ISTIA+1
        DO 48 NR=1,NIA
         NPR=NPR+1
         HDIAG(NPR+NROOT)=
     *   (G(NRR)-HPP+DIA(NRR)*FPP+DPP*F1(NRR))*SXN(NPR)**2
         NRR=NRR+NIA+1
C
C Make this matrix element large for forbidden rotations
C
         IF(NP.LE.NAO.AND.NR.GT.NIO) THEN
          NT=NP
          NU=NR-NIO
          IF(NT.LE.NU) THEN
           HDIAG(NPR+NROOT)=1.D32
           SXN(NPR)=0.0D0
          ELSE
           NTU=ISTZ+ITRI(NT-1)+NU
           IF(IZROT(NTU).NE.0) THEN
             HDIAG(NPR+NROOT)=1.D32
             SXN(NPR)=0.0D0
C             Write(LF,*)'SXHAM. IZROT.eq.1 for NP,NR=',NP,NR
           END IF
          ENDIF
         ENDIF

C
C Make HDIAG large for rotations forbidden by IXSYM input
C
         IF(IXSYM(NR+IX).NE.IXSYM(NP+NIO+IX)) THEN
C           Write(LF,*)'SXHAM. IXSYM forbids NP,NR=',NP,NR
           HDIAG(NPR+NROOT)=1.D32
           SXN(NPR)=0.0D0
         END IF
C
48      CONTINUE
50     CONTINUE
C       Write(LF,*)'SXHAM: The IXSYM array='
C       Write(LF,'(1x,40i2)')(ixsym(i),i=1,ntot)
C
C Test print of all matrices
C
       IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,1000) ISYM
1000    FORMAT(/1X,'Matrices in SXHAM for symmetry block',I2)
        Write(LF,1100) (DDIAG(I),I=1,NO)
1100    FORMAT(/1X,'Diagonal of the full density matrix'/(1X,5G18.6))
        Write(LF,1200) (SXN(I),I=ISTBM+1,ISTBM+NIA*NAE)
1200    FORMAT(/1X,'Normalization constants for SX-states'/(1X,5G18.6))
        Write(LF,1300) (DF(I),I=1,NAO*NO)
1300    FORMAT(/1X,'The matrix sum(D(pr)F(rq)) -DF'/(1X,5G18.6))
        Write(LF,1400) (F1(I),I=ISTIA+1,ISTIA+NIA**2)
1400    FORMAT(/1X,'The occupied part of the Fock matrix'/(1X,5G18.6))
        Write(LF,1500) (F2(I),I=ISTAE+1,ISTAE+NAE**2)
1500    FORMAT(/1X,'The act-ext part of the Fock matrix'/(1X,5G18.6))
        Write(LF,1600) (G(I),I=ISTIA+1,ISTIA+NIA**2)
1600    FORMAT(/1X,'The G-matrix'/(1X,5G18.6))
        Write(LF,1700) (H(I),I=ISTH+1,ISTH+NAO*NAE)
1700    FORMAT(/1X,'The H-matrix'/(1X,5G18.6))
        Write(LF,1800) (HDIAG(I),I=ISTBM+1+NROOT,ISTBM+NAE*NIA+NROOT)
1800    FORMAT(/1X,'Diagonal of the SX Hamiltonian'/(1X,5G18.6))
       ENDIF
C
       IASHI=IASHI+NAO
       ISTD=ISTD+ITRI(NAO+1)
98     ISTIA=ISTIA+NIA**2
       ISTAE=ISTAE+NAE**2
       ISTFP=ISTFP+ITRI(NO+1)
       ISTBM=ISTBM+NIA*NAE
       ISTH=ISTH+NAO*NAE
       ISTZ=ISTZ+(NAO**2-NAO)/2
       IX1=IX1+NBAS(ISYM)
C
100   CONTINUE
C
C Level shift section
C
       XLEV=LVSHFT
       SXSHFT=0.0D0
       HDMIN=1.D32
       DO 152 I=NROOT+1,NROOT+NSXS
        IF(HDIAG(I).LT.HDMIN)  HDMIN=HDIAG(I)
152    CONTINUE
       SXSHFT=MAX(XLEV-HDMIN,0.0D0)
C
C no level shift if input (alpha) is zero
C
       IF(LVSHFT.EQ.0.0D0) SXSHFT=0.0D0
       DO 154 I=NROOT+1,NROOT+NSXS
        HDIAG(I)=HDIAG(I)+SXSHFT
154    CONTINUE
       IF(IPRLEV.GE.DEBUG) THEN
         Write(LF,1101) HDMIN,SXSHFT
1101   FORMAT(1X,'Lowest diagonal element has the value',F12.6,
     *        ' A level shift of',F12.6,' has been used')
       END IF
C
C Add diagonal elements for the reference space (CI-states)
C
      HDIAG(1)=0.0D0
      IF(ICICP.NE.0) THEN
       IROOT1=IROOT(1)
       DO 56 I=1,NROOT
        HDIAG(I)=ENER(I,ITER)-ENER(IROOT1,ITER)
56     CONTINUE
      ENDIF
      RETURN
      END
