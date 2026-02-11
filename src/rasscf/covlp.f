!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SUBROUTINE COVLP(C1IN,C2IN,DIA,PA,SXN,C1,C2,X,OVL)
      use rasscf_global, only: NROOT, NSXS, ITRI
      use PrintLevel, only: DEBUG
      use output_ras, only: LF,IPRLOC
      use general_data, only: NSYM,NASH,NISH,NSSH

      IMPLICIT None
!
! RASSCF program version IBM-3090: SX section
!
! Purpose:Calculation of the overlap between two super-CI
! vectors C1IN and C2IN. The result is given in OVL.
! C1,C2, and X are scratch areas.
!
! ********** IBM-3090 Release 89 01 25 **********
!PAM01 Added: replace correct overlap by adding a diagonal
!PAM01 quantity to the overlap of brillouin states.
!
      REAL*8 C1IN(*),C2IN(*),DIA(*),SXN(*),X(*),C1(*),C2(*),PA(*)
      REAL*8 OVL
      Character(LEN=16), Parameter :: ROUTINE='COVLP   '
      Integer iPrLev
      REAL*8 C1C2, FAC, OVLADD, PRQS, TERM
      REAL*8, External:: DDot_
      Integer I, iAshI, iAshJ, iC1, iC2, ISTBM, ISTC2, iSTIA, ISYM,     &
     &        JSYM, NAE, NAEJ, NAO, NAOJ, NEO, NI, NIA, NIAJ, NIO,      &
     &        NIOJ, NP, NQ, NT, NTT, NTUT, NTUVX, NU, NUT, NV, NVT,     &
     &        NVXT, NX, NXT

      IPRLEV=IPRLOC(4)
      IF(IPRLEV.ge.DEBUG) THEN
       WRITE(LF,*)' Entering ',ROUTINE
      END IF
!PAM02 Note structure of SX-vectors: First NROOT elements are special.
!PAM02 Elements NROOT+1,..,NROOT+NSXS contain the usual SX elements.
!PAM02 NROOT=1 always right now. Part of the code is prepared for using
!PAM02 several roots, so most of the code must use the general case.
      OVL=0.0D0
      DO I=1,NROOT
        OVL=OVL+C1IN(I)*C2IN(I)
      END DO

!PAM01 Adding overlap from small shift of SX overlaps:
      OVL=OVL+(1.0D-6)*DDOT_(NSXS,C1IN(NROOT+1),1,C2IN(NROOT+1),1)
!
! renormalize the C vector (simple element-by-element scaling).
!
      DO I=1,NSXS
       C1(I)=SXN(I)*C1IN(I+NROOT)
       C2(I)=SXN(I)*C2IN(I+NROOT)
      END DO

      ISTIA=0
      ISTBM=0
      IASHI=0
      DO ISYM=1,NSYM
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NIA=NIO+NAO
       NEO=NSSH(ISYM)
       NAE=NAO+NEO
       IF(NIA.EQ.0.OR.NAE.EQ.0) GO TO 97
!
! p is secondary (p = q = a)
!
       IF(NEO.NE.0) THEN
        CALL DGEMM_('N','N',                                            &
     &              NIA,NEO,NIA,                                        &
     &              1.0d0,DIA(ISTIA+1),NIA,                             &
     &              C1(ISTBM+1+NIA*NAO),NIA,                            &
     &              0.0d0,X,NIA)
        OVLADD=DDOT_(NIA*NEO,X,1,C2(ISTBM+1+NIA*NAO),1)
        OVL=OVL+OVLADD
       ENDIF
  97   CONTINUE
       ISTIA=ISTIA+NIA**2
       ISTBM=ISTBM+NIA*NAE
       IASHI=IASHI+NAO
      END DO

! A very long loop over  symmetry
      ISTIA=0
      ISTBM=0
      IASHI=0
      DO ISYM=1,NSYM
       NIO=NISH(ISYM)
       NAO=NASH(ISYM)
       NIA=NIO+NAO
       NEO=NSSH(ISYM)
       NAE=NAO+NEO
!
! r is inactive (r = s = i); p and q are active
!
       IF(NIO.NE.0.AND.NAO.NE.0) THEN
        IC1=ISTBM
        DO NP=NIO+1,NIA
         IC2=ISTBM
         DO NQ=NIO+1,NIA
          C1C2=0.0D0
          DO NI=1,NIO
           C1C2=C1C2+C1(IC1+NI)*C2(IC2+NI)
          END DO
          FAC=-DIA(ISTIA+NIA*(NP-1)+NQ)
          IF(NP.EQ.NQ) FAC=FAC+2.0D0
          OVLADD=C1C2*FAC
          OVL=OVL+OVLADD
          IC2=IC2+NIA
         END DO
         IC1=IC1+NIA
        END DO
       ENDIF
!
! r,s active and p,q active  (p,r=t,u; q,s=v,x)
!
       DO NT=2,NAO
        NTT=NT+IASHI
        DO NU=1,NT-1
         NUT=NU+IASHI
         NTUT=ITRI(NTT)+NUT

         IASHJ=0
         TERM=0.0D0
         ISTC2=0
         DO JSYM=1,NSYM
          NAOJ=NASH(JSYM)
          NIOJ=NISH(JSYM)
          NIAJ=NIOJ+NAOJ
          NAEJ=NAOJ+NSSH(JSYM)
          IF(NAOJ.GT.1) THEN
           IF(JSYM.EQ.ISYM) THEN
!--------
            DO NV=2,NAOJ
             NVT=NV+IASHJ
             DO NX=1,NV-1
              NXT=NX+IASHJ
              NVXT=ITRI(NVT)+NXT
              NTUVX=ITRI(MAX(NTUT,NVXT))+MIN(NTUT,NVXT)
              PRQS=-4.0D0*PA(NTUVX)
              IF(NU.EQ.NX) PRQS=PRQS+DIA(ISTIA+NIA*(NT+NIO-1)+NV+NIO)
              IF(NT.EQ.NV) PRQS=PRQS+DIA(ISTIA+NIA*(NU+NIO-1)+NX+NIO)
              IF(NT.EQ.NX) PRQS=PRQS-DIA(ISTIA+NIA*(NU+NIO-1)+NV+NIO)
              IF(NU.EQ.NV) PRQS=PRQS-DIA(ISTIA+NIA*(NT+NIO-1)+NX+NIO)
              TERM=TERM+PRQS*C2(ISTC2+NIAJ*(NV-1)+NIOJ+NX)
             END DO
            END DO
!--------
           ELSE
            DO NV=2,NAOJ
             NVT=NV+IASHJ
             DO NX=1,NV-1
              NXT=NX+IASHJ
              NVXT=ITRI(NVT)+NXT
              NTUVX=ITRI(MAX(NTUT,NVXT))+MIN(NTUT,NVXT)
              PRQS=-4.0D0*PA(NTUVX)
              TERM=TERM+PRQS*C2(ISTC2+NIAJ*(NV-1)+NIOJ+NX)
             END DO
            END DO
           ENDIF
!--------
          ENDIF
          ISTC2=ISTC2+NIAJ*NAEJ
          IASHJ=IASHJ+NAOJ
         END DO
         OVLADD=C1(ISTBM+NIA*(NT-1)+NIO+NU)*TERM
         OVL=OVL+OVLADD
        END DO
       END DO

       ISTIA=ISTIA+NIA**2
       ISTBM=ISTBM+NIA*NAE
       IASHI=IASHI+NAO

! End of very long loop over  symmetry
      END DO
!
      IF(IPRLEV.GE.DEBUG) THEN
       Write(LF,'(1X,A,F15.9)') ' OVERLAP IN COVLP:',OVL
      END IF

      END SUBROUTINE COVLP
