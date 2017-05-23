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
* Copyright (C) 1987, Bjorn O. Roos                                    *
************************************************************************
*--------------------------------------------*
* 1987  B. O. ROOS                           *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND, SWEDEN                 *
*--------------------------------------------*
      Subroutine tr2nsa3(CMO,X1,nX1,X2,nX2,pqUs,npqUS,
     &                   pqrU,npqrU)
C
C SECOND ORDER TWO-ELECTRON TRANSFORMATION ROUTINE
C
C THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
C (ISP,ISQ,ISR,ISS) WITH ISP.GE.ISQ AND ISR.GE.ISS.
C P,Q,R,S are SO indices.
C A,B are MO indices, counting only non-frozen and non-deleted.
C T,U are occupied MO indices, only non-frozen and non-deleted.
C INTEGRALS (AB/TU) ARE ALWAYS GENERATED
C EXCHANGE INTEGRALS (AT/BU) ARE GENERATED AS FOLLOWS:
C (AT/BU) IF ISP.GE.ISR
C (AT/UB) IF ISP.GT.ISS AND ISP.NE.ISQ
C (TA/BU) IF ISQ.GT.ISR AND ISP.NE.ISQ
C (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ
C
*     This and tr2NsB routines transform non-squared AO integrals. The
*     transformed MO integrals are stored as the same as Tr2Sq
*     subroutine does.
*
      Implicit real*8(a-h,o-z)
CPAM98      COMMON/INTTRA/ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NBPQ,NBRS,IRRST,
CPAM98     &              NOCP,NOCQ,NOCR,NOCS,NPQ,LADX,LRUPQ,LURPQ,LTUPQ,
CPAM98     &              NOP,NOQ,NOR,NOS,LMOP,LMOQ,LMOR,LMOS,LMOP2,LMOQ2,
CPAM98     &              LMOR2,LMOS2,IAD13,ITP,ITQ,ITR,ITS

#include "rasdim.fh"
#include "caspt2.fh"
#include "trafo.fh"
#include "intgrl.fh"

#include "SysDef.fh"
      DIMENSION CMO(NCMO)
      Dimension X1(nX1),X2(nX2)
      Dimension PQRU(nPQRU),PQUS(nPQUS)

      Call Qenter('tr2nsa3')
      NSYMP=NSYM*(NSYM+1)/2
      NOTU=NOCR*NOCS
      IF(ISR.EQ.ISS) NOTU=(NOCR**2+NOCR)/2
      NOUS=NOCR*NBS
      NORU=NBR*NOCS
      icc  =NOP*NOQ*NOCR*NOCS
      icxc1=NOP*NOCQ*NOR*NOCS
      icxc3=NOP*NOCQ*NOCR*NOS
      icxc5=NOCP*NOQ*NOR*NOCS
      icxc7=NOCP*NOQ*NOCR*NOS
c
c Check for in core or out of core transformation
c
C     1. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|US) ON UNIT LUHLF1
      IPQMX1=NBPQ
      IF(NBPQ*NOUS.GT.LURPQ) THEN
       IPQMX1=LURPQ/NOUS
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|US)',IPQMX1
       IAD1S=0
       CALL dDAFILE(LUHLF1,0,PQUS,IPQMX1,IAD1S)
      ENDIF
      IAD1=0
      IOUT1=0
C     2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
      IPQMX2=NBPQ
      IF(NBPQ*NORU.GT.LRUPQ) THEN
       IPQMX2=LRUPQ/NORU
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
       IAD2S=0
       CALL dDAFILE(LUHLF2,0,PQRU,IPQMX2,IAD2S)
      ENDIF
      IAD2=0
      IOUT2=0
*C     3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
*      IPQMX3=NBPQ
*      IF(NBPQ*NOTU.GT.LTUPQ) THEN
*       IPQMX3=LTUPQ/NOTU
*c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
*       IAD3S=0
*       CALL dDAFILE(LUHLF3,0,PQTU,IPQMX3,IAD3S)
*      ENDIF
*      IAD3=0
*      IOUT3=0
c
c-----------------------------------------------------------------------
c Second half transformation Case 3 (AT,UB)
c  Calculated if ISR.ne.ISS. Both type 1 and 2. ISQ.ne.ISR, since
c  equality makes integral to be zero symmetrically.
c  Case 4 (BU,TA) need not be calculated, since always ISP.gt.ISS.
c-----------------------------------------------------------------------
      NOTU=NOCQ*NOCR
      If(ISR.ne.ISS.and.icxc3.ne.0)then
       LAS=(LRUPQ+LTUPQ)/NOTU
       LS=LAS/NOP
       IF(LS.GT.NBS) LS=NBS
       LAS=NOP*LS
       IAD2S=0
       CALL dDAFILE(LUHLF2,0,PQRU,LAS,IAD2S)
       IAD2=0
c Loop over u,s pair
       IS=0
       Do NS=1,NBS
        IS=IS+1
        Do NU=1,NOCR
c  Square if necessary
         IUS=NBS*(NU-1)+NS
         IPQST=1+NBPQ*(IUS-1)
         IF(IPQMX1.LT.NBPQ) THEN
          Call RBuf_tra2(LUHLF1,PQUS,NBPQ,IPQMX1,NOUS,IUS,IPQST,IAD1S)
         ENDIF
         If(ISP.eq.ISQ)then
          Call Square(PQUS(IPQST),X2,1,NBP,NBP)
         Else
          call dcopy_(NBPQ,PQUS(IPQST),1,X2,1)
         Endif
c Always ISQ.ne.ISR  i.e. s(T).ne.s(U)
c  (pq,Us) -> (pT,Us)
          CALL DGEMM_('T','N',
     &                NBP,NOCQ,NBQ,
     &                1.0d0,X2,NBQ,
     &                CMO(LMOQ2),NBQ,
     &                0.0d0,X1,NBP)
c  (pT,Us) -> (AT,Us)
          CALL DGEMM_('T','N',
     &                NOCQ,NOP,NBP,
     &                1.0d0,X1,NBP,
     &                CMO(LMOP),NBP,
     &                0.0d0,X2,NOCQ)
c  Store buffer
        IF(IS.GT.LS) THEN
         IS=1
cvv         DO I=1,NOTU
cvv          CALL dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
cvv         Enddo
        CALL dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
        ENDIF
c  Sort
         NAT=0
         DO NA=1,NOP
          DO NT=1,NOCQ
           ITU=NOCR*(NT-1)+NU-1
           IF(ISQ.LT.ISR) ITU=NOCQ*(NU-1)+NT-1
           NAT=NAT+1
           PQRU(LAS*ITU+NOP*(IS-1)+NA)=X2(NAT)
          Enddo
         Enddo
c End of loop over u,s pair
        Enddo
       Enddo
c Store last buffer
       IF(LS.LT.NBS) THEN
cvv        DO I=1,NOTU
cvv         CALL dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
cvv        Enddo
       CALL dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
       ENDIF
c
c  Trasform last index
c
       If(ISQ.ge.ISR)then
c  Store(type 1)
        ISPQRS=((ISQ**2-ISQ)/2+ISR-1)*NSYMP+(ISP**2-ISP)/2+ISS
        IAD2M(2,ISPQRS)=IAD13
        NTMAX=NOCQ
        NUMAX=NOCR
       Else
c  Store(type 2)
        ISPQRS=((ISR**2-ISR)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISS
        IAD2M(3,ISPQRS)=IAD13
        NTMAX=NOCR
        NUMAX=NOCQ
       Endif
c Loop over t,u pair, always ISQ.ne.ISR
       IST=1-NBS*NOP
       KKTU=0
       Do NT=1,NTMAX
        Do NU=1,NUMAX
         IST=IST+NBS*NOP
         KKTU=KKTU+1
         IF(LS.LT.NBS)THEN
          Call RBuf_tra2(LUHLF2,PQRU,NBS*NOP,LAS,NOTU,KKTU,IST,IAD2S)
         ENDIF
c  (AT,Us) -> (AT,UB)
         CALL DGEMM_('T','T',
     &               NOS,NOP,NBS,
     &               1.0d0,CMO(LMOS2),NBS,
     &               PQRU(IST),NOP,
     &               0.0d0,X2,NOS)
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
         CALL GADSum(X2,NOP*NOS)
         Call dDAFILE(LUINTM,1,X2,NOP*NOS,IAD13)
        Enddo
       Enddo
c End of loop over t,u pair
      Endif
c-----------------------------------------------------------------------
c Case 7 (TA,UB) and 8 (UB,TA)
c  Calculated if ISP.ne.ISQ.and.ISR.ne.ISS
c   Case 7, if ISQ.ge.ISS, always type 1
c   Case 8, if ISQ.lt.ISS.and.ISP.ne.ISR, always type 2
c-----------------------------------------------------------------------
      NOTU=NOCP*NOCR
      If(ISS.gt.ISQ.and.ISP.eq.ISR)goto 200
      If(ISP.ne.ISQ.and.ISR.ne.ISS.and.icxc7.ne.0)then
       LAS=(LRUPQ+LTUPQ)/NOTU
       LS=LAS/NOQ
       IF(LS.GT.NBS) LS=NBS
       LAS=NOQ*LS
       IAD2S=0
       CALL dDAFILE(LUHLF2,0,PQRU,LAS,IAD2S)
       IAD2=0
       IRU=0
c Loop over u,s pair
       IS=0
       Do NS=1,NBS
        IS=IS+1
        Do NU=1,NOCR
         IRU=NBS*(NU-1)+NS
         IPQST=1+NBPQ*(IRU-1)
         IF(IPQMX1.LT.NBPQ) THEN
          Call RBuf_tra2(LUHLF1,PQUS,NBPQ,IPQMX1,NOUS,IRU,IPQST,IAD1S)
         ENDIF
c Always ISP.gt.ISQ
c Square unnecessary
         If(ISP.eq.ISR)then
c (pq,Us) -> (Tq,Us)
          CALL DGEMM_('N','N',
     &                NBQ,NOCP-NU+1,NBP,
     &                1.0d0,PQUS(IPQST),NBQ,
     &                CMO(LMOP2+NBP*(NU-1)),NBP,
     &                0.0d0,X1,NBQ)
c (Tq,Us) -> (TA,Us)
          CALL DGEMM_('T','N',
     &                NOCP-NU+1,NOQ,NBQ,
     &                1.0d0,X1,NBQ,
     &                CMO(LMOQ),NBQ,
     &                0.0d0,X2,NOCP-NU+1)
         Else
c (pq,Us) -> (Tq,Us)
          CALL DGEMM_('N','N',
     &                NBQ,NOCP,NBP,
     &                1.0d0,PQUS(IPQST),NBQ,
     &                CMO(LMOP2),NBP,
     &                0.0d0,X1,NBQ)
c (Tq,Us) -> (TA,Us)
          CALL DGEMM_('T','N',
     &                NOCP,NOQ,NBQ,
     &                1.0d0,X1,NBQ,
     &                CMO(LMOQ),NBQ,
     &                0.0d0,X2,NOCP)
         Endif
c Store buffer
         IF(IS.GT.LS) THEN
          IS=1
cvv          DO I=1,NOTU
cvv           CALL dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
cvv          ENDDO
         CALL dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
         ENDIF
c Sorting
* Note: LAS is supposed to be small enough that (NOCP*NOCR)*LAS is
* at most = (LRUPQ+LTUPQ). NU can be as large as NOCR, so ITU
* can become NOCR*NOCP-1.
* Thus (LAS*ITU+NOQ*(IS-1)+NA can become as large as:
* LAS*(NOCP*NOCR-1)+NOQ*(NBS-1)+NOQ=LAS*(NOCP*NOCR-1)+NOQ*NBS
* which may become nearly as large as (LRUPQ+LTUPQ)+NOQ*NBS.
* But when called, only the size LRUPQ has been reserved for
* use (tractl:LW6=LW5+LRUPQ) by the array PQRU!
* This may be deliberate, if it is intended that PQRU is
* overlaying i.e. extended above the LW6 address which is
* then assumed not to be used any longer....
         NAT=0
         DO NA=1,NOQ
          NTM=1
          IF(ISP.EQ.ISR) NTM=NU
          DO NT=NTM,NOCP
           ITU=NOCR*(NT-1)+NU-1
           If(ISP.eq.ISR)ITU=(NT*NT-NT)/2+NU-1
           NAT=NAT+1
        itst=LAS*ITU+NOQ*(IS-1)+NA
           PQRU(LAS*ITU+NOQ*(IS-1)+NA)=X2(NAT)
          Enddo
         Enddo
c End of Loop over u,s pair
        Enddo
       Enddo
c Store the last buffer
       IF(LS.LT.NBS) THEN
cvv        DO I=1,NOTU
cvv         CALL dDAFILE(LUHLF2,1,PQRU(1+LAS*(I-1)),LAS,IAD2)
cvv        Enddo
       CALL dDAFILE(LUHLF2,1,PQRU,LAS*NOTU,IAD2)
       ENDIF
       If(ISQ.ge.ISS) then
c   Store(Only type1)
        ISPQRS=((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISQ**2-ISQ)/2+ISS
        IAD2M(2,ISPQRS)=IAD13
       Elseif(ISP.ne.ISR.and.ISQ.lt.ISS) then
c   Store(Only type2)
        ISPQRS=((ISP**2-ISP)/2+ISR-1)*NSYMP+(ISS**2-ISS)/2+ISQ
        IAD2M(3,ISPQRS)=IAD13
       Endif
       IST=1-NOQ*NBS
       KKTU=0
c Loop over t,u pair, If(ISP.eq.ISR)then loop should be triangle
       Do NT=1,NOCP
        Num=NOCR
        If(ISP.eq.ISR)Num=NT
        Do NU=1,Num
c  (TA,Us) -> (TA,UB)
         IST=IST+NOQ*NBS
         KKTU=KKTU+1
         IF(LS.LT.NBS)THEN
          Call RBuf_tra2(LUHLF2,PQRU,NBS*NOQ,LAS,NOTU,KKTU,IST,IAD2S)
         ENDIF
         If(ISQ.ge.ISS)then
          CALL DGEMM_('T','T',
     &                NOS,NOQ,NBS,
     &                1.0d0,CMO(LMOS),NBS,
     &                PQRU(IST),NOQ,
     &                0.0d0,X2,NOS)
         Else if(ISP.ne.ISR.and.ISS.gt.ISQ) then
          CALL DGEMM_('N','N',
     &                NOQ,NOS,NBS,
     &                1.0d0,PQRU(IST),NOQ,
     &                CMO(LMOS),NBS,
     &                0.0d0,X2,NOQ)
         Endif
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
         Call GADSum(X2,NOQ*NOS)
         CALL dDAFILE(LUINTM,1,X2,NOQ*NOS,IAD13)
        Enddo
       Enddo
c End of Loop over t,u pair
      Endif
  200 Continue

      CALL QEXIT('tr2nsa3')
      Return
      End
