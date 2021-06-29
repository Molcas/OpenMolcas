!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1987, Bjorn O. Roos                                    *
!***********************************************************************
!--------------------------------------------*
! 1987  B. O. ROOS                           *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND, SWEDEN                 *
!--------------------------------------------*
      Subroutine tr2NsA2(CMO,X1,nX1,X2,nX2,                             &
     &                   pqrU,npqrU, pqTU,npqTU)
!
! SECOND ORDER TWO-ELECTRON TRANSFORMATION ROUTINE
!
! THIS ROUTINE IS CALLED FOR EACH SYMMETRY BLOCK OF INTEGRALS
! (ISP,ISQ,ISR,ISS) WITH ISP.GE.ISQ AND ISR.GE.ISS.
! P,Q,R,S are SO indices.
! A,B are MO indices, counting only non-frozen and non-deleted.
! T,U are occupied MO indices, only non-frozen and non-deleted.
! INTEGRALS (AB/TU) ARE ALWAYS GENERATED
! EXCHANGE INTEGRALS (AT/BU) ARE GENERATED AS FOLLOWS:
! (AT/BU) IF ISP.GE.ISR
! (AT/UB) IF ISP.GT.ISS AND ISP.NE.ISQ
! (TA/BU) IF ISQ.GT.ISR AND ISP.NE.ISQ
! (TA/UB) IF ISQ.GE.ISS AND ISP.NE.ISQ
!
!     This and tr2NsB routines transform non-squared AO integrals. The
!     transformed MO integrals are stored as the same as Tr2Sq
!     subroutine does.
!
      Implicit real*8(a-h,o-z)
!PAM98      COMMON/INTTRA/ISP,ISQ,ISR,ISS,NBP,NBQ,NBR,NBS,NBPQ,NBRS,IRRST,
!PAM98     &              NOCP,NOCQ,NOCR,NOCS,NPQ,LADX,LRUPQ,LURPQ,LTUPQ,
!PAM98     &              NOP,NOQ,NOR,NOS,LMOP,LMOQ,LMOR,LMOS,LMOP2,LMOQ2,
!PAM98     &              LMOR2,LMOS2,IAD13,ITP,ITQ,ITR,ITS
#include "rasdim.fh"
#include "caspt2.fh"
#include "trafo.fh"
#include "intgrl.fh"
#include "SysDef.fh"
      DIMENSION CMO(NCMO)
      Dimension X1(nX1),X2(nX2)
      Dimension PQTU(nPQTU),PQRU(nPQRU)

      NSYMP=NSYM*(NSYM+1)/2
      NOTU=NOCR*NOCS
      IF(ISR.EQ.ISS) NOTU=(NOCR**2+NOCR)/2
      NORU=NBR*NOCS
      icc  =NOP*NOQ*NOCR*NOCS
      icxc1=NOP*NOCQ*NOR*NOCS
      icxc5=NOCP*NOQ*NOR*NOCS
!
! Check for in core or out of core transformation
!
!     2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
      IPQMX2=NBPQ
      IF(NBPQ*NORU.GT.LRUPQ) THEN
       IPQMX2=LRUPQ/NORU
!      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
       IAD2S=0
       CALL dDAFILE(LUHLF2,0,PQRU,IPQMX2,IAD2S)
      ENDIF
!     3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
      IPQMX3=NBPQ
      IF(NBPQ*NOTU.GT.LTUPQ) THEN
       IPQMX3=LTUPQ/NOTU
!      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,PQTU,IPQMX3,IAD3S)
      ENDIF
      IAD3=0
!
!====================================================
! Second half transformation
!====================================================
!-------------------------------------------
! Second half transformation (AB,TU) coulomb
!  Always calculated.
!-------------------------------------------
! Loop over t,u pair,
      If(icc.ne.0) then
       ISPQRS=((ISR**2-ISR)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISQ
       IAD2M(1,ISPQRS)=IAD13
       IPQTU=0
       ITU=0
       Do NT=1,NOCR
        Num=NOCS
        If(ISR.eq.ISS)Num=NT
        Do NU=1,Num
         IPQTU=1+NBPQ*ITU
         ITU=ITU+1
! Read sorted integral from disk
         IF(IPQMX3.LT.NBPQ) THEN
          Call RBuf_tra2(LUHLF3,PQTU,NBPQ,IPQMX3,NOTU,ITU,IPQTU,IAD3S)
         ENDIF
         If(ISP.eq.ISQ) then
!  Square if necessary
          CALL SQUARE(PQTU(IPQTU),X2,1,NBP,NBP)
!  (pq,TU) -> (Aq,TU)
          CALL DGEMM_('N','N',                                          &
     &                NBQ,NOP,NBP,                                      &
     &                1.0d0,X2,NBQ,                                     &
     &                CMO(LMOP),NBP,                                    &
     &                0.0d0,X1,NBQ)
!  (Aq,TU) -> (AB,TU)
          CALL MXMT(X1,       NBQ,1,                                    &
     &              CMO(LMOQ),1,NBQ,                                    &
     &              X2,                                                 &
     &              NOP,NBQ)
          IX2=(NOP+NOP**2)/2
         Else
!  (pq,TU) -> (Aq,TU)
          CALL DGEMM_('N','N',                                          &
     &                NBQ,NOP,NBP,                                      &
     &                1.0d0,PQTU(IPQTU),NBQ,                            &
     &                CMO(LMOP),NBP,                                    &
     &                0.0d0,X1,NBQ)
!  (Aq,TU) -> (AB,TU)
          CALL DGEMM_('T','N',                                          &
     &                NOQ,NOP,NBQ,                                      &
     &                1.0d0,CMO(LMOQ),NBQ,                              &
     &                X1,NBQ,                                           &
     &                0.0d0,X2,NOQ)
          IX2=NOP*NOQ
         Endif
!  Store (AB,TU) of this t,u pair
         Call GADSum(X2,IX2)
         Call dDAFILE(LUINTM,1,X2,IX2,IAD13)
        Enddo
       Enddo
! End of loop over t,u pair
      Endif
!-----------------------------------------------------------------------
! Second half transformation Case 1 (AT,BU)
!  Always calculated. Both type 1 and 2. Case 2 (BU,AT) can be abandoned
!  since always ISP.gt.ISR(equality can be removed by symmetry)
!-----------------------------------------------------------------------
! Loop over r,u pair
      NOTU=NOCQ*NOCS
      IF(ISQ.EQ.ISS) NOTU=(NOCQ**2+NOCQ)/2
      If(icxc1.ne.0)then
       LAR=LTUPQ/NOTU
       LR=LAR/NOP
       IF(LR.GT.NBR) LR=NBR
       LAR=NOP*LR
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,PQTU,LAR,IAD3S)
       IAD3=0
       IR=0
       Do NR=1,NBR
        IR=IR+1
        Do NU=1,NOCS
! Square if necessary
         IRU=NBR*(NU-1)+NR
         IPQST=1+NBPQ*(IRU-1)
         IF(IPQMX2.LT.NBPQ) THEN
          Call RBuf_tra2(LUHLF2,PQRU,NBPQ,IPQMX2,NORU,IRU,IPQST,IAD2S)
         ENDIF
         If(ISP.eq.ISQ)then
          Call Square(PQRU(IPQST),X2,1,NBP,NBP)
         Else
          call dcopy_(NBPQ,PQRU(IPQST),1,X2,1)
         Endif
!  If(ISQ.eq.ISS)then triangular
         If(ISQ.eq.ISS)then
!  (pq,rU) -> (pT,rU)
          CALL DGEMM_('T','N',                                          &
     &                NBP,NOCQ-NU+1,NBQ,                                &
     &                1.0d0,X2,NBQ,                                     &
     &                CMO(LMOQ2+NBQ*(NU-1)),NBQ,                        &
     &                0.0d0,X1,NBP)
!  (pT,rU) -> (AT,rU)
          CALL DGEMM_('T','N',                                          &
     &                NOCQ-NU+1,NOP,NBP,                                &
     &                1.0d0,X1,NBP,                                     &
     &                CMO(LMOP2),NBP,                                   &
     &                0.0d0,X2,NOCQ-NU+1)
         Else
!  (pq,rU) -> (pT,rU)
          CALL DGEMM_('T','N',                                          &
     &                NBP,NOCQ,NBQ,                                     &
     &                1.0d0,X2,NBQ,                                     &
     &                CMO(LMOQ2),NBQ,                                   &
     &                0.0d0,X1,NBP)
!  (pT,rU) -> (AT,rU)
          CALL DGEMM_('T','N',                                          &
     &                NOCQ,NOP,NBP,                                     &
     &                1.0d0,X1,NBP,                                     &
     &                CMO(LMOP),NBP,                                    &
     &                0.0d0,X2,NOCQ)
         Endif
!  Store buffer
         IF(IR.GT.LR) THEN
          IR=1
!vv          DO I=1,NOTU
!vv           CALL dDAFILE(LUHLF3,1,PQTU(1+LAR*(I-1)),LAR,IAD3)
!vv          End do
         CALL dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
         ENDIF
!  Sort
         NAT=0
         DO NA=1,NOP
          NTM=1
          IF(ISQ.EQ.ISS) NTM=NU
          DO NT=NTM,NOCQ
           ITU=NOCS*(NT-1)+NU-1
           IF(ISQ.LT.ISS) ITU=NOCQ*(NU-1)+NT-1
           IF(ISQ.EQ.ISS) ITU=(NT**2-NT)/2+NU-1
           NAT=NAT+1
           PQTU(LAR*ITU+NOP*(IR-1)+NA)=X2(NAT)
!          if(isp.eq.7.and.isq.eq.1.and.isr.eq.6)
!    &     write(6,'(d13.6)')pqtu(1)
          Enddo
         Enddo
! End of loop over r,u pair
        Enddo
       Enddo
! Store last buffer
       IF(LR.LT.NBR) THEN
!vv        DO I=1,NOTU
!vv         CALL dDAFILE(LUHLF3,1,PQTU(1+LAR*(I-1)),LAR,IAD3)
!vv        Enddo
       CALL dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
       ENDIF
! Transform fourth index
       IF(ISQ.GE.ISS) THEN
! Exchange type 1
        ISPQRS=((ISQ**2-ISQ)/2+ISS-1)*NSYMP+(ISP**2-ISP)/2+ISR
        IAD2M(2,ISPQRS)=IAD13
        NTMAX=NOCQ
        NUMAX=NOCS
       ELSE
! Exchange type 2
        ISPQRS=((ISS**2-ISS)/2+ISQ-1)*NSYMP+(ISP**2-ISP)/2+ISR
        IAD2M(3,ISPQRS)=IAD13
        NTMAX=NOCS
        NUMAX=NOCQ
       ENDIF
! Loop over t,u pair, If(ISQ.eq.ISS) loop should be triangle
       IST=1-NOP*NBR
       KKTU=0
       DO NT=1,NTMAX
        NUM=NUMAX
        IF(ISQ.EQ.ISS) NUM=NT
        DO NU=1,NUM
         IST=IST+NOP*NBR
         KKTU=KKTU+1
         IF(LR.LT.NBR)THEN
          Call RBuf_tra2(LUHLF3,PQTU,NBR*NOP,LAR,NOTU,KKTU,IST,IAD3S)
         ENDIF
!  (AT,rU) -> (AT,BU)
         CALL DGEMM_('T','T',                                           &
     &               NOR,NOP,NBR,                                       &
     &               1.0d0,CMO(LMOR2),NBR,                              &
     &               PQTU(IST),NOP,                                     &
     &               0.0d0,X2,NOR)
!
!       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
!
         Call GADSum(X2,NOP*NOR)
         Call dDAFILE(LUINTM,1,X2,NOP*NOR,IAD13)
        Enddo
       Enddo
! End of loop over t,u pair
      Endif
!-----------------------------------------------------------------------
! Case 5 (TA,BU) and 6 (UB,AT)
!  Calculated if ISP.ne.ISQ.or.(ISQ.lt.ISR.and.ISP.eq.ISR)
!  If ISQ.ge.ISR then Case 5, which alwats gives type 1 integrals
!  If ISQ.lt.ISR.and.ISP.ne.ISR
!                then Case 6, which alwats gives type 2 integrals
!-----------------------------------------------------------------------
      NOTU=NOCP*NOCS
! ISP.eq.ISR in Case 6 should be skipped.
      If(ISQ.lt.ISR.and.ISP.eq.ISR)goto 100
      If(ISP.ne.ISQ.and.icxc5.NE.0)then
       LAR=LTUPQ/NOTU
       LR=LAR/NOQ
       IF(LR.GT.NBR) LR=NBR
       LAR=NOQ*LR
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,PQTU,LAR,IAD3S)
       IAD3=0
       IRU=0
       IR=0
! Loop over r,u pair
       Do NR=1,NBR
        IR=IR+1
        Do NU=1,NOCS
!  Square is unnecessary
         IRU=NBR*(NU-1)+NR
         IPQST=1+NBPQ*(IRU-1)
         IF(IPQMX2.LT.NBPQ) THEN
          Call RBuf_tra2(LUHLF2,PQRU,NBPQ,IPQMX2,NORU,IRU,IPQST,IAD2S)
         ENDIF
!  Always ISP.gt.ISS i.e. s(T).gt.s(U)
!  (pq,rU) -> (Tq,rU)
         CALL DGEMM_('N','N',                                           &
     &               NBQ,NOCP,NBP,                                      &
     &               1.0d0,PQRU(IPQST),NBQ,                             &
     &               CMO(LMOP2),NBP,                                    &
     &               0.0d0,X1,NBQ)
!  (Tq,rU) -> (TA,rU)
         CALL DGEMM_('T','N',                                           &
     &               NOCP,NOQ,NBQ,                                      &
     &               1.0d0,X1,NBQ,                                      &
     &               CMO(LMOQ2),NBQ,                                    &
     &               0.0d0,X2,NOCP)
! Store buffer
         IF(IR.GT.LR) THEN
          IR=1
!vv          DO I=1,NOTU
!vv           CALL dDAFILE(LUHLF3,1,PQTU(1+LAR*(I-1)),LAR,IAD3)
!vv          Enddo
         CALL dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
         ENDIF
! Sorting
         NAT=0
         DO NA=1,NOQ
          DO NT=1,NOCP
           ITU=NOCS*(NT-1)+NU-1
           NAT=NAT+1
           PQTU(LAR*ITU+NOQ*(IR-1)+NA)=X2(NAT)
          Enddo
         Enddo
! End of loop over r,u pair
        Enddo
       Enddo
! Store last buffer
       IF(LR.LT.NBR) THEN
!vv        DO I=1,NOTU
!vv         CALL dDAFILE(LUHLF3,1,PQTU(1+LAR*(I-1)),LAR,IAD3)
!vv        Enddo
       CALL dDAFILE(LUHLF3,1,PQTU,LAR*NOTU,IAD3)
       ENDIF
       If(ISQ.ge.ISR)then
!   Store(Only type1)
        ISPQRS=((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISQ**2-ISQ)/2+ISR
        IAD2M(2,ISPQRS)=IAD13
       Elseif(ISP.ne.ISR.and.ISR.gt.ISQ) then
!   Store(Only type2)
        ISPQRS=((ISP**2-ISP)/2+ISS-1)*NSYMP+(ISR**2-ISR)/2+ISQ
        IAD2M(3,ISPQRS)=IAD13
       Endif
! Loop over t,u pair
       IST=1-NOQ*NBR
       KKTU=0
       Do NT=1,NOCP
        Do NU=1,NOCS
         IST=IST+NOQ*NBR
         KKTU=KKTU+1
         IF(LR.LT.NBR)THEN
          Call RBuf_tra2(LUHLF3,PQTU,NBR*NOQ,LAR,NOTU,KKTU,IST,IAD3S)
         ENDIF
!  (TA,rU) -> (TA,BU)
         If(ISQ.ge.ISR)then
          CALL DGEMM_('T','T',                                          &
     &                NOR,NOQ,NBR,                                      &
     &                1.0d0,CMO(LMOR),NBR,                              &
     &                PQTU(IST),NOQ,                                    &
     &                0.0d0,X2,NOR)
         Else if(ISP.ne.ISR.and.ISR.gt.ISQ) then
          CALL DGEMM_('N','N',                                          &
     &                NOQ,NOR,NBR,                                      &
     &                1.0d0,PQTU(IST),NOQ,                              &
     &                CMO(LMOR),NBR,                                    &
     &                0.0d0,X2,NOQ)
         Endif
!
!       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
!
         CALL GADSum(X2,NOQ*NOR)
         Call dDAFILE(LUINTM,1,X2,NOQ*NOR,IAD13)
        Enddo
       Enddo
! End of loop over t,u pair
      Endif
  100 Continue

      Return
      End
