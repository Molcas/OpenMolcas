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
      Subroutine tr2NsB(CMO,X1,X2,pqrs,TUrs,lBuf,MAXRS)
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
C     ********** IBM-3090 RELEASE 87 09 14 **********
C     Replace MXMA with DGEMM P-AA Malmqvist 1992-05-06.
C
*
*     This and tr2NsB routines trasform non-squared AO integrals. The
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
      Dimension X1(*),X2(*)
      Dimension PQRS(*),TURS(*)

      CALL QEnter('Tr2NsB')
      icc=NOCP*NOCQ*NOR*NOS

      If(ISP.gt.ISR)then


      NSYMP=NSYM*(NSYM+1)/2
      NOTU=NOCP*NOCQ
      IF(ISP.EQ.ISQ) NOTU=(NOCP**2+NOCP)/2

C SORT OF PARTIALLY TRANSFORMED INTEGRALS (TU/RS) ON UNIT LUHLF3
      IPQMX3=NBRS
      IF(NBRS*NOTU.GT.LTUPQ) THEN
       IPQMX3=LTUPQ/NOTU
c      WRITE(*,*)'OUT OF CORE SORT FOR INTEGRALS (TU/RS)',IPQMX3,nbrs
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,TURS,IPQMX3,IAD3S)
      ENDIF
      IAD3=0
      IOUT3=0

c MaxRS should be given
      IRS=0
      LRS=0
      NRS=0
      Kread=0
      Nread=NBRS/MaxRS
      Nrest=mod(NBRS,MaxRS)
      If(Nrest.eq.0)then
       Nrest=MaxRS
      Else
       Nread=Nread+1
      Endif

       If(icc.ne.0)then

c Loop over r,s pair
        Do NR=1,NBR
         NumRS=NBS
         If(ISR.eq.ISS)NumRS=NR
         Do NS=1,NumRS
          IRS=IRS+1

c Loop over p,q pair
          If(LRS.eq.NRS)then
           Kread=Kread+1
           IPQ=0
           LPQ=0
           NPQ=0
           iRc=0
           iOpt=1
           IRSST=1-NBRS
           Do NP=1,NBP
            NumPQ=NBQ
            If(ISP.eq.ISQ)NumPQ=NP
            Do NQ=1,NumPQ
             IPQ=IPQ+1

c  Read integrals (pq,rs)
             If(LPQ.eq.NPQ) then
              Call Rdord(iRc,iOpt,ISP,ISQ,ISR,ISS,X1,lBuf,nPQ)
              IF(IRC.GT.1) THEN
                WRITE(6,*)' ERROR RETURN CODE IRC=',IRC
                WRITE(6,*)' FROM RDORD, CALLED FROM TRA2.'
                CALL Abend
              END IF
              iOpt=2
              LPQ=0
              IRSST=1-NBRS
             Endif
             LPQ=LPQ+1
             IRSST=IRSST+NBRS
c Here, we have all the (pq,rs) for this p,q pair
c Copy X1 and construct (pq,rs) for this p,q pair
             Length=MaxRS
             If(Kread.eq.Nread)Length=Nrest
             call dcopy_(Length,X1(IRSST+MaxRS*(Kread-1)),1,
     &                         PQRS(IPQ),NBPQ)
c End of loop over p,q pair
            Enddo
           Enddo
           LRS=0
           NRS=Length
          Endif
          LRS=LRS+1

c Transfer for this r,s pair
          If(ISP.eq.ISQ)then
           Call Square(PQRS(NBPQ*(LRS-1)+1),X2,1,NBP,NBP)
c (pq,rs) -> (pU,rs)
           CALL DGEMM_('T','N',
     &                 NBP,NOCQ,NBQ,
     &                 1.0d0,X2,NBQ,
     &                 CMO(LMOQ2),NBQ,
     &                 0.0d0,X1,NBP)
c (pU,rs) -> (TU,rs)
           CALL MXMT(X1,       NBP,1,
     &               CMO(LMOP2),1,NBP,
     &               X2,
     &               NOCP,NBP)
          Else
           call dcopy_(NBPQ,PQRS(NBPQ*(LRS-1)+1),1,X2,1)
c (pq,rs) -> (pU,rs)
           CALL DGEMM_('T','N',
     &                 NBP,NOCQ,NBQ,
     &                 1.0d0,X2,NBQ,
     &                 CMO(LMOQ2),NBQ,
     &                 0.0d0,X1,NBP)
c (pU,rs) -> (TU,rs)
           CALL DGEMM_('T','N',
     &                 NOCQ,NOCP,NBP,
     &                 1.0d0,X1,NBP,
     &                 CMO(LMOP2),NBP,
     &                 0.0d0,X2,NOCQ)
          Endif
c Store buffer
          IF(IRS.GT.IPQMX3) THEN
           IRS=1
cvv           DO I=1,NOTU
cvv            CALL dDAFILE(LUHLF3,1,TURS(1+IPQMX3*(I-1)),IPQMX3,IAD3)
cvv           Enddo
         CALL dDAFILE(LUHLF3,1,TURS,IPQMX3*NOTU,IAD3)
          ENDIF
c Sorting
          call dcopy_(NOTU,X2,1,TURS(IRS),IPQMX3)
c End of loop over r,s pair
         Enddo
        Enddo
c Store the last buffer
        IF(IPQMX3.LT.NBRS) THEN
cvv         DO I=1,NOTU
cvv          CALL dDAFILE(LUHLF3,1,TURS(1+IPQMX3*(I-1)),IPQMX3,IAD3)
cvv         Enddo
        CALL dDAFILE(LUHLF3,1,TURS,IPQMX3*NOTU,IAD3)
        ENDIF
       Endif

       If(icc.ne.0)then
        ISPQRS=((ISP**2-ISP)/2+ISQ-1)*NSYMP+(ISR**2-ISR)/2+ISS
        IAD2M(1,ISPQRS)=IAD13
        ITU=0
c Loop over t,u pair
        Do NT=1,NOCP
         Num=NOCQ
         If(ISP.eq.ISQ)Num=NT
         Do NU=1,Num
          IPQST=1+NBRS*ITU
          ITU=ITU+1
c Read buffer
          IF(IPQMX3.LT.NBRS) THEN
           Call RBuf_tra2(LUHLF3,TURS,NBRS,IPQMX3,NOTU,ITU,IPQST,IAD3S)
          ENDIF
          If(ISR.eq.ISS)then
c  Square
           Call Square(TURS(IPQST),X2,1,NBR,NBR)
c  (TU,rs) -> (TU,sB)
           CALL DGEMM_('T','N',
     &                 NBR,NOS,NBS,
     &                 1.0d0,X2,NBS,
     &                 CMO(LMOS2),NBS,
     &                 0.0d0,X1,NBR)
c  (TU,sB) -> (TU,AB)
           CALL MXMT(X1,       NBR,1,
     &               CMO(LMOR2),1,NBR,
     &               X2,
     &               NOR,NBR)
           IX2=(NOR*NOR+NOR)/2
          Else
           call dcopy_(NBRS,TURS(IPQST),1,X2,1)
c  (TU,rs) -> (TU,As)
           CALL DGEMM_('T','N',
     &                 NBR,NOS,NBS,
     &                 1.0d0,X2,NBS,
     &                 CMO(LMOS2),NBS,
     &                 0.0d0,X1,NBR)
c  (TU,As) -> (TU,AB)
           CALL DGEMM_('T','N',
     &                 NOS,NOR,NBR,
     &                 1.0d0,X1,NBR,
     &                 CMO(LMOR2),NBR,
     &                 0.0d0,X2,NOS)
           IX2=NOR*NOS
          Endif
c  Store (TU,AB) of this t,u pair
C
C       WRITE THESE BLOCK OF INTEGRALS ON LUINTM
C
          CALL GADSum(X2,IX2)
          CALL dDAFILE(LUINTM,1,X2,IX2,IAD13)
c End of Loop over t,u pair
         Enddo
        Enddo
       Endif
      Endif

      CALL QEXIT('Tr2NsB')
      Return
      End
