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
      Subroutine tr2NsA1(CMO,X1,nX1,X2,nX2,X3,nX3,
     &      pqUS,npqUS, pqRU,npqRU, pqTU,npqTU, lBuf)
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
      Dimension X1(nX1),X2(nX2),X3(nX3)
      Dimension PQTU(nPQTU),pqRU(npqRU),pqUS(npqUS)

      Call Qenter('Tr2NsA1')
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
cvv prevent integer overflow
      IF(1.0d0*NBPQ*NOUS.GT.LURPQ) THEN
       IPQMX1=LURPQ/NOUS
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|US)',IPQMX1
       IAD1S=0
       CALL dDAFILE(LUHLF1,0,pqUS,IPQMX1,IAD1S)
      ENDIF
      IAD1=0
      IOUT1=0
C     2. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|RU) ON UNIT LUHLF2
      IPQMX2=NBPQ
cvv prevent integer overflow
      IF(1.0d0*NBPQ*NORU.GT.LRUPQ) THEN
       IPQMX2=LRUPQ/NORU
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|RU)',IPQMX2
       IAD2S=0
       CALL dDAFILE(LUHLF2,0,pqRU,IPQMX2,IAD2S)
      ENDIF
      IAD2=0
      IOUT2=0
C     3. SORT OF PARTIALLY TRANSFORMED INTEGRALS (PQ|TU) ON UNIT LUHLF3
      IPQMX3=NBPQ
cvv prevent integer overflow
      IF(1.0d0*NBPQ*NOTU.GT.LTUPQ) THEN
       IPQMX3=LTUPQ/NOTU
c      WRITE(*,*) 'OUT OF CORE SORT FOR INTEGRALS (PQ|TU)',IPQMX3
       IAD3S=0
       CALL dDAFILE(LUHLF3,0,PQTU,IPQMX3,IAD3S)
      ENDIF
      IAD3=0
      IOUT3=0
c
c First half transformation
c
      LPQ=0
      NPQ=0
      iRc=0
      iOpt=1
      IRSST=1-NBRS
c Loop over p,q symmetry pair, if(ISP.eq.ISQ) loop should be triangle
      Do IP=1,NBP
       Num=NBQ
       if(ISP.eq.ISQ) Num=IP
       Do IQ=1,Num
        IOUT1=IOUT1+1
        IOUT2=IOUT2+1
        IOUT3=IOUT3+1
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
c  Square if necessary
        If(ISR.eq.ISS)then
         Call Square(X1(IRSST),X2,1,NBS,NBS)
        Else
         call dcopy_(NBRS,X1(IRSST),1,X2,1)
        Endif

c====================================================
c  First half transformation to (pq,Us), if ISR.ne.ISS
c  For exchange case3,4 (AT,UB), case7,8 (TA,UB)
c====================================================
        If(icxc3.ne.0.or.icxc7.ne.0)then
         If(ISR.ne.ISS)then
c  (pq,rs) -> (pq,Us)
          CALL DGEMM_('N','N',
     &                NBS,NOCR,NBR,
     &                1.0d0,X2,NBS,
     &                CMO(LMOR2),NBR,
     &                0.0d0,X3,NBS)
c  (pq,Us) Sorting
          IF(IOUT1.GT.IPQMX1) THEN
           IOUT1=1
cvv           DO I=1,NOUS
cvv            CALL dDAFILE(LUHLF1,1,pqUS(1+IPQMX1*(I-1)),IPQMX1,IAD1)
cvv           Enddo
         CALL dDAFILE(LUHLF1,1,pqUS,IPQMX1*NOUS,IAD1)
          ENDIF
          CALL DCOPY_(NOUS,X3,1,pqUS(IOUT1),IPQMX1)
         Endif
        Endif
c
c  First half transformation to (pq,rU)
c  For coulomb (AB,TU), exchange case1,2 (AT,BU), case5,6 (TA,BU)
c
        If(icc.ne.0.or.icxc1.ne.0.or.icxc5.ne.0)then
c  (pq,rs) -> (pq,rU)
         CALL DGEMM_('T','N',
     &               NBR,NOCS,NBS,
     &               1.0d0,X2,NBS,
     &               CMO(LMOS2),NBS,
     &               0.0d0,X3,NBR)
c  (pq,rU) Sorting
         IF(IOUT2.GT.IPQMX2) THEN
          IOUT2=1
cvv          DO I=1,NORU
cvv           CALL dDAFILE(LUHLF2,1,pqRU(1+IPQMX2*(I-1)),IPQMX2,IAD2)
cvv          Enddo
        CALL dDAFILE(LUHLF2,1,pqRU,IPQMX2*NORU,IAD2)
         ENDIF
         CALL DCOPY_(NORU,X3,1,pqRU(IOUT2),IPQMX2)

c  First half transformation to (pq,TU), if(ISR.eq.ISS)then triangle
c  (pq,rU) -> (pq,TU)
         If(icc.ne.0)then
          If(ISR.eq.ISS)then
           CALL MXMT(X3,        NBR,1,
     &               CMO(LMOR2),1,NBR,
     &               X2,
     &               NOCR,NBR)
          Else
           CALL DGEMM_('T','N',
     &                 NOCS,NOCR,NBR,
     &                 1.0d0,X3,NBR,
     &                 CMO(LMOR2),NBR,
     &                 0.0d0,X2,NOCS)
          Endif
c  (pq,TU) Sorting
          IF(IOUT3.GT.IPQMX3) THEN
           IOUT3=1
cvv           DO I=1,NOTU
cvv            CALL dDAFILE(LUHLF3,1,PQTU(1+IPQMX3*(I-1)),IPQMX3,IAD3)
cvv           Enddo
         CALL dDAFILE(LUHLF3,1,PQTU,IPQMX3*NOTU,IAD3)
          ENDIF
          CALL DCOPY_(NOTU,X2,1,PQTU(IOUT3),IPQMX3)
         Endif
        Endif
       Enddo
      Enddo
c  Store last buffer
      IF(IPQMX1.LT.NBPQ) THEN
cvv       DO I=1,NOUS
cvv        CALL dDAFILE(LUHLF1,1,pqUS(1+IPQMX1*(I-1)),IPQMX1,IAD1)
cvv       Enddo
      CALL dDAFILE(LUHLF1,1,pqUS,IPQMX1*NOUS,IAD1)
      ENDIF
      IF(IPQMX2.LT.NBPQ) THEN
cvv       DO I=1,NORU
cvv        CALL dDAFILE(LUHLF2,1,pqRU(1+IPQMX2*(I-1)),IPQMX2,IAD2)
cvv       Enddo
      CALL dDAFILE(LUHLF2,1,pqRU,IPQMX2*NORU,IAD2)
      ENDIF
      IF(IPQMX3.LT.NBPQ) THEN
cvv       DO I=1,NOTU
cvv        CALL dDAFILE(LUHLF3,1,PQTU(1+IPQMX3*(I-1)),IPQMX3,IAD3)
cvv       Enddo
      CALL dDAFILE(LUHLF3,1,PQTU,IPQMX3*NOTU,IAD3)
      ENDIF

      CALL QEXIT('Tr2NsA1')
      Return
      End
