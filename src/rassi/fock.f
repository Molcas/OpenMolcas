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
      SUBROUTINE FOCK_RASSI(DINAO,FOCKAO)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rassi.fh"
#include "WrkSpc.fh"
#include "symmul.fh"
      DIMENSION FOCKAO(NBSQ),DINAO(NBSQ)
      DIMENSION KEEP(8),NBSX(8)
      DIMENSION NBTRPR(8)
      LOGICAL   ISQARX
* Purpose:
*  ADD THE TWO-ELECTRON PART OF A FOCK MATRIX USING THE
*  INACTIVE DENSITY MATRIX. THE FOCK MATRIX SHOULD CONTAIN THE
*  ONE-ELECTRON HAMILTONIAN MATRIX BEFORE THE CALL. THE MATRICES
*  ARE STORED IN SYMMETRY-BLOCKED SQUARE FORMAT.

C RETRIEVE BASE DATA FROM UNIT LUORD:
      IRC=0
      CALL GETORD(IRC,ISQARX,NSYMX,NBSX,KEEP)
C ALLOCATE WORK SPACE FOR FOLDED DENSITY MATRIX, FOR ONE
C TRIANGULAR SYMMETRY-BLOCK OF THE TWO-ELECTRON CONTRIBUTION,
C FOR AN INTEGRAL BUFFER, AND FOR EXPANDING TRIANGULAR INTEGRAL
C MATRICES INTO SQUARE FORMAT:
      NSQBUF=NBMX**2
      INTBUF=MAX(NSQBUF,256*256)
      CALL GETMEM('PQRS  ','ALLO','REAL',LPQRS,INTBUF)
      CALL GETMEM('DTRI  ','ALLO','REAL',LDTRI,NBTRI)
      CALL GETMEM('SQBUF ','ALLO','REAL',LSQBUF,NSQBUF)
CPAM00 Nr of triangular matrices in previous symmetry blocks
      NBTR=0
      DO ISYM=1,NSYM
        NBTRPR(ISYM)=NBTR
        NB=NBASF(ISYM)
        NBTR=NBTR+(NB*(NB+1))/2
      END DO
      NFTRI=NBTR
      CALL GETMEM('FTRI  ','ALLO','REAL',LFTRI,NFTRI)
      CALL DCOPY_(NFTRI,0.0D0,0,WORK(LFTRI),1)
C FOLD THE D MATRIX:
      IDF=LDTRI
      DO ISYM=1,NSYM
       NB=NBASF(ISYM)
       NPQ=NBSQPR(ISYM)
       DO NP=1,NB
         DO NQ=1,NP
           DF=DINAO(NPQ+NB*(NP-1)+NQ)+DINAO(NPQ+NB*(NQ-1)+NP)
           IF(NP.EQ.NQ) DF=0.5D0*DF
           WORK(IDF)=DF
           IDF=IDF+1
         END DO
       END DO
      END DO
      CONTINUE
C LOOP OVER ISP:
      DO 104 ISP=1,NSYM
        NBP=NBASF(ISP)
        NIP=NISH(ISP)
        KEEPP=KEEP(ISP)
C LOOP OVER ISQ:
        DO 103 ISQ=1,ISP
          NBQ=NBASF(ISQ)
          NIQ=NISH(ISQ)
          NSPQ=MUL(ISP,ISQ)
          KEEPQ=KEEP(ISQ)
C LOOP OVER ISR:
          DO 102 ISR=1,ISP
            NBR=NBASF(ISR)
            NIR=NISH(ISR)
            NSPQR=MUL(NSPQ,ISR)
            KEEPR=KEEP(ISR)
C LOOP OVER ISS:
            ISSMX=ISR
            IF(ISP.EQ.ISR) ISSMX=ISQ
            DO 101 ISS=1,ISSMX
              IF(ISS.NE.NSPQR) GO TO 101
              NBS=NBASF(ISS)
              NIS=NISH(ISS)
              KEEPS=KEEP(ISS)
              KEEPT=KEEPP+KEEPQ+KEEPR+KEEPS
              NBT=NBP*NBQ*NBR*NBS
              IF(KEEPT.NE.0.AND.NBT.NE.0) THEN
                CALL ABEND()
              ENDIF
          IF(NBT.EQ.0) GOTO 101
          INUSE=0
          IF((ISP.EQ.ISQ).AND.(NBP.GT.0).AND.(NIR.GT.0)) INUSE=1
          IF((ISP.EQ.ISQ).AND.(NIP.GT.0).AND.(NBR.GT.0)) INUSE=1
          IF((ISP.EQ.ISR).AND.(NBP.GT.0).AND.(NIQ.GT.0)) INUSE=1
          IF((ISP.EQ.ISR).AND.(NIP.GT.0).AND.(NBQ.GT.0)) INUSE=1
          IF(INUSE.EQ.0) GOTO 101

C SIZES OF INTEGRAL MATRICES:
          NBPQ=NBP*NBQ
          IF(ISP.EQ.ISQ) NBPQ=(NBP**2+NBP)/2
          NBRS=NBR*NBS
          IF(ISR.EQ.ISS) NBRS=(NBR**2+NBR)/2
C READ THE AO INTEGRALS INTO WORK(LPQRS) WHENEVER NEEDED.
C ACCUMULATE CONTRIBUTIONS TO FOCK MATRIX INTO WORK(LFTRI), OR
C DIRECTLY INTO FOCKAO.
          IRC=0
          IOPT=1
          IPQ=0
          LPQ=0
          NPQ=0
          IRSST=LPQRS-NBRS
          DO 40 NP=1,NBP
            NQM=NBQ
            IF(ISP.EQ.ISQ) NQM=NP
            DO 30 NQ=1,NQM
              IPQ=IPQ+1
C READ A NEW BUFFER OF INTEGRALS
              IF(LPQ.EQ.NPQ) THEN
                CALL RDORD(IRC,IOPT,ISP,ISQ,ISR,ISS,
     *                     WORK(LPQRS),INTBUF,NPQ)
                IOPT=2
                LPQ=0
                IRSST=LPQRS-NBRS
C COULOMB CONTRIBUTION:
                IF(ISP.EQ.ISQ) THEN
                 IF(NIR.GT.0) THEN
                  IFPQ=NBTRPR(ISP)+IPQ
                  IDRS=NBTRPR(ISR)+1
                  NPQM=MIN(NPQ,NBPQ-IPQ+1)
                  CALL DGEMV_('T',NBRS,NPQM,1.0D0,WORK(LPQRS),NBRS,
     *                  WORK(LDTRI-1+IDRS),1,1.0D0,WORK(LFTRI-1+IFPQ),1)
                 ENDIF
CPAM00 Added code, to obtain RSPQ contributions from the PQRS integrals
                 IF(ISP.GT.ISR .AND. NIP.GT.0) THEN
                  IFRS=NBTRPR(ISR)+1
                  IDPQ=NBTRPR(ISP)+IPQ
                  NPQM=MIN(NPQ,NBPQ-IPQ+1)
                  CALL DGEMV_('N',NBRS,NPQM,1.0D0,WORK(LPQRS),NBRS,
     *                  WORK(LDTRI-1+IDPQ),1,1.0D0,WORK(LFTRI-1+IFRS),1)
                 END IF
                ENDIF
              ENDIF
              LPQ=LPQ+1
              IRSST=IRSST+NBRS
C EXCHANGE CONTRIBUTIONS:
              IF(ISP.EQ.ISR) THEN
                IF(ISR.EQ.ISS)
     *              CALL SQUARE(WORK(IRSST),WORK(LSQBUF),1,NBR,NBR)
                IF((ISP.NE.ISQ.OR.NP.NE.NQ).AND.(NIR.GT.0)) THEN
                  IF(ISR.EQ.ISS) THEN
                    CALL DGEMV_('N',NBS,NBR,-0.5D0,WORK(LSQBUF),NBS,
     *                   DINAO(NBSQPR(ISP)+NBR*(NP-1)+1),1,
     *                   1.0D0,FOCKAO(NBSQPR(ISQ)+NQ),NBQ)
                  ELSE
                    CALL DGEMV_('N',NBS,NBR,-0.5D0,WORK(IRSST),NBS,
     *                   DINAO(NBSQPR(ISP)+NBR*(NP-1)+1),1,
     *                   1.0D0,FOCKAO(NBSQPR(ISQ)+NQ),NBQ)
                  ENDIF
                ENDIF
                IF(ISP.EQ.ISS) THEN
                  IF(NIR.GT.0) THEN
                    CALL DGEMV_('N',NBS,NBR,-0.5D0,WORK(LSQBUF),NBS,
     *                    DINAO(NBSQPR(ISQ)+NBR*(NQ-1)+1),1,
     *                    1.0D0,FOCKAO(NBSQPR(ISP)+NP),NBP)
                  END IF
                ELSE
                  IF(NIS.GT.0) THEN
                    CALL DGEMV_('T',NBS,NBR,-0.5D0,WORK(IRSST),NBS,
     *                     DINAO(NBSQPR(ISQ)+NBS*(NQ-1)+1),1,
     *                     1.0D0,FOCKAO(NBSQPR(ISP)+NP),NBP)
                   END IF
                 ENDIF
               ENDIF
30          CONTINUE
40        CONTINUE
C
101      CONTINUE
102     CONTINUE
103    CONTINUE
C -- ADD COULOMB CONTRIBUTIONS FROM TRIANGULAR STORAGE TO FOCKAO:
104   CONTINUE
C END OF QUADRUPLE SYMMETRY LOOP.
CPAM00 OK -- Now add Coulomb contributions:
        IFPQ=0
        DO ISP=1,NSYM
          NBP=NBASF(ISP)
          DO NP=1,NBP
            DO NQ=1,NP
              IFPQ=IFPQ+1
              FPQ=WORK(LFTRI-1+IFPQ)
              II=NBSQPR(ISP)+NBP*(NP-1)+NQ
              FOCKAO(II)=FOCKAO(II)+FPQ
              IF(NQ.NE.NP) THEN
                II=NBSQPR(ISP)+NBP*(NQ-1)+NP
                FOCKAO(II)=FOCKAO(II)+FPQ
              END IF
            END DO
          END DO
        END DO
      CALL GETMEM('      ','FREE','REAL',LPQRS,INTBUF)
      CALL GETMEM('      ','FREE','REAL',LDTRI,NBTRI)
      CALL GETMEM('      ','FREE','REAL',LFTRI,NFTRI)
      CALL GETMEM('      ','FREE','REAL',LSQBUF,NSQBUF)
*
      Call GADSum(FOCKAO,NBSQ)
      RETURN
      END
