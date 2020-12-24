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
* Copyright (C) Per Ake Malmqvist                                      *
************************************************************************
      SUBROUTINE TRACHO3(CMO)
      USE CHOVEC_IO
      IMPLICIT NONE
* ----------------------------------------------------------------
#include "rasdim.fh"
#include "warnings.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "chocaspt2.fh"
#include "choptr.fh"
#include "choglob.fh"
#include "WrkSpc.fh"
#include "output.fh"
**********************************************************************
*  Author : P. A. Malmqvist
**********************************************************************
      REAL*8 CMO(NBSQT)

      INTEGER NCES(8),ip_HTVec(8)
      Integer, External :: Cho_IRange
      INTEGER ISTART(8),NUSE(8)
      INTEGER IC,ICASE,IRC,ILOC
      INTEGER JSTART
      INTEGER JRED,JRED1,JRED2,JREDC,JNUM,JV1,JV2
      INTEGER IASTA,IAEND,IISTA,IIEND
      INTEGER NA,NASZ,NI,NISZ,NBUFFY,NPQ
      INTEGER IB,IBATCH,IBATCH_TOT,IBSTA,IBEND,NBATCH
      INTEGER IP_LFT,IP_LHT
      INTEGER ISYM,JSYM,ISYMA,ISYMB,ISYP,ISYQ
      INTEGER N,N1,N2
      INTEGER ip_buffy,ip_chspc,ip_ftspc,ip_htspc,ipnt
      INTEGER NUMV,NVECS_RED,NHTOFF,MUSED

      REAL*8, EXTERNAL :: DDOT_

**********************************************************************
* ======================================================================
* This section deals with density matrices and CMO''s
* Offsets into CMO arrays:
      IC=0
      DO ISYM=1,NSYM
       NCES(ISYM)=IC
       IC=IC+NBAS(ISYM)**2
      END DO

* ======================================================================
      CALL GETMEM('CHSPC','ALLO','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTSPC','ALLO','REAL',IP_HTSPC,NHTSPC)
      CALL GETMEM('FTSPC','ALLO','REAL',IP_FTSPC,NFTSPC)
* ======================================================================

      !IBATCH_TOT=0
* Loop over JSYM
      DO JSYM=1,NSYM
      IBATCH_TOT=NBTCHES(JSYM)

      IF(NUMCHO_PT2(JSYM).EQ.0) GOTO 1000

      ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(jSym-1))
      JRED1=iWork(ipnt)
      JRED2=iWork(ipnt-1+NumCho_PT2(jSym))

* Loop over JRED
      DO JRED=JRED1,JRED2

      CALL Cho_X_nVecRS(JRED,JSYM,JSTART,NVECS_RED)
      IF(NVECS_RED.EQ.0) GOTO 999

      ILOC=3
      CALL CHO_X_SETRED(IRC,ILOC,JRED)
* For a reduced set, the structure is known, including
* the mapping between reduced index and basis set pairs.
* The reduced set is divided into suitable batches.
* First vector is JSTART. Nr of vectors in r.s. is NVECS_RED.

* Determine batch length for this reduced set.
* Make sure to use the same formula as in the creation of disk
* address tables, etc, above:
      NBATCH=1+(NVECS_RED-1)/MXNVC

* Loop over IBATCH
      JV1=JSTART
      DO IBATCH=1,NBATCH
      IBATCH_TOT=IBATCH_TOT+1

      JNUM=NVLOC_CHOBATCH(IBATCH_TOT)
      JV2=JV1+JNUM-1

      JREDC=JRED
* Read a batch of reduced vectors
      CALL CHO_VECRD(WORK(IP_CHSPC),NCHSPC,JV1,JV2,JSYM,
     &                        NUMV,JREDC,MUSED)
      IF(NUMV.ne.JNUM) THEN
        write(6,*)' Rats! CHO_VECRD was called, assuming it to'
        write(6,*)' read JNUM vectors. Instead it returned NUMV'
        write(6,*)' vectors: JNUM, NUMV=',JNUM,NUMV
        write(6,*)' Back to the drawing board?'
        CALL QUIT(_RC_INTERNAL_ERROR_)
      END IF
      IF(JREDC.NE.JRED) THEN
        write(6,*)' Rats! It was assumed that the Cholesky vectors'
        write(6,*)' in HALFTRNSF all belonged to a given reduced'
        write(6,*)' set, but they don''t!'
        write(6,*)' JRED, JREDC:',JRED,JREDC
        write(6,*)' Back to the drawing board?'
        write(6,*)' Let the program continue and see what happens.'
      END IF

* Frozen half-transformation:
      NHTOFF=0
      DO ISYMA=1,NSYM
       ISYMB=MUL(ISYMA,JSYM)
       IP_HTVEC(ISYMA)=IP_HTSPC+NHTOFF
       ISTART(ISYMA)=1
       NUSE(ISYMA)=NFRO(ISYMA)
       NHTOFF=NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      END DO
      CALL HALFTRNSF(IRC,WORK(IP_CHSPC),NCHSPC,1,JV1,JNUM,JNUM,
     &     JSYM,JREDC,CMO,ISTART,NUSE,IP_HTVEC)

* Inactive half-transformation:
* Vectors of type HALF(K,J,B) = Sum(CHO(AB,J)*CMO(A,K) where
* A,B are basis functions of symmetry ISYMA, ISYMB,
* K is inactive of symmetry ISYMA, J is vector number in 1..NUMV
* numbered within the present batch.
* Symmetry block ISYMA,ISYMB is found at WORK(IP_HTVEC(ISYMA)
      NHTOFF=0
      DO ISYMA=1,NSYM
       ISYMB=MUL(ISYMA,JSYM)
       IP_HTVEC(ISYMA)=IP_HTSPC+NHTOFF
       ISTART(ISYMA)=NFRO(ISYMA)+1
       NUSE(ISYMA)=NISH(ISYMA)
       NHTOFF=NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      END DO
      CALL HALFTRNSF(IRC,WORK(IP_CHSPC),NCHSPC,1,JV1,JNUM,JNUM,
     &     JSYM,JREDC,CMO,ISTART,NUSE,IP_HTVEC)

* Loop over ISYQ
      DO ISYQ=1,NSYM
       ISYP=MUL(ISYQ,JSYM)

       N=NBAS(ISYP)
* ---------------------------------------------------
       N1=NASH(ISYP)
       N2=NISH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
       IP_LFT=IP_FTSPC
*   Compute fully transformed TK
       IF(N1*N2.GT.0) THEN
        CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,WORK(IP_LHT),WORK(IP_LFT))
        CALL CHOVEC_SAVE(WORK(IP_LFT),1,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
       N1=NSSH(ISYP)
       N2=NISH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP)+NASH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
       IP_LFT=IP_FTSPC
*   Compute fully transformed AK
       IF(N1*N2.GT.0) THEN

C     CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,WORK(IP_LHT),WORK(IP_LFT))

C =SVC= modified for using boxed ordering of pairs, note that the boxed
C routine is less efficient than the original one (loop over J values)
        NA=N1
        NI=N2
C Allocate memory for small buffer used in FULLTRNSF_BOXED
        NBUFFY=NA*NI
        CALL GETMEM('BUFFY','ALLO','REAL',IP_BUFFY,NBUFFY)
C Loop over boxes
        DO IASTA=1,NA,nSecBX
         IAEND=MIN(IASTA-1+nSecBX,NA)
         NASZ=1+IAEND-IASTA
         DO IISTA=1,NI,nInaBX
          IIEND=MIN(IISTA-1+nInaBX,NI)
          NISZ=1+IIEND-IISTA
C =SVC= note that WITHIN this box, the index of the outer box A (P in
C the FULLTRNSF subroutine, i.e. the secondary orbital index) is the
C fast-running index, as LFT([A],[I],J) = CMO(P,[A])^T * LHT([I],J,P)^T
C with P=1,NB.  So if used in e.g. ADDRHS as BRA(c,l,J), making an inner
C loop over secondary orbital index c is more efficient.
          CALL FULLTRNSF_BOXED (IASTA,IISTA,NASZ,NISZ,NA,NI,
     &                          N,CMO(IC+N*(IASTA-1)),JNUM,
     &                          WORK(IP_LHT),WORK(IP_LFT),
     &                          WORK(IP_BUFFY))
         ENDDO
        ENDDO
        CALL GETMEM('BUFFY','FREE','REAL',IP_BUFFY,NBUFFY)
        CALL CHOVEC_SAVE(WORK(IP_LFT),4,ISYQ,JSYM,IBATCH_TOT)
       END IF
* End loop ISYQ
      END DO
* ---------------------------------------------------

* Active scaled natural orbitals half-transformation:
* Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W) where
* A,B are basis functions of symmetry ISYMA, ISYMB,
* W is active of symmetry ISYMA, J is vector number in 1..NUMV
* numbered within the present batch.
* Symmetry block ISYMA,ISYMB is found at WORK(IP_HTVEC(ISYMA)
      NHTOFF=0
      DO ISYMA=1,NSYM
       ISYMB=MUL(ISYMA,JSYM)
       IP_HTVEC(ISYMA)=IP_HTSPC+NHTOFF
       ISTART(ISYMA)=NFRO(ISYMA)+NISH(ISYMA)+1
       NUSE(ISYMA)=NASH(ISYMA)
       NHTOFF=NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      END DO
* ---------------------------------------------------
* Active half-transformation:
      CALL HALFTRNSF(IRC,WORK(IP_CHSPC),NCHSPC,1,JV1,JNUM,JNUM,
     &    JSYM,JREDC,CMO,ISTART,NUSE,IP_HTVEC)


      DO ISYQ=1,NSYM
       ISYP=MUL(ISYQ,JSYM)

       N=NBAS(ISYP)
* ---------------------------------------------------
* Loop over ISYQ
       N1=NASH(ISYP)
       N2=NASH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
       IP_LFT=IP_FTSPC
* Compute fully transformed TV
       IF(N1*N2.GT.0) THEN
        CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,WORK(IP_LHT),WORK(IP_LFT))
        CALL CHOVEC_SAVE(WORK(IP_LFT),2,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
       N1=NSSH(ISYP)
       N2=NASH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP)+NASH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
       IP_LFT=IP_FTSPC
*   Compute fully transformed AV
       IF(N1*N2.GT.0) THEN
        CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,WORK(IP_LHT),WORK(IP_LFT))
        CALL CHOVEC_SAVE(WORK(IP_LFT),3,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
* End loop ISYQ
      END DO

* 800  CONTINUE
* End loop IBATCH
       JV1=JV1+JNUM
      END DO

* End loop JRED
  999 CONTINUE
      END DO

* End loop JSYM
 1000 CONTINUE
      END DO

      ! if using the RHS on-demand, we need all cholesky vectors on each
      ! process, collect them here
      IF (RHSDIRECT) THEN
        IP_LFT=IP_FTSPC
        DO JSYM=1,NSYM
          IBSTA=NBTCHES(JSYM)+1
          IBEND=NBTCHES(JSYM)+NBTCH(JSYM)
          DO IB=IBSTA,IBEND
            DO ISYQ=1,NSYM
              DO ICASE=1,4
                NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
                IF (NPQ.EQ.0) CYCLE
                CALL CHOVEC_LOAD(WORK(IP_LFT),ICASE,ISYQ,JSYM,IB)
                CALL CHOVEC_COLL(WORK(IP_LFT),ICASE,ISYQ,JSYM,IB)
              END DO
            END DO
          END DO
        END DO
      END IF

      CALL GETMEM('CHSPC','FREE','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTSPC','FREE','REAL',IP_HTSPC,NHTSPC)
      CALL GETMEM('FTSPC','FREE','REAL',IP_FTSPC,NFTSPC)

      RETURN
      END
