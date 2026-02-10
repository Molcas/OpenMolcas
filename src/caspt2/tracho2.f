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
      SUBROUTINE TRACHO2(CMO,NCMO,DREF,NDREF,FFAO,FIAO,FAAO,IF_TRNSF)
      USE CHOVEC_IO, only: NVLOC_CHOBATCH,NPQ_CHOTYPE,chovec_save,
     &                     chovec_load,chovec_coll
      use Cholesky, only: InfVec, nDimRS
      use ChoCASPT2, only: MXCHARR,MXNVC,NCHSPC,NFTSPC,NHTSPC,
     &                     NUMCHO_PT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: nBTri, ECore, nBasT, nBSqT, nInaBx,
     &                         nSecBx, nSym, PotNuc, nBas, nFro,
     &                         nIsh, nAsh, RHSDirect, nBtches, Mul,
     &                         nSsh, nBtch
      IMPLICIT NONE
* ----------------------------------------------------------------
#include "warnings.h"
      INTEGER NCMO, NDREF
      REAL*8 CMO(NCMO),DREF(NDREF),
     &       FFAO(NBTRI),FIAO(NBTRI),FAAO(NBTRI)
      LOGICAL IF_TRNSF

      INTEGER NCES(8),ip_HTVec(8)
      INTEGER ISTART(8),NUSE(8)

      REAL*8 E,ECORE1,ECORE2

      REAL*8 FACTC,FACTXA,FACTXI

      INTEGER I,J,IC,IA,ICASE,IRC,ILOC
      INTEGER JSTART
      INTEGER JRED,JRED1,JRED2,JREDC,JNUM,JV1,JV2
      INTEGER IASTA,IAEND,IISTA,IIEND
      INTEGER NA,NASZ,NI,NISZ,NBUFFY,NF,NK,NW,NPQ,NRS
      INTEGER IB,IBATCH,IBATCH_TOT,IBSTA,IBEND,NB,NBATCH
      INTEGER IDFIJ,IDIIJ,IDAIJ
      INTEGER IP_LHT
      INTEGER LC,LO,LSC,LSO
      INTEGER ISFA,ISFF,ISFI
      INTEGER ISYM,JSYM,ISYMA,ISYMB,ISYMK,ISYMW,ISYP,ISYQ
      INTEGER N,N1,N2
      INTEGER ip_htspc
      INTEGER NUMV,NVECS_RED,NHTOFF,MUSED

      REAL*8 SCL

      REAL*8, EXTERNAL :: DDOT_
      REAL*8, ALLOCATABLE:: OCC(:), CNAT(:), DF(:), DI(:), DA(:)
      REAL*8, ALLOCATABLE:: VEC(:), DF_RED(:), DI_RED(:), DA_RED(:)
      REAL*8, ALLOCATABLE:: FA_RED(:), FF_RED(:), FI_RED(:)
      REAL*8, ALLOCATABLE:: BUFFY(:), CHSPC(:), FTSPC(:), HTSPC(:)

************************************************************************
* ======================================================================
* This section deals with density matrices and CMO''s
* Offsets into CMO arrays:
      IC=0
      DO ISYM=1,NSYM
       NCES(ISYM)=IC
       IC=IC+NBAS(ISYM)**2
      END DO

* Compute natural orbitals for the reference wave function:
      Call mma_allocate(OCC,NBasT,Label='OCC')
      Call mma_allocate(CNAT,NBSQT,Label='CNAT')
*      write(6,*)' Active/Active density matrix, triangular'
*      IF( NASHT.GT.0 ) THEN
*        call TRIPRT(' ',' ',DREF,NASHT)
*      ENDIF
      CALL REF_NATO(DREF,CMO,OCC,CNAT)

c Initialize Fock matrices in AO basis to zero:
      FFAO(:)=0.0D0
      FIAO(:)=0.0D0
      FAAO(:)=0.0D0
* Construct density matrix for frozen orbitals
      Call mma_allocate(DF,NBTRI,Label='DF')
      DO ISYM=1,NSYM
       ISTART(ISYM)=1
       NUSE(ISYM)=NFRO(ISYM)
      END DO
      CALL GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,OCC,DF)
* Construct density matrix for inactive orbitals
      Call mma_allocate(DI,NBTRI,Label='DI')
      DO ISYM=1,NSYM
       ISTART(ISYM)=NFRO(ISYM)+1
       NUSE(ISYM)=NISH(ISYM)
      END DO
      CALL GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,OCC,DI)
* Same, for active density:
      Call mma_allocate(DA ,NBTRI,Label='DA')
      DO ISYM=1,NSYM
       ISTART(ISYM)=NFRO(ISYM)+NISH(ISYM)+1
       NUSE(ISYM)=NASH(ISYM)
      END DO
      CALL GDMAT(NSYM,NBAS,ISTART,NUSE,CNAT,OCC,DA)
* The Cholesky routines want density matrices in a particular storage, and
* also the off-diagonal elements should be doubled. Double them:
      IDFIJ=1
      IDIIJ=1
      IDAIJ=1
      DO ISYM=1,NSYM
       DO I=1,NBAS(ISYM)
        DO J=1,I-1
         DF(IDFIJ)=2.0D0*DF(IDFIJ)
         DI(IDIIJ)=2.0D0*DI(IDIIJ)
         DA(IDAIJ)=2.0D0*DA(IDAIJ)
         IDFIJ=IDFIJ+1
         IDIIJ=IDIIJ+1
         IDAIJ=IDAIJ+1
        END DO
        IDFIJ=IDFIJ+1
        IDIIJ=IDIIJ+1
        IDAIJ=IDAIJ+1
       END DO
      END DO

* Scale natural orbitals by multiplying with square root of half the
* occupation number -- This allows computing the exchange contribution
* to Fock matrices using the same formula as for closed shells.
      LSO=1
      LSC=1
      DO ISYM=1,NSYM
       NF=NFRO(ISYM)
       NI=NISH(ISYM)
       NA=NASH(ISYM)
       NB=NBAS(ISYM)
       LO=LSO+NF+NI
       LC=LSC+NB*(NF+NI)
       DO IA=1,NA
        SCL=SQRT(0.5D0*OCC(LO))
        CALL DSCAL_(NB,SCL,CNAT(LC),1)
        LO=LO+1
        LC=LC+NB
       END DO
       LSO=LSO+NB
       LSC=LSC+NB**2
      END DO
* ======================================================================
      CALL mma_allocate(CHSPC,NCHSPC,LABEL='CHSPC')
      CALL mma_allocate(HTSPC,NHTSPC,LABEL='HTSPC')
      IP_HTSPC=1
      IF (IF_TRNSF) THEN
       CALL mma_allocate(FTSPC,NFTSPC,LABEL='FTSPC')
      END IF
* ======================================================================

      !IBATCH_TOT=0
* Loop over JSYM
      DO JSYM=1,NSYM
      IBATCH_TOT=NBTCHES(JSYM)

*     write(6,*)' Tracho2 JSYM=',JSYM
*     write(6,*)'    NUMCHO_PT2(JSYM)=',NUMCHO_PT2(JSYM)
      IF(NUMCHO_PT2(JSYM).EQ.0) GOTO 1000

      JRED1=InfVec(1,2,jSym)
      JRED2=InfVec(NumCho_PT2(jSym),2,jSym)
*     write(6,*)'tracho2:  JRED1,JRED2:',JRED1,JRED2

      IF(JSYM.EQ.1) THEN
* Allocate space for temporary vector 'Vec' used for Coulomb contrib to
* Fock matrices:
       CALL mma_allocate(VEC,MXCHARR,LABEL='VEC')
* Local density matrices, which will be needed if JSYM=1. At the same time,
* allocate Fock matrices with the same structure and initialize to zero.
       CALL mma_allocate(DF_RED,MXCHARR,LABEL='DF_RED')
       CALL mma_allocate(DI_RED,MXCHARR,LABEL='DI_RED')
       CALL mma_allocate(DA_RED,MXCHARR,LABEL='DA_RED')
       CALL mma_allocate(FF_RED,MXCHARR,LABEL='FF_RED')
       CALL mma_allocate(FI_RED,MXCHARR,LABEL='FI_RED')
       CALL mma_allocate(FA_RED,MXCHARR,LABEL='FA_RED')
      END IF


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
*      JEND=JSTART+NVECS_RED-1
*      write(6,*)'  JRED:  JSTART,JEND:',JRED,JSTART,JEND

      IF(JSYM.EQ.1) THEN
      NRS=NDIMRS(JSYM,JRED)
      CALL DCOPY_(NRS,[0.0D0],0,DF_RED,1)
      CALL full2red(DF,DF_Red)
      CALL DCOPY_(NRS,[0.0D0],0,DI_RED,1)
      CALL full2red(DI,DI_Red)
      CALL DCOPY_(NRS,[0.0D0],0,DA_RED,1)
      CALL full2red(DA,DA_Red)
      CALL DCOPY_(NRS,[0.0D0],0,FF_RED,1)
      CALL DCOPY_(NRS,[0.0D0],0,FI_RED,1)
      CALL DCOPY_(NRS,[0.0D0],0,FA_RED,1)
      END IF

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
      CALL CHO_VECRD(CHSPC,NCHSPC,JV1,JV2,JSYM,NUMV,JREDC,MUSED)
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

      IF (JSYM.EQ.1) THEN
* Coulomb contribution to Fock arrays.
* V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
* Starting at CHSPC is now an array of vectors, conceptually
* L(rs,J), where temporarily we can regard J as ranging 1..JNUM, and
* the layout of pair indices rs is unknown ('reduced storage', a secret
* inside cholesky.) Compute array V(J) at temporary space VEC:
       CALL DGEMV_('T',NRS,JNUM,1.0D0,CHSPC,NRS,
     &            DF_RED,1,0.0D0,VEC,1)
* F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
             FactC=1.0D0
       CALL DGEMV_('N',NRS,JNUM,FactC,CHSPC,NRS,
     &             VEC,1,1.0D0,FF_RED,1)
* The same thing, now for the inactive and active density matrices:
       CALL DGEMV_('T',NRS,JNUM,1.0D0,CHSPC,NRS,
     &            DI_RED,1,0.0D0,VEC,1)
       CALL DGEMV_('N',NRS,JNUM,FactC,CHSPC,NRS,
     &             VEC,1,1.0D0,FI_RED,1)
       CALL DGEMV_('T',NRS,JNUM,1.0D0,CHSPC,NRS,
     &             DA_RED,1,0.0D0,VEC,1)
       CALL DGEMV_('N',NRS,JNUM,FactC,CHSPC,NRS,
     &             VEC,1,1.0D0,FA_RED,1)
*      write(6,*)' Finished Coulomb contributions to Fock matrix.'
*      write(6,*)' Frozen Fock mat at FF_RED'
*      write(6,'(1x,8f10.4)')(FF_RED(i),i=1,nRS)
*      write(6,*)' Inactive Fock mat at FI_RED'
*      write(6,'(1x,8f10.4)')(FI_RED(i),i=1,nRS)
*      write(6,*)' Active Fock matrix at FA_RED.'
*      write(6,'(1x,8f10.4)')(FA_RED(i),i=1,nRS)
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
      CALL HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,
     &     JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)
* Frozen contributions to exchange:
          FactXI=-1.0D0
          ISFF=1
          DO ISYMB=1,NSYM
           iSymk = MUL(jSym,iSymb)
C ---------------------------------------------------------------------
c *** Compute the LT part of the FROZEN exchange matrix ********
C     FF(ab) = FF(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
C ---------------------------------------------------------------------
           NK = NFRO(ISYMK)
           NB = NBAS(ISYMB)
           If (NB*NK.ne.0) Then
            CALL DGEMM_TRI('T','N',NB,NB,NK*JNUM,
     &                  FactXI,HTSPC(ip_HTVec(iSymk)),NK*JNUM,
     &                  HTSPC(ip_HTVec(iSymk)),NK*JNUM,
     &                  1.0D0,FFAO(ISFF),NB)
           EndIf
           ISFF = ISFF+(NB*(NB+1))/2
          END DO

* Inactive half-transformation:
* Vectors of type HALF(K,J,B) = Sum(CHO(AB,J)*CMO(A,K) where
* A,B are basis functions of symmetry ISYMA, ISYMB,
* K is inactive of symmetry ISYMA, J is vector number in 1..NUMV
* numbered within the present batch.
* Symmetry block ISYMA,ISYMB is found at HTSPC(IP_HTVEC(ISYMA)
      NHTOFF=0
      DO ISYMA=1,NSYM
       ISYMB=MUL(ISYMA,JSYM)
       IP_HTVEC(ISYMA)=IP_HTSPC+NHTOFF
       ISTART(ISYMA)=NFRO(ISYMA)+1
       NUSE(ISYMA)=NISH(ISYMA)
       NHTOFF=NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      END DO
      CALL HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,
     &     JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)
* Inactive contributions to exchange:
          FactXI=-1.0D0
          ISFI=1
          DO ISYMB=1,NSYM
           iSymk = MUL(jSym,iSymb)
C ---------------------------------------------------------------------
c *** Compute the LT part of the INACTIVE exchange matrix ********
C     FI(ab) = FI(ab) + FactXI * sum_Jk  LkJ,a * LkJ,b
C ---------------------------------------------------------------------
           NK = NISH(iSymk)
           NB = NBAS(ISYMB)
           If (NB*NK.ne.0) Then
           CALL DGEMM_TRI('T','N',NB,NB,NK*JNUM,
     &                   FactXI,HTSPC(ip_HTVec(iSymk)),NK*JNUM,
     &                   HTSPC(ip_HTVec(iSymk)),NK*JNUM,
     &                   1.0D0,FIAO(ISFI),NB)
           EndIf
           ISFI = ISFI+(NB*(NB+1))/2
          END DO
*      write(6,*)' Inactive Fock mat in FIAO'
*      write(6,'(1x,8f10.4)')(FIAO(i),i=1,nbtri)

      IF (IF_TRNSF) THEN
* Loop over ISYQ
      DO ISYQ=1,NSYM
       ISYP=MUL(ISYQ,JSYM)

       N=NBAS(ISYP)
* ---------------------------------------------------
       N1=NASH(ISYP)
       N2=NISH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
*   Compute fully transformed TK
       IF(N1*N2.GT.0) THEN
        CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)
        CALL CHOVEC_SAVE(FTSPC,1,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
       N1=NSSH(ISYP)
       N2=NISH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP)+NASH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
*   Compute fully transformed AK
       IF(N1*N2.GT.0) THEN

C     CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)

C =SVC= modified for using boxed ordering of pairs, note that the boxed
C routine is less efficient than the original one (loop over J values)
        NA=N1
        NI=N2
C Allocate memory for small buffer used in FULLTRNSF_BOXED
        NBUFFY=NA*NI
        CALL mma_allocate(BUFFY,NBUFFY,Label='BUFFY')
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
     &                          HTSPC(IP_LHT),FTSPC,
     &                          BUFFY)
         ENDDO
        ENDDO
        CALL mma_deallocate(BUFFY)
        CALL CHOVEC_SAVE(FTSPC,4,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
* End loop ISYQ
      END DO
      END IF
* Active scaled natural orbitals half-transformation:
* Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W) where
* A,B are basis functions of symmetry ISYMA, ISYMB,
* W is active of symmetry ISYMA, J is vector number in 1..NUMV
* numbered within the present batch.
* Symmetry block ISYMA,ISYMB is found at HTSPC(IP_HTVEC(ISYMA)
      NHTOFF=0
      DO ISYMA=1,NSYM
       ISYMB=MUL(ISYMA,JSYM)
       IP_HTVEC(ISYMA)=IP_HTSPC+NHTOFF
       ISTART(ISYMA)=NFRO(ISYMA)+NISH(ISYMA)+1
       NUSE(ISYMA)=NASH(ISYMA)
       NHTOFF=NHTOFF+NUSE(ISYMA)*NBAS(ISYMB)*JNUM
      END DO
      CALL HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,
     &    JSYM,JREDC,CNAT,NBSQT,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)
* Active (scaled) contributions to exchange:
      FactXA=-1.0D0
      ISFA=1
      DO ISYMB=1,NSYM
       iSymw = MUL(jSym,iSymb)
C ---------------------------------------------------------------------
c *** Compute the LT part of the ACTIVE exchange matrix ********
C     FA(ab) = FA(ab) + FactXA * sum_Jw  LwJ,a * LwJ,b
C ---------------------------------------------------------------------
       NW = NASH(iSymw)
       NB = NBAS(ISYMB)
       If (NB*NW.ne.0) Then
       CALL DGEMM_TRI('T','N',NB,NB,NW*JNUM,
     &               FactXA,HTSPC(ip_HTVec(iSymw)),NW*JNUM,
     &               HTSPC(ip_HTVec(iSymw)),NW*JNUM,
     &               1.0D0,FAAO(ISFA),NB)
       EndIf
       ISFA = ISFA+(NB*(NB+1))/2
      END DO

*      write(6,*)' Active Fock matrix in FAAO.'
*      write(6,'(1x,8f10.4)')(FAAO(i),i=1,nbtri)

* ---------------------------------------------------
* Active half-transformation:
      CALL HALFTRNSF(IRC,CHSPC,NCHSPC,1,JV1,JNUM,JNUM,
     &    JSYM,JREDC,CMO,NCMO,ISTART,NUSE,IP_HTVEC,HTSPC,NHTSPC)


      IF (IF_TRNSF) THEN
      DO ISYQ=1,NSYM
       ISYP=MUL(ISYQ,JSYM)

       N=NBAS(ISYP)
* ---------------------------------------------------
* Loop over ISYQ
       N1=NASH(ISYP)
       N2=NASH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
* Compute fully transformed TV
       IF(N1*N2.GT.0) THEN
        CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)
        CALL CHOVEC_SAVE(FTSPC,2,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
       N1=NSSH(ISYP)
       N2=NASH(ISYQ)
       IC=1+NCES(ISYP) +(NFRO(ISYP)+NISH(ISYP)+NASH(ISYP))*N
       IP_LHT=IP_HTVEC(ISYQ)
*   Compute fully transformed AV
       IF(N1*N2.GT.0) THEN
        CALL FULLTRNSF(N1,N2,N,CMO(IC),JNUM,HTSPC(IP_LHT),FTSPC)
        CALL CHOVEC_SAVE(FTSPC,3,ISYQ,JSYM,IBATCH_TOT)
       END IF
* ---------------------------------------------------
* End loop ISYQ
      END DO
      END IF

* 800  CONTINUE
* End loop IBATCH
       JV1=JV1+JNUM
      END DO

      IF (jSym.eq.1) THEN
* Add Coulomb contributions in local Fock matrices (in 'reduced storage')
* into the global ones:
        CALL red2full(FFAO,FF_RED)
        CALL red2full(FIAO,FI_RED)
        CALL red2full(FAAO,FA_RED)
      END IF
* End loop JRED
  999 CONTINUE
      END DO

      IF (jSym.eq.1) THEN
* Deallocate local density and fock matrices
        CALL mma_deallocate(VEC)
        CALL mma_deallocate(DF_RED)
        CALL mma_deallocate(DI_RED)
        CALL mma_deallocate(DA_RED)
        CALL mma_deallocate(FF_RED)
        CALL mma_deallocate(FI_RED)
        CALL mma_deallocate(FA_RED)
      END IF
* End loop JSYM
 1000 CONTINUE
      END DO

      ! if using the RHS on-demand, we need all cholesky vectors on each
      ! process, collect them here
      IF (IF_TRNSF.AND.RHSDIRECT) THEN
        DO JSYM=1,NSYM
          IBSTA=NBTCHES(JSYM)+1
          IBEND=NBTCHES(JSYM)+NBTCH(JSYM)
          DO IB=IBSTA,IBEND
            DO ISYQ=1,NSYM
              DO ICASE=1,4
                NPQ=NPQ_CHOTYPE(ICASE,ISYQ,JSYM)
                IF (NPQ.EQ.0) CYCLE
                CALL CHOVEC_LOAD(FTSPC,ICASE,ISYQ,JSYM,IB)
                CALL CHOVEC_COLL(FTSPC,ICASE,ISYQ,JSYM,IB)
              END DO
            END DO
          END DO
        END DO
      END IF

* Synchronize and add the contributions from all nodes into each node:
      CALL GADGOP(FFAO,NBTRI,'+')
      CALL GADGOP(FIAO,NBTRI,'+')
      CALL GADGOP(FAAO,NBTRI,'+')

* Two-electron contribution to the effective core energy
      ECORE2=0.5D0*DDOT_(NBTRI,DF,1,FFAO,1)
c Add OneHam to finalize frozen Fock matrix in AO basis.
c (It is in fact an effective one-electron Hamiltonian).
      CALL ADD1HAM(FFAO)
* The contraction of frozen Fock matrix with frozen density:
      E=DDOT_(NBTRI,DF,1,FFAO,1)
* Correct for double-counting two-electron part:
      E=E-ECORE2
* One-electron part:
      ECORE1=E-ECORE2
* Nuclear repulsion energy:
      ECORE=POTNUC+ECORE1+ECORE2

#ifdef _DEBUGPRINT_
       WRITE(6,'(6X,A,ES20.10)') 'NUCLEAR REPULSION ENERGY:',POTNUC
       WRITE(6,'(6X,A,ES20.10)') 'ONE-ELECTRON CORE ENERGY:',ECORE1
       WRITE(6,'(6X,A,ES20.10)') 'TWO-ELECTRON CORE ENERGY:',ECORE2
       WRITE(6,'(6X,A,ES20.10)') '       TOTAL CORE ENERGY:',ECORE
#endif

      Call mma_deallocate(OCC)
      Call mma_deallocate(CNAT)
      Call mma_deallocate(DF)
      Call mma_deallocate(DI)
      Call mma_deallocate(DA)

      CALL mma_deallocate(CHSPC)
      CALL mma_deallocate(HTSPC)
      IF (IF_TRNSF) THEN
       CALL mma_deallocate(FTSPC)
      END IF

#ifdef _DEBUGPRINT_
        WRITE(6,'(6X,A)')'TEST PRINT FROM TRACHO2.'
        WRITE(6,'(6X,A)')
        write(6,*)' NSYM:',NSYM
        write(6,*)' NBAS:',(NBAS(ISYM),ISYM=1,8)
        WRITE(6,'(6X,A)')
        WRITE(6,'(6X,A)')'***** FROZEN FOCK MATRIX ***** '
        ISFF=1
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          IF( NB.GT.0 ) THEN
            WRITE(6,'(6X,A)')
            WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
            call TRIPRT(' ',' ',FFAO(ISFF),NB)
            ISFI=ISFF+(NB*(NB+1))/2
          ENDIF
        END DO
        WRITE(6,'(6X,A)')
        WRITE(6,'(6X,A)')'***** INACTIVE FOCK MATRIX ***** '
        ISFI=1
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          IF( NB.GT.0 ) THEN
            WRITE(6,'(6X,A)')
            WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
            call TRIPRT(' ',' ',FIAO(ISFI),NB)
            ISFI=ISFI+(NB*(NB+1))/2
          ENDIF
        END DO
        WRITE(6,'(6X,A)')
        WRITE(6,'(6X,A)')'***** ACTIVE FOCK MATRIX ***** '
        ISFA=1
        DO ISYM=1,NSYM
          NB=NBAS(ISYM)
          IF( NB.GT.0 ) THEN
            WRITE(6,'(6X,A)')
            WRITE(6,'(6X,A,I2)')'SYMMETRY SPECIES:',ISYM
            call TRIPRT(' ',' ',FAAO(ISFA),NB)
            ISFA=ISFA+(NB*(NB+1))/2
          ENDIF
        END DO
#endif

      END SUBROUTINE TRACHO2
