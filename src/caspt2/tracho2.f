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
      SUBROUTINE TRACHO2(CMO,DREF,FFAO,FIAO,FAAO,IF_TRNSF)
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
      REAL*8 CMO(NBSQT),DREF(NDREF),
     &       FFAO(NBTRI),FIAO(NBTRI),FAAO(NBTRI)
      LOGICAL IF_TRNSF

      INTEGER NCES(8),ip_HTVec(8)
      Integer, External :: Cho_IRange
      INTEGER ISTART(8),NUSE(8)

      REAL*8 E,ECORE1,ECORE2

      REAL*8 FACTC,FACTXA,FACTXI

      INTEGER I,J,IC,IA,ICASE,IRC,ILOC
      INTEGER JSTART,JEND
      INTEGER JRED,JRED1,JRED2,JREDC,JNUM,JV1,JV2
      INTEGER IASTA,IAEND,IISTA,IIEND
      INTEGER NA,NASZ,NI,NISZ,NBUFFY,NF,NK,NW,NPQ,NRS
      INTEGER IB,IBATCH,IBATCH_TOT,IBSTA,IBEND,NB,NBATCH
      INTEGER IDFIJ,IDIIJ,IDAIJ
      INTEGER IP_LFT,IP_LHT
      INTEGER IPDA,IPDA_RED,IPDF,IPDF_RED,IPDI,IPDI_RED
      INTEGER LC,LCNAT,LO,LOCC,LSC,LSO
      INTEGER LFA_RED,LFF_RED,LFI_RED
      INTEGER ISFA,ISFF,ISFI
      INTEGER ISYM,JSYM,ISYMA,ISYMB,ISYMK,ISYMW,ISYP,ISYQ
      INTEGER N,N1,N2
      INTEGER ip_buffy,ip_chspc,ip_ftspc,ip_htspc,ip_v,ipnt
      INTEGER NUMV,NV,NVECS_RED,NVTOT,NHTOFF,MUSED

      REAL*8 SCL

      REAL*8, EXTERNAL :: DDOT_

**********************************************************************
      Call QEnter('TRACHO2')
* ======================================================================
* This section deals with density matrices and CMO''s
* Offsets into CMO arrays:
      IC=0
      DO ISYM=1,NSYM
       NCES(ISYM)=IC
       IC=IC+NBAS(ISYM)**2
      END DO

* Compute natural orbitals for the reference wave function:
      Call Getmem('OCC','ALLO','REAL',LOCC,NBasT)
      Call Getmem('CNAT','ALLO','REAL',LCNAT,NBSQT)
*      write(6,*)' Active/Active density matrix, triangular'
*      IF( NASHT.GT.0 ) THEN
*        call TRIPRT(' ',' ',DREF,NASHT)
*      ENDIF
      CALL REF_NATO(DREF,CMO,WORK(LOCC),WORK(LCNAT))

c Initialize Fock matrices in AO basis to zero:
      CALL DCOPY_(NBTRI,[0.0D0],0,FFAO,1)
      CALL DCOPY_(NBTRI,[0.0D0],0,FIAO,1)
      CALL DCOPY_(NBTRI,[0.0D0],0,FAAO,1)
* Construct density matrix for frozen orbitals
      Call Getmem('DF','ALLO','REAL',ipDF,NBTRI)
      DO ISYM=1,NSYM
       ISTART(ISYM)=1
       NUSE(ISYM)=NFRO(ISYM)
      END DO
      CALL GDMAT(NSYM,NBAS,ISTART,NUSE,
     &                          WORK(LCNAT),WORK(LOCC),WORK(IPDF))
* Construct density matrix for inactive orbitals
      Call Getmem('DI','ALLO','REAL',ipDI,NBTRI)
      DO ISYM=1,NSYM
       ISTART(ISYM)=NFRO(ISYM)+1
       NUSE(ISYM)=NISH(ISYM)
      END DO
      CALL GDMAT(NSYM,NBAS,ISTART,NUSE,
     &                          WORK(LCNAT),WORK(LOCC),WORK(IPDI))
* Same, for active density:
      Call Getmem('DA ','ALLO','REAL',ipDA ,NBTRI)
      DO ISYM=1,NSYM
       ISTART(ISYM)=NFRO(ISYM)+NISH(ISYM)+1
       NUSE(ISYM)=NASH(ISYM)
      END DO
      CALL GDMAT(NSYM,NBAS,ISTART,NUSE,
     &                          WORK(LCNAT),WORK(LOCC),WORK(IPDA ))
* The Cholesky routines want density matrices in a particular storage, and
* also the off-diagonal elements should be doubled. Double them:
      IDFIJ=IPDF
      IDIIJ=IPDI
      IDAIJ=IPDA
      DO ISYM=1,NSYM
       DO I=1,NBAS(ISYM)
        DO J=1,I-1
         WORK(IDFIJ)=2.0D0*WORK(IDFIJ)
         WORK(IDIIJ)=2.0D0*WORK(IDIIJ)
         WORK(IDAIJ)=2.0D0*WORK(IDAIJ)
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
      LSO=LOCC
      LSC=LCNAT
      DO ISYM=1,NSYM
       NF=NFRO(ISYM)
       NI=NISH(ISYM)
       NA=NASH(ISYM)
       NB=NBAS(ISYM)
       LO=LSO+NF+NI
       LC=LSC+NB*(NF+NI)
       DO IA=1,NA
        SCL=SQRT(0.5D0*WORK(LO))
        CALL DSCAL_(NB,SCL,WORK(LC),1)
        LO=LO+1
        LC=LC+NB
       END DO
       LSO=LSO+NB
       LSC=LSC+NB**2
      END DO
* ======================================================================
      CALL GETMEM('CHSPC','ALLO','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTSPC','ALLO','REAL',IP_HTSPC,NHTSPC)
      IF (IF_TRNSF) THEN
       CALL GETMEM('FTSPC','ALLO','REAL',IP_FTSPC,NFTSPC)
      END IF
* ======================================================================

      !IBATCH_TOT=0
* Loop over JSYM
      DO JSYM=1,NSYM
      IBATCH_TOT=NBTCHES(JSYM)

*      write(6,*)' Tracho2 JSYM=',JSYM
*      write(6,*)'    NUMCHO_PT2(JSYM)=',NUMCHO_PT2(JSYM)
      IF(NUMCHO_PT2(JSYM).EQ.0) GOTO 1000

      ipnt=ip_InfVec+MaxVec_PT2*(1+InfVec_N2_PT2*(jSym-1))
      JRED1=iWork(ipnt)
      JRED2=iWork(ipnt-1+NumCho_PT2(jSym))
*      write(6,*)'  JRED1,JRED2:',JRED1,JRED2

      IF(JSYM.EQ.1) THEN
* Allocate space for temporary vector 'V' used for Coulomb contrib to
* Fock matrices:
       CALL GETMEM('V_VECTOR','ALLO','REAL',IP_V,MXCHARR)
* Local density matrices, which will be needed if JSYM=1. At the same time,
* allocate Fock matrices with the same structure and initialize to zero.
       CALL GETMEM('DF_RED','ALLO','REAL',IPDF_RED,MXCHARR)
       CALL GETMEM('DI_RED','ALLO','REAL',IPDI_RED,MXCHARR)
       CALL GETMEM('DA_RED','ALLO','REAL',IPDA_RED,MXCHARR)
       CALL GETMEM('FF_RED','ALLO','REAL',LFF_RED,MXCHARR)
       CALL GETMEM('FI_RED','ALLO','REAL',LFI_RED,MXCHARR)
       CALL GETMEM('FA_RED','ALLO','REAL',LFA_RED,MXCHARR)
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
      JEND=JSTART+NVECS_RED-1
*      write(6,*)'  JRED:  JSTART,JEND:',JRED,JSTART,JEND

      IF(JSYM.EQ.1) THEN
      NRS=IWORK(IP_NDIMRS-1+NSYM*(JRED-1)+JSYM)
      CALL DCOPY_(NRS,[0.0D0],0,WORK(IPDF_RED),1)
      CALL full2red(Work(ipDF),Work(ipDF_Red))
      CALL DCOPY_(NRS,[0.0D0],0,WORK(IPDI_RED),1)
      CALL full2red(Work(ipDI),Work(ipDI_Red))
      CALL DCOPY_(NRS,[0.0D0],0,WORK(IPDA_RED),1)
      CALL full2red(Work(ipDA),Work(ipDA_Red))
      CALL DCOPY_(NRS,[0.0D0],0,WORK(LFF_RED),1)
      CALL DCOPY_(NRS,[0.0D0],0,WORK(LFI_RED),1)
      CALL DCOPY_(NRS,[0.0D0],0,WORK(LFA_RED ),1)
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

      IF (JSYM.EQ.1) THEN
* Coulomb contribution to Fock arrays.
* V{#J} <- V{#J}  +  sum_rs  L(rs,{#J}) * D(rs)
* Starting at Work(IP_CHSPC) is now an array of vectors, conceptually
* L(rs,J), where temporarily we can regard J as ranging 1..JNUM, and
* the layout of pair indices rs is unknown ('reduced storage', a secret
* inside cholesky.) Compute array V(J) at temporary space ip_V:
       CALL DGEMV_('T',NRS,JNUM,1.0D0,WORK(IP_CHSPC),NRS,
     &            WORK(IPDF_RED),1,0.0D0,WORK(IP_V),1)
* F(rs){#J} <- F(rs){#J} + FactC * sum_J L(rs,{#J})*V{#J}
             FactC=1.0D0
       CALL DGEMV_('N',NRS,JNUM,FactC,WORK(IP_CHSPC),NRS,
     &             WORK(IP_V),1,1.0D0,WORK(LFF_RED),1)
* The same thing, now for the inactive and active density matrices:
       CALL DGEMV_('T',NRS,JNUM,1.0D0,WORK(IP_CHSPC),NRS,
     &            WORK(IPDI_RED),1,0.0D0,WORK(IP_V),1)
       CALL DGEMV_('N',NRS,JNUM,FactC,WORK(IP_CHSPC),NRS,
     &             WORK(IP_V),1,1.0D0,WORK(LFI_RED),1)
       CALL DGEMV_('T',NRS,JNUM,1.0D0,WORK(IP_CHSPC),NRS,
     &             WORK(IPDA_RED),1,0.0D0,WORK(IP_V),1)
       CALL DGEMV_('N',NRS,JNUM,FactC,WORK(IP_CHSPC),NRS,
     &             WORK(IP_V),1,1.0D0,WORK(LFA_RED),1)
*      write(6,*)' Finished Coulomb contributions to Fock matrix.'
*      write(6,*)' Frozen Fock mat at Work(LFF_RED)'
*      write(6,'(1x,8f10.4)')(Work(LFF_RED+i),i=0,nRS-1)
*      write(6,*)' Inactive Fock mat at Work(LFI_RED)'
*      write(6,'(1x,8f10.4)')(Work(LFI_RED+i),i=0,nRS-1)
*      write(6,*)' Active Fock matrix at Work(LFA_RED).'
*      write(6,'(1x,8f10.4)')(Work(LFA_RED+i),i=0,nRS-1)
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
     &                  FactXI,Work(ip_HTVec(iSymk)),NK*JNUM,
     &                  Work(ip_HTVec(iSymk)),NK*JNUM,
     &                  1.0D0,FFAO(ISFF),NB)
           EndIf
           ISFF = ISFF+(NB*(NB+1))/2
          END DO

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
     &                   FactXI,Work(ip_HTVec(iSymk)),NK*JNUM,
     &                   Work(ip_HTVec(iSymk)),NK*JNUM,
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
* ---------------------------------------------------
* End loop ISYQ
      END DO
      END IF
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
      CALL HALFTRNSF(IRC,WORK(IP_CHSPC),NCHSPC,1,JV1,JNUM,JNUM,
     &    JSYM,JREDC,WORK(LCNAT),ISTART,NUSE,IP_HTVEC)
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
     &               FactXA,Work(ip_HTVec(iSymw)),NW*JNUM,
     &               Work(ip_HTVec(iSymw)),NW*JNUM,
     &               1.0D0,FAAO(ISFA),NB)
       EndIf
       ISFA = ISFA+(NB*(NB+1))/2
      END DO

*      write(6,*)' Active Fock matrix in FAAO.'
*      write(6,'(1x,8f10.4)')(FAAO(i),i=1,nbtri)

* ---------------------------------------------------
* Active half-transformation:
      CALL HALFTRNSF(IRC,WORK(IP_CHSPC),NCHSPC,1,JV1,JNUM,JNUM,
     &    JSYM,JREDC,CMO,ISTART,NUSE,IP_HTVEC)


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
      END IF

* 800  CONTINUE
* End loop IBATCH
       JV1=JV1+JNUM
      END DO

      IF (jSym.eq.1) THEN
* Add Coulomb contributions in local Fock matrices (in 'reduced storage')
* into the global ones:
        CALL red2full(FFAO,Work(LFF_RED))
        CALL red2full(FIAO,Work(LFI_RED))
        CALL red2full(FAAO,Work(LFA_RED))
      END IF
* End loop JRED
  999 CONTINUE
      END DO

      IF (jSym.eq.1) THEN
* Deallocate local density and fock matrices
        CALL GETMEM('V_VECTOR','FREE','REAL',IP_V,MXCHARR)
        CALL GETMEM('DF_RED','FREE','REAL',IPDF_RED,MXCHARR)
        CALL GETMEM('DI_RED','FREE','REAL',IPDI_RED,MXCHARR)
        CALL GETMEM('DA_RED','FREE','REAL',IPDA_RED,MXCHARR)
        CALL GETMEM('FF_RED','FREE','REAL',LFF_RED,MXCHARR)
        CALL GETMEM('FI_RED','FREE','REAL',LFI_RED,MXCHARR)
        CALL GETMEM('FA_RED','FREE','REAL',LFA_RED,MXCHARR)
      END IF
* End loop JSYM
 1000 CONTINUE
      END DO

      ! if using the RHS on-demand, we need all cholesky vectors on each
      ! process, collect them here
      IF (IF_TRNSF.AND.RHSDIRECT) THEN
        IP_LFT=IP_FTSPC
        DO JSYM=1,NSYM
          NVTOT=NVTOT_CHOSYM(JSYM)
          IBSTA=NBTCHES(JSYM)+1
          IBEND=NBTCHES(JSYM)+NBTCH(JSYM)
          DO IB=IBSTA,IBEND
            NV=NVLOC_CHOBATCH(IB)
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

* Synchronize and add the contributions from all nodes into each node:
      CALL GADGOP(FFAO,NBTRI,'+')
      CALL GADGOP(FIAO,NBTRI,'+')
      CALL GADGOP(FAAO,NBTRI,'+')

* Two-electron contribution to the effective core energy
      ECORE2=0.5D0*DDOT_(NBTRI,WORK(IPDF),1,FFAO,1)
c Add OneHam to finalize frozen Fock matrix in AO basis.
c (It is in fact an effective one-electron Hamiltonian).
      CALL ADD1HAM(FFAO)
* The contraction of frozen Fock matrix with frozen density:
      E=DDOT_(NBTRI,WORK(IPDF),1,FFAO,1)
* Correct for double-counting two-electron part:
      E=E-ECORE2
* One-electron part:
      ECORE1=E-ECORE2
* Nuclear repulsion energy:
      ECORE=POTNUC+ECORE1+ECORE2

#ifdef _DEBUGPRINT_
       WRITE(6,'(6X,A,E20.10)') 'NUCLEAR REPULSION ENERGY:',POTNUC
       WRITE(6,'(6X,A,E20.10)') 'ONE-ELECTRON CORE ENERGY:',ECORE1
       WRITE(6,'(6X,A,E20.10)') 'TWO-ELECTRON CORE ENERGY:',ECORE2
       WRITE(6,'(6X,A,E20.10)') '       TOTAL CORE ENERGY:',ECORE
#endif

      Call Getmem('OCC','FREE','REAL',LOCC,NBasT)
      Call Getmem('CNAT','FREE','REAL',LCNAT,NBSQT)
      Call Getmem('DF','FREE','REAL',ipDF,NBTRI)
      Call Getmem('DI','FREE','REAL',ipDI,NBTRI)
      Call Getmem('DA ','FREE','REAL',ipDA ,NBTRI)

      CALL GETMEM('CHSPC','FREE','REAL',IP_CHSPC,NCHSPC)
      CALL GETMEM('HTSPC','FREE','REAL',IP_HTSPC,NHTSPC)
      IF (IF_TRNSF) THEN
       CALL GETMEM('FTSPC','FREE','REAL',IP_FTSPC,NFTSPC)
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

      Call QExit('TRACHO2')
      RETURN
      END
