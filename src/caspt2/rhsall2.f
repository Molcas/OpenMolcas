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
      SUBROUTINE RHSALL2(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      USE CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_global, only:iPrGlb, FIMO, PIQK, Buff, idxb
      use PrintLevel, only: VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, NISH, NASH, NSSH, NASHT, NBTCHES,
     &                         NBTCH, NAES
#ifdef _DEBUGPRINT_
      use caspt2_module, only: NASUP, NISUP
#endif
      IMPLICIT None
* ----------------------------------------------------------------
* Code for processing all the cholesky vectors
* in construction of caspt2 right-hand-side array
* Also form the active two-electron integrals 'TUVX'.
* ================================================================
#include "warnings.h"
      integer(kind=iwp), intent(in):: IVEC
*
      integer(kind=iwp), Parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) nSh(8,3)
      integer(kind=iwp), SAVE :: NUMERR=0
      real(kind=wp), allocatable:: TUVX(:), BRA(:), KET(:)

      integer(kind=iwp),allocatable:: BGRP(:,:)
      integer(kind=iwp) IB, IB1, IB2, IBEND, IBGRP, IBSTA, iOffi, iOffK,
     &                  iOffp, iOffQ, ISYI, ISYK, ISYP, ISYQ, JSYM,
     &                  LBRASM, LKETSM, MXBGRP, MXPIQK, NADDBUF, NBGRP,
     &                  nBra, NBRASM, NCHOBUF, NG1, NG2, NI, NK, nKet,
     &                  NKETSM, NP, NPI, NQ, NQK, NTUVX, NV
#ifdef _DEBUGPRINT_
      real(kind=wp) DNRM2
      real(kind=wp), external:: RHS_DDot
      integer(kind=iwp) ISYM, lg_W, NAS, NIS, ICASE
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
*                                                                      *
************************************************************************
*                                                                      *

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'
      END IF

*
* TUVX RHSX        Na^4
* TJVX RHSA        Na^3 Ni
* TJVL RHSB        Na^2 Ni^2
* AJVX RHSD1       Na^2 Ni   Ns
* AJCL RHSH             Ni^2 Ns^2     N^4
* AUVX RHSC        Na^3      Ns
* AUCX RHSF        Na^2      Ns^2
* AUVL RHSD2       Na^2 Ni   Ns
* AUCL RHSG        Na   Ni   Ns^2     N^3
* AJVL RHSE        Na   Ni^2 Ns       N^3
*
*     Initialize RHS array as 'vector' nr IVEC=IRHS on LUSOLV:
*
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate and clear TUVX, two-electron integrals for active
*     orbital indices only: Simple storage, same as for GAMMA2.
*     TUVX is kept allocated until end of subroutine.
      NG1=NASHT**2
      NG2=NG1**2
      NTUVX=NG2
      CALL mma_allocate(TUVX,NTUVX,Label='TUVX')
      TUVX(:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
      DO JSYM=1,NSYM
*
       IB1=NBTCHES(JSYM)+1
       IB2=NBTCHES(JSYM)+NBTCH(JSYM)
*
       MXBGRP=IB2-IB1+1
       IF (MXBGRP.LE.0) CYCLE
       CALL mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
       IBGRP=1
       DO IB=IB1,IB2
        BGRP(1,IBGRP)=IB
        BGRP(2,IBGRP)=IB
        IBGRP=IBGRP+1
       END DO
       NBGRP=MXBGRP

       CALL MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,NCHOBUF,MXPIQK,NADDBUF)
       IF (IPRGLB.GT.VERBOSE) THEN
         WRITE(6,*)
         WRITE(6,'(A,I12)') '  Number of Cholesky batches: ',IB2-IB1+1
         WRITE(6,'(A,I12)') '  Number of batch groups:     ',NBGRP
         WRITE(6,*)
       END IF

* buffers are kept allocated until the end of JSYM loop.
       CALL mma_allocate(PIQK,MXPIQK,Label='PIQK')
       CALL mma_allocate(BUFF,NADDBUF,Label='BUFF')
       CALL mma_allocate(IDXB,NADDBUF,Label='IDXB')

       CALL mma_allocate(BRA,NCHOBUF,Label='BRA')
       CALL mma_allocate(KET,NCHOBUF,Label='KET')
*
*      Loop over groups of batches of Cholesky vectors
*
*      IBSTEP=1
*
*      DO IBSTA=IB1,IB2,IBSTEP
       DO IBGRP=1,NBGRP
        IBSTA=BGRP(1,IBGRP)
        IBEND=BGRP(2,IBGRP)

        NV=0
        DO IB=IBSTA,IBEND
           NV=NV+NVLOC_CHOBATCH(IB)
        END DO

        IF (IPRGLB.GT.VERBOSE) THEN
         WRITE(6,'(A,I12)') '  Cholesky vectors in this group = ', NV
         WRITE(6,*)
        END IF
*                                                                      *
************************************************************************
*                                                                      *
*      Read kets (Cholesky vectors) in the form L(VX), all symmetries:
*
       Call Get_Cholesky_Vectors(Active,Active,JSYM,
     &                           KET,nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble contributions to TUVX integrals
*      Reuse the ket vectors as L(TU) bra vectors
*
       LBRASM=1
       DO ISYI=1,NSYM
        NI=NASH(ISYI)
        iOffi=NAES(iSYI)
        IF(NI.EQ.0) Cycle
        ISYP=Mul(ISYI,JSYM)
        NP=NASH(ISYP)
        iOffp=NAES(iSYP)
        IF(NP.EQ.0) Cycle
        NPI=NP*NI
        NBRASM=NPI*NV
        LKETSM=1

        DO ISYK=1,NSYM
         NK=NASH(ISYK)
         iOffK=NAES(iSYK)
         IF(NK.EQ.0) Cycle
         ISYQ=Mul(ISYK,JSYM)
         NQ=NASH(ISYQ)
         iOffQ=NAES(iSYQ)
         IF(NQ.EQ.0) Cycle
         NQK=NQ*NK
         NKETSM=NQK*NV
*
         IF (NPI*NQK.GT.mxPIQK) THEN
           WRITE(6,*) 'NPIQK larger than mxPIQK in TUVX, bug?'
           Call AbEnd()
         END IF
         CALL DGEMM_('N','T',NPI,NQK,NV,One,KET(LBRASM),NPI,
     &        KET(LKETSM),NQK,Zero,PIQK,NPI)
*
         Call ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,
     &                TUVX,nTUVX,PIQK,NPI*NQK,NUMERR)
*
         LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
       END DO
*                                                                      *
************************************************************************
*                                                                      *
*      Read bra (Cholesky vectors) in the form L(TJ): All symmetries
*
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                           BRA,nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble contributions to TJVX
*      Loop over the bras and kets, form <A|0>
*
      Call Process_RHS_Block(Inactive,Active,Active,Active,
     &                       'A ',
     &                       BRA,nBra,KET,nKet,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
*      TJVL RHSB
*      TJVL: Use TJ buffer as if it was VL, form <B|0>
*
      Call Process_RHS_Block(Inactive,Active,Inactive,Active,
     &                       'B ',
     &                       BRA,nBra,BRA,nBra,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
* We still have L(VX) vectors in core, at KET.
*
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                           BRA,nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVX RHSD1
* Loop over the bra and ket vectors.
*
      Call Process_RHS_Block(Inactive,Virtual,Active,Active,
     &                       'D1',
     &                       BRA,nBra,KET,nKet,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* AJCL RHSH
* AJCL: Use AJ buffer still in core as if it was CL, form <H|0>
*
      Call Process_RHS_Block(Inactive,Virtual,Inactive,Virtual,
     &                       'H ',
     &                       BRA,nBra,BRA,nBra,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read Bra (Cholesky vectors)= L(AU)
*
       Call Get_Cholesky_Vectors(Active,Virtual,JSYM,
     &                           BRA,nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVX RHSC
* AUVX: Loop over the bras and kets
*
      Call Process_RHS_Block(Active,Virtual,Active,Active,
     &                       'C ',
     &                       BRA,nBra,KET,nKet,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* AUCX RHSF
* AUCX: Use AU buffer still in core as if it was CX, form <F|0>
*
      Call Process_RHS_Block(Active,Virtual,Active,Virtual,
     &                       'F ',
     &                       BRA,nBra,BRA,nBra,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(VL), all symmetries:
*
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                           KET,nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVL RHSD2
* Loop over bras and kets, form <D2|0>.
*
      Call Process_RHS_Block(Active,Virtual,Inactive,Active,
     &                       'D2',
     &                       BRA,nBra,KET,nKet,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(CL), all symmetries:
*
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                           KET,nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUCL RHSG
* Loop over bras and kets, form  <G|0>
*
      Call Process_RHS_Block(Active,Virtual,Inactive,Virtual,
     &                       'G ',
     &                       BRA,nBra,KET,nKet,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* Read bra vectors AJ
*
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,
     &                           BRA,nBra,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets in the form L(VL)
*
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                           KET,nKet,
     &                           IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVL RHSE
* AJVL: Loop over bras and kets. Form <E|0>
*
      Call Process_RHS_Block(Inactive,Virtual,Inactive,Active,
     &                       'E ',
     &                       BRA,nBra,KET,nKet,
     &                       nSh,JSYM,IVEC,NV)
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over batches, IB
      END DO
*                                                                      *
************************************************************************
*                                                                      *
      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)
      CALL mma_deallocate(PIQK)
      CALL mma_deallocate(BUFF)
      CALL mma_deallocate(IDXB)
      CALL mma_deallocate(BGRP)
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over JSYM
      END DO
*                                                                      *
************************************************************************
*                                                                      *
* Synchronized add RHS partial arrays from all nodes into each node.

C-SVC: read the DRA's from disk and copy them all to LUSOLV to continue
C      in serial mode.  FIXME: this call has to be removed when we reach
C      full parallel capabilities
*     CALL SYNRHS(IVEC)
C-SVC: at this point, the RHS elements are on disk, both in LUSOLV and
C      as DRAs with the name RHS_XX_XX_XX with XX a number representing
C      the case, symmetry, and rhs vector respectively.

* The RHS elements of Cases A, C, D1  need a correction:
      CALL MODRHS(IVEC,FIMO,SIZE(FIMO))

#ifdef _DEBUGPRINT_
* compute and print RHS fingerprints
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') 'Case','Sym','Fingerprint'
      WRITE(6,'(1X,A4,1X,A3,1X,A18)') '====','===','==========='
      DO ICASE=1,13
        DO ISYM=1,NSYM
          NAS=NASUP(ISYM,ICASE)
          NIS=NISUP(ISYM,ICASE)
          IF (NAS*NIS.NE.0) THEN
            CALL RHS_ALLO (NAS,NIS,lg_W)
            CALL RHS_READ (NAS,NIS,lg_W,iCASE,iSYM,iVEC)
            DNRM2 = RHS_DDOT(NAS,NIS,lg_W,lg_W)
            WRITE(6,'(1X,I4,1X,I3,1X,F18.11)') ICASE,ISYM,DNRM2
          END IF
        END DO
      END DO
#endif

* Synchronized add tuvx partial arrays from all nodes into each node.
      CALL CHO_GADGOP(TUVX,NTUVX,'+')
* Put TUVX on disk for possible later use:
      CALL PT2_PUT(NTUVX,'TUVX',TUVX)
      CALL mma_deallocate(TUVX)
*                                                                      *
************************************************************************
*                                                                      *
      END SUBROUTINE RHSALL2

      Subroutine Get_Cholesky_Vectors(ITK,ITQ,JSYM,
     &                                Array,nArray,
     &                                IBSTA,IBEND)
      use definitions, only: iwp, wp
      USE CHOVEC_IO, only: NPQ_CHOTYPE, NVLOC_CHOBATCH, IDLOC_CHOGROUP
      use caspt2_global, only: LUDRA
      use caspt2_module, only: NSYM
      IMPLICIT None
      integer(kind=iwp), Intent(in):: ITK,ITQ,JSYM,IBSTA,IBEND
      integer(kind=iwp), Intent(Out):: nArray
      real(kind=wp), intent(Out):: Array(*)

      integer(kind=iwp) ICASE, LKETSM, ISYK, NQK, IB, NV, NKETSM, IDISK

      ! ugly hack to convert separate k/q orbital types into a specific
      ! case
      ICASE=ITK*ITQ
      IF (ICASE.EQ.3) THEN
        ICASE=4
      ELSE
        ICASE=ICASE/2
      END IF

      LKETSM=1
      DO ISYK=1,NSYM
        NQK=NPQ_CHOTYPE(ICASE,ISYK,JSYM)
        IF(NQK.EQ.0) CYCLE
        DO IB=IBSTA,IBEND
          NV=NVLOC_CHOBATCH(IB)
          NKETSM=NQK*NV
          IDISK=IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
          CALL DDAFILE(LUDRA,2,Array(LKETSM),NKETSM,IDISK)
          LKETSM=LKETSM+NKETSM
        END DO
      END DO
      nArray=LKETSM-1
*
      End Subroutine Get_Cholesky_Vectors

      Subroutine Process_RHS_Block(ITI,ITP,ITK,ITQ,
     &                             Case,
     &                             Cho_Bra,nBra,Cho_Ket,nKet,
     &                             nSh,JSYM,IVEC,NV)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use caspt2_global, only: iPrGlb, PIQK, BUFF, idxb
      use PrintLevel, only: DEBUG
      use caspt2_module, only: NSYM
      IMPLICIT None
      integer(kind=iwp), Intent(in):: ITI,ITP,ITK,ITQ
      Character(LEN=2), intent(in)::  Case
      integer(kind=iwp), intent(in):: nBra, nKet
      real(kind=wp), intent(in):: Cho_Bra(nBra), Cho_Ket(nKet)
      integer(kind=iwp), intent(in):: nSh(8,3), JSYM, iVec, nV

      integer(kind=iwp) ISYI, ISYK, ISYP, ISYQ, KPI, KQK, LBRASM,
     &                  LKETSM, NBRASM, NI, NK, NKETSM, NP, NPI, NPIQK,
     &                  NQ, NQK
      integer(kind=iwp) mxPIQK, nBuff
      mxPIQK=Size(PIQK)
      nBuff=Size(BUFF)
*
*
      IF (iPrGlb.GE.DEBUG) THEN
        WRITE(6,*) 'Processing RHS block '//Case
      END IF

      LBRASM=1
      DO ISYI=1,NSYM
         NI=NSH(ISYI,ITI)
         IF(NI.EQ.0) Cycle
         ISYP=Mul(ISYI,JSYM)
         NP=NSH(ISYP,ITP)
         IF(NP.EQ.0) Cycle
         NPI=NP*NI
         NBRASM=NPI*NV
*
         LKETSM=1
         DO ISYK=1,NSYM
            NK=NSH(ISYK,ITK)
            IF(NK.EQ.0) Cycle
            ISYQ=Mul(ISYK,JSYM)
            NQ=NSH(ISYQ,ITQ)
            IF(NQ.EQ.0) Cycle
            NQK=NQ*NK
            NKETSM=NQK*NV
*
C SVC: we need an NPI*NQK to store the 2-electron integrals, and 2
C buffers (values+indices) for sorting them.  Later, we can try to get
C rid of the buffer that stores the values and only use an index buffer
C and the two-electron integrals for the scatter operation.  For the
C buffer, any size can be taken, but assuming there is enough memory
C available, it's set to the size of the two-electron integrals unless
C larger than some predefined maximum buffer size.
            NPIQK=NPI*NQK
            IF (NPIQK.GT.MXPIQK) THEN
              IF (Case.eq.'H') THEN
                KPI=MXPIQK/NQK
                NPIQK=KPI*NQK
              ELSE IF (Case.eq.'G') THEN
                KQK=MXPIQK/NPI
                NPIQK=NPI*KQK
              ELSE
                WRITE(6,*) ' NPIQK > MXPIQK and case != G or H'
                WRITE(6,'(A,A2)')  ' CASE =   ', Case
                WRITE(6,'(A,I12)') ' NPIQK =  ', NPIQK
                WRITE(6,'(A,I12)') ' MXPIQK = ', MXPIQK
                WRITE(6,*) ' This should not happen, please report.'
                CALL AbEnd()
              END IF
            END IF

C-SVC: sanity check
            IF (NPIQK.LE.0) THEN
              WRITE(6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
              CALL AbEnd()
            END IF
*
            If (Case.eq.'A ') Then
               CALL ADDRHSA(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'B ') Then
               CALL ADDRHSB(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'D1') Then
               CALL ADDRHSD1(IVEC,JSYM,ISYI,ISYK,
     &                       NP,NI,NQ,NK,PIQK,
     &                       nBuff,Buff,idxb,
     &                       Cho_Bra(LBRASM),
     &                       Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'H ') Then
               CALL ADDRHSH(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,NPIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'C ') Then
               CALL ADDRHSC(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'F ') Then
               CALL ADDRHSF(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'D2') Then
               CALL ADDRHSD2(IVEC,JSYM,ISYI,ISYK,
     &                       NP,NI,NQ,NK,PIQK,
     &                       nBuff,Buff,idxb,
     &                       Cho_Bra(LBRASM),
     &                       Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'G ') Then
               CALL ADDRHSG(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,NPIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else If (Case.eq.'E ') Then
               CALL ADDRHSE(IVEC,JSYM,ISYI,ISYK,
     &                      NP,NI,NQ,NK,PIQK,
     &                      nBuff,Buff,idxb,
     &                      Cho_Bra(LBRASM),
     &                      Cho_Ket(LKETSM),NV)
            Else
               Call Abend()
            End If
*
         LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
      END DO
*
      End Subroutine Process_RHS_Block

      Subroutine ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,
     &                   TUVX,nTUVX,PIQK,nPIQK,
     &                   NUMERR)
      use definitions, only: iwp, wp
#ifndef _DEBUGPRINT_
      use Constants, only: One
#endif
      Implicit None
      integer(kind=iwp), intent(in):: NP,NI,NQ,NK,NASHT,iOffP,iOffI,
     &                                iOffQ,iOffK
      integer(kind=iwp), intent(in):: nTUVX, nPIQK
      real(kind=wp), intent(in)::  PIQK(nPIQK)
      real(kind=wp), intent(inout):: TUVX(nTUVX)
      integer(kind=iwp), intent(inout):: NUMERR

      integer(kind=iwp) iU, iUVX1, iUVX2, iV, iVX1, iVX2, iX, iX1, iX2
#ifndef _DEBUGPRINT_
#include "macros.fh"
      unused_var(NUMERR)
#else
      integer(kind=iwp) iPIQK, iT, iTUVX
#endif
*
* Add into correct positions in TUVX:
*
      DO iX=0,NK-1
         iX1=NASHT*(iX+iOffK)
         iX2=NQ   * iX
         DO iV=0,NQ-1
            iVX1=NASHT*(iV+iOffQ+iX1)
            iVX2=NI   *(iV      +iX2)
            DO iU=0,NI-1
               iUVX1=NASHT*(iU+iOffI+iVX1)
               iUVX2=NP   *(iU+      iVX2)
#ifdef _DEBUGPRINT_
               DO iT=1,NP
                  iTUVX=iT+iOffP+iUVX1
                  iPIQK=iT      +iUVX2
* Temporary test statements -- remove after debug!
                  IF(ITUVX.LT.1 .or. ITUVX.gt.NTUVX+1) THEN
                     ITUVX=NTUVX+1
                     NUMERR=NUMERR+1
                     IF (NUMERR.GT.100) THEN
                        WRITE(6,*)' THIS IS TOO MUCH -- STOP.'
                        CALL QUIT(_RC_INTERNAL_ERROR_)
                     END IF
                  END IF
* End of temporary test statements
                  TUVX(iTUVX)=TUVX(iTUVX)+PIQK(iPIQK)
               END DO
#else
               Call DaXpY_(nP,One,PIQK(1+      iUVX2),1,
     &                             TUVX(1+iOffP+iUVX1),1)
#endif
            END DO
         END DO
      END DO
*
      End Subroutine ADDTUVX

      SUBROUTINE MEMORY_ESTIMATE(JSYM,LBGRP,NBGRP,NCHOBUF,NPIQK,NADDBUF)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp
      USE CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_global, only: iParRHS,iPrGlb,iStpGrd
      use PrintLevel, only: VERBOSE
      use stdalloc, only: mma_MaxDBLE
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use caspt2_module, only: NSYM, NASHT, NISUP, NISH, NASH,
     &                         NSSH, NTU, NTUV, NASH, NIGEJ, NIGTJ,
     &                         NAGEB, NAGTB, NTGEU, NTGTU, NBTCHES,
     &                         NBTCH

      IMPLICIT None
      integer(kind=iwp), intent(in):: JSYM
      integer(kind=iwp), intent(out):: NBGRP
      integer(kind=iwp), intent(out):: LBGRP(2,*)
      integer(kind=iwp), intent(out) :: NCHOBUF,NPIQK,NADDBUF

      integer(kind=iwp) IB, IB1, IB2, IBGRP, ICASE, ISYI, ISYK, ISYP,
     &                  ISYQ, MAXBUFF, MAXCHOL, MAXPIQK, MINBUFF,
     &                  MINCHOL, MINGOOD, MINNICE, MINPIQK, MINSLOW,
     &                  MXAVAIL, MXBATCH, MXCHOVEC, MXNPITOT, MXRHS,
     &                  NCHOVEC, NCHUNK, NI, NK, NP, NPI, NQ, NQK, NV,
     &                  NVECTOT
      integer(kind=iwp), external:: iPARDIV
      integer(kind=iwp), Parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) nSh(8,3)
      Logical(kind=iwp) :: call_from_grad
      integer(kind=iwp) :: ITYPE(4,9)=reshape([
     &                                Inactive,Active,Active,Active,
     &                                Inactive,Active,Inactive,Active,
     &                                Inactive,Virtual,Active,Active,
     &                                Inactive,Virtual,Inactive,Virtual,
     &                                Active,Virtual,Active,Active,
     &                                Active,Virtual,Active,Virtual,
     &                                Active,Virtual,Inactive,Active,
     &                                Active,Virtual,Inactive,Virtual,
     &                                Inactive,Virtual,Inactive,Active]
     &                               ,[4,9])
      integer(kind=iwp) ISYM

      call_from_grad = .false.
      if (iStpGrd == -1) call_from_grad = .true.

      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)

CSVC: compute the maximum RHS size that will be allocated per node, this
C     will be later allocated when addrhs routines are called. If global
C     arrays are used, this may actually be done outside of GETMEM, but
C     let's just reserve it anyway.
      ISYM=JSYM
      MXRHS=2*NTU(ISYM)*NISUP(ISYM,5)
      DO ISYI=1,NSYM
        ISYM=ISYI
C case A
        MXRHS=Max(MXRHS,NTUV(ISYM)*NISH(ISYM))
C case GP,GM
        MXRHS=Max(MXRHS,NASH(ISYM)*NISUP(ISYM,10)
     &                 +NASH(ISYM)*NISUP(ISYM,11))

        ISYM=Mul(JSYM,ISYI)
C case C
        MXRHS=Max(MXRHS,NTUV(ISYM)*NSSH(ISYM))
C case EP,EM
        MXRHS=Max(MXRHS,NASH(ISYM)*NISUP(ISYM,6)
     &                 +NASH(ISYM)*NISUP(ISYM,7))
        DO ISYK=1,NSYM
          ISYM=Mul(ISYI,ISYK)
C case BP,BM
          MXRHS=Max(MXRHS,NTGEU(ISYM)*NIGEJ(ISYM))
          MXRHS=Max(MXRHS,NTGTU(ISYM)*NIGTJ(ISYM))
C case HP,HM
          MXRHS=Max(MXRHS,NAGEB(ISYM)*NIGEJ(ISYM))
          MXRHS=Max(MXRHS,NAGTB(ISYM)*NIGTJ(ISYM))
C case F
          MXRHS=Max(MXRHS,NTGEU(ISYM)*NAGEB(ISYM))
          MXRHS=Max(MXRHS,NTGTU(ISYM)*NAGTB(ISYM))

          ISYM=Mul(ISYI,Mul(JSYM,ISYK))
C case D1,D2
          MXRHS=Max(MXRHS,2*NTU(ISYM)*NISUP(ISYM,5))
        End Do
      End Do
      MXRHS = iPARDIV(MXRHS,0)

CSVC: determine maximum pair index size per symmetry and in total.
*     MXBFSZ=0
      MXNPITOT=0
      DO ISYI=1,NSYM
        ISYP=Mul(ISYI,JSYM)
        NI=MAX(NISH(ISYI),NASH(ISYI))
        NP=MAX(NASH(ISYP),NSSH(ISYP))
        MXNPITOT=MXNPITOT+NP*NI
*       MXBFSZ=MXBFSZ+NP*NI*NJSCT
      END DO

CSVC: determine maximum and minimum size to hold integral matrix.
C     cases here correspond to A, B, D1, H, C, F, D2, G, E.
C     with number icase =      1, 2,  3, 4, 5, 6,  7, 8, 9.
C     both case G and H can use blocking, currently this is done over
C     the number of inactive orbitals in one of the pair indices.
C     for case G, the minimum needed is NAU*NC, because of blocking NL.
C     for case H, the minimum needed is NA*NCL, because of blocking NJ.
C     to start, reserve space for TUVX integrals (NASHT**4)
      MAXPIQK=NASHT**4
      MINPIQK=NASHT**4
      DO ICASE=1,9
        DO ISYI=1,NSYM
          NI=NSH(ISYI,ITYPE(1,ICASE))
          ISYP=Mul(ISYI,JSYM)
          NP=NSH(ISYP,ITYPE(2,ICASE))
          NPI=NP*NI
          DO ISYK=1,NSYM
            NK=NSH(ISYK,ITYPE(3,ICASE))
            ISYQ=Mul(ISYK,JSYM)
            NQ=NSH(ISYQ,ITYPE(4,ICASE))
            NQK=NQ*NK
            MAXPIQK=MAX(MAXPIQK,NPI*NQK)
            IF (ICASE.EQ.4) THEN
              MINPIQK=MAX(MINPIQK,NP*NQK)
            ELSE IF (ICASE.EQ.8) THEN
              MINPIQK=MAX(MINPIQK,NPI*NQ)
            ELSE
              MINPIQK=MAX(MINPIQK,NPI*NQK)
            END IF
          END DO
        END DO
      END DO
      MAXBUFF=NINT(SQRT(DBLE(MAXPIQK)))
      MINBUFF=NINT(SQRT(DBLE(MINPIQK)))

      ! In some cases, we do not use the buffer arrays
      if (call_from_grad .or. iParRHS == 2) then
        MAXBUFF = 1
        MINBUFF = 1
      end if

CSVC: total number of cholesky vectors
      MXBATCH=0
      NVECTOT=0
      IB1=NBTCHES(JSYM)+1
      IB2=NBTCHES(JSYM)+NBTCH(JSYM)
      DO IB=IB1,IB2
        NV=NVLOC_CHOBATCH(IB)
        MXBATCH=MAX(MXBATCH,NV)
        NVECTOT=NVECTOT+NV
      END DO
      MAXCHOL=MXNPITOT*NVECTOT
      MINCHOL=MXNPITOT*MXBATCH

CSVC: can we fit this all in memory?
      CALL mma_MaxDBLE(MXAVAIL)

      MINNICE=MXRHS+MAXPIQK+2*MAXBUFF+2*MAXCHOL
      MINGOOD=MXRHS+MINPIQK+2*MINBUFF+2*MAXCHOL
      MINSLOW=MXRHS+MINPIQK+2*MINBUFF+2*MINCHOL

      if (call_from_grad) then
#ifdef _MOLCAS_MPP_
        if (is_real_par()) then
          !! One more vector is needed for buffer
          MAXCHOL = MAXCHOL + MXNPITOT
          MINNICE=2*MXRHS+MAXPIQK+2*MAXBUFF+4*MAXCHOL
          MINGOOD=2*MXRHS+MINPIQK+2*MINBUFF+4*MAXCHOL
          MINSLOW=2*MXRHS+MINPIQK+2*MINBUFF+4*MINCHOL
        else
#endif
          MINNICE=  MXRHS+MAXPIQK+2*MAXBUFF+4*MAXCHOL
          MINGOOD=  MXRHS+MINPIQK+2*MINBUFF+4*MAXCHOL
          MINSLOW=  MXRHS+MINPIQK+2*MINBUFF+4*MINCHOL
#ifdef _MOLCAS_MPP_
        end if
#endif
      end if

      IF (IPRGLB.GT.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(A,I1)') '  Memory estimates in RHSALL, SYM ', JSYM
        WRITE(6,'(A,2X,I16)') '   allocatable:    ', MXAVAIL
        WRITE(6,'(A,2X,I16)') '   recommended:    ', MINNICE
        WRITE(6,'(A,2X,I16)') '   convenient:     ', MINGOOD
        WRITE(6,'(A,2X,I16)') '   minimum:        ', MINSLOW
        WRITE(6,*)
        IF (MXAVAIL.GE.MINNICE) THEN
          WRITE(6,*) ' I can use all cholesky vectors at once'
          WRITE(6,*) ' as well as the whole integral matrixi.'
        ELSE IF (MXAVAIL.GE.MINGOOD) THEN
          WRITE(6,*) ' I will group batches of cholesky vectors'
          WRITE(6,*) ' and then maximize use of the integral matrix.'
        ELSE IF (MXAVAIL.GE.MINSLOW) THEN
          WRITE(6,*) ' I will at least try to group batches.'
        ELSE
          WRITE(6,*) ' Do you see my problem?'
        END IF
      END IF

      IF (MXAVAIL.GE.MINNICE) THEN
C group all batches, take the maximum needed for integrals and try
C to max out the buffer size, check it is larger than the minimum.
        NCHOBUF=MAXCHOL
        NBGRP=1
        LBGRP(1,1)=IB1
        LBGRP(2,1)=IB2
        NADDBUF=MAXBUFF
        NPIQK=MAXPIQK
      ELSE IF (MXAVAIL.GE.MINGOOD) THEN
C group all batches, take smaller buffer size, and try to max out
C integrals, and check they are lager than minimum needed
        NCHOBUF=MAXCHOL
        NBGRP=1
        LBGRP(1,1)=IB1
        LBGRP(2,1)=IB2
        NADDBUF=MINBUFF
        NCHUNK=(MXAVAIL-MXRHS-2*NADDBUF-2*NCHOBUF)/MINPIQK
        if (call_from_grad)
     *    NCHUNK=(MXAVAIL-2*MXRHS-2*NADDBUF-4*NCHOBUF)/MINPIQK
        NPIQK=MINPIQK*NCHUNK
      ELSE IF (MXAVAIL.GE.MINSLOW) THEN
        NADDBUF=MINBUFF
        NPIQK=MINPIQK
        NCHOBUF=(MXAVAIL-MXRHS-NPIQK-2*NADDBUF)/2
        if (call_from_grad)
     *    NCHOBUF=(MXAVAIL-2*MXRHS-NPIQK-4*NADDBUF)/2
C create batch groups that have at most MXCHOVEC cholesky vectors
        MXCHOVEC=MAX(NCHOBUF/MXNPITOT,1)
        NCHOVEC=0
        IBGRP=1
        LBGRP(1,IBGRP)=IB1
        DO IB=IB1,IB2
          NV=NVLOC_CHOBATCH(IB)
          NCHOVEC=NCHOVEC+NV
          IF (NCHOVEC.GT.MXCHOVEC) THEN
            LBGRP(2,IBGRP)=IB-1
            IBGRP=IBGRP+1
            LBGRP(1,IBGRP)=IB
            NCHOVEC=NV
          END IF
        END DO
        LBGRP(2,IBGRP)=IB2
        NBGRP=IBGRP
      ELSE
        WRITE(6,*)
        WRITE(6,*) '  Not enough memory in RHSLL2...'
        WRITE(6,'(A,I16)') '   allocatable:    ', MXAVAIL
        WRITE(6,'(A,I16)') '   minimum:        ', MINSLOW
        Call AbEnd()
      END IF

CSVC: sanity check, should not happen.
      IF (MXRHS.GT.MXAVAIL-2*NCHOBUF-NPIQK-2*NADDBUF) THEN
        WRITE(6,*)
        WRITE(6,*) 'RHSALL2: RHS allocation will starve.'
        WRITE(6,*) 'Possible bug in memory estimate.'
        WRITE(6,*) 'This should not happen, please report.'
        WRITE(6,*)
        WRITE(6,'(2X,A8,2X,I14)') 'MXAVAIL ', MXAVAIL
        WRITE(6,'(2X,A8,2X,I14)') 'MXRHS   ', MXRHS
        WRITE(6,'(2X,A8,2X,I14)') 'NCHOBUF ', NCHOBUF
        WRITE(6,'(2X,A8,2X,I14)') 'NPIQK   ', NPIQK
        WRITE(6,'(2X,A8,2X,I14)') 'NADDBUF ', NADDBUF
        CALL AbEnd()
      END IF

      IF (IPRGLB.GT.VERBOSE) THEN
        WRITE(6,*)
        WRITE(6,'(A16,2A16)') '  Buffer sizes:',
     &                        '           used',
     &                        '          ideal'
        WRITE(6,'(A16,2I16)') '  ChoVecs:  ', 2*NCHOBUF, 2*MAXCHOL
        WRITE(6,'(A16,2I16)') '  Integral: ',   NPIQK  ,   MAXPIQK
        WRITE(6,'(A16,2I16)') '  Scatter:  ', 2*NADDBUF, 2*MAXBUFF
        WRITE(6,*)
      END IF

      END SUBROUTINE MEMORY_ESTIMATE
