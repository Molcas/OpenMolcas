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
      SUBROUTINE RHSALL2(IVEC)
      use Symmetry_Info, only: Mul
      use definitions, only: iwp, wp
      use constants, only: Zero, One
      USE CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_global, only:iPrGlb, FIMO, PIQK, Buff, idxb
      use PrintLevel, only: VERBOSE
      use stdalloc, only: mma_allocate, mma_deallocate
      use caspt2_module, only: NSYM, NISH, NASH, NSSH, NASHT, NBTCHES,  &
     &                         NBTCH, NAES
#ifdef _DEBUGPRINT_
      use caspt2_module, only: NASUP, NISUP
#endif
      IMPLICIT None
! ----------------------------------------------------------------
! Code for processing all the cholesky vectors
! in construction of caspt2 right-hand-side array
! Also form the active two-electron integrals 'TUVX'.
! ================================================================
#include "warnings.h"
      integer(kind=iwp), intent(in):: IVEC
!
      integer(kind=iwp), Parameter :: Inactive=1, Active=2, Virtual=3
      integer(kind=iwp) nSh(8,3)
      integer(kind=iwp), SAVE :: NUMERR=0
      real(kind=wp), allocatable:: TUVX(:), BRA(:), KET(:)

      integer(kind=iwp),allocatable:: BGRP(:,:)
      integer(kind=iwp) IB, IB1, IB2, IBEND, IBGRP, IBSTA, iOffi, iOffK,&
     &                  iOffp, iOffQ, ISYI, ISYK, ISYP, ISYQ, JSYM,     &
     &                  LBRASM, LKETSM, MXBGRP, MXPIQK, NADDBUF, NBGRP, &
     &                  nBra, NBRASM, NCHOBUF, NG1, NG2, NI, NK, nKet,  &
     &                  NKETSM, NP, NPI, NQ, NQK, NTUVX, NV
#ifdef _DEBUGPRINT_
      real(kind=wp) DNRM2
      real(kind=wp), external:: RHS_DDot
      integer(kind=iwp) ISYM, lg_W, NAS, NIS, ICASE
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      Call ICopy(NSYM,NISH,1,nSh(1,Inactive),1)
      Call ICopy(NSYM,NASH,1,nSh(1,Active  ),1)
      Call ICopy(NSYM,NSSH,1,nSh(1,Virtual ),1)
!                                                                      *
!***********************************************************************
!                                                                      *

      IF (IPRGLB.GE.VERBOSE) THEN
        WRITE(6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'
      END IF

!
! TUVX RHSX        Na^4
! TJVX RHSA        Na^3 Ni
! TJVL RHSB        Na^2 Ni^2
! AJVX RHSD1       Na^2 Ni   Ns
! AJCL RHSH             Ni^2 Ns^2     N^4
! AUVX RHSC        Na^3      Ns
! AUCX RHSF        Na^2      Ns^2
! AUVL RHSD2       Na^2 Ni   Ns
! AUCL RHSG        Na   Ni   Ns^2     N^3
! AJVL RHSE        Na   Ni^2 Ns       N^3
!
!     Initialize RHS array as 'vector' nr IVEC=IRHS on LUSOLV:
!
!***********************************************************************
!                                                                      *
!                                                                      *
!***********************************************************************
!                                                                      *
!     Allocate and clear TUVX, two-electron integrals for active
!     orbital indices only: Simple storage, same as for GAMMA2.
!     TUVX is kept allocated until end of subroutine.
      NG1=NASHT**2
      NG2=NG1**2
      NTUVX=NG2
      CALL mma_allocate(TUVX,NTUVX,Label='TUVX')
      TUVX(:)=Zero
!                                                                      *
!***********************************************************************
!                                                                      *
      DO JSYM=1,NSYM
!
       IB1=NBTCHES(JSYM)+1
       IB2=NBTCHES(JSYM)+NBTCH(JSYM)
!
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

! buffers are kept allocated until the end of JSYM loop.
       CALL mma_allocate(PIQK,MXPIQK,Label='PIQK')
       CALL mma_allocate(BUFF,NADDBUF,Label='BUFF')
       CALL mma_allocate(IDXB,NADDBUF,Label='IDXB')

       CALL mma_allocate(BRA,NCHOBUF,Label='BRA')
       CALL mma_allocate(KET,NCHOBUF,Label='KET')
!
!      Loop over groups of batches of Cholesky vectors
!
!      IBSTEP=1
!
!      DO IBSTA=IB1,IB2,IBSTEP
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
!                                                                      *
!***********************************************************************
!                                                                      *
!      Read kets (Cholesky vectors) in the form L(VX), all symmetries:
!
       Call Get_Cholesky_Vectors(Active,Active,JSYM,                    &
     &                           KET,SIZE(KET),nKet,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
!      Assemble contributions to TUVX integrals
!      Reuse the ket vectors as L(TU) bra vectors
!
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
!
         IF (NPI*NQK.GT.mxPIQK) THEN
           WRITE(6,*) 'NPIQK larger than mxPIQK in TUVX, bug?'
           Call AbEnd()
         END IF
         CALL DGEMM_('N','T',NPI,NQK,NV,One,KET(LBRASM),NPI,            &
     &        KET(LKETSM),NQK,Zero,PIQK,NPI)
!
         Call ADDTUVX(NP,NI,NQ,NK,NASHT,iOffP,iOffI,iOffQ,iOffK,        &
     &                TUVX,nTUVX,PIQK,NPI*NQK,NUMERR)
!
         LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
       END DO
!                                                                      *
!***********************************************************************
!                                                                      *
!      Read bra (Cholesky vectors) in the form L(TJ): All symmetries
!
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,                  &
     &                           BRA,SIZE(BRA),nBra,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
!      Assemble contributions to TJVX
!      Loop over the bras and kets, form <A|0>
!
      Call Process_RHS_Block(Inactive,Active,Active,Active,             &
     &                       'A ',                                      &
     &                       BRA,nBra,KET,nKet,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
!      TJVL RHSB
!      TJVL: Use TJ buffer as if it was VL, form <B|0>
!
      Call Process_RHS_Block(Inactive,Active,Inactive,Active,           &
     &                       'B ',                                      &
     &                       BRA,nBra,BRA,nBra,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
! We still have L(VX) vectors in core, at KET.
!
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,                 &
     &                           BRA,SIZE(BRA),nBra,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
! AJVX RHSD1
! Loop over the bra and ket vectors.
!
      Call Process_RHS_Block(Inactive,Virtual,Active,Active,            &
     &                       'D1',                                      &
     &                       BRA,nBra,KET,nKet,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! AJCL RHSH
! AJCL: Use AJ buffer still in core as if it was CL, form <H|0>
!
      Call Process_RHS_Block(Inactive,Virtual,Inactive,Virtual,         &
     &                       'H ',                                      &
     &                       BRA,nBra,BRA,nBra,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read Bra (Cholesky vectors)= L(AU)
!
       Call Get_Cholesky_Vectors(Active,Virtual,JSYM,                   &
     &                           BRA,SIZE(BRA),nBra,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
! AUVX RHSC
! AUVX: Loop over the bras and kets
!
      Call Process_RHS_Block(Active,Virtual,Active,Active,              &
     &                       'C ',                                      &
     &                       BRA,nBra,KET,nKet,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! AUCX RHSF
! AUCX: Use AU buffer still in core as if it was CX, form <F|0>
!
      Call Process_RHS_Block(Active,Virtual,Active,Virtual,             &
     &                       'F ',                                      &
     &                       BRA,nBra,BRA,nBra,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read kets (Cholesky vectors) in the form L(VL), all symmetries:
!
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,                  &
     &                           KET,SIZE(KET),nKet,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
! AUVL RHSD2
! Loop over bras and kets, form <D2|0>.
!
      Call Process_RHS_Block(Active,Virtual,Inactive,Active,            &
     &                       'D2',                                      &
     &                       BRA,nBra,KET,nKet,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read kets (Cholesky vectors) in the form L(CL), all symmetries:
!
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,                 &
     &                           KET,SIZE(KET),nKet,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
! AUCL RHSG
! Loop over bras and kets, form  <G|0>
!
      Call Process_RHS_Block(Active,Virtual,Inactive,Virtual,           &
     &                       'G ',                                      &
     &                       BRA,nBra,KET,nKet,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read bra vectors AJ
!
       Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,                 &
     &                           BRA,SIZE(BRA),nBra,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
! Read kets in the form L(VL)
!
       Call Get_Cholesky_Vectors(Inactive,Active,JSYM,                  &
     &                           KET,SIZE(KET),nKet,                    &
     &                           IBSTA,IBEND)
!                                                                      *
!***********************************************************************
!                                                                      *
! AJVL RHSE
! AJVL: Loop over bras and kets. Form <E|0>
!
      Call Process_RHS_Block(Inactive,Virtual,Inactive,Active,          &
     &                       'E ',                                      &
     &                       BRA,nBra,KET,nKet,                         &
     &                       nSh,JSYM,IVEC,NV)
!                                                                      *
!***********************************************************************
!                                                                      *
! End of loop over batches, IB
      END DO
!                                                                      *
!***********************************************************************
!                                                                      *
      CALL mma_deallocate(BRA)
      CALL mma_deallocate(KET)
      CALL mma_deallocate(PIQK)
      CALL mma_deallocate(BUFF)
      CALL mma_deallocate(IDXB)
      CALL mma_deallocate(BGRP)
!                                                                      *
!***********************************************************************
!                                                                      *
! End of loop over JSYM
      END DO
!                                                                      *
!***********************************************************************
!                                                                      *
! Synchronized add RHS partial arrays from all nodes into each node.

!-SVC: read the DRA's from disk and copy them all to LUSOLV to continue
!      in serial mode.  FIXME: this call has to be removed when we reach
!      full parallel capabilities
!     CALL SYNRHS(IVEC)
!-SVC: at this point, the RHS elements are on disk, both in LUSOLV and
!      as DRAs with the name RHS_XX_XX_XX with XX a number representing
!      the case, symmetry, and rhs vector respectively.

! The RHS elements of Cases A, C, D1  need a correction:
      CALL MODRHS(IVEC,FIMO,SIZE(FIMO))

#ifdef _DEBUGPRINT_
! compute and print RHS fingerprints
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

! Synchronized add tuvx partial arrays from all nodes into each node.
      CALL CHO_GADGOP(TUVX,NTUVX,'+')
! Put TUVX on disk for possible later use:
      CALL PT2_PUT(NTUVX,'TUVX',TUVX)
      CALL mma_deallocate(TUVX)
!                                                                      *
!***********************************************************************
!                                                                      *
      END SUBROUTINE RHSALL2
