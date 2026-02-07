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
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
!     based on rhsall2.f
!
!     In principle, For E = T_{ij}^{ab}*(ia|jb).
!     With p and q for general orbitals,
!     L_{pq} = (pa|jb)*T_{qj}^{ab} + (ip|jb)*T_{ij}^{qb}
!            + (ia|pb)*T_{iq}^{ab} + (ia|jp)*T_{ij}^{aq}
!     For the first term, with P and Q for auxiliary orbitals
!     L_{pq}(1) = (pa|jb)*T_{qj}^{ab}
!               = (pa|P)*(P|jb) * T_{qj}^{ab} (2c-2e is omitted?)
!               = (pa|P) * tilde{T}_{qa}^P
!               = C_{mu p}*C_{nu a}*(mu nu|P) * tilde{T}_{qa}^P
!               = C_{mu p} * (mu nu|P) * V_{nu q}^P
!     where
!     tilde{T}_{ia}^P = T_{ij}^{ab} * (P|jb)
!     V_{mu p}^P      = tilde{T}_{pa}^P * C_{mu a}
!
!     Dimension of tilde{T} will be the same as that of (ia|P) in MO,
!     where (i,a) = (inact,active), (inact,virtual), (active,virtual)
!
!     tilde{T} is constructed in this file.
!     tilde{T} -> V_{mu p}^P transformations, contraction with 3c-2e,
!     and construction of the orbital Lagrangian (L_{pq}) is elsewhere
!     ... maybe in OLagVVVO
!
!-----------------------------------------------------------------------
!
!     For the ERI derivative calculation,
!     d(mu nu|rho sigma)/da
!     = d/da (mu nu|P) (P|Q)^-1 (Q|rho sigma)
!     = d(mu nu|P)/da (P|Q)^-1 (Q|rho sigma)
!       + (mu nu|P) d(P|Q)^-1/da (Q|rho sigma)
!       + (mu nu|P) (P|Q)^-1 d(Q|rho sigma)/da
!     = d(mu nu|P)/da (P|Q)^-1 (Q|rho sigma)
!       - (mu nu|P) (P|R)^-1 d(R|S)/da (S|Q)^-1 (Q|rho sigma)
!       + (mu nu|P) (P|Q)^-1 d(Q|rho sigma)/da
!     = d(mu nu|P)/da (tP|rho sigma)
!       - (mu nu|tP) d(P|Q)/da (tQ|rho sigma)
!       + (mu nu|tP) d(P|rho sigma)/da
!     where (mu nu|tP) = (mu nu|Q)*(Q|P)^-1
!
!     D_{mu nu rho sigma}*d(mu nu|rho sigma)/da
!     = d(mu nu|P)/da (tP|rho sigma) * D_{mu nu rho sigma}
!       - (mu nu|tP) d(P|Q)/da (tQ|rho sigma) * D_{mu nu rho sigma}
!       + (mu nu|tP) d(P|rho sigma)/da * D_{mu nu rho sigma}
!     = d(mu nu|P)/da tD_{mu nu}^tP
!       - (mu nu|tP) d(P|Q)/da tD_{mu nu}^tQ
!       + tD_{rho sigma}^tP d(P|rho sigma)/da
!     where tD_{mu nu}^tP = D_{mu nu rho sigma} * (rho sigma|tP)
!     In practice, tD_{pq}^tP is constructed and saved in disk.
!     This will be read when 3c-2e ERI derivatives are concerned,
!     and MO->AO transformations will be done on-the-fly.
!     Note that MO coefficients of CASPT2 have to be used.
!
!     For 2c-2e ERI derivatives,
!     D(tP,tQ) = tD_{pq} * C_{mu p} C_{nu q} * (mu nu|tP)
!     then saved.
!
      Subroutine OLagNS_RI(iSym0,DPT2C,DPT2Canti,A_PT2)

      Use CHOVEC_IO, only: NVLOC_CHOBATCH
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: do_csf, iStpGrd
      use PrintLevel, only: verbose
      use EQSOLV, only: IVECC2
      use ChoCASPT2, only: MaxVec_PT2
      use stdalloc, only: mma_allocate, mma_deallocate
      use definitions, only: iwp,wp,u6
      use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
      use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array
      USE Para_Info, ONLY: Is_Real_Par, King
#endif
      use caspt2_module, only: NACTEL, NSYM, NFRO, NISH, NIES, NASH,
     &                         NAES, NSSH, NSES, NORBT, NINABX,
     &                         NSECBX, NBSQT, MUL
      use caspt2_module, only: NTUV, NTU, NTGEU, NTGTU, NIGEJ, NIGTJ,
     &                         NAGEB, NAGTB, NTUVES, NTUES, NTGEUES,
     &                         NTGTUES, NIGEJES, NIGTJES, NAGEBES,
     &                         NAGTBES, NASUP, NISUP, NINDEP, NBTCH,
     &                         NBTCHES
      use Constants, only: Zero, One, Quart, Half, Two, Three, OneHalf

      implicit none

#include "warnings.h"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

      integer(kind=iwp), intent(in) :: iSym0
      real(kind=wp), intent(inout) :: DPT2C(*), DPT2Canti(*),
     &  A_PT2(MaxVec_PT2,MaxVec_PT2)

      integer(kind=iwp), parameter :: Inactive=1, Active=2, Virtual=3

      integer(kind=iwp),allocatable :: BGRP(:,:)
      real(kind=wp),allocatable :: BRA(:),KET(:),BRAD(:),KETD(:),
     &                             PIQK(:)
#ifdef _MOLCAS_MPP_
      integer(kind=iwp), allocatable :: map2(:)
      integer(kind=iwp) :: myRank, NPROCS, i, lg_V, lg_V1, ndim2
#endif

      integer(kind=iwp) :: nSh(8,3), iSym, JSYM, IB1, IB2, MXBGRP,
     &  IBGRP, IB, NBGRP, iStpGrd_sav, NCHOBUF, MXPIQK, NADDBUF, IOFFCV,
     &  IBSTA, IBEND, NV, nBra, nKet, NVI, JOFFCV, JBGRP, JBSTA, JBEND,
     &  NVJ, JB
      real(kind=wp) :: SCLNEL

      integer(kind=iwp) :: ICASE, ISYT, ISYU, ISYV, ISYW, ISYX, ISYJ,
     &  ISYL, ISYA, ISYC, ISYJL, ISYAC, ISYII, ISIJ, ISAB, ISI, ISA,
     &  ISW,
     &  NAS, NIS, NIN, NAS1, NASP, NISP, NINP, NASM, NISM, NINM, NWA,
     &  NWBP, NWBM, NWC, NWD,
     &  NWEP, NWEM, NWFP, NWFM, NWGP, NWGM, NWHP, NWHM, NW, NVEC,
     &  ipT, ipTP, ipTM, ipTanti,
     &  IT, ITABS, iTtot, IU, IUABS, iUtot, IV, IVABS, iVtot, IVMAX,
     &  IX, IXABS, iXtot, IXMAX, IJ, IJABS, iJtot, IL, ILABS,
     &  IA, IAABS, iAtot, IC, ICABS,
     &  IW1, IW2, IW, nOrbA, IO, IO1, IO2,
     &  NBXSZA, NBXSZC, NBXSZJ, NBXSZL,
     &  IASTA, IAEND, NASZ, ICSTA, ICEND, NCSZ, IJSTA, IJEND, NJSZ,
     &  ILSTA, ILEND, ILMAX, IAJSTA, IAJ, KAJ, ICLSTA, ICL, KCL,
     &  IAGEC, IAGTC, IJGEL, IJGTL, JGEL, JGTL
      real(kind=wp) :: SCL, SCL1
      real(kind=wp), parameter :: SQ2 = SQRT(Two), SQ3 = SQRT(Three),
     &  SQ05 = SQRT(Half), SQ32 = SQRT(OneHalf)

      nSh(1:nSym,Inactive) = NISH(1:nSym)
      nSh(1:nSym,Active  ) = NASH(1:nSym)
      nSh(1:nSym,Virtual ) = NSSH(1:nSym)

      IF (IPRGLB >= VERBOSE) THEN
        WRITE(u6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'
      END IF

      ! iVec = iVecX
      iSym = iSym0
      SCLNEL = One/real(MAX(1,NACTEL),kind=wp)
*                                                                      *
************************************************************************
*                                                                      *
      DO JSYM=1,NSYM
*
      IB1=NBTCHES(JSYM)+1
      IB2=NBTCHES(JSYM)+NBTCH(JSYM)
*
      MXBGRP=IB2-IB1+1
      IF (MXBGRP <= 0) CYCLE
      call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
      IBGRP=1
      DO IB=IB1,IB2
       BGRP(1,IBGRP) = IB
       BGRP(2,IBGRP) = IB
       IBGRP=IBGRP+1
      END DO
      NBGRP=MXBGRP

      !! With iStpGrd = -1, we try to allocate 4 large arrays
      iStpGrd_sav = iStpGrd
      iStpGrd = -1
      CALL MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,
     &                     NCHOBUF,MXPIQK,NADDBUF)
      iStpGrd = iStpGrd_sav
      IF (IPRGLB > VERBOSE) THEN
        WRITE(u6,*)
        WRITE(u6,'(A,I12)') '  Number of Cholesky batches: ',IB2-IB1+1
        WRITE(u6,'(A,I12)') '  Number of batch groups:     ',NBGRP
        WRITE(u6,*)
      END IF
* buffers are kept allocated until the end of JSYM loop.
      call mma_allocate(PIQK,MXPIQK,Label='PIQK')
      call mma_allocate(BRA,NCHOBUF,Label='BRABUF')
      call mma_allocate(KET,NCHOBUF,Label='KETBUF')
      call mma_allocate(BRAD,NCHOBUF,Label='BRAD')
      call mma_allocate(KETD,NCHOBUF,Label='KETD')
!
!     Loop over groups of batches of Cholesky vectors
!
      IOFFCV = 1
      DO IBGRP=1,NBGRP

      IBSTA=BGRP(1,IBGRP)
      IBEND=BGRP(2,IBGRP)

      NV=0
      DO IB=IBSTA,IBEND
        NV=NV+NVLOC_CHOBATCH(IB)
      END DO

      IF (IPRGLB > VERBOSE) THEN
        WRITE(u6,'(A,I12)') '  Cholesky vectors in this group = ', NV
        WRITE(u6,*)
      END IF

#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        myRank = GA_NodeID()
        NPROCS = GA_nNodes()

        call mma_allocate(MAP2,NPROCS,Label='MAP2')
        MAP2(:) = 0
        MAP2(myRank+1) = NV
        call GAIGOP(MAP2,NPROCS,'+')
        ndim2 = sum(map2)

        do i = nprocs, 2, -1
          map2(i) = sum(map2(1:i-1))+1
        end do
        map2(1) = 1
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Read kets (Cholesky vectors) in the form L(VX), all symmetries:
*
      Call Get_Cholesky_Vectors(Active,Active,JSYM,KET,nKet,
     &                          IBSTA,IBEND)
      KETD(1:nKet) = Zero
*                                                                      *
************************************************************************
*                                                                      *
*       Read bra (Cholesky vectors) in the form L(TJ): All symmetries
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,BRA,nBra,
     &                          IBSTA,IBEND)
      BRAD(1:nBra) = Zero
*                                                                      *
************************************************************************
*                                                                      *
*      Assemble contributions to TJVX
*      Loop over the bras and kets, form <A|0>
*
      Call OLagNS_RI2(Inactive,Active,Active,Active,
     &                'A ',BRA,KET,BRAD,KETD)
*                                                                      *
************************************************************************
*                                                                      *
*      TJVL RHSB
*      TJVL: Use TJ buffer as if it was VL, form <B|0>
*
      nKet = nBra
      Call OLagNS_RI2(Inactive,Active,Inactive,Active,
     &                'B ',BRA,BRA,BRAD,BRAD)
*                                                                      *
************************************************************************
*                                                                      *
* Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
* We still have L(VX) vectors in core, at KETS.
*
      Call Cholesky_Vectors(1,Inactive,Active,JSYM,BRAD,nBra,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,BRA,nBra,
     &                          IBSTA,IBEND)
      BRAD(1:nBra) = Zero
*                                                                      *
************************************************************************
*                                                                      *
* AJVX RHSD1
* Loop over the bra and ket vectors.
*
      Call OLagNS_RI2(Inactive,Virtual,Active,Active,
     &                'D1',BRA,KET,BRAD,KETD)
*                                                                      *
************************************************************************
*                                                                      *
* AJCL RHSH
* AJCL: Use AJ buffer still in core as if it was CL, form <H|0>
*
      nKet = nBra
      Call OLagNS_RI2(Inactive,Virtual,Inactive,Virtual,
     &                'H ',BRA,BRA,BRAD,BRAD)
*                                                                      *
************************************************************************
*                                                                      *
* Read Bra (Cholesky vectors)= L(AU)
*
      Call Cholesky_Vectors(1,Inactive,Virtual,JSYM,BRAD,nBra,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Active,Virtual,JSYM,BRA,nBra,
     &                          IBSTA,IBEND)
      BRAD(1:nBra) = Zero
!                                                                      *
************************************************************************
*                                                                      *
* AUVX RHSC
* AUVX: Loop over the bras and kets
*
      Call OLagNS_RI2(Active,Virtual,Active,Active,
     &                'C ',BRA,KET,BRAD,KETD)
*                                                                      *
************************************************************************
*                                                                      *
* AUCX RHSF
* AUCX: Use AU buffer still in core as if it was CX, form <F|0>
*
      nKet = nBra
      Call OLagNS_RI2(Active,Virtual,Active,Virtual,
     &                'F ',BRA,BRA,BRAD,BRAD)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(VL), all symmetries:
*
      Call Cholesky_Vectors(1,Active,Active,JSYM,KETD,nKet,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,KET,nKet,
     &                          IBSTA,IBEND)
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,KETD,nKet,
     &                      IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUVL RHSD2
* Loop over bras and kets, form <D2|0>.
*
      Call OLagNS_RI2(Active,Virtual,Inactive,Active,
     &                'D2',BRA,KET,BRAD,KETD)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets (Cholesky vectors) in the form L(CL), all symmetries:
*
      Call Cholesky_Vectors(1,Inactive,Active,JSYM,KETD,nKet,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,KET,nKet,
     &                          IBSTA,IBEND)
      Call Cholesky_Vectors(2,Inactive,Virtual,JSYM,KETD,nKet,
     &                      IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AUCL RHSG
* Loop over bras and kets, form  <G|0>
*
      Call OLagNS_RI2(Active,Virtual,Inactive,Virtual,
     &                'G ',BRA,KET,BRAD,KETD)
*                                                                      *
************************************************************************
*                                                                      *
* Read bra vectors AJ
*
      Call Cholesky_Vectors(1,Active,Virtual,JSYM,BRAD,nBra,
     &                      IBSTA,IBEND)
      Call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,BRA,nBra,
     &                          IBSTA,IBEND)
      BRAD(1:nBra) = KETD(1:nBra)
*                                                                      *
************************************************************************
*                                                                      *
* Read kets in the form L(VL)
*
      Call Get_Cholesky_Vectors(Inactive,Active,JSYM,
     &                          KET,nKet,
     &                          IBSTA,IBEND)
      Call Cholesky_Vectors(2,Inactive,Active,JSYM,KETD,nKet,
     &                      IBSTA,IBEND)
*                                                                      *
************************************************************************
*                                                                      *
* AJVL RHSE
* AJVL: Loop over bras and kets. Form <E|0>
*
      Call OLagNS_RI2(Inactive,Virtual,Inactive,Active,
     &                'E ',BRA,KET,BRAD,KETD)
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over batches, IB
      Call Cholesky_Vectors(1,Inactive,Virtual,JSYM,BRAD,nBra,
     &                      IBSTA,IBEND)
      Call Cholesky_Vectors(1,Inactive,Active,JSYM,KETD,nKet,
     &                      IBSTA,IBEND)
      !! Construct A_PT2
      NVI = NV
      JOFFCV = 1
      DO JBGRP=1,NBGRP

        JBSTA=BGRP(1,JBGRP)
        JBEND=BGRP(2,JBGRP)

        NVJ=0
        DO JB=JBSTA,JBEND
          NVJ=NVJ+NVLOC_CHOBATCH(JB)
        END DO

        !! BraAI
        Call Cnst_A_PT2(Inactive,Active)

        !! BraSI
        Call Cnst_A_PT2(Inactive,Virtual)

        !! BraSA
        Call Cnst_A_PT2(Active,Virtual)

        !! BraAA
        Call Cnst_A_PT2(Active,Active)
        JOFFCV = JOFFCV + NVJ
      END DO !! end of JBGRP loop
      IOFFCV = IOFFCV + NVI
      END DO !! end of IBGRP loop
*                                                                      *
************************************************************************
*                                                                      *
      call mma_deallocate(BRA)
      call mma_deallocate(KET)
      call mma_deallocate(BRAD)
      call mma_deallocate(KETD)
      call mma_deallocate(PIQK)
      call mma_deallocate(BGRP)
*                                                                      *
************************************************************************
*                                                                      *
* End of loop over JSYM
      END DO
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) call mma_deallocate(map2)
#endif
*                                                                      *
************************************************************************
*                                                                      *
* Synchronized add RHS partial arrays from all nodes into each node.

!-SVC: read the DRA's from disk and copy them all to LUSOLV to continue
!      in serial mode.  FIXME: this call has to be removed when we reach
!      full parallel capabilities
*     CALL SYNRHS(IVEC)
!-SVC: at this point, the RHS elements are on disk, both in LUSOLV and
!      as DRAs with the name RHS_XX_XX_XX with XX a number representing
!      the case, symmetry, and rhs vector respectively.
!

      A_PT2(:,:) = Two*A_PT2(:,:)

      If (NBGRP /= 0) SCLNEL = SCLNEL/real(NBGRP,kind=wp)
      DPT2C(1:NBSQT) = DPT2C(1:NBSQT)*SCLNEL
      If (do_csf) DPT2Canti(1:NBSQT) = DPT2Canti(1:NBSQT)*SCLNEL

#ifdef _MOLCAS_MPP_
      If (is_real_par()) then
        CALL GADSUM (A_PT2,MaxVec_PT2**2)
      end if
#endif

      Return

      Contains
!
!-----------------------------------------------------------------------
!
      Subroutine OLagNS_RI2(ITI,ITP,ITK,ITQ,Case,Cho_Bra,Cho_Ket,
     &                      Cho_BraD,Cho_KetD)

      use caspt2_global, only:iPrGlb
      use PrintLevel, only: debug

      implicit none

      integer(kind=iwp), intent(in) :: ITI, ITP, ITK, ITQ
      character(len=2), intent(in) :: Case
      real(kind=wp), intent(in) :: Cho_Bra(*), Cho_Ket(*)
      real(kind=wp), intent(inout) :: Cho_BraD(*), Cho_KetD(*)

      integer(kind=iwp) :: LBRASM, ISYI, NI, ISYP, NP, NPI, NBRASM,
     &  LKETSM, ISYK, NK, ISYQ, NQ, NQK, NKETSM, NPIQK, KPI, KQK

      real(kind=wp) :: TotCPU0, TotWall0, TotCPU1, TotWall1

      IF (iPrGlb >= DEBUG) THEN
        WRITE(u6,*) 'Processing RHS block '//Case
      END IF

      LBRASM=1
      CALL CWTime(TotCPU0,TotWall0)
      DO ISYI=1,NSYM
        NI=NSH(ISYI,ITI)
        IF(NI == 0) CYCLE
        ISYP=MUL(ISYI,JSYM)
        NP=NSH(ISYP,ITP)
        IF(NP == 0) CYCLE
        NPI=NP*NI
        NBRASM=NPI*NV

        LKETSM=1
        DO ISYK=1,NSYM
          NK=NSH(ISYK,ITK)
          IF(NK == 0) CYCLE
          ISYQ=MUL(ISYK,JSYM)
          NQ=NSH(ISYQ,ITQ)
          IF(NQ == 0) CYCLE
          NQK=NQ*NK
          NKETSM=NQK*NV
*
! SVC: we need an NPI*NQK to store the 2-electron integrals, and 2
! buffers (values+indices) for sorting them.  Later, we can try to get
! rid of the buffer that stores the values and only use an index buffer
! and the two-electron integrals for the scatter operation.  For the
! buffer, any size can be taken, but assuming there is enough memory
! available, it's set to the size of the two-electron integrals unless
! larger than some predefined maximum buffer size.
          NPIQK=NPI*NQK
          IF (NPIQK > MXPIQK) THEN
            IF (Case == 'H') THEN
              KPI=MXPIQK/NQK
              NPIQK=KPI*NQK
            ELSE IF (Case == 'G') THEN
              KQK=MXPIQK/NPI
              NPIQK=NPI*KQK
            ELSE
              WRITE(u6,*) ' NPIQK > MXPIQK and case != G or H'
              WRITE(u6,'(A,A2)')  ' CASE =   ', Case
              WRITE(u6,'(A,I12)') ' NPIQK =  ', NPIQK
              WRITE(u6,'(A,I12)') ' MXPIQK = ', MXPIQK
              WRITE(u6,*) ' This should not happen, please report.'
              CALL AbEnd()
            END IF
          END IF

          IF (NPIQK <= 0) THEN
            WRITE(u6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
            CALL AbEnd()
          END IF
*
          !! NBUFF(=nAddBuf) is removed
          If (Case == 'A ') Then
             CALL OLagNS_RI_A(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'B ') Then
             CALL OLagNS_RI_B(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'D1') Then
             CALL OLagNS_RI_D1(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'H ') Then
             CALL OLagNS_RI_H(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'C ') Then
             CALL OLagNS_RI_C(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'F ') Then
             CALL OLagNS_RI_F(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'D2') Then
             CALL OLagNS_RI_D2(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'G ') Then
             CALL OLagNS_RI_G(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else If (Case == 'E ') Then
             CALL OLagNS_RI_E(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,
     &                        Cho_Bra(LBRASM),Cho_Ket(LKETSM),
     &                        Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
          Else
             Call Abend()
          End If

          LKETSM=LKETSM+NKETSM
        END DO
        LBRASM=LBRASM+NBRASM
      END DO
      CALL CWTime(TotCPU1,TotWall1)
      IF (IPRGLB >= VERBOSE) THEN
        write(u6,'(" CPU/Wall Time (Case ",A2,"):",2f10.2)')
     &    Case,totcpu1-totcpu0,totwall1-totwall0
      END IF

      Return

      End Subroutine OLagNS_RI2
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_A(ISYI,ISYK,NT,NJ,NV,NX,TJVX,NTJVX,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KTUV
      use caspt2_global, only: do_csf

      implicit none

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NT, NJ, NV, NX,
     &  NTJVX, NCHO
      real(kind=wp), intent(out) :: TJVX(NT,NJ,NV,NX)
      real(kind=wp), intent(in) :: Cho_Bra(NT,NJ,NCHO),
     &  Cho_Ket(NV,NX,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NT,NJ,NCHO),
     &  Cho_KetD(NV,NX,NCHO)

      integer(kind=iwp) :: IT, IJ, IV, IX

      ISYJ = ISYI
      ISYX = ISYK

      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYJ
      IF(NINDEP(ISYM,1) == 0) RETURN
      NAS=NTUV(ISYM)
      NIS=NISH(ISYM)
      NWA=NAS*NIS
      IF(NWA == 0) RETURN
!
!     ---- A
!
      !! Read the T-amplitude
      ICASE=1
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN /= 0) Then
        nIS = nISup(iSym,iCase)
        nVec = nAS*nIS
        If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            ! copy global array to local buffer
            Call RHS_ALLO(nAS,nIS,lg_V)
            CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
            ipT = Allocate_GA_Array(nAS*nIS,'ipT')
            CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A(1),nAS)
            If (do_csf) Then
              CALL RHS_READ_C(lg_V,iCase,iSym,7)
              ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
              CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A(1),nAS)
            End If
            CALL RHS_FREE(lg_V)
            CALL GASYNC()
          ELSE
#endif
            Call RHS_ALLO(nAS,nIS,ipT)
            CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
            If (do_csf) Then
              Call RHS_ALLO(nAS,nIS,ipTanti)
              CALL RHS_READ_C(ipTanti,iCase,iSym,7)
            End If
#ifdef _MOLCAS_MPP_
          END IF
#endif
        End If
      End If

      Call DCopy_(NTJVX,[Zero],0,TJVX,1)

      nOrbT = nFro(iSyT)+nIsh(iSyT)+nAsh(iSyT)+nSsh(iSyT)
      DO IT=1,NT
        ITABS=IT+NAES(ISYT)
        iTtot = iT + nFro(iSyT) + nIsh(iSyT)
        DO IJ=1,NJ
          ! IJABS=IJ+NIES(ISYJ)
          iJtot = iJ + nFro(iSyJ)

          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            IF (ISYV == ISYX) THEN !! not sure
              !! ONEADD contributions
              IW1=KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
              IW2=IJ
              IW=IW1+NAS*(IW2-1)

              DPT2C(iTtot+nOrbT*(iJtot-1))
     &          = DPT2C(iTtot+nOrbT*(iJtot-1))
     &          + Two*GA_Arrays(ipT)%A(IW)
              If (do_csf) Then
                DPT2Canti(iTtot+nOrbT*(iJtot-1))
     &            = DPT2Canti(iTtot+nOrbT*(iJtot-1))
     &            + Two*GA_Arrays(ipTanti)%A(IW)
              End If
            END IF
            DO IX=1,NX
              IXABS=IX+NAES(ISYX)
              IW1=KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
              IW2=IJ
              IW=IW1+NAS*(IW2-1)
              TJVX(IT,IJ,IV,IX) = GA_Arrays(ipT)%A(IW)
            END DO
          END DO
        END DO
      END DO

      Call DGEMM_('T','N',NV*NX,NCHO,NT*NJ,
     &            One,TJVX(1,1,1,1),NT*NJ,Cho_Bra(1,1,1),NT*NJ,
     &            One,Cho_KetD(1,1,1),NV*NX)
      Call DGEMM_('N','N',NT*NJ,NCHO,NV*NX,
     &            One,TJVX(1,1,1,1),NT*NJ,Cho_Ket(1,1,1),NV*NX,
     &            One,Cho_BraD(1,1,1),NT*NJ)

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipT)
          if (do_CSF) call deallocate_GA_array(ipTanti)
        ELSE
#endif
          CALL RHS_FREE(ipT)
          If (do_csf) CALL RHS_FREE(ipTanti)

#ifdef _MOLCAS_MPP_
        END IF
#endif

      RETURN

      End Subroutine OLagNS_RI_A
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_B(ISYI,ISYK,NT,NJ,NV,NL,TJVL,NTJVL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KIGTJ, KIGEJ, KTGTU, KTGEU

      implicit none

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NT, NJ, NV, NL,
     &  NTJVL, NCHO
      real(kind=wp), intent(out) :: TJVL(NT,NJ,NV,NL)
      real(kind=wp), intent(in) :: Cho_Bra(NT,NJ,NCHO),
     &  Cho_Ket(NV,NL,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NT,NJ,NCHO),
     &  Cho_KetD(NV,NL,NCHO)

      integer(kind=iwp) :: IT, IV, IJ, IL

      ISYJ = ISYI
      ISYL = ISYK

      ISYT=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      IF(ISYT < ISYV) RETURN
      ISYM=MUL(ISYJ,ISYL) !!

      IF(NINDEP(ISYM,2) > 0) THEN
* The plus combination:
       ICASE=2
       NASP=NTGEU(ISYM)
       NISP=NIGEJ(ISYM)
       NWBP=NASP*NISP
      ELSE
       NWBP=0
      ENDIF
      IF(NINDEP(ISYM,3) > 0) THEN
* The minus combination:
       ICASE=3
       NASM=NTGTU(ISYM)
       NISM=NIGTJ(ISYM)
       NWBM=NASM*NISM
      ELSE
       NWBM=0
      ENDIF
      If (Max(NWBP,NWBM) <= 0) RETURN

      Call DCopy_(NTJVL,[Zero],0,TJVL,1)

      IF(NWBP > 0.AND.NINDEP(ISYM,2) > 0) THEN
        !! Read the T-amplitude
        ICASE=2
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP /= 0) Then
          nISP = nISup(iSym,iCase)
          nVec = nASP*nISP
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASP,nISP,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
              CALL GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A(1),nASP)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASP,nISP,ipTP)
              CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            END IF
#endif
          End If
        End If

        DO IT=1,NT
          ITABS=IT+NAES(ISYT)
          IVMAX=NV
          IF(ISYV == ISYT) IVMAX=IT
          DO IV=1,IVMAX
            IVABS=IV+NAES(ISYV)
            SCL1=Half
            IW1=KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
            IF(ITABS == IVABS) SCL1=Quart
            DO IJ=1,NJ
              IJABS=IJ+NIES(ISYJ)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                SCL=SCL1
                IF(IJABS >= ILABS) THEN
                 IW2=KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
                 IF(IJABS == ILABS) SCL=SQ2*SCL1
                ELSE
                 IW2=KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
                END IF
                IW=IW1+NASP*(IW2-1)
                TJVL(IT,IJ,IV,IL) = SCL*GA_Arrays(ipTP)%A(IW)
              END DO
            END DO
          END DO
        END DO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTP)
        else
#endif
          CALL RHS_FREE(ipTP)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF

      IF(NINDEP(ISYM,3) > 0) THEN
        !! Read the T-amplitude
        ICASE=3
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM /= 0) Then
          nISM = nISup(iSym,iCase)
          nVec = nASM*nISM
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASM,nISM,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
              CALL GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A(1),nASM)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASM,nISM,ipTM)
              CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            End If
#endif
          End If
        End If

        DO IT=1,NT
          ITABS=IT+NAES(ISYT)
          IVMAX=NV
          IF(ISYV == ISYT) IVMAX=IT-1
          DO IV=1,IVMAX
            IVABS=IV+NAES(ISYV)
            IW1=KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
            DO IJ=1,NJ
              IJABS=IJ+NIES(ISYJ)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                IF(IJABS > ILABS) THEN
                  IW2=KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
                  SCL =  Half
                ELSE IF (IJABS < ILABS) THEN
                  IW2=KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
                  SCL = -Half
                ELSE
                  CYCLE
                END IF
                IW=IW1+NASM*(IW2-1)
                TJVL(IT,IJ,IV,IL) = TJVL(IT,IJ,IV,IL)
     &            + SCL*GA_Arrays(ipTM)%A(IW)
              END DO
            END DO
          END DO
        END DO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTM)
        else
#endif
          CALL RHS_FREE(ipTM)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF

      Call DGEMM_('T','N',NV*NL,NCHO,NT*NJ,
     &            One,TJVL(1,1,1,1),NT*NJ,Cho_Bra(1,1,1),NT*NJ,
     &            One,Cho_KetD(1,1,1),NV*NL)
      Call DGEMM_('N','N',NT*NJ,NCHO,NV*NL,
     &            One,TJVL(1,1,1,1),NT*NJ,Cho_Ket(1,1,1),NV*NL,
     &            One,Cho_BraD(1,1,1),NT*NJ)

      RETURN

      End Subroutine OLagNS_RI_B
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_C(ISYI,ISYK,NA,NU,NV,NX,AUVX,NAUVX,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX
      use caspt2_global, only: do_csf

      implicit none

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NV, NX,
     &  NAUVX, NCHO
      real(kind=wp), intent(out) :: AUVX(NA,NU,NV,NX)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO),
     &  Cho_Ket(NV,NX,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO),
     &  Cho_KetD(NV,NX,NCHO)

      real(kind=wp) :: ValCF
      integer(kind=iwp) :: IA, IU, IV, IX

      ISYU = ISYI
      ISYX = ISYK

      ISYA=MUL(JSYM,ISYU)
      ISYV=MUL(JSYM,ISYX)
      ISYM=ISYA !!
      IF(NINDEP(ISYM,4) == 0) RETURN
      NAS=NTUV(ISYM)
      NIS=NSSH(ISYM)
      NWC=NAS*NIS
      IF(NWC == 0) RETURN
!
!     ---- C
!
      !! Read the T-amplitude
      ICASE=4
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN /= 0) Then
        nIS = nISup(iSym,iCase)
        nVec = nIN*nIS
        If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            Call RHS_ALLO(nAS,nIS,lg_V)
            CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
            ipT = Allocate_GA_Array(nAS*nIS,'ipT')
            CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A(1),nAS)
            If (do_csf) Then
              CALL RHS_READ_C(lg_V,iCase,iSym,7)
              ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
              CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A(1),nAS)
            End If
            Call RHS_FREE(lg_V)
            CALL GASYNC()
          ELSE
#endif
            Call RHS_ALLO(nAS,nIS,ipT)
            CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
            If (do_csf) Then
              Call RHS_ALLO(nAS,nIS,ipTanti)
              CALL RHS_READ_C(ipTanti,iCase,iSym,7)
            End If
#ifdef _MOLCAS_MPP_
          END IF
#endif
        End If
      End If

      Call DCopy_(NAUVX,[Zero],0,AUVX,1)

      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      DO IA=1,NA
        iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          iUtot = iU + nFro(iSyU) + nIsh(iSyU)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            ValCF = Zero
            IF (ISYV == ISYX) THEN !! not sure
              !! ONEADD contributions
              IW1=KTUV(IUABS,IVABS,IVABS)-NTUVES(ISYM)
              IW2=IA
              IW=IW1+NAS*(IW2-1)

              ValCF = GA_Arrays(ipT)%A(IW)*Two
              DPT2C(iAtot+nOrbA*(iUtot-1))
     &          = DPT2C(iAtot+nOrbA*(iUtot-1)) + ValCF
              If (do_csf) Then
                DPT2Canti(iAtot+nOrbA*(iUtot-1))
     &            = DPT2Canti(iAtot+nOrbA*(iUtot-1))
     &            + Two*GA_Arrays(ipTanti)%A(IW)
              End If
              ValCF = ValCF*SCLNEL
            END IF
            DO IX=1,NX
              IXABS=IX+NAES(ISYX)
              IW1=KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
              IW2=IA
              IW=IW1+NAS*(IW2-1)

              AUVX(IA,IU,IV,IX) = AUVX(IA,IU,IV,IX)
     &          + GA_Arrays(ipT)%A(IW)
              AUVX(IA,IX,IU,IX) = AUVX(IA,IX,IU,IX) - ValCF*Half
            END DO
          END DO
        END DO
      END DO

      Call DGEMM_('T','N',NV*NX,NCHO,NA*NU,
     &            One,AUVX(1,1,1,1),NA*NU,Cho_Bra(1,1,1),NA*NU,
     &            One,Cho_KetD(1,1,1),NV*NX)
      Call DGEMM_('N','N',NA*NU,NCHO,NV*NX,
     &            One,AUVX(1,1,1,1),NA*NU,Cho_Ket(1,1,1),NV*NX,
     &            One,Cho_BraD(1,1,1),NA*NU)

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipT)
          if (do_CSF) call deallocate_GA_array(ipTanti)
        ELSE
#endif
          CALL RHS_FREE(ipT)
          If (do_csf) CALL RHS_FREE(ipTanti)
#ifdef _MOLCAS_MPP_
        END IF
#endif

      RETURN

      End Subroutine OLagNS_RI_C
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_D1(ISYI,ISYK,NA,NJ,NV,NX,AJVX,NAJVX,
     &                        Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KTU
      use caspt2_global, only: do_csf

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NJ, NV, NX,
     &  NAJVX, NCHO
      real(kind=wp), intent(_OUT_) :: AJVX(NV,NX,*)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NJ,NCHO),
     &  Cho_Ket(NV,NX,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NJ,NCHO),
     &  Cho_KetD(NV,NX,NCHO)
*      Logical Incore
      integer(kind=iwp) :: IOFFD(8,8)
      integer(kind=iwp) :: ISW, ISA, IASTA, IJSTA, IJ, IA, IX, IV

      ISYJ = ISYI
      ISYX = ISYK

      DO ISW=1,NSYM
       IO=0
       DO ISA=1,NSYM
        IOFFD(ISA,ISW)=IO
        ISI=MUL(ISA,ISW)
        IO=IO+NSSH(ISA)*NISH(ISI)
       END DO
      END DO

      ISYA=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYX)
      ISYM=JSYM !!!
      IF(NINDEP(ISYM,5) == 0) RETURN
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD == 0) RETURN
!
!     ---- D1
!
      !! Read the T-amplitude
      ICASE=5
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN /= 0) Then
        nIS = nISup(iSym,iCase)
        nVec = nAS*nIS
        If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            Call RHS_ALLO(nAS,nIS,lg_V)
            CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
            ipT = Allocate_GA_Array(nAS*nIS,'ipT')
            CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A(1),nAS)
            If (do_csf) Then
              CALL RHS_READ_C(lg_V,iCase,iSym,7)
              ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
              CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A(1),nAS)
            End If
            Call RHS_FREE(lg_V)
            CALL GASYNC()
          ELSE
#endif
            Call RHS_ALLO(nAS,nIS,ipT)
            CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
            If (do_csf) Then
              Call RHS_ALLO(nAS,nIS,ipTanti)
              CALL RHS_READ_C(ipTanti,iCase,iSym,7)
            End If
#ifdef _MOLCAS_MPP_
          END IF
#endif
        End If
      End If

      NBXSZA=NSECBX
      NBXSZJ=NINABX

      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        NASZ=IAEND-IASTA+1
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          NJSZ=IJEND-IJSTA+1

      IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      Call DCopy_(NAJVX,[Zero],0,AJVX,1)

      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      IAJ=0
      DO IJ=IJSTA,IJEND
        iJtot = iJ + nFro(iSyJ)
        DO IA=IASTA,IAEND
          iAtot = iA + nFro(iSyA) + nIsh(iSyA) + nAsh(iSyA)
          IAJ = IAJ + 1
          DO IX=1,NX
            IXABS=IX+NAES(ISYX)
            iXtot = iX + nFro(iSyX) + nIsh(iSyX)
            DO IV=1,NV
              IVABS=IV+NAES(ISYV)
              iVtot = iV + nFro(iSyV) + nIsh(iSyV)
              IW1=KTU(IVABS,IXABS)-NTUES(ISYM)
              IW2=IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
              IW=IW1+NAS*(IW2-1)

              If (iVtot == iXtot) Then
                DPT2C(iAtot+nOrbA*(iJtot-1))
     &            = DPT2C(iAtot+nOrbA*(iJtot-1))
     &            + Two*GA_Arrays(ipT)%A(IW)
                If (do_csf) Then
                  DPT2Canti(iAtot+nOrbA*(iJtot-1))
     &              = DPT2Canti(iAtot+nOrbA*(iJtot-1))
     &              + Two*GA_Arrays(ipTanti)%A(IW)
                End If
              End If
              AJVX(IV,IX,IAJ) = GA_Arrays(ipT)%A(IW)
            END DO
          END DO
        END DO
      END DO

      Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NX,
     &            One,AJVX(1,1,1),NV*NX,Cho_Ket(1,1,1),NV*NX,
     &            One,Cho_BraD(IAJSTA,1,1),NA*NJ)
      Call DGEMM_('N','N',NV*NX,NCHO,NASZ*NJSZ,
     &            One,AJVX(1,1,1),NV*NX,Cho_Bra(IAJSTA,1,1),NA*NJ,
     &            One,Cho_KetD(1,1,1),NV*NX)

        ENDDO
      ENDDO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipT)
          if (do_CSF) call deallocate_GA_array(ipTanti)
        ELSE
#endif
          CALL RHS_FREE(ipT)
          If (do_csf) CALL RHS_FREE(ipTanti)
#ifdef _MOLCAS_MPP_
        END IF
#endif

      RETURN

      End Subroutine OLagNS_RI_D1
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_D2(ISYI,ISYK,NA,NU,NV,NL,AUVL,NAUVL,
     &                        Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KTU

      implicit none

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NV, NL,
     &  NAUVL, NCHO
      real(kind=wp), intent(out) :: AUVL(NA,NU,NV,NL)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO),
     &  Cho_Ket(NV,NL,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO),
     &  Cho_KetD(NV,NL,NCHO)
*      Logical Incore
      integer(kind=iwp) :: IOFFD(8,8)
      integer(kind=iwp) :: ISYW, ISYA, IA, IU, IV, IL

      ISYU = ISYI
      ISYL = ISYK

      DO ISYW=1,NSYM
       IO=0
       DO ISYA=1,NSYM
        IOFFD(ISYA,ISYW)=IO
        ISYII=MUL(ISYA,ISYW)
        IO=IO+NSSH(ISYA)*NISH(ISYII)
       END DO
      END DO

      ISYA=MUL(JSYM,ISYU)
      ISYV=MUL(JSYM,ISYL)
      ISYM=MUL(ISYU,ISYV)
      IF(NINDEP(ISYM,5) == 0) RETURN
      NAS1=NTU(ISYM)
      NAS=2*NAS1
      NIS=NISUP(ISYM,5)
      NWD=NAS*NIS
      IF(NWD == 0) RETURN
!
!     ---- D2
!
      !! Read the T-amplitude
      ICASE=5
      nIN = nINDEP(iSym,iCase)
      nAS = nASup(iSym,iCase)
      If (nIN /= 0) Then
        nIS = nISup(iSym,iCase)
        nVec = nAS*nIS
        If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            Call RHS_ALLO(nAS,nIS,lg_V)
            CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
            ipT = Allocate_GA_Array(nAS*nIS,'ipT')
            CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A(1),nAS)
            If (do_csf) Then
              CALL RHS_READ_C(lg_V,iCase,iSym,7)
              ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
              CALL GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A(1),nAS)
            End If
            Call RHS_FREE(lg_V)
            CALL GASYNC()
          ELSE
#endif
            Call RHS_ALLO(nAS,nIS,ipT)
            CALL RHS_READ_C(ipT,iCase,iSym,iVecC2)
            If (do_csf) Then
              Call RHS_ALLO(nAS,nIS,ipTanti)
              CALL RHS_READ_C(ipTanti,iCase,iSym,7)
            End If
#ifdef _MOLCAS_MPP_
          END IF
#endif
        End If
      End If

      Call DCopy_(NAUVL,[Zero],0,AUVL,1)

      DO IA=1,NA
        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          DO IV=1,NV
            IVABS=IV+NAES(ISYV)
            DO IL=1,NL
              IW1=NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
              IW2=IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
              IW=IW1+NAS*(IW2-1)
              AUVL(IA,IU,IV,IL) = GA_Arrays(ipT)%A(IW)
            END DO
          END DO
        END DO
      END DO

      Call DGEMM_('T','N',NV*NL,NCHO,NA*NU,
     &            One,AUVL,NA*NU,Cho_Bra,NA*NU,
     &            One,Cho_KetD,NV*NL)
      Call DGEMM_('N','N',NA*NU,NCHO,NV*NL,
     &            One,AUVL,NA*NU,Cho_Ket,NV*NL,
     &            One,Cho_BraD,NA*NU)

      if (nIN /= 0 .and. nVec /= 0) then
#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipT)
          if (do_CSF) call deallocate_GA_array(ipTanti)
        else
#endif
          CALL RHS_FREE(ipT)
          If (do_csf) CALL RHS_FREE(ipTanti)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF

      RETURN

      End Subroutine OLagNS_RI_D2
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_E(ISYI,ISYK,NA,NJ,NV,NL,AJVL,NAJVL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KIGTJ, KIGEJ

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NJ, NV, NL,
     &  NAJVL, NCHO
      real(kind=wp), intent(_OUT_) :: AJVL(NV,NL,*)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NJ,NCHO),
     &  Cho_Ket(NV,NL,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NJ,NCHO),
     &  Cho_KetD(NV,NL,NCHO)
*      Logical Incore
      integer(kind=iwp) :: IOFF1(8), IOFF2(8)
      integer(kind=iwp) :: ISA, IASTA, IJSTA, IJ, IA, IV, IL

      ISYJ = ISYI
      ISYL = ISYK

      ISYA=MUL(JSYM,ISYJ)
      ISYV=MUL(JSYM,ISYL)
      ISYM=ISYV
      ISYJL=MUL(ISYJ,ISYL)

! Set up offset table:
      IO1=0
      IO2=0
      DO ISA=1,NSYM
        IOFF1(ISA)=IO1
        IOFF2(ISA)=IO2
        ISIJ=MUL(ISA,ISYM)
        IO1=IO1+NSSH(ISA)*NIGEJ(ISIJ)
        IO2=IO2+NSSH(ISA)*NIGTJ(ISIJ)
      END DO

      NAS=NASH(ISYM)
      NISP=NISUP(ISYM,6)
      NISM=NISUP(ISYM,7)
      ! NIS=NISP+NISM
      NWEP=NAS*NISP
      NWEM=NAS*NISM
      NW=NWEP+NWEM
      If (NW == 0) RETURN
!
!     ---- EP
!
      IF (NWEP > 0) THEN
        !! Read the T-amplitude
        ICASE=6
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP /= 0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASP,nISP,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
              CALL GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A(1),nASP)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASP,nISP,ipTP)
              CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            END IF
#endif
          End If
        End If

        NBXSZA=NSECBX
        NBXSZJ=NINABX

        DO IASTA=1,NA,NBXSZA
          IAEND=MIN(IASTA-1+NBXSZA,NA)
          NASZ=IAEND-IASTA+1
          DO IJSTA=1,NJ,NBXSZJ
            IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
            NJSZ=IJEND-IJSTA+1

        IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        Call DCopy_(NAJVL,[Zero],0,AJVL,1)

        IAJ=0
        DO IJ=IJSTA,IJEND
          IJABS=IJ+NIES(ISYJ)
          DO IA=IASTA,IAEND
            ! IAABS=IA+NSES(ISYA)
            IAJ=IAJ+1

            DO IV=1,NV
              ! IVABS=IV+NAES(ISYV)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                SCL=SQ05
                IF(IJABS >= ILABS) THEN
                  JGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                  IF(IJABS == ILABS) SCL=One
                ELSE
                  JGEL=KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
                END IF
                IW1=IV
                IW2=IA+NA*(JGEL-1)+IOFF1(ISYA)
                IW=IW1+NAS*(IW2-1)

                AJVL(IV,IL,IAJ) = SCL*GA_Arrays(ipTP)%A(IW)
              END DO
            END DO
          END DO
        END DO

        Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NL,
     &              One,AJVL(1,1,1),NV*NL,Cho_Ket(1,1,1),NV*NL,
     &              One,Cho_BraD(IAJSTA,1,1),NA*NJ)
        Call DGEMM_('N','N',NV*NL,NCHO,NASZ*NJSZ,
     &              One,AJVL(1,1,1),NV*NL,Cho_Bra(IAJSTA,1,1),NA*NJ,
     &              One,Cho_KetD(1,1,1),NV*NL)

          ENDDO
        ENDDO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTP)
        else
#endif
          CALL RHS_FREE(ipTP)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF
!
!     ---- EM
!
      IF (NWEM > 0) THEN
        !! Read the T-amplitude
        ICASE=7
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM /= 0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASM,nISM,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
              CALL GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A(1),nASM)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASM,nISM,ipTM)
              CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            End If
#endif
          End If
        End If

        NBXSZA=NSECBX
        NBXSZJ=NINABX

        DO IASTA=1,NA,NBXSZA
          IAEND=MIN(IASTA-1+NBXSZA,NA)
          NASZ=IAEND-IASTA+1
          DO IJSTA=1,NJ,NBXSZJ
            IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
            NJSZ=IJEND-IJSTA+1

        IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        Call DCopy_(NAJVL,[Zero],0,AJVL,1)

        IAJ=0
        DO IJ=IJSTA,IJEND
          IJABS=IJ+NIES(ISYJ)
          DO IA=IASTA,IAEND
            ! IAABS=IA+NSES(ISYA)
            IAJ=IAJ+1

            DO IV=1,NV
              ! IVABS=IV+NAES(ISYV)
              DO IL=1,NL
                ILABS=IL+NIES(ISYL)
                IF(IJABS /= ILABS) THEN
                  IF(IJABS > ILABS) THEN
                    SCL=SQ32
                    JGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                  ELSE
                    SCL=-SQ32
                    JGTL=KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
                  END IF
                  IW1=IV
                  IW2=IA+NA*(JGTL-1)+IOFF2(ISYA)
                  IW=IW1+NAS*(IW2-1)

                  AJVL(IV,IL,IAJ) = SCL*GA_Arrays(ipTM)%A(IW)
                END IF
              END DO
            END DO
          END DO
        END DO

        Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NL,
     &              One,AJVL(1,1,1),NV*NL,Cho_Ket(1,1,1),NV*NL,
     &              One,Cho_BraD(IAJSTA,1,1),NA*NJ)
        Call DGEMM_('N','N',NV*NL,NCHO,NASZ*NJSZ,
     &              One,AJVL(1,1,1),NV*NL,Cho_Bra(IAJSTA,1,1),NA*NJ,
     &              One,Cho_KetD(1,1,1),NV*NL)

          ENDDO
        ENDDO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTM)
        else
#endif
          CALL RHS_FREE(ipTM)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF

      RETURN

      End Subroutine OLagNS_RI_E
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_F(ISYI,ISYK,NA,NU,NC,NX,AUCX,NAUCX,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KTGTU, KTGEU, KAGTB, KAGEB

      implicit none

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NC, NX,
     &  NAUCX, NCHO
      real(kind=wp), intent(out) :: AUCX(NA,NU,NC,NX)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO),
     &  Cho_Ket(NC,NX,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO),
     &  Cho_KetD(NC,NX,NCHO)

      integer(kind=iwp) :: IU, IX, IA, IC

      ISYU = ISYI
      ISYX = ISYK

      IF(ISYU < ISYX) RETURN

      ISYA=MUL(JSYM,ISYU)
      ISYC=MUL(JSYM,ISYX)
      ISYM=MUL(ISYU,ISYX) !!

      IF(NINDEP(ISYM,8) > 0) THEN
* The plus combination:
       NASP=NTGEU(ISYM)
       NISP=NAGEB(ISYM)
       NWFP=NASP*NISP
      ELSE
       NWFP=0
      ENDIF
      IF(NINDEP(ISYM,9) > 0) THEN
       ICASE=9
* The minus combination:
       NASM=NTGTU(ISYM)
       NISM=NAGTB(ISYM)
       NWFM=NASM*NISM
      ELSE
       NWFM=0
      ENDIF
      If (NWFP+NWFM <= 0) RETURN

      Call DCopy_(NAUCX,[Zero],0,AUCX,1)
!
!     ---- FP
!
      IF (NWFP > 0.AND.NINDEP(ISYM,8) > 0) THEN
        !! Read the T-amplitude
        ICASE=8
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP /= 0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASP,nISP,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
              CALL GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A(1),nASP)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASP,nISP,ipTP)
              CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            END IF
#endif
          End If
        End If

        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          IXMAX=NX
          IF(ISYU == ISYX) IXMAX=IU
          DO IX=1,IXMAX
            IXABS=IX+NAES(ISYX)
            SCL1=Half
            IF(IUABS == IXABS) SCL1=Quart
            IW1=KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              DO IC=1,NC
                ICABS=IC+NSES(ISYC)
                SCL=SCL1
                IF(IAABS >= ICABS) THEN
                  IW2=KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
                  IF(IAABS == ICABS) SCL=SQ2*SCL1
                ELSE
                  IW2=KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
                END IF
                IW=IW1+NASP*(IW2-1)
                AUCX(IA,IU,IC,IX) = SCL*GA_Arrays(ipTP)%A(IW)
              END DO
            END DO
          END DO
        END DO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTP)
        else
#endif
          CALL RHS_FREE(ipTP)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF
!
!     ---- FM
!
      IF (NWFM > 0.AND.NINDEP(ISYM,9) > 0) THEN
        !! Read the T-amplitude
        ICASE=9
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM /= 0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASM,nISM,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
              CALL GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A(1),nASM)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASM,nISM,ipTM)
              CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            End If
#endif
          End If
        End If

        DO IU=1,NU
          IUABS=IU+NAES(ISYU)
          IXMAX=NX
          IF(ISYU == ISYX) IXMAX=IU-1
          DO IX=1,IXMAX
            IXABS=IX+NAES(ISYX)
            IW1=KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              DO IC=1,NC
                ICABS=IC+NSES(ISYC)
                IF(IAABS > ICABS) THEN
                  IW2=KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
                  SCL = -Half
                ELSE IF(IAABS < ICABS) THEN
                  IW2=KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
                  SCL =  Half
                ELSE
                  CYCLE
                END IF
                IW=IW1+NASM*(IW2-1)
                AUCX(IA,IU,IC,IX) = AUCX(IA,IU,IC,IX)
     &            + SCL*GA_Arrays(ipTM)%A(IW)
              END DO
            END DO
          END DO
        END DO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTM)
        else
#endif
        CALL RHS_FREE(ipTM)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF

      Call DGEMM_('T','N',NC*NX,NCHO,NA*NU,
     &            One,AUCX,NA*NU,Cho_Bra,NA*NU,
     &            One,Cho_KetD,NC*NX)
      Call DGEMM_('N','N',NA*NU,NCHO,NC*NX,
     &            One,AUCX,NA*NU,Cho_Ket,NC*NX,
     &            One,Cho_BraD,NA*NU)

      RETURN

      End Subroutine OLagNS_RI_F
!
!-----------------------------------------------------------------------
!
      SUBROUTINE OLagNS_RI_G(ISYI,ISYK,NA,NU,NC,NL,AUCL,NAUCL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KAGEB, KAGTB

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NC, NL,
     &  NAUCL, NCHO
      real(kind=wp), intent(_OUT_) :: AUCL(NA,NU,*)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO),
     &  Cho_Ket(NC,NL,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO),
     &  Cho_KetD(NC,NL,NCHO)
*      Logical Incore
      integer(kind=iwp) :: IOFF1(8), IOFF2(8)
      integer(kind=iwp) :: ISI, ICSTA, ILSTA, IL, IC, IA, IU

      ISYU = ISYI
      ISYL = ISYK

      ISYA=MUL(JSYM,ISYU)
      ISYC=MUL(JSYM,ISYL)
      ISYM=ISYU
      ISYAC=MUL(ISYA,ISYC)
! Set up offset table:
      IO1=0
      IO2=0
      DO ISI=1,NSYM
        IOFF1(ISI)=IO1
        IOFF2(ISI)=IO2
        ISAB=MUL(ISI,ISYM)
        IO1=IO1+NISH(ISI)*NAGEB(ISAB)
        IO2=IO2+NISH(ISI)*NAGTB(ISAB)
      END DO

!   Allocate W with parts WP,WM
      NAS=NASH(ISYM)
      NISP=NISUP(ISYM,10)
      NISM=NISUP(ISYM,11)
      ! NIS=NISP+NISM
      NWGP=NAS*NISP
      NWGM=NAS*NISM
      ! NWG=NWGP+NWGM
!
!     LDGP=NAS
!     LDGM=NAS
!
!     ---- GP
!
      IF (NWGP > 0) THEN
        !! Read the T-amplitude
        ICASE=10
        nINP = nINDEP(iSym,iCase)
        nASP = nASup(iSym,iCase)
        If (nINP /= 0) Then
          nISP = nISup(iSym,iCase)
          nVec = nINP*nISP
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASP,nISP,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
              CALL GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A(1),nASP)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASP,nISP,ipTP)
              CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            END IF
#endif
          End If
        End If

        NBXSZC=NSECBX
        KCL=NAUCL/(NA*NU)
        NBXSZL=KCL/NC
        IF (NBXSZL <= 0) THEN
          Write (u6,*) 'Not enough memory in ADDRHSG, I give up'
          CALL Abend()
        ENDIF

        DO ICSTA=1,NC,NBXSZC
          ICEND=MIN(ICSTA-1+NBXSZC,NC)
          NCSZ=ICEND-ICSTA+1
          DO ILSTA=1,NL,NBXSZL
            ILEND=MIN(ILSTA-1+NBXSZL,NL)
            ! NLSZ=ILEND-ILSTA+1

            ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
            Call DCopy_(NAUCL,[Zero],0,AUCL,1)

        ICL=0
        ! IBUF=0
        DO IL=ILSTA,ILEND
          ! ILABS=IL+NIES(ISYL)
          DO IC=ICSTA,ICEND
            ICABS=IC+NSES(ISYC)
            ICL=ICL+1

            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              SCL=SQ05
              IF(IAABS >= ICABS) THEN
               IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
               IF(IAABS == ICABS) SCL=One
              ELSE
               IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
              END IF
              DO IU=1,NU
                ! IUABS=IU+NAES(ISYU)
                IW1=IU
                IW2=IL+NL*(IAGEC-1)+IOFF1(ISYL)
                IW=IW1+NAS*(IW2-1)

                AUCL(IA,IU,ICL) = SCL*GA_Arrays(ipTP)%A(IW)
              END DO
            END DO
          END DO
        END DO

        Call DGEMM_('N','N',NA*NU,NCHO,ICL,
     &              One,AUCL,NA*NU,Cho_Ket(ICLSTA,1,1),NC*NL,
     &              One,Cho_BraD(1,1,1),NA*NU)
        Call DGEMM_('T','N',ICL,NCHO,NA*NU,
     &              One,AUCL,NA*NU,Cho_Bra(1,1,1),NA*NU,
     &              One,Cho_KetD(ICLSTA,1,1),NC*NL)

          ENDDO
        ENDDO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTP)
        else
#endif
          CALL RHS_FREE(ipTP)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF
!
!     ---- GM
!
      IF (NWGM > 0) THEN
        !! Read the T-amplitude
        ICASE=11
        nINM = nINDEP(iSym,iCase)
        nASM = nASup(iSym,iCase)
        If (nINM /= 0) Then
          nISM = nISup(iSym,iCase)
          nVec = nINM*nISM
          If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
            IF (Is_Real_Par()) THEN
              ! copy global array to local buffer
              Call RHS_ALLO(nASM,nISM,lg_V)
              CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
              ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
              CALL GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A(1),nASM)
              CALL RHS_FREE(lg_V)
              CALL GASYNC()
            ELSE
#endif
              Call RHS_ALLO(nASM,nISM,ipTM)
              CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
            End If
#endif
          End If
        End If

        NBXSZC=NSECBX
        KCL=NAUCL/(NA*NU)
        NBXSZL=KCL/NC
        IF (NBXSZL <= 0) THEN
          Write (u6,*) 'Not enough memory in ADDRHSG, I give up'
          CALL Abend()
        ENDIF


        DO ICSTA=1,NC,NBXSZC
          ICEND=MIN(ICSTA-1+NBXSZC,NC)
          NCSZ=ICEND-ICSTA+1
          DO ILSTA=1,NL,NBXSZL
            ILEND=MIN(ILSTA-1+NBXSZL,NL)
            ! NLSZ=ILEND-ILSTA+1

            ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
            Call DCopy_(NAUCL,[Zero],0,AUCL,1)

        ICL=0
        ! IBUF=0
        DO IL=ILSTA,ILEND
          ! ILABS=IL+NIES(ISYL)
          DO IC=ICSTA,ICEND
            ICABS=IC+NSES(ISYC)
            ICL=ICL+1

            DO IA=1,NA
              IAABS=IA+NSES(ISYA)
              IF(IAABS > ICABS) THEN
                IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL=SQ32
              ELSE IF (IAABS < ICABS) Then
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQ32
              ELSE
                CYCLE
              End If
              DO IU=1,NU
                ! IUABS=IU+NAES(ISYU)
                IW1=IU
                IW2=IL+NL*(IAGTC-1)+IOFF2(ISYL)
                IW=IW1+NAS*(IW2-1)

                AUCL(IA,IU,ICL) = SCL*GA_Arrays(ipTM)%A(IW)
              END DO
            END DO
          END DO
        END DO

        Call DGEMM_('N','N',NA*NU,NCHO,ICL,
     &              One,AUCL,NA*NU,Cho_Ket(ICLSTA,1,1),NC*NL,
     &              One,Cho_BraD(1,1,1),NA*NU)
        Call DGEMM_('T','N',ICL,NCHO,NA*NU,
     &              One,AUCL,NA*NU,Cho_Bra(1,1,1),NA*NU,
     &              One,Cho_KetD(ICLSTA,1,1),NC*NL)

          ENDDO
        ENDDO

#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTM)
        else
#endif
          CALL RHS_FREE(ipTM)
#ifdef _MOLCAS_MPP_
        end if
#endif
      END IF

      RETURN

      End Subroutine OLagNS_RI_G
!
!-----------------------------------------------------------------------
!
      !! ADDRHSH
      Subroutine OLagNS_RI_H(ISYI,ISYK,NA,NJ,NC,NL,AJCL,NAJCL,
     &                       Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

      USE SUPERINDEX, only: KAGEB, KAGTB, KIGEJ, KIGTJ

      implicit none

#include "intent.fh"
#include "macros.fh"

      integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NJ, NC, NL,
     &  NAJCL, NCHO
      real(kind=wp), intent(_OUT_) :: AJCL(NC*NL,*)
      real(kind=wp), intent(in) :: Cho_Bra(NA,NJ,NCHO),
     &  Cho_Ket(NC,NL,NCHO)
      real(kind=wp), intent(inout) :: Cho_BraD(NA,NJ,NCHO),
     &  Cho_KetD(NC,NL,NCHO)

      integer(kind=iwp) :: IASTA, IJSTA, ICSTA, ILSTA, IJ, IA, IL, IC

      ISYJ = ISYI
      ISYL = ISYK

      IF(ISYJ < ISYL) Return
      ISYA=MUL(JSYM,ISYJ)
      ISYC=MUL(JSYM,ISYL)
      ISYM=MUL(ISYA,ISYC)
      ISYAC=ISYM
      ISYJL=ISYM

      NASP=NAGEB(ISYM)
      NISP=NIGEJ(ISYM)
      NWHP=NASP*NISP
      IF(NWHP == 0) Return
!     if (nwhp == 0) write (u6,*) cho_bra(1,1,1) !! avoid unused tenta
!     if (nwhp == 0) write (u6,*) cho_ketd(1,1,1) !! avoid unused tenta
      NASM=NAGTB(ISYM)
      NISM=NIGTJ(ISYM)
      NWHM=NASM*NISM
!
!     LDHP=NASP
!     LDHM=NASM
!
!     ---- HP
!
      !! Read the T-amplitude
      ICASE=12
      nINP = nINDEP(iSym,iCase)
      nASP = nASup(iSym,iCase)
      nVec = 0
      If (nINP /= 0) Then
        nISP = nISup(iSym,iCase)
        nVec = nINP*nISP
        If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            ! copy global array to local buffer
            Call RHS_ALLO(nASP,nISP,lg_V)
            CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
            ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
            CALL GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A(1),nASP)
            CALL RHS_FREE(lg_V)
            CALL GASYNC()
          ELSE
#endif
            Call RHS_ALLO(nASP,nISP,ipTP)
            CALL RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
          END IF
#endif
        End If
      End If

      NBXSZA=NSECBX
      KAJ=NAJCL/(NC*NL)
      NBXSZJ=KAJ/NA
      IF (NBXSZJ <= 0) THEN
        Write (u6,*) 'Not enough memory in ADDRHSH, I give up'
        CALL Abend()
      ENDIF

      NBXSZC=NSECBX
      NBXSZL=NINABX

      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        NASZ=IAEND-IASTA+1
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          NJSZ=IJEND-IJSTA+1

           IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
           Call DCopy_(NAJCL,[Zero],0,AJCL,1)

          DO ICSTA=1,NC,NBXSZC
            ICEND=MIN(ICSTA-1+NBXSZC,NC)
            NCSZ=ICEND-ICSTA+1
            DO ILSTA=1,NL,NBXSZL
              ILEND=MIN(ILSTA-1+NBXSZL,NL)
              ! NLSZ=ILEND-ILSTA+1

              ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

      IAJ=0
      ! IBUF=0
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ == ISYL) ILMAX=IJ
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
          IAJ=IAJ+1

          ICL=0
          DO IL=ILSTA,MIN(ILEND,ILMAX)
            ILABS=IL+NIES(ISYL)
            SCL1=One
            IJGEL=KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
            IF(IJABS == ILABS) SCL1=SQ05
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
              ICL=ICL+1

              SCL=SCL1
              IF(IAABS >= ICABS) THEN
                IAGEC=KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                IF(IAABS == ICABS) SCL=SQ2*SCL1
              ELSE
                IAGEC=KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
              END IF
              IW=IAGEC+NAGEB(ISYM)*(IJGEL-1)

              AJCL(ICLSTA+ICL-1,IAJ) = SCL*GA_Arrays(ipTP)%A(IW)
            END DO
          END DO
        END DO
      END DO

            ENDDO
          ENDDO
          Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NC*NL,
     &                One,AJCL(1,1),NC*NL,Cho_Ket(1,1,1),NC*NL,
     &                One,Cho_BraD(IAJSTA,1,1),NA*NJ)
          Call DGEMM_('N','N',NC*NL,NCHO,NASZ*NJSZ,
     &                One,AJCL(1,1),NC*NL,Cho_Ket(IAJSTA,1,1),NC*NL,
     &                One,Cho_BraD(1,1,1),NA*NJ)

        ENDDO
      ENDDO

      if (nVec /= 0) then
#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTP)
        else
#endif
          CALL RHS_FREE(ipTP)
#ifdef _MOLCAS_MPP_
        end if
#endif
      end if
!
!     ---- HM
!
      !! Read the T-amplitude
      ICASE=13
      nINM = nINDEP(iSym,iCase)
      nASM = nASup(iSym,iCase)
      nVec = 0
      If (nINM /= 0) Then
        nISM = nISup(iSym,iCase)
        nVec = nASM*nISM
        If (nVec /= 0) Then
#ifdef _MOLCAS_MPP_
          IF (Is_Real_Par()) THEN
            ! copy global array to local buffer
            Call RHS_ALLO(nASM,nISM,lg_V)
            CALL RHS_READ_C(lg_V,iCase,iSym,iVecC2)
            ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
            CALL GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A(1),nASM)
            CALL RHS_FREE(lg_V)
            CALL GASYNC()
          ELSE
#endif
            Call RHS_ALLO(nASM,nISM,ipTM)
            CALL RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#ifdef _MOLCAS_MPP_
          End If
#endif
        End If
      End If

      NBXSZA=NSECBX
      KAJ=NAJCL/(NC*NL)
      NBXSZJ=KAJ/NA
      IF (NBXSZJ <= 0) THEN
        Write (u6,*) 'Not enough memory in ADDRHSH, I give up'
        CALL Abend()
      ENDIF
      NBXSZC=NSECBX
      NBXSZL=NINABX

      DO IASTA=1,NA,NBXSZA
        IAEND=MIN(IASTA-1+NBXSZA,NA)
        NASZ=IAEND-IASTA+1
        DO IJSTA=1,NJ,NBXSZJ
          IJEND=MIN(IJSTA-1+NBXSZJ,NJ)
          NJSZ=IJEND-IJSTA+1

          IAJSTA=1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
          Call DCopy_(NAJCL,[Zero],0,AJCL,1)

          DO ICSTA=1,NC,NBXSZC
            ICEND=MIN(ICSTA-1+NBXSZC,NC)
            NCSZ=ICEND-ICSTA+1
            DO ILSTA=1,NL,NBXSZL
              ILEND=MIN(ILSTA-1+NBXSZL,NL)
              ! NLSZ=ILEND-ILSTA+1

              ICLSTA=1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

      IAJ=0
      ! IBUF=0
      DO IJ=IJSTA,IJEND
        IJABS=IJ+NIES(ISYJ)
        ILMAX=NL
        IF(ISYJ == ISYL) ILMAX=IJ-1
        DO IA=IASTA,IAEND
          IAABS=IA+NSES(ISYA)
          IAJ=IAJ+1

          ICL=0
          DO IL=ILSTA,MIN(ILMAX,ILEND)
            ILABS=IL+NIES(ISYL)
            IJGTL=KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
            DO IC=ICSTA,ICEND
              ICABS=IC+NSES(ISYC)
              ICL=ICL+1

              IF (IAABS > ICABS) THEN
                IAGTC=KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL= SQ3
              ELSE IF(IAABS < ICABS) THEN
                IAGTC=KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL=-SQ3
              ELSE
                Cycle
              ENDIF
              IW=IAGTC+NAGTB(ISYM)*(IJGTL-1)

              AJCL(ICLSTA+ICL-1,IAJ) = SCL*GA_Arrays(ipTM)%A(IW)
            END DO
          END DO
        END DO
      END DO

            ENDDO
          ENDDO
          Call DGEMM_('T','N',NASZ*NJSZ,NCHO,NC*NL,
     &                One,AJCL(1,1),NC*NL,Cho_Ket(1,1,1),NC*NL,
     &                One,Cho_BraD(IAJSTA,1,1),NA*NJ)
          Call DGEMM_('N','N',NC*NL,NCHO,NASZ*NJSZ,
     &                One,AJCL(1,1),NC*NL,Cho_Ket(IAJSTA,1,1),NC*NL,
     &                One,Cho_BraD(1,1,1),NA*NJ)
        ENDDO
      ENDDO

      if (nVec /= 0) then
#ifdef _MOLCAS_MPP_
        IF (Is_Real_Par()) THEN
          call deallocate_GA_array(ipTM)
        else
#endif
        CALL RHS_FREE(ipTM)
#ifdef _MOLCAS_MPP_
        end if
#endif
      end if

      unused_var(cho_bra)
      unused_var(cho_ketd)

      Return

      End Subroutine OLagNS_RI_H
!
!-----------------------------------------------------------------------
!
      Subroutine Cnst_A_PT2(block1,block2)

      implicit none

      integer(kind=iwp), intent(in) :: block1, block2

      integer(kind=iwp) :: ndim1
#ifdef _MOLCAS_MPP_
      logical(kind=iwp) :: bstat
      integer(kind=iwp) :: ILOV1, IHIV1, JLOV1, JHIV1, MV1, iRank
#endif

      ndim1 = nSh(iSym0,block1)*nSh(iSym0,block2)
      if (ndim1 == 0) return

#ifdef _MOLCAS_MPP_
      if (is_real_par()) then
        !! The BRAD vector is first put from local BRA(:) to the
        !! global LG_V1 array. The data is shared among all processes.
        !! The bare Cholesky vector (KET(:)) is local.

        !! NV(I): Number of local Cholesky vectors for BRA(:)
        !! NDIM2: number of total Cholesky vectors for this batch
        !!        (usually NV*NPROCS)
        !! NVJ: number of local Cholesky vectors for KET(:)
        bStat = GA_CREATE_IRREG(MT_DBL,NDIM1,NDIM2,'BRAD',
     &                          1,1,MAP2,NPROCS,LG_V1)
        Call Cholesky_Vectors(2,block1,block2,JSYM,BRA,
     &                        nBra,IBSTA,IBEND)
        !! finds out the range of the global array g_a that process
        !! iproc owns.
        CALL GA_DISTRIBUTION(LG_V1,MYRANK,ILOV1,IHIV1,JLOV1,JHIV1)
        !! provides access to local data in the specified patch of
        !! the array owned by the calling process
        CALL GA_ACCESS(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,MV1,NDIM1)
        CALL DCOPY_(NDIM1*(JHIV1-JLOV1+1),BRA,1,DBL_MB(MV1),1)
        !! releases access to a global array, when the data was
        !! accessed for writing
        CALL GA_RELEASE_UPDATE(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1)
        CALL GA_SYNC()

        !! ket is ndim1*NVJ dimension
        Call Get_Cholesky_Vectors(block1,block2,JSYM,KET,nKet,
     &                            JBSTA,JBEND)

        do iRank = 0, NPROCS-1
          CALL GA_DISTRIBUTION(LG_V1,iRank,ILOV1,IHIV1,JLOV1,JHIV1)
          CALL GA_GET(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,BRA,NDIM1)
          Call DGEMM_('T','N',JHIV1-JLOV1+1,NVJ,ndim1,
     &                One,BRA,ndim1,KET,ndim1,
     &    One,A_PT2(IOFFCV+JLOV1-1,JOFFCV+MAP2(myRank+1)-1),MaxVec_PT2)
        end do
        bStat =  GA_DESTROY(LG_V1)
      else
#endif
        Call Cholesky_Vectors(2,block1,block2,JSYM,BRA,nBra,
     &                        IBSTA,IBEND)
        Call Get_Cholesky_Vectors(block1,block2,JSYM,KET,nKet,
     &                            JBSTA,JBEND)
        Call DGEMM_('T','N',NVI,NVJ,ndim1,
     &              One,BRA,ndim1,KET,ndim1,
     &              One,A_PT2(IOFFCV,JOFFCV),MaxVec_PT2)
#ifdef _MOLCAS_MPP_
      end if
#endif

      End Subroutine Cnst_A_PT2

      End Subroutine OLagNS_RI
!
!-----------------------------------------------------------------------
!
      Subroutine Cholesky_Vectors(MODE,ITK,ITQ,JSYM,Array,nArray,
     &                            IBSTA,IBEND)

      USE CHOVEC_IO, only: NVLOC_CHOBATCH, IDLOC_CHOGROUP, NPQ_CHOTYPE
      use caspt2_module, only: NSYM
      use definitions, only: wp, iwp

      implicit none

#include "intent.fh"

      integer(kind=iwp), intent(in) :: MODE, ITK, ITQ, JSYM, IBSTA,
     &  IBEND
      integer(kind=iwp), intent(out) :: nArray
      real(kind=wp), intent(_OUT_) :: Array(*)

      integer(kind=iwp) :: ICASE, LKETSM, LUCDER, ISYK, NQK, IB, NV,
     &  NKETSM, IDISK

      ! ugly hack to convert separate k/q orbital types into a specific
      ! case
      ICASE=ITK*ITQ
      IF (ICASE == 3) THEN
        ICASE=4
      ELSE
        ICASE=ICASE/2
      END IF

      LKETSM=1
      LUCDER = 63 ! tentative
      DO ISYK=1,NSYM
        NQK=NPQ_CHOTYPE(ICASE,ISYK,JSYM)
        IF(NQK == 0) CYCLE
        DO IB=IBSTA,IBEND
          NV=NVLOC_CHOBATCH(IB)
          NKETSM=NQK*NV
          IDISK=IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
          !! MODE = 1: write
          !! MODE = 2: read
          CALL DDAFILE(LUCDER,MODE,Array(LKETSM),NKETSM,IDISK)
          LKETSM=LKETSM+NKETSM
        END DO
      END DO
      nArray=LKETSM-1
*
      End Subroutine Cholesky_Vectors
