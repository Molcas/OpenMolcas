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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
! based on rhsall2
!
! In principle, For E = T_{ij}^{ab}*(ia|jb).
! With p and q for general orbitals,
! L_{pq} = (pa|jb)*T_{qj}^{ab} + (ip|jb)*T_{ij}^{qb}
!        + (ia|pb)*T_{iq}^{ab} + (ia|jp)*T_{ij}^{aq}
! For the first term, with P and Q for auxiliary orbitals
! L_{pq}(1) = (pa|jb)*T_{qj}^{ab}
!           = (pa|P)*(P|jb) * T_{qj}^{ab} (2c-2e is omitted?)
!           = (pa|P) * tilde{T}_{qa}^P
!           = C_{mu p}*C_{nu a}*(mu nu|P) * tilde{T}_{qa}^P
!           = C_{mu p} * (mu nu|P) * V_{nu q}^P
! where
! tilde{T}_{ia}^P = T_{ij}^{ab} * (P|jb)
! V_{mu p}^P      = tilde{T}_{pa}^P * C_{mu a}
!
! Dimension of tilde{T} will be the same as that of (ia|P) in MO,
! where (i,a) = (inact,active), (inact,virtual), (active,virtual)
!
! tilde{T} is constructed in this file.
! tilde{T} -> V_{mu p}^P transformations, contraction with 3c-2e,
! and construction of the orbital Lagrangian (L_{pq}) is elsewhere
! ... maybe in OLagVVVO
!
!-----------------------------------------------------------------------
!
! For the ERI derivative calculation,
! d(mu nu|rho sigma)/da
! = d/da (mu nu|P) (P|Q)^-1 (Q|rho sigma)
! = d(mu nu|P)/da (P|Q)^-1 (Q|rho sigma)
!   + (mu nu|P) d(P|Q)^-1/da (Q|rho sigma)
!   + (mu nu|P) (P|Q)^-1 d(Q|rho sigma)/da
! = d(mu nu|P)/da (P|Q)^-1 (Q|rho sigma)
!   - (mu nu|P) (P|R)^-1 d(R|S)/da (S|Q)^-1 (Q|rho sigma)
!   + (mu nu|P) (P|Q)^-1 d(Q|rho sigma)/da
! = d(mu nu|P)/da (tP|rho sigma)
!   - (mu nu|tP) d(P|Q)/da (tQ|rho sigma)
!   + (mu nu|tP) d(P|rho sigma)/da
! where (mu nu|tP) = (mu nu|Q)*(Q|P)^-1
!
! D_{mu nu rho sigma}*d(mu nu|rho sigma)/da
! = d(mu nu|P)/da (tP|rho sigma) * D_{mu nu rho sigma}
!   - (mu nu|tP) d(P|Q)/da (tQ|rho sigma) * D_{mu nu rho sigma}
!   + (mu nu|tP) d(P|rho sigma)/da * D_{mu nu rho sigma}
! = d(mu nu|P)/da tD_{mu nu}^tP
!   - (mu nu|tP) d(P|Q)/da tD_{mu nu}^tQ
!   + tD_{rho sigma}^tP d(P|rho sigma)/da
! where tD_{mu nu}^tP = D_{mu nu rho sigma} * (rho sigma|tP)
! In practice, tD_{pq}^tP is constructed and saved in disk.
! This will be read when 3c-2e ERI derivatives are concerned,
! and MO->AO transformations will be done on-the-fly.
! Note that MO coefficients of CASPT2 have to be used.
!
! For 2c-2e ERI derivatives,
! D(tP,tQ) = tD_{pq} * C_{mu p} C_{nu q} * (mu nu|tP)
! then saved.

subroutine OLagNS_RI(iSym0,NBSQT,MaxVec_PT2,DPT2C,DPT2Canti,A_PT2)

use Symmetry_Info, only: Mul
use CHOVEC_IO, only: NVLOC_CHOBATCH
use caspt2_global, only: iPrGlb
use caspt2_global, only: do_csf, iStpGrd
use PrintLevel, only: VERBOSE
use EQSOLV, only: IVECC2
use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: iwp, wp, u6
use fake_GA, only: GA_Arrays
#ifdef _MOLCAS_MPP_
use fake_GA, only: Allocate_GA_Array, Deallocate_GA_Array
use Para_Info, only: Is_Real_Par
#endif
use caspt2_module, only: NACTEL, NSYM, NFRO, NISH, NIES, NASH, NAES, NSSH, NSES, NORBT, NINABX, NSECBX
use caspt2_module, only: NTUV, NTU, NTGEU, NTGTU, NIGEJ, NIGTJ, NAGEB, NAGTB, NTUVES, NTUES, NTGEUES, NTGTUES, NIGEJES, NIGTJES, &
                         NAGEBES, NAGTBES, NASUP, NISUP, NINDEP, NBTCH, NBTCHES
use Constants, only: Zero, One, Quart, Half, Two, Three, OneHalf

#include "intent.fh"
#include "macros.fh"

implicit none
#include "warnings.h"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif
integer(kind=iwp), intent(in) :: iSym0, NBSQT, MaxVec_PT2
real(kind=wp), intent(inout) :: DPT2C(NBSQT), DPT2Canti(NBSQT), A_PT2(MaxVec_PT2,MaxVec_PT2)
integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3
integer(kind=iwp), allocatable :: BGRP(:,:)
real(kind=wp), allocatable :: BRA(:), KET(:), BRAD(:), KETD(:), PIQK(:)
#ifdef _MOLCAS_MPP_
integer(kind=iwp), allocatable :: map2(:)
integer(kind=iwp) :: myRank, NPROCS, i, lg_V, lg_V1, ndim2
#endif
integer(kind=iwp) :: nSh(8,3), iSym, JSYM, IB1, IB2, MXBGRP, IBGRP, IB, NBGRP, iStpGrd_sav, NCHOBUF, MXPIQK, NADDBUF, IOFFCV, &
                     IBSTA, IBEND, NV, nBra, nKet, NVI, JOFFCV, JBGRP, JBSTA, JBEND, NVJ, JB
real(kind=wp) :: SCLNEL
integer(kind=iwp) :: ICASE, ISYT, ISYU, ISYV, ISYX, ISYJ, ISYL, ISYA, ISYC, ISYJL, ISYAC, ISYII, ISIJ, ISAB, ISI, NAS, NIS, NIN, &
                     NAS1, NASP, NISP, NINP, NASM, NISM, NINM, NWA, NWBP, NWBM, NWC, NWD, NWEP, NWEM, NWFP, NWFM, NWGP, NWGM, &
                     NWHP, NWHM, NW, NVEC, ipT, ipTP, ipTM, ipTanti, ITABS, iTtot, IUABS, iUtot, IVABS, iVtot, IVMAX, IXABS, &
                     iXtot, IXMAX, IJABS, iJtot, ILABS, IAABS, iAtot, ICABS, IW1, IW2, IW, nOrbA, IO, IO1, IO2, NBXSZA, NBXSZC, &
                     NBXSZJ, NBXSZL, IAEND, NASZ, ICEND, NCSZ, IJEND, NJSZ, ILEND, ILMAX, IAJSTA, IAJ, KAJ, ICLSTA, ICL, KCL, &
                     IAGEC, IAGTC, IJGEL, IJGTL, JGEL, JGTL
real(kind=wp) :: SCL, SCL1
real(kind=wp), parameter :: SQ2 = sqrt(Two), SQ3 = sqrt(Three), SQ05 = sqrt(Half), SQ32 = sqrt(OneHalf)

nSh(1:nSym,Inactive) = NISH(1:nSym)
nSh(1:nSym,Active) = NASH(1:nSym)
nSh(1:nSym,Virtual) = NSSH(1:nSym)

if (IPRGLB >= VERBOSE) write(u6,'(1X,A)') ' Using RHSALL2+ADDRHS algorithm'

!iVec = iVecX
iSym = iSym0
SCLNEL = One/real(max(1,NACTEL),kind=wp)
!                                                                      *
!***********************************************************************
!                                                                      *
do JSYM=1,NSYM

  IB1 = NBTCHES(JSYM)+1
  IB2 = NBTCHES(JSYM)+NBTCH(JSYM)

  MXBGRP = IB2-IB1+1
  if (MXBGRP <= 0) cycle
  call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
  IBGRP = 1
  do IB=IB1,IB2
    BGRP(1,IBGRP) = IB
    BGRP(2,IBGRP) = IB
    IBGRP = IBGRP+1
  end do
  NBGRP = MXBGRP

  !! With iStpGrd = -1, we try to allocate 4 large arrays
  iStpGrd_sav = iStpGrd
  iStpGrd = -1
  call MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,NCHOBUF,MXPIQK,NADDBUF)
  iStpGrd = iStpGrd_sav
  if (IPRGLB > VERBOSE) then
    write(u6,*)
    write(u6,'(A,I12)') '  Number of Cholesky batches: ',IB2-IB1+1
    write(u6,'(A,I12)') '  Number of batch groups:     ',NBGRP
    write(u6,*)
  end if
  ! buffers are kept allocated until the end of JSYM loop.
  call mma_allocate(PIQK,MXPIQK,Label='PIQK')
  call mma_allocate(BRA,NCHOBUF,Label='BRABUF')
  call mma_allocate(KET,NCHOBUF,Label='KETBUF')
  call mma_allocate(BRAD,NCHOBUF,Label='BRAD')
  call mma_allocate(KETD,NCHOBUF,Label='KETD')

  ! Loop over groups of batches of Cholesky vectors

  IOFFCV = 1
  do IBGRP=1,NBGRP

    IBSTA = BGRP(1,IBGRP)
    IBEND = BGRP(2,IBGRP)

    NV = 0
    do IB=IBSTA,IBEND
      NV = NV+NVLOC_CHOBATCH(IB)
    end do

    if (IPRGLB > VERBOSE) then
      write(u6,'(A,I12)') '  Cholesky vectors in this group = ',NV
      write(u6,*)
    end if

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      myRank = GA_NodeID()
      NPROCS = GA_nNodes()

      call mma_allocate(MAP2,NPROCS,Label='MAP2')
      MAP2(:) = 0
      MAP2(myRank+1) = NV
      call GAIGOP(MAP2,NPROCS,'+')
      ndim2 = sum(map2)

      do i=nprocs,2,-1
        map2(i) = sum(map2(1:i-1))+1
      end do
      map2(1) = 1
    end if
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets (Cholesky vectors) in the form L(VX), all symmetries:

    call Get_Cholesky_Vectors(Active,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    KETD(1:nKet) = Zero
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read bra (Cholesky vectors) in the form L(TJ): All symmetries

    call Get_Cholesky_Vectors(Inactive,Active,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    BRAD(1:nBra) = Zero
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Assemble contributions to TJVX
    ! Loop over the bras and kets, form <A|0>

    call OLagNS_RI2(Inactive,Active,Active,Active,'A ',size(BRA),size(KET),BRA,KET,BRAD,KETD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! TJVL RHSB
    ! TJVL: Use TJ buffer as if it was VL, form <B|0>

    !nKet = nBra
    call OLagNS_RI2(Inactive,Active,Inactive,Active,'B ',size(BRA),size(BRA),BRA,BRA,BRAD,BRAD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read bra (Cholesky vectors) in the form L(AJ), form <D1|0>
    ! We still have L(VX) vectors in core, at KETS.

    call Cholesky_Vectors(1,Inactive,Active,JSYM,BRAD,size(BRAD),nBra,IBSTA,IBEND)
    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    BRAD(1:nBra) = Zero
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJVX RHSD1
    ! Loop over the bra and ket vectors.

    call OLagNS_RI2(Inactive,Virtual,Active,Active,'D1',size(BRA),size(KET),BRA,KET,BRAD,KETD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJCL RHSH
    ! AJCL: Use AJ buffer still in core as if it was CL, form <H|0>

    !nKet = nBra
    call OLagNS_RI2(Inactive,Virtual,Inactive,Virtual,'H ',size(BRA),size(BRA),BRA,BRA,BRAD,BRAD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read Bra (Cholesky vectors)= L(AU)

    call Cholesky_Vectors(1,Inactive,Virtual,JSYM,BRAD,size(BRAD),nBra,IBSTA,IBEND)
    call Get_Cholesky_Vectors(Active,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    BRAD(1:nBra) = Zero
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUVX RHSC
    ! AUVX: Loop over the bras and kets

    call OLagNS_RI2(Active,Virtual,Active,Active,'C ',size(BRA),size(KET),BRA,KET,BRAD,KETD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUCX RHSF
    ! AUCX: Use AU buffer still in core as if it was CX, form <F|0>

    !     nKet = nBra
    call OLagNS_RI2(Active,Virtual,Active,Virtual,'F ',size(BRA),size(BRA),BRA,BRA,BRAD,BRAD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets (Cholesky vectors) in the form L(VL), all symmetries:

    call Cholesky_Vectors(1,Active,Active,JSYM,KETD,size(KETD),nKet,IBSTA,IBEND)
    call Get_Cholesky_Vectors(Inactive,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    call Cholesky_Vectors(2,Inactive,Active,JSYM,KETD,size(KETD),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUVL RHSD2
    ! Loop over bras and kets, form <D2|0>.

    call OLagNS_RI2(Active,Virtual,Inactive,Active,'D2',size(BRA),size(KET),BRA,KET,BRAD,KETD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets (Cholesky vectors) in the form L(CL), all symmetries:

    call Cholesky_Vectors(1,Inactive,Active,JSYM,KETD,size(KETD),nKet,IBSTA,IBEND)
    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    call Cholesky_Vectors(2,Inactive,Virtual,JSYM,KETD,size(KETD),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AUCL RHSG
    ! Loop over bras and kets, form  <G|0>

    call OLagNS_RI2(Active,Virtual,Inactive,Virtual,'G ',size(BRA),size(KET),BRA,KET,BRAD,KETD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read bra vectors AJ

    call Cholesky_Vectors(1,Active,Virtual,JSYM,BRAD,size(BRAD),nBra,IBSTA,IBEND)
    call Get_Cholesky_Vectors(Inactive,Virtual,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    BRAD(1:nBra) = KETD(1:nBra)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Read kets in the form L(VL)

    call Get_Cholesky_Vectors(Inactive,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)
    call Cholesky_Vectors(2,Inactive,Active,JSYM,KETD,size(KETD),nKet,IBSTA,IBEND)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! AJVL RHSE
    ! AJVL: Loop over bras and kets. Form <E|0>

    call OLagNS_RI2(Inactive,Virtual,Inactive,Active,'E ',size(BRA),size(KET),BRA,KET,BRAD,KETD)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! End of loop over batches, IB
    call Cholesky_Vectors(1,Inactive,Virtual,JSYM,BRAD,size(BRAD),nBra,IBSTA,IBEND)
    call Cholesky_Vectors(1,Inactive,Active,JSYM,KETD,size(KETD),nKet,IBSTA,IBEND)
      !! Construct A_PT2
    NVI = NV
    JOFFCV = 1
    do JBGRP=1,NBGRP

      JBSTA = BGRP(1,JBGRP)
      JBEND = BGRP(2,JBGRP)

      NVJ = 0
      do JB=JBSTA,JBEND
        NVJ = NVJ+NVLOC_CHOBATCH(JB)
      end do

        !! BraAI
      call Cnst_A_PT2(Inactive,Active)

        !! BraSI
      call Cnst_A_PT2(Inactive,Virtual)

        !! BraSA
      call Cnst_A_PT2(Active,Virtual)

        !! BraAA
      call Cnst_A_PT2(Active,Active)
      JOFFCV = JOFFCV+NVJ
    end do !! end of JBGRP loop
    IOFFCV = IOFFCV+NVI
  end do !! end of IBGRP loop
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call mma_deallocate(BRA)
  call mma_deallocate(KET)
  call mma_deallocate(BRAD)
  call mma_deallocate(KETD)
  call mma_deallocate(PIQK)
  call mma_deallocate(BGRP)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! End of loop over JSYM
end do
#ifdef _MOLCAS_MPP_
if (Is_Real_Par()) call mma_deallocate(map2)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Synchronized add RHS partial arrays from all nodes into each node.

!-SVC: read the DRA's from disk and copy them all to LUSOLV to continue
!      in serial mode.  FIXME: this call has to be removed when we reach
!      full parallel capabilities
!call SYNRHS(IVEC)
!-SVC: at this point, the RHS elements are on disk, both in LUSOLV and
!      as DRAs with the name RHS_XX_XX_XX with XX a number representing
!      the case, symmetry, and rhs vector respectively.

A_PT2(:,:) = Two*A_PT2(:,:)

if (NBGRP /= 0) SCLNEL = SCLNEL/real(NBGRP,kind=wp)
DPT2C(1:NBSQT) = DPT2C(1:NBSQT)*SCLNEL
if (do_csf) DPT2Canti(1:NBSQT) = DPT2Canti(1:NBSQT)*SCLNEL

#ifdef _MOLCAS_MPP_
if (is_real_par()) call GADGOP(A_PT2,MaxVec_PT2**2,'+')
#endif

return

contains

subroutine OLagNS_RI2(ITI,ITP,ITK,ITQ,case,nBra,nKet,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD)

  use caspt2_global, only: iPrGlb
  use PrintLevel, only: DEBUG

  integer(kind=iwp), intent(in) :: ITI, ITP, ITK, ITQ, nBra, nKet
  character(len=2), intent(in) :: case
  real(kind=wp), intent(in) :: Cho_Bra(nBra), Cho_Ket(nKet)
  real(kind=wp), intent(inout) :: Cho_BraD(nBra), Cho_KetD(nKet)
  integer(kind=iwp) :: LBRASM, ISYI, NI, ISYP, NP, NPI, NBRASM,LKETSM, ISYK, NK, ISYQ, NQ, NQK, NKETSM, NPIQK, KPI, KQK
  real(kind=wp) :: TotCPU0, TotWall0, TotCPU1, TotWall1

  if (iPrGlb >= DEBUG) write(u6,*) 'Processing RHS block '//case

  LBRASM = 1
  call CWTime(TotCPU0,TotWall0)
  do ISYI=1,NSYM
    NI = NSH(ISYI,ITI)
    if (NI == 0) cycle
    ISYP = Mul(ISYI,JSYM)
    NP = NSH(ISYP,ITP)
    if (NP == 0) cycle
    NPI = NP*NI
    NBRASM = NPI*NV

    LKETSM = 1
    do ISYK=1,NSYM
      NK = NSH(ISYK,ITK)
      if (NK == 0) cycle
      ISYQ = Mul(ISYK,JSYM)
      NQ = NSH(ISYQ,ITQ)
      if (NQ == 0) cycle
      NQK = NQ*NK
      NKETSM = NQK*NV

      ! SVC: we need an NPI*NQK to store the 2-electron integrals, and 2
      ! buffers (values+indices) for sorting them.  Later, we can try to get
      ! rid of the buffer that stores the values and only use an index buffer
      ! and the two-electron integrals for the scatter operation.  For the
      ! buffer, any size can be taken, but assuming there is enough memory
      ! available, it's set to the size of the two-electron integrals unless
      ! larger than some predefined maximum buffer size.
      NPIQK = NPI*NQK
      if (NPIQK > MXPIQK) then
        if (case == 'H') then
          KPI = MXPIQK/NQK
          NPIQK = KPI*NQK
        else if (case == 'G') then
          KQK = MXPIQK/NPI
          NPIQK = NPI*KQK
        else
          write(u6,*) ' NPIQK > MXPIQK and case != G or H'
          write(u6,'(A,A2)') ' CASE =   ',case
          write(u6,'(A,I12)') ' NPIQK =  ',NPIQK
          write(u6,'(A,I12)') ' MXPIQK = ',MXPIQK
          write(u6,*) ' This should not happen, please report.'
          call AbEnd()
        end if
      end if

      if (NPIQK <= 0) then
        write(u6,'(1X,A)') ' ADDRHS: zero-sized NPIQK'
        call AbEnd()
      end if

      !! NBUFF(=nAddBuf) is removed
      select case (case)
        case ('A ')
          call OLagNS_RI_A(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('B ')
          call OLagNS_RI_B(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('D1')
          call OLagNS_RI_D1(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('H ')
          call OLagNS_RI_H(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('C ')
          call OLagNS_RI_C(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('F ')
          call OLagNS_RI_F(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('D2')
          call OLagNS_RI_D2(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('G ')
          call OLagNS_RI_G(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case ('E ')
          call OLagNS_RI_E(ISYI,ISYK,NP,NI,NQ,NK,PIQK,NPIQK,Cho_Bra(LBRASM),Cho_Ket(LKETSM),Cho_BraD(LBRASM),Cho_KetD(LKETSM),NV)
        case default
          call Abend()
      end select

      LKETSM = LKETSM+NKETSM
    end do
    LBRASM = LBRASM+NBRASM
  end do
  call CWTime(TotCPU1,TotWall1)
  if (IPRGLB >= VERBOSE) write(u6,'(" CPU/Wall Time (Case ",A2,"):",2f10.2)') case,totcpu1-totcpu0,totwall1-totwall0

  return

end subroutine OLagNS_RI2

subroutine OLagNS_RI_A(ISYI,ISYK,NT,NJ,NV,NX,TJVX,NTJVX,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KTUV
  use caspt2_global, only: do_csf

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NT, NJ, NV, NX, NTJVX, NCHO
  real(kind=wp), intent(out) :: TJVX(NT,NJ,NV,NX)
  real(kind=wp), intent(in) :: Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NT,NJ,NCHO), Cho_KetD(NV,NX,NCHO)
  integer(kind=iwp) :: IT, IJ, IV, IX

  ISYJ = ISYI
  ISYX = ISYK

  ISYT = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYX)
  ISYM = ISYJ
  if (NINDEP(ISYM,1) == 0) return
  NAS = NTUV(ISYM)
  NIS = NISH(ISYM)
  NWA = NAS*NIS
  if (NWA == 0) return

  ! ---- A

  !! Read the T-amplitude
  ICASE = 1
  nIN = nINDEP(iSym,iCase)
  nAS = nASup(iSym,iCase)
  if (nIN /= 0) then
    nIS = nISup(iSym,iCase)
    nVec = nAS*nIS
    if (nVec /= 0) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        ! copy global array to local buffer
        call RHS_ALLO(nAS,nIS,lg_V)
        call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
        ipT = Allocate_GA_Array(nAS*nIS,'ipT')
        call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A,nAS)
        if (do_csf) then
          call RHS_READ_C(lg_V,iCase,iSym,7)
          ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
          call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A,nAS)
        end if
        call RHS_FREE(lg_V)
        call GASYNC()
      else
#     endif
        call RHS_ALLO(nAS,nIS,ipT)
        call RHS_READ_C(ipT,iCase,iSym,iVecC2)
        if (do_csf) then
          call RHS_ALLO(nAS,nIS,ipTanti)
          call RHS_READ_C(ipTanti,iCase,iSym,7)
        end if
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end if
  end if

  call DCopy_(NTJVX,[Zero],0,TJVX,1)

  nOrbT = nFro(iSyT)+nIsh(iSyT)+nAsh(iSyT)+nSsh(iSyT)
  do IT=1,NT
    ITABS = IT+NAES(ISYT)
    iTtot = iT+nFro(iSyT)+nIsh(iSyT)
    do IJ=1,NJ
      !IJABS = IJ+NIES(ISYJ)
      iJtot = iJ+nFro(iSyJ)

      do IV=1,NV
        IVABS = IV+NAES(ISYV)
        if (ISYV == ISYX) then !! not sure
          !! ONEADD contributions
          IW1 = KTUV(ITABS,IVABS,IVABS)-NTUVES(ISYM)
          IW2 = IJ
          IW = IW1+NAS*(IW2-1)

          DPT2C(iTtot+nOrbT*(iJtot-1)) = DPT2C(iTtot+nOrbT*(iJtot-1))+Two*GA_Arrays(ipT)%A(IW)
          if (do_csf) DPT2Canti(iTtot+nOrbT*(iJtot-1)) = DPT2Canti(iTtot+nOrbT*(iJtot-1))+Two*GA_Arrays(ipTanti)%A(IW)
        end if
        do IX=1,NX
          IXABS = IX+NAES(ISYX)
          IW1 = KTUV(ITABS,IVABS,IXABS)-NTUVES(ISYM)
          IW2 = IJ
          IW = IW1+NAS*(IW2-1)
          TJVX(IT,IJ,IV,IX) = GA_Arrays(ipT)%A(IW)
        end do
      end do
    end do
  end do

  call DGEMM_('T','N',NV*NX,NCHO,NT*NJ,One,TJVX,NT*NJ,Cho_Bra,NT*NJ,One,Cho_KetD,NV*NX)
  call DGEMM_('N','N',NT*NJ,NCHO,NV*NX,One,TJVX,NT*NJ,Cho_Ket,NV*NX,One,Cho_BraD,NT*NJ)

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    call deallocate_GA_array(ipT)
    if (do_CSF) call deallocate_GA_array(ipTanti)
  else
# endif
    call RHS_FREE(ipT)
    if (do_csf) call RHS_FREE(ipTanti)

# ifdef _MOLCAS_MPP_
  end if
# endif

  return

end subroutine OLagNS_RI_A

subroutine OLagNS_RI_B(ISYI,ISYK,NT,NJ,NV,NL,TJVL,NTJVL,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KIGTJ, KIGEJ, KTGTU, KTGEU

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NT, NJ, NV, NL, NTJVL, NCHO
  real(kind=wp), intent(out) :: TJVL(NT,NJ,NV,NL)
  real(kind=wp), intent(in) :: Cho_Bra(NT,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NT,NJ,NCHO), Cho_KetD(NV,NL,NCHO)
  integer(kind=iwp) :: IT, IV, IJ, IL

  ISYJ = ISYI
  ISYL = ISYK

  ISYT = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYL)
  if (ISYT < ISYV) return
  ISYM = Mul(ISYJ,ISYL) !!

  if (NINDEP(ISYM,2) > 0) then
    ! The plus combination:
    ICASE = 2
    NASP = NTGEU(ISYM)
    NISP = NIGEJ(ISYM)
    NWBP = NASP*NISP
  else
    NWBP = 0
  end if
  if (NINDEP(ISYM,3) > 0) then
    ! The minus combination:
    ICASE = 3
    NASM = NTGTU(ISYM)
    NISM = NIGTJ(ISYM)
    NWBM = NASM*NISM
  else
    NWBM = 0
  end if
  if (max(NWBP,NWBM) <= 0) return

  call DCopy_(NTJVL,[Zero],0,TJVL,1)

  if ((NWBP > 0) .and. (NINDEP(ISYM,2) > 0)) then
    !! Read the T-amplitude
    ICASE = 2
    nINP = nINDEP(iSym,iCase)
    nASP = nASup(iSym,iCase)
    if (nINP /= 0) then
      nISP = nISup(iSym,iCase)
      nVec = nASP*nISP
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASP,nISP,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
          call GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A,nASP)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASP,nISP,ipTP)
          call RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    do IT=1,NT
      ITABS = IT+NAES(ISYT)
      IVMAX = NV
      if (ISYV == ISYT) IVMAX = IT
      do IV=1,IVMAX
        IVABS = IV+NAES(ISYV)
        SCL1 = Half
        IW1 = KTGEU(ITABS,IVABS)-NTGEUES(ISYM)
        if (ITABS == IVABS) SCL1 = Quart
        do IJ=1,NJ
          IJABS = IJ+NIES(ISYJ)
          do IL=1,NL
            ILABS = IL+NIES(ISYL)
            SCL = SCL1
            if (IJABS >= ILABS) then
              IW2 = KIGEJ(IJABS,ILABS)-NIGEJES(ISYM)
              if (IJABS == ILABS) SCL = SQ2*SCL1
            else
              IW2 = KIGEJ(ILABS,IJABS)-NIGEJES(ISYM)
            end if
            IW = IW1+NASP*(IW2-1)
            TJVL(IT,IJ,IV,IL) = SCL*GA_Arrays(ipTP)%A(IW)
          end do
        end do
      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTP)
    else
#   endif
      call RHS_FREE(ipTP)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  if (NINDEP(ISYM,3) > 0) then
    !! Read the T-amplitude
    ICASE = 3
    nINM = nINDEP(iSym,iCase)
    nASM = nASup(iSym,iCase)
    if (nINM /= 0) then
      nISM = nISup(iSym,iCase)
      nVec = nASM*nISM
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASM,nISM,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
          call GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A,nASM)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASM,nISM,ipTM)
          call RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    do IT=1,NT
      ITABS = IT+NAES(ISYT)
      IVMAX = NV
      if (ISYV == ISYT) IVMAX = IT-1
      do IV=1,IVMAX
        IVABS = IV+NAES(ISYV)
        IW1 = KTGTU(ITABS,IVABS)-NTGTUES(ISYM)
        do IJ=1,NJ
          IJABS = IJ+NIES(ISYJ)
          do IL=1,NL
            ILABS = IL+NIES(ISYL)
            if (IJABS > ILABS) then
              IW2 = KIGTJ(IJABS,ILABS)-NIGTJES(ISYM)
              SCL = Half
            else if (IJABS < ILABS) then
              IW2 = KIGTJ(ILABS,IJABS)-NIGTJES(ISYM)
              SCL = -Half
            else
              cycle
            end if
            IW = IW1+NASM*(IW2-1)
            TJVL(IT,IJ,IV,IL) = TJVL(IT,IJ,IV,IL)+SCL*GA_Arrays(ipTM)%A(IW)
          end do
        end do
      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTM)
    else
#   endif
      call RHS_FREE(ipTM)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  call DGEMM_('T','N',NV*NL,NCHO,NT*NJ,One,TJVL,NT*NJ,Cho_Bra,NT*NJ,One,Cho_KetD,NV*NL)
  call DGEMM_('N','N',NT*NJ,NCHO,NV*NL,One,TJVL,NT*NJ,Cho_Ket,NV*NL,One,Cho_BraD,NT*NJ)

  return

end subroutine OLagNS_RI_B

subroutine OLagNS_RI_C(ISYI,ISYK,NA,NU,NV,NX,AUVX,NAUVX,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX
  use caspt2_global, only: do_csf

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NV, NX, NAUVX, NCHO
  real(kind=wp), intent(out) :: AUVX(NA,NU,NV,NX)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NX,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO), Cho_KetD(NV,NX,NCHO)
  real(kind=wp) :: ValCF
  integer(kind=iwp) :: IA, IU, IV, IX

  ISYU = ISYI
  ISYX = ISYK

  ISYA = Mul(JSYM,ISYU)
  ISYV = Mul(JSYM,ISYX)
  ISYM = ISYA !!
  if (NINDEP(ISYM,4) == 0) return
  NAS = NTUV(ISYM)
  NIS = NSSH(ISYM)
  NWC = NAS*NIS
  if (NWC == 0) return

  ! ---- C

  !! Read the T-amplitude
  ICASE = 4
  nIN = nINDEP(iSym,iCase)
  nAS = nASup(iSym,iCase)
  if (nIN /= 0) then
    nIS = nISup(iSym,iCase)
    nVec = nIN*nIS
    if (nVec /= 0) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call RHS_ALLO(nAS,nIS,lg_V)
        call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
        ipT = Allocate_GA_Array(nAS*nIS,'ipT')
        call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A,nAS)
        if (do_csf) then
          call RHS_READ_C(lg_V,iCase,iSym,7)
          ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
          call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A,nAS)
        end if
        call RHS_FREE(lg_V)
        call GASYNC()
      else
#     endif
        call RHS_ALLO(nAS,nIS,ipT)
        call RHS_READ_C(ipT,iCase,iSym,iVecC2)
        if (do_csf) then
          call RHS_ALLO(nAS,nIS,ipTanti)
          call RHS_READ_C(ipTanti,iCase,iSym,7)
        end if
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end if
  end if

  call DCopy_(NAUVX,[Zero],0,AUVX,1)

  nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
  do IA=1,NA
    iAtot = iA+nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)
    do IU=1,NU
      IUABS = IU+NAES(ISYU)
      iUtot = iU+nFro(iSyU)+nIsh(iSyU)
      do IV=1,NV
        IVABS = IV+NAES(ISYV)
        ValCF = Zero
        if (ISYV == ISYX) then !! not sure
          !! ONEADD contributions
          IW1 = KTUV(IUABS,IVABS,IVABS)-NTUVES(ISYM)
          IW2 = IA
          IW = IW1+NAS*(IW2-1)

          ValCF = GA_Arrays(ipT)%A(IW)*Two
          DPT2C(iAtot+nOrbA*(iUtot-1)) = DPT2C(iAtot+nOrbA*(iUtot-1))+ValCF
          if (do_csf) DPT2Canti(iAtot+nOrbA*(iUtot-1)) = DPT2Canti(iAtot+nOrbA*(iUtot-1))+Two*GA_Arrays(ipTanti)%A(IW)
          ValCF = ValCF*SCLNEL
        end if
        do IX=1,NX
          IXABS = IX+NAES(ISYX)
          IW1 = KTUV(IUABS,IVABS,IXABS)-NTUVES(ISYM)
          IW2 = IA
          IW = IW1+NAS*(IW2-1)

          AUVX(IA,IU,IV,IX) = AUVX(IA,IU,IV,IX)+GA_Arrays(ipT)%A(IW)
          AUVX(IA,IX,IU,IX) = AUVX(IA,IX,IU,IX)-ValCF*Half
        end do
      end do
    end do
  end do

  call DGEMM_('T','N',NV*NX,NCHO,NA*NU,One,AUVX,NA*NU,Cho_Bra,NA*NU,One,Cho_KetD,NV*NX)
  call DGEMM_('N','N',NA*NU,NCHO,NV*NX,One,AUVX,NA*NU,Cho_Ket,NV*NX,One,Cho_BraD,NA*NU)

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    call deallocate_GA_array(ipT)
    if (do_CSF) call deallocate_GA_array(ipTanti)
  else
# endif
    call RHS_FREE(ipT)
    if (do_csf) call RHS_FREE(ipTanti)
# ifdef _MOLCAS_MPP_
  end if
# endif

  return

end subroutine OLagNS_RI_C

subroutine OLagNS_RI_D1(ISYI,ISYK,NA,NJ,NV,NX,AJVX,NAJVX,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KTU
  use caspt2_global, only: do_csf

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NJ, NV, NX, NAJVX, NCHO
  real(kind=wp), intent(_OUT_) :: AJVX(NV,NX,NA*NJ)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NJ,NCHO), Cho_Ket(NV,NX,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NJ,NCHO), Cho_KetD(NV,NX,NCHO)
  !logical Incore
  integer(kind=iwp) :: IOFFD(8,8)
  integer(kind=iwp) :: ISW, ISA, IASTA, IJSTA, IJ, IA, IX, IV

  ISYJ = ISYI
  ISYX = ISYK

  do ISW=1,NSYM
    IO = 0
    do ISA=1,NSYM
      IOFFD(ISA,ISW) = IO
      ISI = Mul(ISA,ISW)
      IO = IO+NSSH(ISA)*NISH(ISI)
    end do
  end do

  ISYA = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYX)
  ISYM = JSYM !!!
  if (NINDEP(ISYM,5) == 0) return
  NAS1 = NTU(ISYM)
  NAS = 2*NAS1
  NIS = NISUP(ISYM,5)
  NWD = NAS*NIS
  if (NWD == 0) return

  ! ---- D1

  !! Read the T-amplitude
  ICASE = 5
  nIN = nINDEP(iSym,iCase)
  nAS = nASup(iSym,iCase)
  if (nIN /= 0) then
    nIS = nISup(iSym,iCase)
    nVec = nAS*nIS
    if (nVec /= 0) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call RHS_ALLO(nAS,nIS,lg_V)
        call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
        ipT = Allocate_GA_Array(nAS*nIS,'ipT')
        call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A,nAS)
        if (do_csf) then
          call RHS_READ_C(lg_V,iCase,iSym,7)
          ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
          call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A,nAS)
        end if
        call RHS_FREE(lg_V)
        call GASYNC()
      else
#     endif
        call RHS_ALLO(nAS,nIS,ipT)
        call RHS_READ_C(ipT,iCase,iSym,iVecC2)
        if (do_csf) then
          call RHS_ALLO(nAS,nIS,ipTanti)
          call RHS_READ_C(ipTanti,iCase,iSym,7)
        end if
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end if
  end if

  NBXSZA = NSECBX
  NBXSZJ = NINABX

  do IASTA=1,NA,NBXSZA
    IAEND = min(IASTA-1+NBXSZA,NA)
    NASZ = IAEND-IASTA+1
    do IJSTA=1,NJ,NBXSZJ
      IJEND = min(IJSTA-1+NBXSZJ,NJ)
      NJSZ = IJEND-IJSTA+1

      IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      call DCopy_(NAJVX,[Zero],0,AJVX,1)

      nOrbA = nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)+nSsh(iSyA)
      IAJ = 0
      do IJ=IJSTA,IJEND
        iJtot = iJ+nFro(iSyJ)
        do IA=IASTA,IAEND
          iAtot = iA+nFro(iSyA)+nIsh(iSyA)+nAsh(iSyA)
          IAJ = IAJ+1
          do IX=1,NX
            IXABS = IX+NAES(ISYX)
            iXtot = iX+nFro(iSyX)+nIsh(iSyX)
            do IV=1,NV
              IVABS = IV+NAES(ISYV)
              iVtot = iV+nFro(iSyV)+nIsh(iSyV)
              IW1 = KTU(IVABS,IXABS)-NTUES(ISYM)
              IW2 = IOFFD(ISYA,ISYM)+IJ+NJ*(IA-1)
              IW = IW1+NAS*(IW2-1)

              if (iVtot == iXtot) then
                DPT2C(iAtot+nOrbA*(iJtot-1)) = DPT2C(iAtot+nOrbA*(iJtot-1))+Two*GA_Arrays(ipT)%A(IW)
                if (do_csf) DPT2Canti(iAtot+nOrbA*(iJtot-1)) = DPT2Canti(iAtot+nOrbA*(iJtot-1))+Two*GA_Arrays(ipTanti)%A(IW)
              end if
              AJVX(IV,IX,IAJ) = GA_Arrays(ipT)%A(IW)
            end do
          end do
        end do
      end do

      call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NX,One,AJVX,NV*NX,Cho_Ket,NV*NX,One,Cho_BraD(IAJSTA,1,1),NA*NJ)
      call DGEMM_('N','N',NV*NX,NCHO,NASZ*NJSZ,One,AJVX,NV*NX,Cho_Bra(IAJSTA,1,1),NA*NJ,One,Cho_KetD,NV*NX)

    end do
  end do

# ifdef _MOLCAS_MPP_
  if (Is_Real_Par()) then
    call deallocate_GA_array(ipT)
    if (do_CSF) call deallocate_GA_array(ipTanti)
  else
# endif
    call RHS_FREE(ipT)
    if (do_csf) call RHS_FREE(ipTanti)
# ifdef _MOLCAS_MPP_
  end if
# endif

  return

end subroutine OLagNS_RI_D1

subroutine OLagNS_RI_D2(ISYI,ISYK,NA,NU,NV,NL,AUVL,NAUVL,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KTU

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NV, NL, NAUVL, NCHO
  real(kind=wp), intent(out) :: AUVL(NA,NU,NV,NL)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NV,NL,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO), Cho_KetD(NV,NL,NCHO)
  !logical Incore
  integer(kind=iwp) :: IOFFD(8,8)
  integer(kind=iwp) :: ISYW, ISYA, IA, IU, IV, IL

  ISYU = ISYI
  ISYL = ISYK

  do ISYW=1,NSYM
    IO = 0
    do ISYA=1,NSYM
      IOFFD(ISYA,ISYW) = IO
      ISYII = Mul(ISYA,ISYW)
      IO = IO+NSSH(ISYA)*NISH(ISYII)
    end do
  end do

  ISYA = Mul(JSYM,ISYU)
  ISYV = Mul(JSYM,ISYL)
  ISYM = Mul(ISYU,ISYV)
  if (NINDEP(ISYM,5) == 0) return
  NAS1 = NTU(ISYM)
  NAS = 2*NAS1
  NIS = NISUP(ISYM,5)
  NWD = NAS*NIS
  if (NWD == 0) return

  ! ---- D2

  !! Read the T-amplitude
  ICASE = 5
  nIN = nINDEP(iSym,iCase)
  nAS = nASup(iSym,iCase)
  if (nIN /= 0) then
    nIS = nISup(iSym,iCase)
    nVec = nAS*nIS
    if (nVec /= 0) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        call RHS_ALLO(nAS,nIS,lg_V)
        call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
        ipT = Allocate_GA_Array(nAS*nIS,'ipT')
        call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipT)%A,nAS)
        if (do_csf) then
          call RHS_READ_C(lg_V,iCase,iSym,7)
          ipTanti = Allocate_GA_Array(nAS*nIS,'ipTanti')
          call GA_GET(lg_V,1,nAS,1,nIS,GA_Arrays(ipTanti)%A,nAS)
        end if
        call RHS_FREE(lg_V)
        call GASYNC()
      else
#     endif
        call RHS_ALLO(nAS,nIS,ipT)
        call RHS_READ_C(ipT,iCase,iSym,iVecC2)
        if (do_csf) then
          call RHS_ALLO(nAS,nIS,ipTanti)
          call RHS_READ_C(ipTanti,iCase,iSym,7)
        end if
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end if
  end if

  call DCopy_(NAUVL,[Zero],0,AUVL,1)

  do IA=1,NA
    do IU=1,NU
      IUABS = IU+NAES(ISYU)
      do IV=1,NV
        IVABS = IV+NAES(ISYV)
        do IL=1,NL
          IW1 = NAS1+KTU(IVABS,IUABS)-NTUES(ISYM)
          IW2 = IOFFD(ISYA,ISYM)+IL+NL*(IA-1)
          IW = IW1+NAS*(IW2-1)
          AUVL(IA,IU,IV,IL) = GA_Arrays(ipT)%A(IW)
        end do
      end do
    end do
  end do

  call DGEMM_('T','N',NV*NL,NCHO,NA*NU,One,AUVL,NA*NU,Cho_Bra,NA*NU,One,Cho_KetD,NV*NL)
  call DGEMM_('N','N',NA*NU,NCHO,NV*NL,One,AUVL,NA*NU,Cho_Ket,NV*NL,One,Cho_BraD,NA*NU)

  if ((nIN /= 0) .and. (nVec /= 0)) then
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipT)
      if (do_CSF) call deallocate_GA_array(ipTanti)
    else
#   endif
      call RHS_FREE(ipT)
      if (do_csf) call RHS_FREE(ipTanti)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  return

end subroutine OLagNS_RI_D2

subroutine OLagNS_RI_E(ISYI,ISYK,NA,NJ,NV,NL,AJVL,NAJVL,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KIGTJ, KIGEJ

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NJ, NV, NL, NAJVL, NCHO
  real(kind=wp), intent(_OUT_) :: AJVL(NV,NL,NA*NJ)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NJ,NCHO), Cho_Ket(NV,NL,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NJ,NCHO), Cho_KetD(NV,NL,NCHO)
  !logical Incore
  integer(kind=iwp) :: IOFF1(8), IOFF2(8)
  integer(kind=iwp) :: ISA, IASTA, IJSTA, IJ, IA, IV, IL

  ISYJ = ISYI
  ISYL = ISYK

  ISYA = Mul(JSYM,ISYJ)
  ISYV = Mul(JSYM,ISYL)
  ISYM = ISYV
  ISYJL = Mul(ISYJ,ISYL)

  ! Set up offset table:
  IO1 = 0
  IO2 = 0
  do ISA=1,NSYM
    IOFF1(ISA) = IO1
    IOFF2(ISA) = IO2
    ISIJ = Mul(ISA,ISYM)
    IO1 = IO1+NSSH(ISA)*NIGEJ(ISIJ)
    IO2 = IO2+NSSH(ISA)*NIGTJ(ISIJ)
  end do

  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,6)
  NISM = NISUP(ISYM,7)
  !NIS = NISP+NISM
  NWEP = NAS*NISP
  NWEM = NAS*NISM
  NW = NWEP+NWEM
  if (NW == 0) return

  ! ---- EP

  if (NWEP > 0) then
    !! Read the T-amplitude
    ICASE = 6
    nINP = nINDEP(iSym,iCase)
    nASP = nASup(iSym,iCase)
    if (nINP /= 0) then
      nISP = nISup(iSym,iCase)
      nVec = nINP*nISP
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASP,nISP,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
          call GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A,nASP)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASP,nISP,ipTP)
          call RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    NBXSZA = NSECBX
    NBXSZJ = NINABX

    do IASTA=1,NA,NBXSZA
      IAEND = min(IASTA-1+NBXSZA,NA)
      NASZ = IAEND-IASTA+1
      do IJSTA=1,NJ,NBXSZJ
        IJEND = min(IJSTA-1+NBXSZJ,NJ)
        NJSZ = IJEND-IJSTA+1

        IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        call DCopy_(NAJVL,[Zero],0,AJVL,1)

        IAJ = 0
        do IJ=IJSTA,IJEND
          IJABS = IJ+NIES(ISYJ)
          do IA=IASTA,IAEND
            !IAABS = IA+NSES(ISYA)
            IAJ = IAJ+1

            do IV=1,NV
              !IVABS = IV+NAES(ISYV)
              do IL=1,NL
                ILABS = IL+NIES(ISYL)
                SCL = SQ05
                if (IJABS >= ILABS) then
                  JGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                  if (IJABS == ILABS) SCL = One
                else
                  JGEL = KIGEJ(ILABS,IJABS)-NIGEJES(ISYJL)
                end if
                IW1 = IV
                IW2 = IA+NA*(JGEL-1)+IOFF1(ISYA)
                IW = IW1+NAS*(IW2-1)

                AJVL(IV,IL,IAJ) = SCL*GA_Arrays(ipTP)%A(IW)
              end do
            end do
          end do
        end do

        call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NL,One,AJVL,NV*NL,Cho_Ket,NV*NL,One,Cho_BraD(IAJSTA,1,1),NA*NJ)
        call DGEMM_('N','N',NV*NL,NCHO,NASZ*NJSZ,One,AJVL,NV*NL,Cho_Bra(IAJSTA,1,1),NA*NJ,One,Cho_KetD,NV*NL)

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTP)
    else
#   endif
      call RHS_FREE(ipTP)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  ! ---- EM

  if (NWEM > 0) then
    !! Read the T-amplitude
    ICASE = 7
    nINM = nINDEP(iSym,iCase)
    nASM = nASup(iSym,iCase)
    if (nINM /= 0) then
      nISM = nISup(iSym,iCase)
      nVec = nINM*nISM
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASM,nISM,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
          call GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A,nASM)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASM,nISM,ipTM)
          call RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    NBXSZA = NSECBX
    NBXSZJ = NINABX

    do IASTA=1,NA,NBXSZA
      IAEND = min(IASTA-1+NBXSZA,NA)
      NASZ = IAEND-IASTA+1
      do IJSTA=1,NJ,NBXSZJ
        IJEND = min(IJSTA-1+NBXSZJ,NJ)
        NJSZ = IJEND-IJSTA+1

        IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
        call DCopy_(NAJVL,[Zero],0,AJVL,1)

        IAJ = 0
        do IJ=IJSTA,IJEND
          IJABS = IJ+NIES(ISYJ)
          do IA=IASTA,IAEND
            !IAABS = IA+NSES(ISYA)
            IAJ = IAJ+1

            do IV=1,NV
              !IVABS = IV+NAES(ISYV)
              do IL=1,NL
                ILABS = IL+NIES(ISYL)
                if (IJABS /= ILABS) then
                  if (IJABS > ILABS) then
                    SCL = SQ32
                    JGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                  else
                    SCL = -SQ32
                    JGTL = KIGTJ(ILABS,IJABS)-NIGTJES(ISYJL)
                  end if
                  IW1 = IV
                  IW2 = IA+NA*(JGTL-1)+IOFF2(ISYA)
                  IW = IW1+NAS*(IW2-1)

                  AJVL(IV,IL,IAJ) = SCL*GA_Arrays(ipTM)%A(IW)
                end if
              end do
            end do
          end do
        end do

        call DGEMM_('T','N',NASZ*NJSZ,NCHO,NV*NL,One,AJVL,NV*NL,Cho_Ket,NV*NL,One,Cho_BraD(IAJSTA,1,1),NA*NJ)
        call DGEMM_('N','N',NV*NL,NCHO,NASZ*NJSZ,One,AJVL,NV*NL,Cho_Bra(IAJSTA,1,1),NA*NJ,One,Cho_KetD,NV*NL)

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTM)
    else
#   endif
      call RHS_FREE(ipTM)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  return

end subroutine OLagNS_RI_E

subroutine OLagNS_RI_F(ISYI,ISYK,NA,NU,NC,NX,AUCX,NAUCX,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KTGTU, KTGEU, KAGTB, KAGEB

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NC, NX, NAUCX, NCHO
  real(kind=wp), intent(out) :: AUCX(NA,NU,NC,NX)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NX,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO), Cho_KetD(NC,NX,NCHO)
  integer(kind=iwp) :: IU, IX, IA, IC

  ISYU = ISYI
  ISYX = ISYK

  if (ISYU < ISYX) return

  ISYA = Mul(JSYM,ISYU)
  ISYC = Mul(JSYM,ISYX)
  ISYM = Mul(ISYU,ISYX) !!

  if (NINDEP(ISYM,8) > 0) then
    ! The plus combination:
    NASP = NTGEU(ISYM)
    NISP = NAGEB(ISYM)
    NWFP = NASP*NISP
  else
    NWFP = 0
  end if
  if (NINDEP(ISYM,9) > 0) then
    ICASE = 9
    ! The minus combination:
    NASM = NTGTU(ISYM)
    NISM = NAGTB(ISYM)
    NWFM = NASM*NISM
  else
    NWFM = 0
  end if
  if (NWFP+NWFM <= 0) return

  call DCopy_(NAUCX,[Zero],0,AUCX,1)

  ! ---- FP

  if ((NWFP > 0) .and. (NINDEP(ISYM,8) > 0)) then
    !! Read the T-amplitude
    ICASE = 8
    nINP = nINDEP(iSym,iCase)
    nASP = nASup(iSym,iCase)
    if (nINP /= 0) then
      nISP = nISup(iSym,iCase)
      nVec = nINP*nISP
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASP,nISP,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
          call GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A,nASP)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASP,nISP,ipTP)
          call RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    do IU=1,NU
      IUABS = IU+NAES(ISYU)
      IXMAX = NX
      if (ISYU == ISYX) IXMAX = IU
      do IX=1,IXMAX
        IXABS = IX+NAES(ISYX)
        SCL1 = Half
        if (IUABS == IXABS) SCL1 = Quart
        IW1 = KTGEU(IUABS,IXABS)-NTGEUES(ISYM)
        do IA=1,NA
          IAABS = IA+NSES(ISYA)
          do IC=1,NC
            ICABS = IC+NSES(ISYC)
            SCL = SCL1
            if (IAABS >= ICABS) then
              IW2 = KAGEB(IAABS,ICABS)-NAGEBES(ISYM)
              if (IAABS == ICABS) SCL = SQ2*SCL1
            else
              IW2 = KAGEB(ICABS,IAABS)-NAGEBES(ISYM)
            end if
            IW = IW1+NASP*(IW2-1)
            AUCX(IA,IU,IC,IX) = SCL*GA_Arrays(ipTP)%A(IW)
          end do
        end do
      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTP)
    else
#   endif
      call RHS_FREE(ipTP)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  ! ---- FM

  if ((NWFM > 0) .and. (NINDEP(ISYM,9) > 0)) then
    !! Read the T-amplitude
    ICASE = 9
    nINM = nINDEP(iSym,iCase)
    nASM = nASup(iSym,iCase)
    if (nINM /= 0) then
      nISM = nISup(iSym,iCase)
      nVec = nINM*nISM
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASM,nISM,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
          call GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A,nASM)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASM,nISM,ipTM)
          call RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    do IU=1,NU
      IUABS = IU+NAES(ISYU)
      IXMAX = NX
      if (ISYU == ISYX) IXMAX = IU-1
      do IX=1,IXMAX
        IXABS = IX+NAES(ISYX)
        IW1 = KTGTU(IUABS,IXABS)-NTGTUES(ISYM)
        do IA=1,NA
          IAABS = IA+NSES(ISYA)
          do IC=1,NC
            ICABS = IC+NSES(ISYC)
            if (IAABS > ICABS) then
              IW2 = KAGTB(IAABS,ICABS)-NAGTBES(ISYM)
              SCL = -Half
            else if (IAABS < ICABS) then
              IW2 = KAGTB(ICABS,IAABS)-NAGTBES(ISYM)
              SCL = Half
            else
              cycle
            end if
            IW = IW1+NASM*(IW2-1)
            AUCX(IA,IU,IC,IX) = AUCX(IA,IU,IC,IX)+SCL*GA_Arrays(ipTM)%A(IW)
          end do
        end do
      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTM)
    else
#   endif
      call RHS_FREE(ipTM)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  call DGEMM_('T','N',NC*NX,NCHO,NA*NU,One,AUCX,NA*NU,Cho_Bra,NA*NU,One,Cho_KetD,NC*NX)
  call DGEMM_('N','N',NA*NU,NCHO,NC*NX,One,AUCX,NA*NU,Cho_Ket,NC*NX,One,Cho_BraD,NA*NU)

  return

end subroutine OLagNS_RI_F

subroutine OLagNS_RI_G(ISYI,ISYK,NA,NU,NC,NL,AUCL,NAUCL,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NU, NC, NL, NAUCL, NCHO
  real(kind=wp), intent(_OUT_) :: AUCL(NA,NU,NC*NL)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NU,NCHO), Cho_Ket(NC,NL,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NU,NCHO), Cho_KetD(NC,NL,NCHO)
  !logical Incore
  integer(kind=iwp) :: IOFF1(8), IOFF2(8)
  integer(kind=iwp) :: ISI, ICSTA, ILSTA, IL, IC, IA, IU

  ISYU = ISYI
  ISYL = ISYK

  ISYA = Mul(JSYM,ISYU)
  ISYC = Mul(JSYM,ISYL)
  ISYM = ISYU
  ISYAC = Mul(ISYA,ISYC)
  ! Set up offset table:
  IO1 = 0
  IO2 = 0
  do ISI=1,NSYM
    IOFF1(ISI) = IO1
    IOFF2(ISI) = IO2
    ISAB = Mul(ISI,ISYM)
    IO1 = IO1+NISH(ISI)*NAGEB(ISAB)
    IO2 = IO2+NISH(ISI)*NAGTB(ISAB)
  end do

  ! Allocate W with parts WP,WM
  NAS = NASH(ISYM)
  NISP = NISUP(ISYM,10)
  NISM = NISUP(ISYM,11)
  !NIS = NISP+NISM
  NWGP = NAS*NISP
  NWGM = NAS*NISM
  !NWG = NWGP+NWGM

  !LDGP = NAS
  !LDGM = NAS

  ! ---- GP

  if (NWGP > 0) then
    !! Read the T-amplitude
    ICASE = 10
    nINP = nINDEP(iSym,iCase)
    nASP = nASup(iSym,iCase)
    if (nINP /= 0) then
      nISP = nISup(iSym,iCase)
      nVec = nINP*nISP
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASP,nISP,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
          call GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A,nASP)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASP,nISP,ipTP)
          call RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    NBXSZC = NSECBX
    KCL = NAUCL/(NA*NU)
    NBXSZL = KCL/NC
    if (NBXSZL <= 0) then
      write(u6,*) 'Not enough memory in ADDRHSG, I give up'
      call Abend()
    end if

    do ICSTA=1,NC,NBXSZC
      ICEND = min(ICSTA-1+NBXSZC,NC)
      NCSZ = ICEND-ICSTA+1
      do ILSTA=1,NL,NBXSZL
        ILEND = min(ILSTA-1+NBXSZL,NL)
        !NLSZ = ILEND-ILSTA+1

        ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
        call DCopy_(NAUCL,[Zero],0,AUCL,1)

        ICL = 0
        !IBUF = 0
        do IL=ILSTA,ILEND
          !ILABS = IL+NIES(ISYL)
          do IC=ICSTA,ICEND
            ICABS = IC+NSES(ISYC)
            ICL = ICL+1

            do IA=1,NA
              IAABS = IA+NSES(ISYA)
              SCL = SQ05
              if (IAABS >= ICABS) then
                IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                if (IAABS == ICABS) SCL = One
              else
                IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
              end if
              do IU=1,NU
                !IUABS = IU+NAES(ISYU)
                IW1 = IU
                IW2 = IL+NL*(IAGEC-1)+IOFF1(ISYL)
                IW = IW1+NAS*(IW2-1)

                AUCL(IA,IU,ICL) = SCL*GA_Arrays(ipTP)%A(IW)
              end do
            end do
          end do
        end do

        call DGEMM_('N','N',NA*NU,NCHO,ICL,One,AUCL,NA*NU,Cho_Ket(ICLSTA,1,1),NC*NL,One,Cho_BraD,NA*NU)
        call DGEMM_('T','N',ICL,NCHO,NA*NU,One,AUCL,NA*NU,Cho_Bra,NA*NU,One,Cho_KetD(ICLSTA,1,1),NC*NL)

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTP)
    else
#   endif
      call RHS_FREE(ipTP)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  ! ---- GM

  if (NWGM > 0) then
    !! Read the T-amplitude
    ICASE = 11
    nINM = nINDEP(iSym,iCase)
    nASM = nASup(iSym,iCase)
    if (nINM /= 0) then
      nISM = nISup(iSym,iCase)
      nVec = nINM*nISM
      if (nVec /= 0) then
#       ifdef _MOLCAS_MPP_
        if (Is_Real_Par()) then
          ! copy global array to local buffer
          call RHS_ALLO(nASM,nISM,lg_V)
          call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
          ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
          call GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A,nASM)
          call RHS_FREE(lg_V)
          call GASYNC()
        else
#       endif
          call RHS_ALLO(nASM,nISM,ipTM)
          call RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#       ifdef _MOLCAS_MPP_
        end if
#       endif
      end if
    end if

    NBXSZC = NSECBX
    KCL = NAUCL/(NA*NU)
    NBXSZL = KCL/NC
    if (NBXSZL <= 0) then
      write(u6,*) 'Not enough memory in ADDRHSG, I give up'
      call Abend()
    end if

    do ICSTA=1,NC,NBXSZC
      ICEND = min(ICSTA-1+NBXSZC,NC)
      NCSZ = ICEND-ICSTA+1
      do ILSTA=1,NL,NBXSZL
        ILEND = min(ILSTA-1+NBXSZL,NL)
        !NLSZ = ILEND-ILSTA+1

        ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)
        call DCopy_(NAUCL,[Zero],0,AUCL,1)

        ICL = 0
        !IBUF = 0
        do IL=ILSTA,ILEND
          !ILABS = IL+NIES(ISYL)
          do IC=ICSTA,ICEND
            ICABS = IC+NSES(ISYC)
            ICL = ICL+1

            do IA=1,NA
              IAABS = IA+NSES(ISYA)
              if (IAABS > ICABS) then
                IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                SCL = SQ32
              else if (IAABS < ICABS) then
                IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                SCL = -SQ32
              else
                cycle
              end if
              do IU=1,NU
                !IUABS = IU+NAES(ISYU)
                IW1 = IU
                IW2 = IL+NL*(IAGTC-1)+IOFF2(ISYL)
                IW = IW1+NAS*(IW2-1)

                AUCL(IA,IU,ICL) = SCL*GA_Arrays(ipTM)%A(IW)
              end do
            end do
          end do
        end do

        call DGEMM_('N','N',NA*NU,NCHO,ICL,One,AUCL,NA*NU,Cho_Ket(ICLSTA,1,1),NC*NL,One,Cho_BraD,NA*NU)
        call DGEMM_('T','N',ICL,NCHO,NA*NU,One,AUCL,NA*NU,Cho_Bra,NA*NU,One,Cho_KetD(ICLSTA,1,1),NC*NL)

      end do
    end do

#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTM)
    else
#   endif
      call RHS_FREE(ipTM)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  return

end subroutine OLagNS_RI_G

!! ADDRHSH
subroutine OLagNS_RI_H(ISYI,ISYK,NA,NJ,NC,NL,AJCL,NAJCL,Cho_Bra,Cho_Ket,Cho_BraD,Cho_KetD,NCHO)

  use SUPERINDEX, only: KAGEB, KAGTB, KIGEJ, KIGTJ

  integer(kind=iwp), intent(in) :: ISYI, ISYK, NA, NJ, NC, NL, NAJCL, NCHO
  real(kind=wp), intent(_OUT_) :: AJCL(NC*NL,NA*NJ)
  real(kind=wp), intent(in) :: Cho_Bra(NA,NJ,NCHO), Cho_Ket(NC,NL,NCHO)
  real(kind=wp), intent(inout) :: Cho_BraD(NA,NJ,NCHO), Cho_KetD(NC,NL,NCHO)
  integer(kind=iwp) :: IASTA, IJSTA, ICSTA, ILSTA, IJ, IA, IL, IC

  ISYJ = ISYI
  ISYL = ISYK

  if (ISYJ < ISYL) return
  ISYA = Mul(JSYM,ISYJ)
  ISYC = Mul(JSYM,ISYL)
  ISYM = Mul(ISYA,ISYC)
  ISYAC = ISYM
  ISYJL = ISYM

  NASP = NAGEB(ISYM)
  NISP = NIGEJ(ISYM)
  NWHP = NASP*NISP
  if (NWHP == 0) return
  NASM = NAGTB(ISYM)
  NISM = NIGTJ(ISYM)
  NWHM = NASM*NISM

  !LDHP = NASP
  !LDHM = NASM

  ! ---- HP

  !! Read the T-amplitude
  ICASE = 12
  nINP = nINDEP(iSym,iCase)
  nASP = nASup(iSym,iCase)
  nVec = 0
  if (nINP /= 0) then
    nISP = nISup(iSym,iCase)
    nVec = nINP*nISP
    if (nVec /= 0) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        ! copy global array to local buffer
        call RHS_ALLO(nASP,nISP,lg_V)
        call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
        ipTP = Allocate_GA_Array(nASP*nISP,'ipTP')
        call GA_GET(lg_V,1,nASP,1,nISP,GA_Arrays(ipTP)%A,nASP)
        call RHS_FREE(lg_V)
        call GASYNC()
      else
#     endif
        call RHS_ALLO(nASP,nISP,ipTP)
        call RHS_READ_C(ipTP,iCase,iSym,iVecC2)
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end if
  end if

  NBXSZA = NSECBX
  KAJ = NAJCL/(NC*NL)
  NBXSZJ = KAJ/NA
  if (NBXSZJ <= 0) then
    write(u6,*) 'Not enough memory in ADDRHSH, I give up'
    call Abend()
  end if

  NBXSZC = NSECBX
  NBXSZL = NINABX

  do IASTA=1,NA,NBXSZA
    IAEND = min(IASTA-1+NBXSZA,NA)
    NASZ = IAEND-IASTA+1
    do IJSTA=1,NJ,NBXSZJ
      IJEND = min(IJSTA-1+NBXSZJ,NJ)
      NJSZ = IJEND-IJSTA+1

      IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      call DCopy_(NAJCL,[Zero],0,AJCL,1)

      do ICSTA=1,NC,NBXSZC
        ICEND = min(ICSTA-1+NBXSZC,NC)
        NCSZ = ICEND-ICSTA+1
        do ILSTA=1,NL,NBXSZL
          ILEND = min(ILSTA-1+NBXSZL,NL)
          !NLSZ = ILEND-ILSTA+1

          ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

          IAJ = 0
          !IBUF = 0
          do IJ=IJSTA,IJEND
            IJABS = IJ+NIES(ISYJ)
            ILMAX = NL
            if (ISYJ == ISYL) ILMAX = IJ
            do IA=IASTA,IAEND
              IAABS = IA+NSES(ISYA)
              IAJ = IAJ+1

              ICL = 0
              do IL=ILSTA,min(ILEND,ILMAX)
                ILABS = IL+NIES(ISYL)
                SCL1 = One
                IJGEL = KIGEJ(IJABS,ILABS)-NIGEJES(ISYJL)
                if (IJABS == ILABS) SCL1 = SQ05
                do IC=ICSTA,ICEND
                  ICABS = IC+NSES(ISYC)
                  ICL = ICL+1

                  SCL = SCL1
                  if (IAABS >= ICABS) then
                    IAGEC = KAGEB(IAABS,ICABS)-NAGEBES(ISYAC)
                    if (IAABS == ICABS) SCL = SQ2*SCL1
                  else
                    IAGEC = KAGEB(ICABS,IAABS)-NAGEBES(ISYAC)
                  end if
                  IW = IAGEC+NAGEB(ISYM)*(IJGEL-1)

                  AJCL(ICLSTA+ICL-1,IAJ) = SCL*GA_Arrays(ipTP)%A(IW)
                end do
              end do
            end do
          end do

        end do
      end do
      call DGEMM_('T','N',NASZ*NJSZ,NCHO,NC*NL,One,AJCL,NC*NL,Cho_Ket,NC*NL,One,Cho_BraD(IAJSTA,1,1),NA*NJ)
      call DGEMM_('N','N',NC*NL,NCHO,NASZ*NJSZ,One,AJCL,NC*NL,Cho_Ket(IAJSTA,1,1),NC*NL,One,Cho_BraD,NA*NJ)

    end do
  end do

  if (nVec /= 0) then
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTP)
    else
#   endif
      call RHS_FREE(ipTP)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  ! ---- HM

  !! Read the T-amplitude
  ICASE = 13
  nINM = nINDEP(iSym,iCase)
  nASM = nASup(iSym,iCase)
  nVec = 0
  if (nINM /= 0) then
    nISM = nISup(iSym,iCase)
    nVec = nASM*nISM
    if (nVec /= 0) then
#     ifdef _MOLCAS_MPP_
      if (Is_Real_Par()) then
        ! copy global array to local buffer
        call RHS_ALLO(nASM,nISM,lg_V)
        call RHS_READ_C(lg_V,iCase,iSym,iVecC2)
        ipTM = Allocate_GA_Array(nASM*nISM,'ipTM')
        call GA_GET(lg_V,1,nASM,1,nISM,GA_Arrays(ipTM)%A,nASM)
        call RHS_FREE(lg_V)
        call GASYNC()
      else
#     endif
        call RHS_ALLO(nASM,nISM,ipTM)
        call RHS_READ_C(ipTM,iCase,iSym,iVecC2)
#     ifdef _MOLCAS_MPP_
      end if
#     endif
    end if
  end if

  NBXSZA = NSECBX
  KAJ = NAJCL/(NC*NL)
  NBXSZJ = KAJ/NA
  if (NBXSZJ <= 0) then
    write(u6,*) 'Not enough memory in ADDRHSH, I give up'
    call Abend()
  end if
  NBXSZC = NSECBX
  NBXSZL = NINABX

  do IASTA=1,NA,NBXSZA
    IAEND = min(IASTA-1+NBXSZA,NA)
    NASZ = IAEND-IASTA+1
    do IJSTA=1,NJ,NBXSZJ
      IJEND = min(IJSTA-1+NBXSZJ,NJ)
      NJSZ = IJEND-IJSTA+1

      IAJSTA = 1+NJ*(IASTA-1)+NASZ*(IJSTA-1)
      call DCopy_(NAJCL,[Zero],0,AJCL,1)

      do ICSTA=1,NC,NBXSZC
        ICEND = min(ICSTA-1+NBXSZC,NC)
        NCSZ = ICEND-ICSTA+1
        do ILSTA=1,NL,NBXSZL
          ILEND = min(ILSTA-1+NBXSZL,NL)
          !NLSZ = ILEND-ILSTA+1

          ICLSTA = 1+NL*(ICSTA-1)+NCSZ*(ILSTA-1)

          IAJ = 0
          !IBUF = 0
          do IJ=IJSTA,IJEND
            IJABS = IJ+NIES(ISYJ)
            ILMAX = NL
            if (ISYJ == ISYL) ILMAX = IJ-1
            do IA=IASTA,IAEND
              IAABS = IA+NSES(ISYA)
              IAJ = IAJ+1

              ICL = 0
              do IL=ILSTA,min(ILMAX,ILEND)
                ILABS = IL+NIES(ISYL)
                IJGTL = KIGTJ(IJABS,ILABS)-NIGTJES(ISYJL)
                do IC=ICSTA,ICEND
                  ICABS = IC+NSES(ISYC)
                  ICL = ICL+1

                  if (IAABS > ICABS) then
                    IAGTC = KAGTB(IAABS,ICABS)-NAGTBES(ISYAC)
                    SCL = SQ3
                  else if (IAABS < ICABS) then
                    IAGTC = KAGTB(ICABS,IAABS)-NAGTBES(ISYAC)
                    SCL = -SQ3
                  else
                    cycle
                  end if
                  IW = IAGTC+NAGTB(ISYM)*(IJGTL-1)

                  AJCL(ICLSTA+ICL-1,IAJ) = SCL*GA_Arrays(ipTM)%A(IW)
                end do
              end do
            end do
          end do

        end do
      end do
      call DGEMM_('T','N',NASZ*NJSZ,NCHO,NC*NL,One,AJCL,NC*NL,Cho_Ket,NC*NL,One,Cho_BraD(IAJSTA,1,1),NA*NJ)
      call DGEMM_('N','N',NC*NL,NCHO,NASZ*NJSZ,One,AJCL,NC*NL,Cho_Ket(IAJSTA,1,1),NC*NL,One,Cho_BraD,NA*NJ)
    end do
  end do

  if (nVec /= 0) then
#   ifdef _MOLCAS_MPP_
    if (Is_Real_Par()) then
      call deallocate_GA_array(ipTM)
    else
#   endif
      call RHS_FREE(ipTM)
#   ifdef _MOLCAS_MPP_
    end if
#   endif
  end if

  unused_var(cho_bra)
  unused_var(cho_ketd)

  return

end subroutine OLagNS_RI_H

subroutine Cnst_A_PT2(block1,block2)

  integer(kind=iwp), intent(in) :: block1, block2
  integer(kind=iwp) :: ndim1
# ifdef _MOLCAS_MPP_
  logical(kind=iwp) :: bstat
  integer(kind=iwp) :: ILOV1, IHIV1, JLOV1, JHIV1, MV1, iRank
# endif

  ndim1 = nSh(iSym0,block1)*nSh(iSym0,block2)
  if (ndim1 == 0) return

# ifdef _MOLCAS_MPP_
  if (is_real_par()) then
    !! The BRAD vector is first put from local BRA(:) to the
    !! global LG_V1 array. The data is shared among all processes.
    !! The bare Cholesky vector (KET(:)) is local.

    !! NV(I): Number of local Cholesky vectors for BRA(:)
    !! NDIM2: number of total Cholesky vectors for this batch
    !!        (usually NV*NPROCS)
    !! NVJ: number of local Cholesky vectors for KET(:)
    bStat = GA_CREATE_IRREG(MT_DBL,NDIM1,NDIM2,'BRAD',1,1,MAP2,NPROCS,LG_V1)
    call Cholesky_Vectors(2,block1,block2,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    !! finds out the range of the global array g_a that process
    !! iproc owns.
    call GA_DISTRIBUTION(LG_V1,MYRANK,ILOV1,IHIV1,JLOV1,JHIV1)
    !! provides access to local data in the specified patch of
    !! the array owned by the calling process
    call GA_ACCESS(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,MV1,NDIM1)
    call DCOPY_(NDIM1*(JHIV1-JLOV1+1),BRA,1,DBL_MB(MV1),1)
    !! releases access to a global array, when the data was
    !! accessed for writing
    call GA_RELEASE_UPDATE(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1)
    call GA_SYNC()

    !! ket is ndim1*NVJ dimension
    call Get_Cholesky_Vectors(block1,block2,JSYM,KET,size(KET),nKet,JBSTA,JBEND)

    do iRank=0,NPROCS-1
      call GA_DISTRIBUTION(LG_V1,iRank,ILOV1,IHIV1,JLOV1,JHIV1)
      call GA_GET(LG_V1,ILOV1,IHIV1,JLOV1,JHIV1,BRA,NDIM1)
      call DGEMM_('T','N',JHIV1-JLOV1+1,NVJ,ndim1,One,BRA,ndim1,KET,ndim1,One,A_PT2(IOFFCV+JLOV1-1,JOFFCV+MAP2(myRank+1)-1), &
                  MaxVec_PT2)
    end do
    bStat = GA_DESTROY(LG_V1)
    unused_var(bStat)
  else
# endif
    call Cholesky_Vectors(2,block1,block2,JSYM,BRA,size(BRA),nBra,IBSTA,IBEND)
    call Get_Cholesky_Vectors(block1,block2,JSYM,KET,size(KET),nKet,JBSTA,JBEND)
    call DGEMM_('T','N',NVI,NVJ,ndim1,One,BRA,ndim1,KET,ndim1,One,A_PT2(IOFFCV,JOFFCV),MaxVec_PT2)
# ifdef _MOLCAS_MPP_
  end if
# endif

end subroutine Cnst_A_PT2

end subroutine OLagNS_RI
