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

subroutine DEPSAOffC(NCONF,NSTATE,NASHT,NBAST,CLag,DEPSA,FIFA,FIMO,WRK1,WRK2,U0)

use Symmetry_Info, only: Mul
use caspt2_global, only: IPrGlb
use PrintLevel, only: VERBOSE
use caspt2_global, only: ConvInvar, SLag
use sguga, only: SGS, CIS
use caspt2_global, only: LUCIEX, IDCIEX, IDTCEX
use stdalloc, only: mma_allocate, mma_deallocate
use definitions, only: iwp, wp, u6
use Constants, only: Zero, One, Half
use caspt2_module, only: IFXMS, IFRMS, NSYM, STSYM, NFRO, NISH, NASH, NORB, NBAS, ISCF, NBTCH, NBTCHES, NROOTS
!use caspt2_module, only: NSSH
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif

implicit none
integer(kind=iwp), intent(in) :: NCONF, NSTATE, NASHT, NBAST
real(kind=wp), intent(inout) :: CLag(nConf,nState), DEPSA(nAshT,nAshT), WRK1(nBasT,nBasT), WRK2(nBasT**2)
real(kind=wp), intent(in) :: FIFA(nBasT**2), FIMO(nBasT**2), U0(nState,nState)
real(kind=wp), allocatable :: VecST(:,:), VecS1(:,:), VecS2(:,:), VecCID(:,:), VecPre(:), VecFancy(:), VecCIT(:,:), INT1(:), &
                              INT2(:), G2(:)
real(kind=wp), allocatable :: Eact(:)
real(kind=wp) :: Thres, DeltaC, Delta, Delta0, AlphaC, Alpha, ResCI, Beta, Res
real(kind=wp), external :: DDot_
integer(kind=iwp) :: nLev, nMidV, iSym, ID, iState, isyci, MaxIter, Iter
real(kind=wp) :: CPUT, WALLT, CPE, CPTF0, TIOE, TIOTF0
real(kind=wp) :: CPTF1, CPTF2, TIOTF1, TIOTF2

nLev = SGS%nLev
nMidV = CIS%nMidV

Thres = ConvInvar !! 1.0d-07

if (IPRGLB >= VERBOSE) then
  write(u6,*)
  write(u6,'(3X,"Linear Equation for Non-Invariant CASPT2 (threshold =",ES9.2,")")') Thres
  write(u6,*)
  call TIMING(CPTF0,CPE,TIOTF0,TIOE)
end if

! If CASPT2 energy is not invariant with respect to rotations within
! active space (with IPEA shift and/or with RAS reference), the
! active density obtained in constructing CI derivative is no longer
! correct... well, it may be correct, but orbital rotations in the
! active space cannot be parametrized in Z-vector, so analytic
! derivatives cannot be computed with the existing module. So,
! the active density is computed in a differnt way.

! See J. Chem. Phys. 2023, 158, 174112. for details, in particular,
! Section II C 4 "Non-invariance with respect to active MOs"
! To be more specific, this subroutine solves the linear equation
! (Eq. (71)) and computes the second term in Eq. (70) later.
! CLag corresponds to the RHS in Eq. (71).

!! Some post-processing of CI derivative
!! Somehow, this has to be done in the XMS basis
call CLagFinalOffC(nState,SLag)

! ----- Solve the linear equation -----
! A_{IS,JR}*X_{JR} = CLag_{IS}, where A_{IS,JR} is the CI-CI Hessian
! which may be seen in Z-vector

call mma_allocate(VecST,nConf,nState,Label='VecST')
call mma_allocate(VecS1,nConf,nState,Label='VecS1')
call mma_allocate(VecS2,nConf,nState,Label='VecS2')
call mma_allocate(VecCID,nConf,nState,Label='VecCID')
!call mma_allocate(VecS,nState*(nState-1)/2,Label='VecS')
call mma_allocate(VecPre,nConf,Label='VecPre')
call mma_allocate(VecFancy,nState**3,Label='VecFancy')

call mma_allocate(VecCIT,nConf,nState,Label='VecCIT')
call mma_allocate(INT1,nAshT**2,Label='INT1')
call mma_allocate(INT2,nAshT**4,Label='INT2')

call mma_allocate(Eact,nState,Label='Eact')

!! We do not have Cholesky vectors for frozen orbitals,
!! so may be it is not possible to get inactive energies?
!! It can be computed with TimesE2
iSym = 1
call CnstInt(0,nAshT,INT1,INT2)
do iState=1,nState
  ID = IDTCEX(iState)
  if (ISCF == 0) then
    !! quasi-canonical, XMS
    call DDaFile(LUCIEX,2,VecCIT(1,iState),nConf,ID)
  else
    VecCIT(1,iState) = One
  end if
  !! The second term should be removed
  Eact(iState) = Zero
end do
if (ifxms .or. ifrms) then
  !! Transform the CLag and CI vector from XMS to SCF basis
  !! Maybe, in order to define Eact
  call DGEMM_('N','T',nConf,nState,nState,One,CLag,nConf,U0,nState,Zero,VecST,nConf)
  CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)
  call DGEMM_('N','T',nConf,nState,nState,One,VecCIT,nConf,U0,nState,Zero,VecST,nConf)
  VecCIT(1:nConf,1:nState) = VecST(1:nConf,1:nState)
end if
call TimesE2(0,nConf,nState,nAshT,VecCIT,VecS1,INT1,INT2)
do iState=1,nState
  !! scaling with nState is due to the division in TimesE2
  Eact(iState) = -Half*nState*DDot_(nConf,VecS1(1,iState),1,VecCIT(1,iState),1)
end do
isyci = 1

!! Precondition
call CnstInt(2,nAshT,INT1,INT2)
call CnstPrec(ISYCI,nConf,nRoots,NLEV,nMidV,VecPre,VecCIT,INT1,INT2,VecFancy)
call CnstInt(0,nAshT,INT1,INT2)

!! Begin!
VecST(1:nConf,1:nState) = CLag(1:nConf,1:nState)

!! z0 = M^{-1}*r0
VecS2(1:nConf,1:nState) = VecST(1:nConf,1:nState)
call DoPrec(nConf,nRoots,VecST,VecS2,VecS1,VecPre,VecFancy)
!! p0 = z0
VecCId(1:nConf,1:nState) = VecS2(1:nConf,1:nState)
MaxIter = 100
Iter = 1
iSym = 1
! jspin   = 0
! r^T dot z
! r (residue) = ipST
! z (prec. r) = ipS2
! p (...)     = ipCId
! x (solution)= ipCIT
! Ap          = ipS1
! r_{k}z_{k}  = ipST*ipS2 = deltaC
DeltaC = DDot_(nConf*nState,VecST,1,VecS2,1)
Delta = DeltaC
Delta0 = Delta

if (IPRGLB >= VERBOSE) write(u6,*) ' Iteration       Delta           Res(CI)          DeltaC'
VecCIT(1:nConf,1:nState) = Zero
if (Delta0 > abs(Thres)) then
  do Iter=1,MaxIter
    if (nConf == 1) then
      do iState=1,nState
        VecCIT(1,iState) = One
      end do
      exit
    end if
    !! Compute Ap
    !! ipS2 is used as a workind array
    call TimesE2(1,nConf,nState,nAshT,VecCId,VecS1,INT1,INT2)

    !! AlphaC = p^T*A*p
    AlphaC = DDot_(nConf*nState,VecS1,1,VecCId,1)
    !! Alpha = r^T*z / AlphaC
    Alpha = Delta/(AlphaC)
    ! new x of CI
    VecCIT(1:nConf,1:nState) = VecCIT(1:nConf,1:nState)+Alpha*VecCId(1:nConf,1:nState)
    ! new r of CI
    VecST(1:nConf,1:nState) = VecST(1:nConf,1:nState)-Alpha*VecS1(1:nConf,1:nState)
    ResCI = sqrt(DDot_(nConf*nState,VecST,1,VecST,1))
    !! z = M^{-1}*r
    VecS2(1:nConf,1:nState) = VecST(1:nConf,1:nState)
    call DoPrec(nConf,nRoots,VecST,VecS2,VecS1,VecPre,VecFancy)

    !! Append new vectors
    DeltaC = Ddot_(nConf*nState,VecST,1,VecS2,1)
    Beta = DeltaC/Delta
    Delta = DeltaC
    VecCId(1:nConf,1:nState) = Beta*VecCId(1:nConf,1:nState)+VecS2(1:nConf,1:nState)

    if (IPRGLB >= VERBOSE) write(u6,'(I7,4X,ES17.9,ES17.9,ES17.9)') iter,delta/delta0,resci,deltac

    Res = ResCI
    if (Res <= abs(Thres)) exit
  end do

  if (Iter == MaxIter+1) then
    write(u6,*) 'CI iteration for non-invariant CASPT2 did not converge...'
    call abend()
  end if
end if

if (IPRGLB >= VERBOSE) then
  call TIMING(CPTF1,CPE,TIOTF1,TIOE)
  CPUT = CPTF1-CPTF0
  WALLT = TIOTF1-TIOTF0
  write(u6,*)
  write(u6,'(3X,"Linear equation converged in ",I3," steps")') iter-1
  write(u6,'(3X,"CPU and wall time (in s) = ",2F8.2)') CPUT,WALLT
  write(u6,*)
end if

if (IFXMS .or. IFRMS) then
  !! Transform back the CLag from CAS to XMS
  call DGEMM_('N','N',NConf,nState,nState,One,CLag,nConf,U0,nState,Zero,VecST,nConf)
  CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)
end if

call mma_deallocate(VecS1)
call mma_deallocate(VecS2)
call mma_deallocate(VecCId)
!call mma_deallocate(VecS)
call mma_deallocate(VecPre)
call mma_deallocate(VecFancy)

! ----- Construct (a part of) the true active density -----
! Compute the second term in Eq. (70) = Eq. (72)
! The SCF, not XMS, basis is used

do iState=1,nState
  ID = IDCIEX(iState) !! idtcex?
  if (ISCF == 0) then
    if (IFXMS .or. IFRMS) then
      !! Use unrotated (SCF) CI vector
      call LoadCI_XMS('C',1,nConf,nState,VecST(1,iState),iState,U0)
    else
      call DDaFile(LUCIEX,2,VecST(1,iState),nConf,ID)
    end if
  else
    VecST(1,iState) = One
  end if
end do
call mma_allocate(G2,nAshT**4,Label='G2')
call CnstInt(1,nAshT,INT1,INT2)
call CnstDEPSA(nConf,nState,nAshT,VecST,VecCIT,INT1,G2,INT2)
call mma_deallocate(G2)

if (IPRGLB >= VERBOSE) then
  call TIMING(CPTF2,CPE,TIOTF2,TIOE)
  CPUT = CPTF2-CPTF1
  WALLT = TIOTF2-TIOTF1
  write(u6,'(3X,A)') 'Off-diagonal density is constructed'
  write(u6,'(3X,"CPU and wall time (in s) = ",2F8.2)') CPUT,WALLT
  write(u6,*)
end if

call mma_deallocate(VecST)
call mma_deallocate(VecCIT)
call mma_deallocate(INT1)
call mma_deallocate(INT2)
call mma_deallocate(Eact)

contains

subroutine CLagFinalOffC(nState,SLag)

  use caspt2_module, only: REFENE

  integer(kind=iwp), intent(in) :: nState
  real(kind=wp), intent(inout) :: SLag(nState**2)
  real(kind=wp), allocatable :: CI1(:), CI2(:)
  integer(kind=iwp) :: ijst, ilStat, jlStat
  real(kind=wp) :: Scal, Ovl
  real(kind=wp), external :: DDot_

  ! Orthogonalize the partial derivative with respect to CI coeff

  call mma_allocate(VecST,nConf,nState,Label='VecST')
  call DGEMM_('N','T',nConf,nState,nState,One,CLag,nConf,U0,nState,Zero,VecST,nConf)
  CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)

  call mma_allocate(CI1,nConf,Label='CI1')
  call mma_allocate(CI2,nConf,Label='CI2')

  !! Construct SLag
  ijst = 0
  do ilStat=1,nState
    if (ISCF == 0) then
      call LoadCI_XMS('C',1,nConf,nState,CI1,ilStat,U0)
    else
      CI1(1) = One
    end if
    do jlStat=1,ilStat !! -1
      ijst = ilStat+nState*(jlStat-1)
      if (ilStat == jlStat) cycle
      if (ISCF == 0) then
        call LoadCI_XMS('C',1,nConf,nState,CI2,jlStat,U0)
      else
        CI2(1) = One
      end if
      Scal = DDOT_(nConf,CI1,1,CLag(1,jlStat),1)-DDOT_(nConf,CI2,1,CLag(1,ilStat),1)
      Scal = Scal/(REFENE(jlStat)-REFENE(ilStat))
      SLag(ijst) = SLag(ijst)+Scal
      if (IPRGLB >= VERBOSE) write(u6,'(1x,"SLag for State ",i1,"-",i1," = ",f20.10)') ilstat,jlstat,slag(ijst)
    end do
  end do

  !! Projection
  do ilStat=1,nState
    CI1(1:nConf) = CLag(1:nConf,ilStat)
    do jlStat=1,nState
      if (ISCF == 0) then
        call LoadCI_XMS('C',1,nConf,nState,CI2,jlStat,U0)
      else
        CI2(1) = One
      end if
      Ovl = DDot_(nConf,CI1,1,CI2,1)
      CLag(1:nConf,ilStat) = CLag(1:nConf,ilStat)-Ovl*CI2(1:nConf)
    end do
  end do

  call DGEMM_('N','N',nConf,nState,nState,One,CLag,nConf,U0,nState,Zero,VecST,nConf)
  CLag(1:nConf,1:nState) = VecST(1:nConf,1:nState)

  call mma_deallocate(CI1)
  call mma_deallocate(CI2)
  call mma_deallocate(VecST)

end subroutine CLagFinalOffC

subroutine CnstInt(Mode,nAshT,INT1,INT2)

  use CHOVEC_IO, only: NVLOC_CHOBATCH
  use caspt2_module, only: IfChol, NFRO, NBAS

  integer(kind=iwp), intent(in) :: Mode, nAshT
  real(kind=wp), intent(inout) :: INT1(nAshT,nAshT), INT2(nAshT,nAshT,nAshT,nAshT)
  integer(kind=iwp), allocatable :: BGRP(:,:)
  real(kind=wp), allocatable :: KET(:)
  integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3
  integer(kind=iwp) :: nFroI, nIshI, nCorI, nBasI, iAshI, jAshI, iSymA, iSymI, iSymB, iSymJ, JSYM, IB, IB1, IB2, MXBGRP, IBGRP, &
                       NBGRP, NCHOBUF, MXPIQK, NADDBUF, IBSTA, IBEND, NV, nKET, kAshI, lAshI, iT, iU, iTU, iV, iX, iVX, iOrb, jOrb
  !integer(kind=iwp) :: nSh(8,3)
  real(kind=wp) :: Val

  INT1(:,:) = Zero
  Int2(:,:,:,:) = Zero

  nFroI = nFro(iSym)
  nIshI = nIsh(iSym)
  nCorI = nFroI+nIshI
  nBasI = nBas(iSym)

  ! --- One-Electron Integral

  !! Read H_{\mu \nu}
  !IRC = -1
  !IOPT = 6
  !ICOMP = 1
  !ISYLBL = 1
  !call RDONE(IRC,IOPT,'OneHam  ',ICOMP,WRK2,ISYLBL)
  !! triangular -> square transformation
  !call Square(WRK2,WRK1,1,nBasT,nBasT)
  !! AO -> MO transformation
  !call DGemm_('T','N',nBasT,nBasT,nBasT,One,CMOPT2,nBasT,WRK1,nBasT,Zero,WRK2,nBasT)
  !call DGemm_('N','N',nBasT,nBasT,nBasT,One,WRK2,nBasT,CMOPT2,nBasT,Zero,WRK1,nBasT)
  !! Inactive energy
  !do iCorI=1,nFro(iSym)+nIsh(iSym)
  !  RIn_Ene = RIn_Ene+Two*WRK1(iCorI,iCorI)
  !end do
  !! Put in INT1
  !do iAshI=1,nAsh(iSym)
  !  do jAshI=1,nAsh(iSym)
  !    Val = WRK1(nCorI+iAshI,nCorI+jAshI)
  !    INT1(iAshI,jAshI) = INT1(iAshI,jAshI)+Val
  !  end do
  !end do
  do iAshI=1,nAsh(iSym)
    do jAshI=1,nAsh(iSym)
      Val = FIMO(nCorI+iAshI+nBasI*(nCorI+jAshI-1))
      INT1(iAshI,jAshI) = INT1(iAshI,jAshI)+Val
    end do
  end do

  ! --- Two-Electron Integral

  iSymA = 1
  iSymI = 1
  iSymB = 1
  iSymJ = 1
  !if (.not. IfChol) then
  !  do iCorI=1,nFro(iSym)+nIsh(iSym)
  !    iOrb = iCorI
  !    jOrb = iCorI
  !    call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
  !    do jCorI=1,nFro(iSym)+nIsh(iSym)
  !      RIn_Ene = RIn_Ene+Two*WRK1(jCorI,jCorI)
  !    end do
  !    call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
  !    do jCorI=1,nFro(iSym)+nIsh(iSym)
  !      RIn_Ene = RIn_Ene-WRK1(jCorI,jCorI)
  !    end do
  !  end do
  !end if

  if (IfChol) then
    !nSh(1:nSym,Inactive) = NISH(1:nSym)
    !nSh(1:nSym,Active) = NASH(1:nSym)
    !nSh(1:nSym,Virtual) = NSSH(1:nSym)
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

      call MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,NCHOBUF,MXPIQK,NADDBUF)
      call mma_allocate(KET,NCHOBUF,Label='KETBUF')
      !write(u6,*) 'nchobuf= ',nchobuf
      !write(u6,*) 'nbgrp= ',nbgrp
      !write(u6,*) 'nbtch= ',nbtch(jsym)
      do IBGRP=1,NBGRP

        IBSTA = BGRP(1,IBGRP)
        IBEND = BGRP(2,IBGRP)
        !write(u6,*) ibsta,ibend

        NV = 0
        do IB=IBSTA,IBEND
          NV = NV+NVLOC_CHOBATCH(IB)
        end do

        !! int2(tuvx) = (tu|vx)/2
        !! This can be computed without frozen orbitals
        call Get_Cholesky_Vectors(Active,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)

        call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,Half,KET,NASH(JSYM)**2,KET,NASH(JSYM)**2,Zero,INT2,NASH(JSYM)**2)
      end do
      call mma_deallocate(KET)
      call mma_deallocate(BGRP)
    end do
  else
    do iAshI=1,nAsh(iSym)
      iOrb = nCorI+iAshI
      do jAshI=1,nAsh(iSym)
        jOrb = nCorI+jAshI

        call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
        !! Put in INT1
        !do iCorI=1,nFro(iSym)+nIsh(iSym)
        !  INT1(iAshI,jAshI) = INT1(iAshI,jAshI)+Two*WRK1(iCorI,iCorI)
        !end do
        !! Put in INT2
        do kAshI=1,nAsh(iSym)
          do lAshI=1,nAsh(iSym)
            INT2(iAshI,jAshI,kAshI,lAshI) = INT2(iAshI,jAshI,kAshI,lAshI)+WRK1(nCorI+kAshI,nCorI+lAshI)*Half
          end do
        end do

        !call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
        !! Put in INT1
        !do iCorI=1,nFro(iSym)+nIsh(iSym)
        !  INT1(iAshI,jAshI) = INT1(iAshI,jAshI)-WRK1(iCorI,iCorI)
        !end do
      end do
    end do
  end if
# ifdef _MOLCAS_MPP_
  call GADGOP(INT2,nAshT**4,'+')
# endif
  !write(u6,*) 'int2'
  !call sqprt(int2,25)
  !call sqprt(int1,5)
  !call sqprt(fimo,12)
  if (Mode == 0) then
    do IT=1,nAshT
      do iU=1,nAshT
        iTU = iT+nAshT*(iU-1)
        do iV=1,nAshT
          do iX=1,nAshT
            iVX = iV+nAshT*(iX-1)
            if (iVX > iTU) then
              INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX)+INT2(iV,iX,iT,iU)
              INT2(iV,iX,iT,iU) = Zero
            end if
          end do
        end do
      end do
    end do
  end if

  if ((mode == 0) .or. (mode == 1)) then
    do IT=1,nAshT
      do iU=1,nAshT
        do iX=1,nAshT
          INT1(IT,IU) = INT1(IT,IU)-INT2(IT,IX,IX,IU)
        end do
      end do
    end do
  end if

  return

end subroutine CnstInt

!! dens2_rpt2
subroutine TimesE2(Mode,nConf,nState,nAshT,CIin,CIout,INT1,INT2)

  use sguga, only: L2ACT
  use Task_Manager, only: Init_Tsk, Free_Tsk, Rsv_Tsk
  use Constants, only: Two

  integer(kind=iwp), intent(in) :: Mode, nConf, nState, nAshT
  real(kind=wp), intent(in) :: CIin(nConf,nState), INT1(nAshT,nAshT), INT2(nAshT,nAshT,nAshT,nAshT)
  real(kind=wp), intent(out) :: CIout(nConf,nState)
  real(kind=wp), allocatable :: SGM1(:), SGM2(:)
  integer(kind=iwp), allocatable :: TASK(:,:)
  integer(kind=iwp) :: nLev, nTasks, iTask, LT, LU, kState, IST, IT, ISU, IU, ISTU, ISSG, NSGM, LVX, LV, ISV, IV, LX, ISX, ISVX, &
                       IX, ilStat, jlStat
  real(kind=wp) :: Ovl
  !logical tras,uras,vras,xras

  nLev = SGS%nLev

  !  --- H_{IJ}*P_J
  ! <CI1|EtuEvx|CI2>=<CI1|Evx

  nTasks = nLev**2
  call mma_allocate(Task,nTasks,2,Label='TASK')

  iTask = 0
  do LT=1,nLev
    do LU=1,nLev
      iTask = iTask+1
      TASK(iTask,1) = LT
      TASK(iTask,2) = LU
    end do
  end do
  if (iTask /= nTasks) write(u6,*) 'ERROR nTasks'

  call mma_allocate(SGM1,nConf,Label='SGM1')
  call mma_allocate(SGM2,nConf,Label='SGM2')

  CIout(1:nConf,1:nState) = Zero
  do kState=1,nState
    !! Start the actual part of dens2_rpt2
    call Init_Tsk(ID,nTasks)

    do while (Rsv_Tsk(ID,iTask))
      LT = TASK(iTask,1)
      !tras = lt <= nras1(1)
      IST = SGS%ISM(LT)
      IT = L2ACT(LT)
      LU = Task(iTask,2)
      !uras = lu > nras1(1)+nras2(1)
      !if (tras .and. uras) cycle
      !LTU = iTask
      ISU = SGS%ISM(LU)
      IU = L2ACT(LU)
      ISTU = Mul(IST,ISU)
      ISSG = Mul(ISTU,STSYM)
      NSGM = CIS%NCSF(ISSG)
      if (NSGM == 0) cycle
      !! <CIin|Etu
      call GETSGM2(LU,LT,STSYM,CIin(1,kState),nConf,SGM1,NSGM)
      !! <CIin|Etu|CIout>*I1tu
      if (ISTU == 1) CIout(1:NSGM,kState) = CIout(1:NSGM,kState)+INT1(IT,IU)*SGM1(1:NSGM)
      LVX = 0
      do LV=1,NLEV
        ISV = SGS%ISM(LV)
        IV = L2ACT(LV)
        !vras = lv <= nras1(1)
        do LX=1,NLEV
          LVX = LVX+1
          ISX = SGS%ISM(LX)
          ISVX = Mul(ISV,ISX)
          !xras = lx > nras1(1)+nras2(1)
          !if (vras .and. xras) cycle
          if (ISVX /= ISTU) cycle
          IX = L2ACT(LX)
          call GETSGM2(LX,LV,ISSG,SGM1,NSGM,SGM2,NSGM)
          CIout(1:NSGM,kState) = CIout(1:NSGM,kState)+INT2(IT,IU,IV,IX)*SGM2(1:NSGM)
        end do
      end do
    end do
    call Free_Tsk(ID)
    !! End the actual part of dens2_rpt2
  end do

  call mma_deallocate(Task)

# ifdef _MOLCAS_MPP_
  call GAdGOP(CIout,nConf*nState,'+')
# endif

  ! --- -E_{S}*CJ + zL_{KL}

  do kState=1,nState
    CIout(1:nConf,kState) = CIout(1:nConf,kState)+Eact(kState)*CIin(1:nConf,kState)
  end do

  !! Project out the reference vector, just in case
  if (Mode == 1) then
    do ilStat=1,nState
      SGM1(1:nConf) = CIout(1:nConf,ilStat)
      do jlStat=1,nState
        call LoadCI_XMS('C',1,nConf,nState,SGM2,jlStat,U0)
        Ovl = DDot_(nConf,SGM1,1,SGM2,1)
        CIout(1:nConf,ilStat) = CIout(1:nConf,ilStat)-Ovl*SGM2(1:nConf)
      end do
    end do
  end if

  call mma_deallocate(SGM1)
  call mma_deallocate(SGM2)

  CIout(1:nConf,1:nState) = Two/nState*CIout(1:nConf,1:nState)

  return

end subroutine TimesE2

subroutine CnstDEPSA(nConf,nState,nAshT,CI,CIT,G1,G2,INT2)

  use caspt2_module, only: MXCI, NG1, NG2

  integer(kind=iwp), intent(in) :: nConf, nState, nAshT
  real(kind=wp), intent(in) :: CI(nConf,nState), CIT(nConf,nState), INT2(nAshT,nAshT,nAshT,nAshT)
  real(kind=wp), intent(out) :: G1(nAshT,nAshT), G2(nAshT,nAshT,nAshT,nAshT)
  real(kind=wp), allocatable :: SGM1(:), SGM2(:), G1T(:), G2T(:), Fock(:), FockOut(:)
  integer(kind=iwp) :: nLev, kState, ilState, jlState, iS, jS, iA, jA, ip1, ip2, ipS, ijS, kS, lS, kAsh, kAA, lAsh, lAA, iAsh, &
                       ipQ, jAsh, ipM, imo, iOrb, jOrb
  real(kind=wp) :: Wgt, vSLag, rd, EigI, EigJ, OLagIJ, Tmp

  nLev = SGS%nLev

  ! This subroutine computes the second term in Eq. (70) or the RHS of
  ! Eq. (72) in the CASPT2-IPEA gradient paper
  ! CIT corresponds to \overline{Q}, if I remember correctly

  G1(:,:) = Zero
  G2(:,:,:,:) = Zero

  !! Construct transition(?) density matrix
  !! (<CI|Etu|CIT>+<CIT|Etu|CI>)/2, where CIT is the solution
  call mma_allocate(SGM1,MXCI,Label='SGM1')
  call mma_allocate(SGM2,MXCI,Label='SGM2')
  call mma_allocate(G1T,NG1,Label='GT1')
  call mma_allocate(G2T,NG2,Label='GT2')

   !! This is for CASSCF orbital Lagrangian, but this may not contribute
   !call Dens2T_RPT2(NLEV,NCONF,MXCI,CI(1,jState),CI(1,jState),SGM1,SGM2,G1T,G2T)
   !call DaXpY_(NG1,-Half,G1T,1,G1,1)
   !call DaXpY_(NG2,-Half,G2T,1,G2,1)

  do kState=1,nState
    !Wgt = DWgt(iState,iState)
    Wgt = One/nState

    !! <CI|Etu|CIT>+<CIT|Etu|CI> and the t+ u+ x v variant
    call Dens2T_RPT2(NLEV,NCONF,MXCI,CI(1,kState),CIT(1,kState),SGM1,SGM2,G1T,G2T)
    call DaXpY_(NG1,WGT,G1T,1,G1,1)
    call DaXpY_(NG2,WGT,G2T,1,G2,1)

    !! For the orbital contribution of CASSCF Lagrangian
    !! Just add the SLag rotation contributions
    ilState = kState
    do jlState=1,ilState-1
      vSLag = -Half*SLag(ilState,jlState)
      if (abs(vSLag) <= 1.0e-08_wp) cycle
      call Dens2T_RPT2(NLEV,NCONF,MXCI,CI(1,ilState),CI(1,jlState),SGM1,SGM2,G1T,G2T)
      call DaXpY_(NG1,vSLag,G1T,1,G1,1)
      call DaXpY_(NG2,vSLag,G2T,1,G2,1)
    end do
  end do

  call mma_deallocate(SGM1)
  call mma_deallocate(SGM2)
  call mma_deallocate(G1T)
  call mma_deallocate(G2T)

  !! Finally, construct the Fock matrix only for active-active
  !! Should be equivalent to FockGen in MCLR
  call mma_allocate(Fock,nAshT**2,Label='Fock')
  Fock(:) = Zero

  !! 1) FIMO term
  do iS=1,nSym
    if (nBas(iS) > 0) then
      jS = ieor(is-1,iSym-1)+1
      do iA=1,nAsh(is)
        do jA=1,nAsh(js)
          !rd = rDens1(iA+nA(iS),jA+nA(js))
          !ip1 = nBas(iS)*(nIsh(is)+iA-1)+ipCM(is)-1
          !ip2 = nBas(iS)*(nIsh(js)+jA-1)+ipmat(is,js)
          rd = G1(iA,jA)
          ip1 = 1+nFro(jS)+nIsh(jS)+nBas(iS)*(nFro(iS)+nIsh(iS)+iA-1)
          ip2 = 1+nAsh(iS)*(jA-1)
          call DaXpY_(nAsh(iS),Rd,FIMO(ip1),1,Fock(ip2),1)
        end do
      end do
    end if
  end do
  !write(u6,*) 'after 1'
  !call sqprt(fock,nasht)

  !! 2) two-electron term (only CreQADD part)
  do iS=1,nSym
    ipS = ieor(is-1,isym-1)+1
    if (norb(ips) /= 0) then
      do jS=1,nsym
        ijS = ieor(is-1,js-1)+1
        do kS=1,nSym
          ls = ieor(ijs-1,ieor(ks-1,isym-1))+1
          !                                                            *
          !*************************************************************
          !                                                            *
          do kAsh=1,nAsh(kS)
            kAA = kAsh+nFro(kS)+nIsh(kS)
            do lAsh=1,nAsh(lS)
              lAA = lAsh+nFro(lS)+nIsh(lS)

              ! Pick up (pj|kl)

              call Coul(ipS,jS,kS,lS,kAA,lAA,WRK1,WRK2)

              do iAsh=1,nAsh(iS)
                ipQ = nAsh(ipS)*(iAsh-1)
                do jAsh=1,nAsh(jS)
                  ipM = nFro(ipS)+nIsh(ipS)+(nFro(jS)+nIsh(jS)+jAsh-1)*nBas(ipS)
                  call DaXpY_(nAsh(ipS),G2(iAsh,jAsh,kAsh,lAsh)*2,INT2(1,jAsh,kAsh,lAsh),1,Fock(1+ipQ),1)
                  ipM = ipM+nOrb(ipS)

                end do
              end do

            end do
          end do
          !                                                            *
          !*************************************************************
          !                                                            *
        end do ! kS
      end do   ! jS
    end if
  end do       ! iS

  !! 3) anti-symmetrize
  !! 4) Divide by the difference of orbital energies
  call mma_allocate(FockOut,nAshT**2,Label='FockOut')
  do iS=1,nSym
    jS = ieor(iS-1,iSym-1)+1
    if (nAsh(is)*nAsh(jS) /= 0) then
      !! Anti-symmetrize
      call DGeSub(Fock,nAsh(iS),'N',Fock,nAsh(jS),'T',FockOut,nAsh(iS),nAsh(iS),nAsh(jS))

      !! Divide
      imo = 1
      do iAsh=1,nAsh(iSym)
        iOrb = iAsh+nFro(iSym)+nIsh(iSym)
        EigI = FIFA(iMO+iOrb-1+nBas(iSym)*(iOrb-1))
        do jAsh=1,iAsh-1
          jOrb = jAsh+nFro(iSym)+nIsh(iSym)
          EigJ = FIFA(iMO+jOrb-1+nBas(iSym)*(jOrb-1))
          OLagIJ = FockOut(iAsh+nAsh(iSym)*(jAsh-1))
          Tmp = OLagIJ/(EigI-EigJ)
          DEPSA(iAsh,jAsh) = DEPSA(iAsh,jAsh)+Tmp
          DEPSA(jAsh,iAsh) = DEPSA(jAsh,iAsh)+Tmp
        end do
      end do
    end if
  end do

  call mma_deallocate(FockOut)
  call mma_deallocate(Fock)

  return

end subroutine CnstDEPSA

!! PRWF1_CP2
subroutine CnstPrec(ISYCI,NCONF,NROOTS,NLEV,nMidV,PRE,CI,INT1,INT2,Fancy)

  use molcas, only: MXLEV
  use Constants, only: Two, Four

# include "intent.fh"

  integer(kind=iwp), intent(in) :: ISYCI, NCONF, NROOTS, nLev, nMidV
  real(kind=wp), intent(in) :: CI(nConf*nRoots), INT1(NLEV,NLEV), INT2(NLEV,NLEV,NLEV,NLEV)
  real(kind=wp), intent(_OUT_) :: PRE(nConf)
  real(kind=wp), intent(out) :: Fancy(nRoots,nRoots,nRoots)
  real(kind=wp) :: ICS(MXLEV), val, val2, Ene, dnum
  integer(kind=iwp) :: nIpWlk, LENCSF, ISY, LEV, MV, ISYUP, NCI, NUP, ISYDWN, NDWN, ICONF, IUW0, IDW0, IDWN, IUP, ICDPOS, ICDWN, &
                       NNN, IC1, ICUP, K, IDWNSV, ICUPOS, LEV2, iSt, jSt, kSt

  nIpWlk = CIS%nIpWlk

  ! Construct (approximate?) preconditioner for the active linear
  ! equation that should be solved for non-invariant CASPT2 methods
  ! (with IPEA shift)

  ! -- NOTE: THIS PRWF ROUTINE USES THE CONVENTION THAT CI BLOCKS
  ! -- ARE MATRICES CI(I,J), WHERE THE   F I R S T   INDEX I REFERS TO
  ! -- THE   U P P E R   PART OF THE WALK.

  ! SVC: set up a CSF string length as LENCSF
  !LINE = ''
  LENCSF = 0
  ISY = 0
  do LEV=1,NLEV
    if (ISY /= SGS%ISM(LEV)) then
      ISY = SGS%ISM(LEV)
      LENCSF = LENCSF+1
    end if
    LENCSF = LENCSF+1
  end do
  LENCSF = min(LENCSF,256)
  LENCSF = max(LENCSF,10)

  !LINE = ''

  ! -- THE MAIN LOOP IS OVER BLOCKS OF THE ARRAY CI
  !    WITH SPECIFIED MIDVERTEX MV, AND UPPERWALK SYMMETRY ISYUP.
  do MV=1,NMIDV
    do ISYUP=1,NSYM
      NCI = CIS%NOCSF(ISYUP,MV,ISYCI)
      if (NCI == 0) cycle
      NUP = CIS%NOW(1,ISYUP,MV)
      ISYDWN = Mul(ISYUP,ISYCI)
      NDWN = CIS%NOW(2,ISYDWN,MV)
      ICONF = CIS%IOCSF(ISYUP,MV,ISYCI)
      IUW0 = 1-NIPWLK+CIS%IOW(1,ISYUP,MV)
      IDW0 = 1-NIPWLK+CIS%IOW(2,ISYDWN,MV)
      IDWNSV = 0
      do IDWN=1,NDWN
        do IUP=1,NUP
          ICONF = ICONF+1
          !COEF = CI(ICONF)
          ! -- SKIP OR PRINT IT OUT?
          !if (abs(COEF) < THR) cycle
          if (IDWNSV /= IDWN) then
            ICDPOS = IDW0+IDWN*NIPWLK
            ICDWN = CIS%ICASE(ICDPOS)
            ! -- UNPACK LOWER WALK.
            NNN = 0
            do LEV=1,SGS%MIDLEV
              NNN = NNN+1
              if (NNN == 16) then
                NNN = 1
                ICDPOS = ICDPOS+1
                ICDWN = CIS%ICASE(ICDPOS)
              end if
              IC1 = ICDWN/4
              ICS(LEV) = ICDWN-4*IC1
              ICDWN = IC1
            end do
            IDWNSV = IDWN
          end if
          ICUPOS = IUW0+NIPWLK*IUP
          ICUP = CIS%ICASE(ICUPOS)
          ! -- UNPACK UPPER WALK:
          NNN = 0
          do LEV=SGS%MIDLEV+1,NLEV
            NNN = NNN+1
            if (NNN == 16) then
              NNN = 1
              ICUPOS = ICUPOS+1
              ICUP = CIS%ICASE(ICUPOS)
            end if
            IC1 = ICUP/4
            ICS(LEV) = ICUP-4*IC1
            ICUP = IC1
          end do
          ! -- PRINT IT!
          K = 0
          ISY = 0
          PRE(ICONF) = Zero
          do LEV=1,NLEV
            if (ISY /= SGS%ISM(LEV)) then
              ISY = SGS%ISM(LEV)
              K = K+1
              !LINE(K:K) = ' '
            end if
            K = K+1
            !LINE(K:K) = CODE(ICS(LEV))
            if (ICS(LEV) == 0) then
              VAL = Zero
            else if (ICS(LEV) == 3) then
              VAL = Two*INT1(LEV,LEV)
              !L = 0
              !JSY = 0
              do LEV2=1,NLEV
                if (ICS(LEV2) == 0) then
                else if (((LEV == LEV2) .and. (ICS(LEV2) == 3)) .or. ((LEV /= LEV2) .and. (ICS(LEV2) == 1)) .or. &
                         ((LEV /= LEV2) .and. (ICS(LEV2) == 2))) then
                  val2 = Four*int2(lev,lev,lev2,lev2)-Two*int2(lev,lev2,lev,lev2)
                  val = val+val2
                else if ((LEV /= LEV2) .and. (ICS(LEV2) == 3)) then
                  val2 = Four*int2(lev,lev,lev2,lev2)-Two*int2(lev,lev2,lev,lev2)
                  val = val+val2
                end if
              end do
            else
              VAL = INT1(LEV,LEV)
              !L = 0
              !JSY = 0
              do LEV2=1,NLEV
                if ((ICS(LEV2) == 0) .or. (LEV == LEV2)) then
                else if (ICS(LEV2) == 3) then
                  !val2 = Two*int2(lev,lev,lev2,lev2)-One*int2(lev,lev2,lev,lev2)
                  !val2 = Zero
                  !val = val+val2*Half
                else
                  val2 = int2(lev,lev,lev2,lev2)+int2(lev,lev2,lev,lev2)
                  if (ics(lev) == ics(lev2)) val2 = int2(lev,lev,lev2,lev2)-int2(lev,lev2,lev,lev2)
                  val = val+val2
                end if
              end do
            end if
            PRE(ICONF) = PRE(ICONF)+VAL
          end do
        end do
      end do
    end do
  end do

  !! mclr/sa_prec
  !! Prepare so-called fancy preconditioner
  do iSt=1,nRoots
    Ene = Eact(iSt)
    do jSt=1,nRoots
      do kSt=1,nRoots
        Fancy(jSt,kSt,iSt) = Zero
        do iConf=1,nConf
          dnum = PRE(iConf)+Ene
          dnum = sign(max(abs(dnum),1.0e-16_wp),dnum)
          Fancy(jSt,kSt,iSt) = Fancy(jSt,kSt,iSt)+CI(iConf+nConf*(jSt-1))*CI(iConf+nConf*(kSt-1))/dnum
        end do
      end do
    end do
    call MatInvert(Fancy(1,1,iSt),nRoots)
  end do

  return

end subroutine CnstPrec

!! mclr/dminvci_sa
subroutine DoPrec(nConf,nRoots,VecIN,VecOUT,CI,Pre,Fancy)
  ! Apply precondition to CI vectors, taken from the MCLR module

  integer(kind=iwp), intent(in) :: nConf, nRoots
  real(kind=wp), intent(in) :: VecIN(nConf,nRoots), Pre(nConf), Fancy(nRoots,nRoots,nRoots)
  real(kind=wp), intent(out) :: VecOUT(nConf,nRoots), CI(nConf,nRoots)
  integer(kind=iwp) :: iRoots, iConf, jRoots, kRoots
  real(kind=wp) :: rcoeff(nRoots), alpha(nRoots)

  !! Standard inverse of the diagonal elements
  do iRoots=1,nRoots
    do iConf=1,nConf
      VecOUT(iConf,iRoots) = VecIN(iConf,iRoots)/(Pre(iConf)+Eact(iRoots))
    end do
  end do

  !! Construct reference CI vectors
  do iRoots=1,nRoots
    call LoadCI_XMS('C',1,nConf,nState,CI(1,iRoots),iRoots,U0)
  end do

  !! The so-called fancy precondioner
  do iRoots=1,nRoots
    do jRoots=1,nRoots
      rcoeff(jRoots) = DDot_(nconf,VecOUT(1,iRoots),1,CI(1,jRoots),1)
    end do

    do jRoots=1,nRoots
      alpha(jRoots) = Zero
      do kRoots=1,nRoots
        alpha(jRoots) = alpha(jRoots)+Fancy(jRoots,kRoots,iRoots)*rcoeff(kRoots)
      end do
    end do

    do jRoots=1,nRoots
      do iConf=1,nConf
        VecOUT(iConf,iRoots) = VecOUT(iConf,iRoots)-CI(iConf,jRoots)*alpha(jRoots)/(Pre(iConf)+Eact(iRoots))
      end do
    end do
  end do

end subroutine DoPrec

end subroutine DEPSAOffC
