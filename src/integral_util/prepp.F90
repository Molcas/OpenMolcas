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
! Copyright (C) 1992, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine PrepP()
!***********************************************************************
!                                                                      *
! Object: to set up the handling of the 2nd order density matrix.      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             January '92                                              *
!***********************************************************************

use setup, only: mSkal, nSOs
use pso_stuff, only: lSA, Gamma_On, Gamma_MRCISD, lPSO, Case_2C, Case_3C, Case_MP2, nDens, mCMO, LuGam, FnGam, lBin, LuGamma, &
                     mDens, D0, DS, iD0Lbl, DSVar, KCMO, nG1, mG1, G1, nG2, mG2, G2, CMO, DVar, SO2CI, G_ToC, Bin
use iSD_data, only: iSO2Sh
use Basis_Info, only: nBas
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, Half, One
use stdalloc, only: mma_allocate, mma_deallocate
use etwas, only: nCMO, ExFac, CoulFac, nDSO, mIrrep, mBas, nAsh, nIsh
use NAC, only: IsNAC
use mspdft_grad, only: DoGradMSPD

implicit none
#include "dmrginfo_mclr.fh"
!#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
integer nFro(0:7)
integer Columbus
character(len=8) Method
character(len=80) KSDFT
logical DoCholesky
real*8 CoefX, CoefR
real*8, allocatable :: D1ao(:), D1AV(:), Tmp(:,:)
logical Do_Hybrid
real*8 WF_Ratio, PDFT_Ratio
#ifdef _DEBUGPRINT_
character(len=8) RlxLbl
character(len=60) Fmt
integer :: iComp = 1, ipTmp1
#endif
integer iIrrep, iSpin, iSeed, nShell, nPair, nQUad, LgToC, iDisk, n, nAct, iGo, nSA, ij, iBas, jBas, nTsT, nDim1, i, nDim0, nDim2
real*8, external :: Get_ExFac
integer, external :: IsFreeUnit

! Prologue

#ifdef _CD_TIMING_
call CWTIME(PreppCPU1,PreppWall1)
#endif

call StatusLine(' Alaska:',' Prepare the 2-particle matrix')

iD0Lbl = 1

lsa = .false.
Gamma_On = .false.
Gamma_mrcisd = .false.
lPSO = .false.
Case_2C = .false.
Case_3C = .false.
Case_mp2 = .false.

nDens = 0
do iIrrep=0,nIrrep-1
  nDens = nDens+nBas(iIrrep)*(nBas(iIrrep)+1)/2
end do

! Get the method label
call Get_cArray('Relax Method',Method,8)
call Get_iScalar('Columbus',columbus)
nCMo = S%n2Tot
mCMo = S%n2Tot
if ((Method == 'KS-DFT') .or. (Method == 'MCPDFT') .or. (Method == 'MSPDFT') .or. (Method == 'CASDFT')) then
  call Get_iScalar('Multiplicity',iSpin)
  call Get_cArray('DFT functional',KSDFT,80)
  call Get_dScalar('DFT exch coeff',CoefX)
  call Get_dScalar('DFT corr coeff',CoefR)
  ExFac = Get_ExFac(KSDFT)
  CoulFac = One
else
  iSpin = 0
  ExFac = One
  CoulFac = One
end if

! Check the wave function type

!                                                                      *
!***********************************************************************
!                                                                      *
if ((Method == 'RHF-SCF') .or. (Method == 'UHF-SCF') .or. (Method == 'IVO-SCF') .or. (Method == 'MBPT2') .or. &
    (Method == 'KS-DFT') .or. (Method == 'ROHF')) then
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,'(2A)') ' Wavefunction type: ',Method
  if (Method == 'KS-DFT  ') then
    write(6,'(2A)') ' Functional type:   ',KSDFT
    Fmt = '(1X,A26,20X,F18.6)'
    write(6,Fmt) 'Exchange scaling factor',CoefX
    write(6,Fmt) 'Correlation scaling factor',CoefR
  end if
  write(6,*)
# endif
  if (Method == 'MBPT2   ') then
    Case_mp2 = .true.
    call DecideOnCholesky(DoCholesky)
    if (.not. DoCholesky) then
      iSeed = 10
      LuGam = IsFreeUnit(iSeed)
      write(FnGam,'(A6)') 'LuGam'
      call DaName_MF_WA(LuGam,FnGam)
      Gamma_on = .true.
      call Aces_Gamma()
    end if
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if (Method == 'Corr. WF') then
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) ' Wavefunction type: an Aces 2 correlated wavefunction'
  write(6,*)
# endif
  Gamma_On = .true.
  call Aces_Gamma()
else if ((Method(1:7) == 'MR-CISD') .and. (Columbus == 1)) then
  !*********** columbus interface ****************************************
  !do not reconstruct the two-particle density from the one-particle
  !density or partial two-particle densities but simply read them from
  !file

  nShell = mSkal
  nPair = nShell*(nShell+1)/2
  nQuad = nPair*(nPair+1)/2

  ! Allocate Table Of Content for half sorted gammas.

  call mma_Allocate(G_Toc,nQuad+2,Label='G_Toc')

  ! find free unit number
  lgtoc = 61
  lgtoc = isFreeUnit(lgtoc)
  idisk = 0
  ! read table of contents of half-sorted gamma file
  call DaName(lgtoc,'gtoc')
  call ddafile(lgtoc,2,G_Toc,nQuad+2,idisk)
  call Daclos(lgtoc)
  n = int(G_Toc(nQuad+1))
  lbin = int(G_Toc(nQuad+2))
  if (n /= nQuad) then
    call WarningMessage(2,'n /= nQuad')
    write(6,*) 'n,nQuad=',n,nQuad
    call Abend()
  end if

  Gamma_On = .true.
  Gamma_mrcisd = .true.
  ! open gamma file
  LuGamma = 60
  LuGamma = isfreeunit(LuGamma)
  ! closed in closep
  call DaName_MF(LuGamma,'GAMMA')
  ! allocate space for bins
  call mma_Allocate(Bin,2,lBin,Label='Bin')
  ! compute SO2cI array
  call mma_Allocate(SO2cI,2,nSOs,Label='SO2cI')
  call Mk_SO2cI(SO2cI,iSO2Sh,nsos)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((Method == 'RASSCF') .or. (Method == 'CASSCF') .or. (Method == 'GASSCF') .or. (Method == 'MCPDFT') .or. &
         (Method == 'MSPDFT') .or. (Method == 'DMRGSCF') .or. (Method == 'CASDFT')) then

  call Get_iArray('nAsh',nAsh,nIrrep)
  nAct = 0
  do iIrrep=0,nIrrep-1
    nAct = nAct+nAsh(iIrrep)
  end do
  if (nAct > 0) lPSO = .true.

  nDSO = nDens
  mIrrep = nIrrep
  call ICopy(nIrrep,nBas,1,mBas,1)
# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,'(2A)') ' Wavefunction type: ',Method
  if ((Method == 'CASDFT') .or. (Method == 'MCPDFT')) write(6,'(2A)') ' Functional type:   ',KSDFT
  if (Method == 'MSPDFT') write(6,'(2A)') ' MS-PDFT Functional type:   ',KSDFT
  write(6,*)
# endif
  if (Method == 'MCPDFT') lSA = .true.
  if (Method == 'MSPDFT') then
    lSA = .true.
    dogradmspd = .true.
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else if ((Method == 'CASSCFSA') .or. (Method == 'DMRGSCFS') .or. (Method == 'GASSCFSA') .or. (Method == 'RASSCFSA') .or. &
         (Method == 'CASPT2')) then
  call Get_iArray('nAsh',nAsh,nIrrep)
  nAct = 0
  do iIrrep=0,nIrrep-1
    nAct = nAct+nAsh(iIrrep)
  end do
  if (nAct > 0) lPSO = .true.
  nDSO = nDens
  call Get_iScalar('SA ready',iGo)
  if (iGO == 1) lSA = .true.
  mIrrep = nIrrep
  call ICopy(nIrrep,nBas,1,mBas,1)
# ifdef _DEBUGPRINT_
  write(6,*)
  if (lSA) then
    write(6,'(2A)') ' Wavefunction type: State average ',Method(1:6)
  else
    write(6,'(2A)') ' Wavefunction type: ',Method
  end if
  write(6,*)
# endif
  if (Method == 'CASPT2  ') then
    call DecideOnCholesky(DoCholesky)
    !if (.not. DoCholesky) then
    Gamma_On = .true.
    ! It is opened, but not used actually. I just want to
    ! use the Gamma_On flag.
    ! Just to avoid error termination.
    ! Actual working arrys are allocated in drvg1.f
    LuGamma = 65
    call DaName_MF_WA(LuGamma,'GAMMA')
    if (DoCholesky) call mma_allocate(G_Toc,1,Label='G_Toc')
    call mma_allocate(SO2cI,1,1,Label='SO2cI')
    call mma_allocate(Bin,1,1,Label='Bin')
    !end if
  end if
  Method = 'RASSCF  '
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  call WarningMessage(2,'Alaska: Unknown wavefuntion type')
  write(6,*) 'Wavefunction type:',Method
  write(6,*) 'Illegal type of wave function!'
  write(6,*) 'ALASKA cannot continue.'
  call Quit_OnUserError()
end if

! Read the (non) variational 1st order density matrix
! density matrix in AO/SO basis
nsa = 1
if (lsa) nsa = 5
if (Method == 'MCPDFT  ') nsa = 5
if (Method == 'MSPDFT  ') nsa = 5
!AMS modification: add a fifth density slot
mDens = nsa+1
call mma_allocate(D0,nDens,mDens,Label='D0')
D0(:,:) = Zero
call mma_allocate(DVar,nDens,nsa,Label='DVar')
if (.not. gamma_mrcisd) call Get_dArray_chk('D1ao',D0(1,1),nDens)

call Get_D1ao_Var(DVar,nDens)

call mma_Allocate(DS,nDens,Label='DS')
call mma_Allocate(DSVar,nDens,Label='DSVar')
if ((Method == 'UHF-SCF') .or. (Method == 'ROHF') .or. ((Method == 'KS-DFT') .and. (iSpin /= 1)) .or. (Method == 'Corr. WF')) then
  call Get_dArray_chk('D1sao',DS,nDens)
  call Get_D1sao_Var(DSVar,nDens)
else
  DS(:) = Zero
  DSVar(:) = Zero
end if

! This is necessary for the ci-lag

! Unfold density matrix

!*********** columbus interface ****************************************
!do not modify the effective density matrices and fock matrices
if (.not. gamma_mrcisd) then
  ij = -1
  do iIrrep=0,nIrrep-1
    do iBas=1,nBas(iIrrep)
      do jBas=1,iBas-1
        ij = ij+1
        D0(1+ij,1) = Half*D0(1+ij,1)
        DVar(1+ij,1) = Half*DVar(1+ij,1)
        DS(1+ij) = Half*DS(1+ij)
        DSVar(1+ij) = Half*DSVar(1+ij)
      end do
      ij = ij+1
    end do
  end do
end if

#ifdef _DEBUGPRINT_
RlxLbl = 'D1AO    '
call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],D0)
RlxLbl = 'D1AO-Var'
call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DVar)
RlxLbl = 'DSAO    '
call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DS)
RlxLbl = 'DSAO-Var'
call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DSVar)
#endif

! Get the MO-coefficients
!*********** columbus interface ****************************************
! not that the columbus mcscf MO coefficients have been written
! to the RUNFILE !

if ((Method == 'UHF-SCF') .or. (Method == 'ROHF') .or. ((Method == 'KS-DFT') .and. (iSpin /= 1)) .or. (Method == 'Corr. WF')) then
  nsa = 2
else
  nsa = 1
  if (lsa) nsa = 2
end if
kCMO = nsa
call mma_allocate(CMO,mCMO,kCMO,Label='CMO')
call Get_dArray_chk('Last orbitals',CMO(:,1),mCMO)
#ifdef _DEBUGPRINT_
ipTmp1 = 1
do iIrrep=0,nIrrep-1
  call RecPrt(' CMO''s',' ',CMO(ipTmp1,1),nBas(iIrrep),nBas(iIrrep))
  ipTmp1 = ipTmp1+nBas(iIrrep)**2
end do
#endif

! Get additional information in the case of a RASSCF wave function
! Get the number of inactive, active and frozen orbitals
!*********** columbus interface ****************************************
! no need for MRCI gradient
if ((.not. lpso) .or. gamma_mrcisd) goto 1000
call Get_iScalar('nSym',i)
call Get_iArray('nIsh',nIsh,i)
call Get_iArray('nAsh',nAsh,i)
call Get_iArray('nFro',nFro,i)
#ifdef _DEBUGPRINT_
write(6,*) ' nISh=',nISh
write(6,*) ' nASh=',nASh
write(6,*) ' nFro=',nFro
#endif
nAct = 0
nTst = 0
do iIrrep=0,nIrrep-1
  !write(6,*)"nAsh(iIrrep)",nAsh(iIrrep)  ! yma
  nAct = nAct+nAsh(iIrrep)
  nTst = nTst+nFro(iIrrep)
end do
if (nTst /= 0) then
  call WarningMessage(2,'; No frozen orbitals are allowed!; ALASKA cannot continue;')
  call Quit_OnUserError()
end if

! Get the one body density for the active orbitals
! (not needed for SA-CASSCF)
nG1 = nAct*(nAct+1)/2
nsa = 1
if (lsa) nsa = 0
mG1 = nsa
call mma_allocate(G1,nG1,mG1,Label='G1')
if (nsa > 0) then
  call Get_dArray_chk('D1mo',G1(:,1),nG1)
#ifdef _DEBUGPRINT_
  call TriPrt(' G1',' ',G1(:,1),nAct)
#endif
end if

! Get the two body density for the active orbitals
nG2 = nG1*(nG1+1)/2
nsa = 1
if (lsa) nsa = 2
mG2 = nsa
call mma_allocate(G2,nG2,mG2,Label='G2')
!write(6,*) 'got the 2rdm, Ithink.'
if ((Method == 'MCPDFT') .or. (Method == 'MSPDFT')) then
  call Get_dArray_chk('P2MOt',G2,nG2)!PDFT-modified 2-RDM
else
  call Get_dArray_chk('P2mo',G2,nG2)
end if
#ifdef _DEBUGPRINT_
call TriPrt(' G2',' ',G2(1,1),nG1)
#endif
if (lsa) then

  ! CMO1 Ordinary CMOs

  ! CMO2 CMO*Kappa

  call Get_dArray_chk('LCMO',CMO(:,2),mCMO)
# ifdef _DEBUGPRINT_
  ipTmp1 = 1
  do iIrrep=0,nIrrep-1
    call RecPrt('LCMO''s',' ',CMO(ipTmp1,2),nBas(iIrrep),nBas(iIrrep))
    ipTmp1 = ipTmp1+nBas(iIrrep)**2
  end do
# endif

  ! P are stored as
  !                            _                     _
  !   P1=<i|e_pqrs|i> + sum_i <i|e_pqrs|i>+<i|e_pqrs|i>
  !   P2=sum_i <i|e_pqrs|i>

  call Get_dArray_chk('PLMO',G2(:,2),nG2)
  ndim1 = 0
  if (doDMRG) then
    ndim0 = 0  !yma
    do i=1,8
      ndim0 = ndim0+LRras2(i)
    end do
    ndim1 = (ndim0+1)*ndim0/2
    ndim2 = (ndim1+1)*ndim1/2
    do i=1,ng2
      if (i > ndim2) G2(i,2) = 0.0d0
    end do
  end if
  call Daxpy_(ng2,One,G2(:,2),1,G2(:,1),1)
# ifdef _DEBUGPRINT_
  call TriPrt(' G2L',' ',G2(:,2),nG1)
  call TriPrt(' G2T',' ',G2(:,1),nG1)
# endif

  call Get_dArray_chk('D2av',G2(:,2),nG2)
# ifdef _DEBUGPRINT_
  call TriPrt('G2A',' ',G2(:,2),nG1)
# endif

  ! Densities are stored as:
  !
  !  ipd0 AO:
  !
  !  D1 = inactive diagonal density matrix
  !                                _                 _
  !  D2 = <i|E_pq|i> + sum_i <i|E_pq|i>+<i|E_pq|i> + sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> - 1/2 D2
  !
  !  D3 = sum_i <i|E_pq|i> (active)
  !
  !  D4 = sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> (inactive)
  !
  !  G1 = <i|e_ab|i>
  !  G2 = sum i <i|e_ab|i>
  !
  !************************
  !RlxLbl = 'D1AO    '
  !call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0)

  call mma_allocate(Tmp,nDens,2,Label='Tmp')
  call Get_D1I(CMO(1,1),D0(1,1),Tmp,nIsh,nBas,nIrrep)
  call mma_deallocate(Tmp)

  !************************
  !RlxLbl = 'D1AO    '
  !call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0)

  call dcopy_(ndens,DVar,1,D0(1,2),1)
  if (.not. isNAC) call daxpy_(ndens,-Half,D0(1,1),1,D0(1,2),1)
  !RlxLbl = 'D1COMBO  '
  !call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0(1,2))

  ! This is necessary for the kap-lag

  nG1 = nAct*(nAct+1)/2
  call mma_allocate(D1AV,nG1,Label='D1AV')
  call Get_dArray_chk('D1av',D1AV,nG1)
  call Get_D1A(CMO(1,1),D1AV,D0(1,3),nIrrep,nbas,nish,nash,ndens)
  call mma_deallocate(D1AV)
  !************************
  !RlxLbl = 'D1AOA   '
  !call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0(1,3))

  call Get_dArray_chk('DLAO',D0(:,4),nDens)

  ! Getting conditions for hybrid MC-PDFT
  Do_Hybrid = .false.
  call qpg_DScalar('R_WF_HMC',Do_Hybrid)
  if (Do_Hybrid) then
    call Get_DScalar('R_WF_HMC',WF_Ratio)
    PDFT_Ratio = 1.0d0-WF_Ratio
  end if
  !ANDREW - modify D2: should contain only the correction pieces
  if (Method == 'MCPDFT  ') then
    ! Get the D_theta piece
    call mma_allocate(D1ao,nDens)
    call Get_dArray_chk('D1ao',D1ao,ndens)
    ij = 0
    do iIrrep=0,nIrrep-1
      do iBas=1,nBas(iIrrep)
        do jBas=1,iBas-1
          ij = ij+1
          D1ao(ij) = Half*D1ao(ij)
        end do
        ij = ij+1
      end do
    end do
    call daxpy_(ndens,-1d0,D0(1,1),1,D1ao,1)
    !write(6,*) 'do they match?'
    !do i=1,ndens
    !  write(6,*) d1ao(i),DO(i,3)
    !end do

    call daxpy_(ndens,-Half,D0(1,1),1,D0(1,2),1)
    call daxpy_(ndens,-1.0d0,D1ao,1,D0(1,2),1)
    !ANDREW -   Generate new D5 piece:
    D0(:,5) = Zero
    call daxpy_(ndens,0.5d0,D0(1,1),1,D0(1,5),1)
    call daxpy_(ndens,1.0d0,D1ao,1,D0(1,5),1)

    if (do_hybrid) then
      ! add back the wave function parts that are subtracted
      ! this might be inefficient, but should have a clear logic
      call daxpy_(ndens,Half*WF_Ratio,D0(1,1),1,D0(1,2),1)
      call daxpy_(ndens,WF_Ratio,D1ao,1,D0(1,2),1)
      ! scale the pdft part
      call dscal_(ndens,PDFT_Ratio,D0(1,5),1)
    end if
    call mma_deallocate(D1ao)
  else if (Method == 'MSPDFT  ') then
    call Get_DArray('MSPDFTD5        ',D0(1,5),nDens)
    call Get_DArray('MSPDFTD6        ',D0(1,6),nDens)
    call daxpy_(ndens,-1.0d0,D0(1,5),1,D0(1,2),1)
  end if

  !call dcopy_(ndens*5,0.0d0,0,D0,1)
  !call dcopy_(nG2,0.0d0,0,G2,1)

  !************************
  !call dscal_(Length,0.5d0,D0(1,4),1)
  !call dscal_(Length,0.0d0,D0(1,4),1)

  !RlxLbl = 'DLAO    '
  !call PrMtrx(RlxLbl,iD0Lbl,iComp,[1],D0(1,4))
  ! DMRG with the reduced AS
  !if (doDMRG) length = ndim1  !yma
end if
#ifdef _DEBUGPRINT_
call TriPrt(' G2',' ',G2(1,1),nG1)
#endif

! Close 'RELAX' file
1000 continue

! Epilogue, end
#ifdef _CD_TIMING_
call CWTIME(PreppCPU2,PreppWall2)
Prepp_CPU = PreppCPU2-PreppCPU1
Prepp_Wall = PreppWall2-PreppWall1
#endif

return

end subroutine PrepP
