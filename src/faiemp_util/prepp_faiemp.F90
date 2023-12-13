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
! Copyright (C) Ben Swerts                                             *
!***********************************************************************

subroutine PrepP_FAIEMP(nBas_Valence,nBT,nBVT)
!***********************************************************************
!                                                                      *
! Object: to set up the handling of the 2nd order density matrix for   *
!         the calculation of the 2-electron FAIEMP derivatives         *
!                                                                      *
!     Author: Ben Swerts                                               *
!                                                                      *
! Based on PrepP                                                       *
!***********************************************************************

use pso_stuff, only: CMO, D0, DVar, DS, DSVar, G1, G2, Gamma_On, id0Lbl, kCMO, lPSO, lsa, mCMo, mDens, mG1, mG2, nDens, nG1, nG2
use Basis_Info, only: nBas
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Etwas, only: CoulFac, ExFac, mBas, mIrrep, nASh, nCMO, nDSO, nISh
use NAC, only: isNAC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas_Valence(0:7), nBT, nBVT
#include "print.fh"
integer(kind=iwp) :: nFro(0:7), i, iBas, iGo, iIrrep, ij, ipTmp1, iSpin, jBas, nAct, nDens_Valence, nsa, nTst, iRout, iPrint, iComp
logical(kind=iwp) :: lPrint
real(kind=wp) :: CoefX, CoefR
character(len=8) :: RlxLbl, Method
character(len=80) :: KSDFT
real(kind=wp), allocatable :: D1AV(:), Tmp(:)
real(kind=wp), external :: Get_ExFac

!...  Prologue
iRout = 205
iPrint = nPrint(iRout)
lPrint = .true.
iD0Lbl = 1
iComp = 1

nDens = nBT
nDens_Valence = nBVT

lsa = .false.
Gamma_On = .false.
lPSO = .false.

!...  Get the method label
call Get_cArray('Relax Method',Method,8)
nCMo = S%n2Tot
mCMo = S%n2Tot
if (Method == 'KS-DFT  ' .or. Method == 'CASDFT  ') then
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

!... Check the wave function type

!                                                                      *
!***********************************************************************
!                                                                      *
if (Method == 'RHF-SCF ' .or. Method == 'UHF-SCF ' .or. (Method == 'KS-DFT  ' .and. iSpin == 1) .or. Method == 'ROHF    ') then
  if (lPrint) then
    write(u6,*)
    write(u6,'(2A)') ' Wavefunction type: ',Method
    if (Method == 'KS-DFT  ') then
      write(u6,'(2A)') ' Functional type:   ',KSDFT
      write(u6,100) 'Exchange scaling factor',CoefX
      write(u6,100) 'Correlation scaling factor',CoefR
    end if
    write(u6,*)
  end if
!                                                                      *
!***********************************************************************
!                                                                      *
else if (Method == 'Corr. WF') then
  if (lPrint) then
    write(u6,*)
    write(u6,*) ' Wavefunction type: an Aces 2 correlated wavefunction'
    write(u6,*)
  end if
  Gamma_On = .true.
  call Aces_Gamma()
!                                                                      *
!***********************************************************************
!                                                                      *
else if (Method == 'RASSCF  ' .or. Method == 'CASSCF  ' .or. Method == 'CASDFT  ') then

  call Get_iArray('nAsh',nAsh,nIrrep)
  nAct = 0
  do iIrrep=0,nIrrep-1
    nAct = nAct+nAsh(iIrrep)
  end do
  if (nAct > 0) lPSO = .true.

  nDSO = nDens
  mIrrep = nIrrep
  mBas(0:nIrrep-1) = nBas(0:nIrrep-1)
  if (lPrint) then
    write(u6,*)
    write(u6,'(2A)') ' Wavefunction type: ',Method
    if (Method == 'CASDFT  ') write(u6,'(2A)') ' Functional type:   ',KSDFT
    write(u6,*)
  end if
!                                                                      *
!***********************************************************************
!                                                                      *
else if (Method == 'CASSCFSA' .or. Method == 'RASSCFSA') then
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
  mBas(0:nIrrep-1) = nBas(0:nIrrep-1)
  if (lPrint .and. lSA) then
    write(u6,*)
    write(u6,'(2A)') ' Wavefunction type: State average ',Method(1:6)
    write(u6,*)
  else if (lPrint) then
    write(u6,*)
    write(u6,'(2A)') ' Wavefunction type: ',Method
  end if
  Method = 'RASSCF  '
!                                                                      *
!***********************************************************************
!                                                                      *
else
  write(u6,*)
  write(u6,*) ' Wavefunction type:',Method
  write(u6,*) ' Illegal type of wave function!'
  write(u6,*) ' ALASKA cannot continue'
  write(u6,*)
  call Quit_OnUserError()
end if

!...  Read the (non) variational 1st order density matrix
!...  density matrix in AO/SO basis
nsa = 1
if (lsa) nsa = 4
mDens = nsa
call mma_allocate(D0,nDens,mDens,Label='D0')
call mma_allocate(DVar,nDens,mDens,Label='DVar')
D0(:,:) = Zero
DVar(:,:) = Zero
call Get_dArray_chk('D1ao',D0,nDens)
call Get_D1ao_Var(DVar,nDens)

call ReIndexFrag(D0,nDens,nDens_Valence,nBas,nBas_Valence,nIrrep)
call ReIndexFrag(DVar,nDens,nDens_Valence,nBas,nBas_Valence,nIrrep)
call AddFragDens(D0,nDens,nBas_Valence)
call AddFragDens(DVar,nDens,nBas_Valence)

call mma_allocate(DS,nDens,Label='DS')
call mma_allocate(DSVar,nDens,Label='DSVar')
DS(:) = Zero
DSVar(:) = Zero
if (Method == 'UHF-SCF ' .or. Method == 'ROHF    ' .or. Method == 'Corr. WF') then
  call Get_dArray_chk('D1sao',DS,nDens)
  call Get_D1sao_Var(DSVar,nDens)
end if

! Unfold density matrix

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
if (iPrint >= 99) then
  RlxLbl = 'D1AO    '
  call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],D0)
  RlxLbl = 'D1AO-Var'
  call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DVar)
  RlxLbl = 'DSAO    '
  call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DS)
  RlxLbl = 'DSAO-Var'
  call PrMtrx(RlxLbl,[iD0Lbl],iComp,[1],DSVar)
end if

!...  Get the MO-coefficients
if (Method == 'UHF-SCF ' .or. Method == 'ROHF    ' .or. Method == 'Corr. WF') then
  nsa = 2
else
  nsa = 1
  if (lsa) nsa = 2
end if
kCMO = nsa
call mma_allocate(CMO,mCMO,kCMO,Label='CMO')
call Get_dArray_chk('Last orbitals',CMO(:,1),mCMO)
if (iPrint >= 99) then
  ipTmp1 = 1
  do iIrrep=0,nIrrep-1
    call RecPrt(' CMO''s',' ',CMO(ipTmp1,1),nBas_Valence(iIrrep),nBas_Valence(iIrrep))
    ipTmp1 = ipTmp1+nBas_Valence(iIrrep)**2
  end do
end if

!...  Get additional information in the case of a RASSCF wave function
!...  Get the number of inactive, active and frozen orbitals
if (lpso) then
  call Get_iScalar('nSym',i)
  call Get_iArray('nIsh',nIsh,i)
  call Get_iArray('nAsh',nAsh,i)
  call Get_iArray('nFro',nFro,i)
  if (iPrint >= 99) then
    write(u6,*) ' nISh=',nISh
    write(u6,*) ' nASh=',nASh
    write(u6,*) ' nFro=',nFro
  end if
  nAct = 0
  nTst = 0
  do iIrrep=0,nIrrep-1
    nAct = nAct+nAsh(iIrrep)
    nTst = nTst+nFro(iIrrep)
  end do
  if (nTst /= 0) then
    write(u6,*)
    write(u6,*) ' No frozen orbitals are allowed!'
    write(u6,*) ' ALASKA cannot continue'
    write(u6,*)
    call Quit_OnUserError()
  end if

  !...  Get the one body density for the active orbitals
  !     (not needed for SA-CASSCF)
  nG1 = nAct*(nAct+1)/2
  nsa = 1
  if (lsa) nsa = 0
  mG1 = nsa
  call mma_allocate(G1,nG1,mG1,Label='G1')
  if (nsa > 0) then
    call Get_dArray_chk('D1mo',G1(:,1),nG1)
    if (iPrint >= 99) call TriPrt(' G1',' ',G1(:,1),nAct)
  end if

  !...  Get the two body density for the active orbitals
  nG2 = nG1*(nG1+1)/2
  nsa = 1
  if (lsa) nsa = 2
  mG2 = nsa
  call mma_allocate(G2,nG2,mG2,Label='G2')
  call Get_dArray_chk('P2MO',G2(:,1),nG2)
  if (iPrint >= 99) call TriPrt(' G2',' ',G2(1,1),nG1)
  if (lsa) then

    ! CMO1 Ordinary CMO's
    !
    ! CMO2 CMO*Kappa

    call Get_dArray_chk('LCMO',CMO(:,2),mCMO)
    if (iPrint >= 99) then
      ipTmp1 = 1
      do iIrrep=0,nIrrep-1
        call RecPrt('LCMO''s',' ',CMO(ipTmp1,2),nBas_Valence(iIrrep),nBas_Valence(iIrrep))
        ipTmp1 = ipTmp1+nBas_Valence(iIrrep)**2
      end do
    end if

    ! P are stored as
    !                              _                     _
    !   P1 = <i|e_pqrs|i> + sum_i <i|e_pqrs|i>+<i|e_pqrs|i>
    !   P2 = sum_i <i|e_pqrs|i>

    call Get_dArray_chk('PLMO',G2(:,2),nG2)
    call Daxpy_(nG2,One,G2(:,2),1,G2(:,1),1)
    if (iPrint >= 99) call TriPrt(' G2L',' ',G2(:,2),nG1)
    if (iPrint >= 99) call TriPrt(' G2T',' ',G2(:,1),nG1)

    call Get_dArray_chk('D2av',G2(:,2),nG2)
    if (iPrint >= 99) call TriPrt('G2A',' ',G2(:,2),nG2)

    ! Densities are stored as:
    !
    !     ipd0 AO:
    !
    !     D1 = inactive diagonal density matrix
    !                              _                 _
    !     D2 = <i|E_pq|i> + sum_i <i|E_pq|i>+<i|E_pq|i> + sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> - 1/2 D1
    !
    !     D3 = sum_i <i|E_pq|i> (active)
    !
    !     D4 = sum_i sum_o k_po <i|E_oq|i> +k_oq <i|E_po|i> (inactive)
    !
    !     G1 = <i|e_ab|i>
    !     G2 = sum i <i|e_ab|i>
    call mma_allocate(Tmp,2*ndens,label='Tmp')
    call Get_D1I(CMO(1,1),D0(1,1),Tmp,nish,nBas_Valence,nIrrep)
    call mma_deallocate(Tmp)

    call dcopy_(nDens_Valence,DVar,1,D0(1,2),1)
    if (.not. isNAC) call daxpy_(ndens,-Half,D0(1,1),1,D0(1,2),1)
    if (iprint > 90) call PrMtrx('D0',[iD0Lbl],iComp,[1],D0)

    ! This is necessary for the kap-lag

    nG1 = nAct*(NAct+1)/2
    call mma_Allocate(D1AV,nG1,Label='D1AV')
    call Get_dArray_chk('D1av',D1AV,nG1)
    call Get_D1A(CMO(1,1),D1AV,D0(1,3),nIrrep,nBas_Valence,nish,nash,nDens_Valence)
    call mma_deallocate(D1AV)

    call Get_dArray_chk('DLAO',D0(1,4),nDens)
  end if
  if (iPrint >= 99) call TriPrt(' G2',' ',G2(1,1),nG1)

end if

return

100 format(1X,A26,20X,F18.6)

end subroutine PrepP_FAIEMP
