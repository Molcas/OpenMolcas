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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

!#define _DEBUGPRINT_
subroutine PMat_SCF(FstItr,XCf,nXCF,nD)
!***********************************************************************
!                                                                      *
!     purpose: Compute two-electron part of the Fock matrix            *
!                                                                      *
!     output:                                                          *
!       TwoHam  : two-electron part of the Fock matrix constructed by  *
!                 contraction of proper density matrix difference with *
!                 two-electron integrals i) in conventional way ii) in *
!                 direct way                                           *
!                                                                      *
!***********************************************************************

use OFembed, only: Do_OFemb
use InfSCF, only: Algo, DDnOff, Dens, dfkMat, DoFMM, DSCF, EDFT, exFac, FMOMax, FockAO, iDisk, iDummy_Run, ipsLst, Iter, Klockan, &
                  KSDFT, MapDns, MaxBas, MiniDn, MxConstr, nBas, nBB, nBT, nCore, nDens, nDIsc, nIter, nIterP, nOcc, NoExchange, &
                  nSkip, nSym, OneHam, pMTime, PotNuc, PreSch, RFPert, Thize, TimFld, tNorm, Tot_Charge, TwoHam, Vxc
use RICD_Info, only: Do_DCCD
use Int_Options, only: Exfac_Int => ExFac, FckNoClmb, PreSch_Int => PreSch, Thize_Int => Thize
use rctfld_module, only: lRF
use Integral_interfaces, only: Drv2El_dscf
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(inout) :: FstItr
integer(kind=iwp), intent(in) :: nXCf, nD
real(kind=wp), intent(in) :: XCf(nXCf,nD)
integer(kind=iwp) :: Algo_Save, iCharge, iD, iDumm, iM, iMat, iSpin, nT
real(kind=wp) :: Backup_ExFac, CPU1, CPU2, Dummy(1), ERFSelf, TCF2, TCF2_1, Tim1, Tim2, Tim3, Tmp, TWF2, TWF2_1, XCPM, XCPM1, &
                 XCPM2, XWPM, XwPM1, XwPM2
logical(kind=iwp) :: Do_DFT, Do_ESPF, First = .true., Found, ltmp1, ltmp2, NonEq
real(kind=wp), allocatable :: D(:), DnsS(:,:), RFfld(:), Saved(:,:), tVxc(:,:)
real(kind=wp), allocatable, target :: Aux(:,:), Temp(:,:)
real(kind=wp), pointer :: pTwoHam(:,:)
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: EFP_On

if (PmTime) call CWTime(xCPM1,xWPM1)
call Timing(Cpu1,Tim1,Tim2,Tim3)
#ifdef _DEBUGPRINT_
call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat: Enter','T in nDens')
call NrmClc(Vxc(1,1,nDens),nBT*nD,'PMat: Enter','T in nDens')
call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat: Enter','T in iPsLst')
call NrmClc(Vxc(1,1,nDens),nBT*nD,'PMat: Enter','T in iPsLst')
#endif

! Copy the (abs.) value of the Max Offdiag Fmat to a Common Block
! Used in the LK Cholesky algorithm
dFKmat = abs(FMOmax)

! Add contribution due to external potential

TwoHam(:,:,iPsLst) = Zero
iSpin = 1
if (nD == 2) iSpin = 2
call Put_iScalar('Multiplicity',iSpin)
!                                                                      *
!***********************************************************************
!                                                                      *
call DecideOnESPF(Do_ESPF)
if (Do_ESPF .or. lRF .or. (KSDFT /= 'SCF') .or. Do_OFemb .or. EFP_On()) then

  ! Observe that this call always has to be prior to the calls
  ! to Drv2El_dScf and/or FTwoa. This since DrvXV will redefine
  ! ExFac!!!!!

  ! Note, the linear (Oneham) and bilinear (TwoHam) contributions
  ! can be computed with partial densities, as supplied with the
  ! arguments to the routine. For the DFT contributions, however,
  ! not being linear or bilinears, the total density is read from
  ! the runfile (as put there by dmat).

  iCharge = int(Tot_Charge)
  NonEq = .false.
  Do_DFT = .true.
  iDumm = 1
  ltmp1 = iter == 1
  ltmp2 = iter /= 1
  if (nD == 1) then
    call DrvXV(OneHam,TwoHam(:,1,iPsLst),Dens(:,1,iPsLst),PotNuc,nBT,ltmp1,ltmp2,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,'SCF ',Do_DFT)
  else
    call mma_allocate(D,nBT,Label='D')
    D(:) = Dens(:,1,iPsLst)+Dens(:,2,iPsLst)
    call DrvXV(OneHam,TwoHam(:,1,iPsLst),D,PotNuc,nBT,ltmp1,ltmp2,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,'SCF ',Do_DFT)
    call mma_deallocate(D)
    TwoHam(:,2,iPsLst) = TwoHam(:,1,iPsLst)
    if ((MxConstr > 0) .and. (klockan == 1)) then
      call SetUp_iSD()
      call Get_Enondyn_dft(nBT,Dummy,iDumm,'SCF ')
      call Free_iSD()
      klockan = 24
    end if
  end if

  ! Pick up the integrated energy contribution of the external
  ! potential to the total energy.

  call Peek_dScalar('KSDFT energy',EDFT(iter))

  ! Pick up the contribution to the Fock matrix due to the
  ! external field. Note that for some external field the
  ! potential is neither linear nor bi-linear.

  if (KSDFT /= 'SCF') then
    call Get_dArray_chk('dExcdRa',Vxc(:,:,iPsLst),nBT*nD)
  else
    Vxc(:,:,iPsLst) = Zero
  end if

  if (Do_OFemb) then
    call NameRun('AUXRFIL') ! switch the RUNFILE name
    call mma_allocate(tVxc,nBT,nD,Label='tVxc')
    call Get_dArray_chk('dExcdRa',tVxc,nBT*nD)
    Vxc(:,:,iPsLst) = Vxc(:,:,iPsLst)+tVxc(:,:)
    call mma_deallocate(tVxc)
    call NameRun('#Pop')    ! switch back RUNFILE name
  end if
# ifdef _DEBUGPRINT_
  call NrmClc(Vxc(:,:,iPsLst),nBT*nD,'PMat','Optimal V ')
# endif

else if (RFpert .and. First) then

  if (nD == 2) then
    write(u6,*) ' UHF+RF: Not implemented'
    call Abend()
  end if
  call mma_allocate(RFfld,nBT,Label='RFfld')
  call f_Inquire('RUNOLD',Found)
  if (Found) call NameRun('RUNOLD')
  call Get_dScalar('RF Self Energy',ERFself)
  call Get_dArray('Reaction field',RFfld,nBT)
  if (Found) call NameRun('#Pop')
  PotNuc = PotNuc+ERFself
  OneHam(:) = OneHam(:)+RFfld(:)
  do iD=1,nD
    FockAO(:,iD) = OneHam(:)
  end do
  call mma_deallocate(RFfld)

else

  Vxc(:,:,iPsLst) = Zero

end if
!                                                                      *
!***********************************************************************
!                                                                      *
First = .false.

! Compute the two-electron contribution to the Fock matrix

Backup_ExFac = ExFac
if (NoExchange) ExFac = Zero

nT = 1+(nD-1)*2

call mma_allocate(Temp,nBT,nT,Label='Temp')
Temp(:,:) = Zero

if (PmTime) call CWTime(tCF2,tWF2)

if (DSCF .and. (.not. Do_DCCD)) then

  call Drv2El_dscf_Front_End(Temp)

else   ! RICD/Cholesky option

  ! Allocate memory for squared density matrix
  call mma_allocate(DnsS,nBB,nD,Label='DnsS')

  ! Expand the 1-density matrix
  do iD=1,nD
    call Unfold(Dens(1,iD,iPsLst),nBT,DnsS(1,iD),nBB,nSym,nBas)
  end do

  call FockTwo_Drv_scf(nSym,nBas,nBas,nSkip,Dens(:,:,iPsLst),DnsS,Temp,nBT,ExFac,nBB,MaxBas,nD,nOcc,size(nOcc,1),iDummy_run)

  if (Do_DCCD) then
    call mma_Allocate(Saved,size(Temp,1),size(Temp,2),Label='Saved')

    Saved(:,:) = Zero
    call Drv2El_dscf_Front_End(Saved)
    Temp(:,:) = Temp(:,:)+Saved(:,:)

    Algo_save = Algo
    Algo = 0
    Saved(:,:) = Zero
    call FockTwo_Drv_scf(nSym,nBas,nBas,nSkip,Dens(:,:,iPsLst),DnsS,Saved,nBT,ExFac,nBB,MaxBas,nD,nOcc,size(nOcc,1),iDummy_run)
    Temp(:,:) = Temp(:,:)-Saved(:,:)
    Algo = Algo_save
    call mma_deAllocate(Saved)
  end if

  ! Deallocate memory for squared density matrix
  call mma_deallocate(DnsS)

end if

if (PmTime) then
  call CWTime(tCF2_1,tWF2_1)
  tCF2 = tCF2_1-tCF2
  tWF2 = tWF2_1-tWF2
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Add on FMM contributions to Fock matrix, only works for HF

if (DoFMM) call FMMFck(Dens(1,1,iPsLst),Temp(1,1),nDens)
!                                                                      *
!***********************************************************************
!                                                                      *
TwoHam(:,:,iPsLst) = TwoHam(:,:,iPsLst)+Temp(:,1:nD)
#ifdef _DEBUGPRINT_
call NrmClc(Temp,nBT*nD,'PMat_SCF','Temp')
call NrmClc(TwoHam(:,:,iPsLst),nBT*nD,'PMat_SCF','T in iPsLst')
#endif
call mma_deallocate(Temp)

!                                                                      *
!***********************************************************************
!                                                                      *
! Now compute the total two-electron contribution

! Generate the two-electron contribution corresponding to the total
! density.

if (MiniDn .and. (max(0,nIter(nIterP)-1) > 0)) then

  ! Minimized density option
  !
  ! D(k+1) = Sum_i_k C_i D_i + delta(k+1)
  !
  ! G(D(k+1)) = G(delta(k+1)) + Sum_i_k C_i G(D_i)

  call mma_allocate(Aux,nBT,nD,Label='Aux')
  do iMat=1,iter-1

    tmp = sum(abs(XCf(iMat,:)))
    if (tmp == Zero) cycle

    iM = MapDns(iMat)
    if (iM < 0) then
      call RWDTG(-iM,Aux,nBT*nD,'R','TWOHAM',iDisk,size(iDisk,1))
      pTwoHam => Aux
    else
      pTwoHam => TwoHam(:,:,iM)
    end if

    do iD=1,nD
      if (Xcf(iMat,iD) == Zero) cycle
      TwoHam(:,iD,iPsLst) = TwoHam(:,iD,iPsLst)+Xcf(iMat,iD)*pTwoHam(:,iD)
    end do

    nullify(pTwoHam)

  end do
  call mma_deallocate(Aux)

else if (.not. DDnOFF) then

  ! Normal density difference
  !
  ! G(D(k+1)) = G(D(k)) + G(D(k+1)-D(k))

  TwoHam(:,:,iPsLst) = TwoHam(:,:,iPsLst)+TwoHam(:,:,nDens)

end if

! Restore the total density in position iPsLst

Dens(:,:,iPsLst) = Dens(:,:,nDens)
TwoHam(:,:,nDens) = TwoHam(:,:,iPsLst)
Vxc(:,:,nDens) = Vxc(:,:,iPsLst)

! Restore ExFac (if it was changed)
if (NoExchange) ExFac = Backup_ExFac
!                                                                      *
!***********************************************************************
!                                                                      *
TNorm = DDot_(nBT*nD,TwoHam(1,1,iPsLst),1,TwoHam(1,1,iPsLst),1)/real(nD,kind=wp)

#ifdef _DEBUGPRINT_
call NrmClc(Dens(1,1,iPsLst),nBT*nD,'PMat  ','D iPsLst  ')
call NrmClc(Dens(1,1,nDens),nBT*nD,'PMat  ','D nDens   ')
call NrmClc(TwoHam(1,1,iPsLst),nBT*nD,'PMat  ','T iPsLst  ')
call NrmClc(TwoHam(1,1,nDens),nBT*nD,'PMat  ','T nDens   ')
call NrmClc(Vxc(1,1,iPsLst),nBT*nD,'PMat  ','V iPsLst  ')
call NrmClc(Vxc(1,1,nDens),nBT*nD,'PMat  ','V nDens   ')

#endif
call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(5) = TimFld(5)+(Cpu2-Cpu1)
if (PmTime) then
  call CWTime(xCPM2,xWPM2)
  xCPM = xCPM2-xCPM1
  xWPM = xWPM2-xWPM1
  write(u6,'(1X,A,F15.2,A,F15.2,A,/,1X,A,F15.2,A,F15.2,A)') '>>> PMat_SCF: CPU  time:',xCPM,' seconds  (2-el contributions: ', &
                                                            tCF2,' seconds) <<<','>>> PMat_SCF: Wall time:',xWPM, &
                                                            ' seconds  (2-el contributions: ',tWF2,' seconds) <<<'
  call xFlush(u6)
end if

contains

subroutine Drv2El_dscf_Front_End(Tmp)

  real(kind=wp) :: Tmp(:,:)

  ! while the Drv2El_dscf can't handle UHF in a trivial way this interface has
  ! to be used.

  ExFac_Int = ExFac
  Thize_Int = Thize
  PreSch_Int = PreSch
  FckNoClmb = .false.
  if (size(Tmp,2) == 1) then
    call Drv2El_dscf(Dens(:,1,iPsLst),Tmp(:,1),nBT,0,FstItr)
  else

    ! Compute the Coulomb potential for the total density and
    ! exchange of alpha and beta, respectively. Add together
    ! to get the correct contributions to the alpha and beta
    ! Fock matrices.

    ! Set exchange factor to zero and compute only Coulomb
    ! for the total electron density.

    Tmp(:,2) = Dens(:,1,iPsLst)+Dens(:,2,iPsLst)

    ExFac_Int = Zero
    call Drv2El_dscf(Tmp(:,2),Tmp(:,3),nBT,max(nDisc*1024,nCore),FstItr)

    ! alpha exchange
    FckNoClmb = .true.
    Tmp(:,2) = Zero
    ExFac_Int = ExFac
    call Drv2El_dscf(Dens(:,1,iPsLst),Tmp(:,1),nBT,max(nDisc*1024,nCore),FstItr)
    Tmp(:,1) = Two*Tmp(:,1)

    ! beta exchange
    ExFac_Int = ExFac
    call Drv2El_dscf(Dens(:,2,iPsLst),Tmp(:,2),nBT,max(nDisc*1024,nCore),FstItr)
    Tmp(:,2) = Two*Tmp(:,2)

    ! Add together J and K contributions to form the correct
    ! alpha and beta Fock matrices.

    Tmp(:,1) = Tmp(:,1)+Tmp(:,3)
    Tmp(:,2) = Tmp(:,2)+Tmp(:,3)
  end if

end subroutine Drv2El_dscf_Front_End

end subroutine PMat_SCF
