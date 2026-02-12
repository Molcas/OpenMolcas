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
! Copyright (C) 1990, Markus P. Fuelscher                              *
!***********************************************************************

!>  @brief
!>    Generate the Fock-Matrix for frozen and inactive orbitals in
!>      the basis of active orbitals.
!>
!>  @author
!>    Markus P. Fuelscher
!>
!>  @details
!>  Generate the Fock-matrix for the frozen and inactive orbitals.
!>  Compute also the core energy and write to global variable EMY.
!>  Finally, transform the generated Fock-matrix
!>  into the basis of the active orbitals.
!>  Look into chapters 10.8.3 and 10.8.4 of \cite purple_book.
!>  The one body density matrices are required for e.g. reaction field
!>  or DFT calculations. In this case they are used to create a modified
!>  Fock Matrix.
!>
!>  @param[in] CMO The MO-coefficients
!>  @param[out] F The inactive Fock matrix in the basis of the active MO
!>  @param[in,out] FI The inactive Fock matrix in AO-space
!>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu})\f]
!>    In output FI contains also the core energy added to
!>    the diagonal elements.
!>    \f[\sum_{\sigma\rho} D^I_{\sigma\rho}(g_{\mu\nu\sigma\rho} - \frac{1}{2} g_{\mu\sigma\rho\nu}) + \frac{E^{(0)}}{n_{el}} \delta_{\mu\nu} \f]
!>  @param[in] D1I The inactive one-body density matrix in AO-space
!>    \f[D^{\text{AO}, I} = 2 C (C^I)^\dagger \f]
!>    See ::get_D1I_rasscf.
!>  @param[in] D1A The active one-body density matrix in AO-space
!>    \f[ D^{\text{AO}, A} = C^A D^A (C^A)^\dagger \f]
!>    See ::get_D1A_rasscf.
!>  @param[in] D1S The active spin density matrix in AO-space
!>    \f[ D^{\text{AO}, A}_S = C^A (D^A_\alpha - D^A_\beta) (C^A)^\dagger \f]
subroutine SGFCIN(CMO,F,FI,D1I,D1A,D1S)

use RunFile_procedures, only: Get_dExcdRa
use fcidump, only: DumpOnly
use fciqmc, only: DoNECI
use CC_CI_mod, only: Do_CC_CI
use timers, only: TimeDens
use lucia_data, only: INT1, INT1O
use rasscf_global, only: EMY, KSDFT, dftfock, exfac, nac, nacpar, noneq, potnuc, rfpert, tot_charge, tot_el_charge, &
                         tot_nuc_charge, doBlockDMRG, doDMRG
use OneDat, only: sNoNuc, sNoOri
use general_data, only: iSpin, nActEl, nSym, nTot1, nBas, nIsh, nAsh, nFro
use OFEmbed, only: Do_OFemb, OFE_first, FMaux, Rep_EN
use rctfld_module, only: lRF
use PrintLevel, only: DEBUG
use output_ras, only: LF, IPRLOC
#ifdef _DMRG_
use qcmaquis_interface_cfg
#endif
use Constants, only: Zero, One
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
character(len=16), parameter :: ROUTINE = 'SGFCIN  '
real*8, intent(in) :: CMO(*), D1I(*), D1A(*)
real*8, intent(inout) :: FI(*), D1S(*), F(*)
character(len=8) Label
logical First, Dff, Do_DFT, Found
logical Do_ESPF
real*8 :: CASDFT_Funct, dum1, dum2, dum3, dumm(1), Emyn, Eone, Erf1, Erf2, Erfx, Etwo, potnuc_ref, Time(2)
integer :: i, iadd, ibas, icharge, iComp, ioff, iopt, iprlev, ntmpfck, irc, iSyLbl, iSym, iTu, j, mxna, mxnb, nAt, nst, nt, ntu, &
           nu, nvxc
real*8, allocatable :: TmpFckI(:), Tmpx(:)
real*8, allocatable :: Tmp0(:), Tmp1(:), Tmp2(:), Tmp3(:), Tmp4(:), Tmp5(:), Tmp6(:), Tmp7(:), Tmpz(:), X0(:), X1(:), X2(:), X3(:)
real*8, external :: dDot_

! Local print level (if any)
IPRLEV = IPRLOC(3)
IPRLEV = 0000
if (IPRLEV >= DEBUG) write(LF,*) ' Entering ',ROUTINE

! Generate molecular charges
call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
iRc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
Tot_Nuc_Charge = Tmp0(nTot1+4)
if (iRc /= 0) then
  write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
  write(LF,*) 'Label = ',Label
  write(LF,*) 'iRc = ',iRc
  call Abend()
end if
call mma_deallocate(Tmp0)
Tot_El_Charge = Zero
do iSym=1,nSym
  Tot_El_Charge = Tot_El_Charge-2.0d0*dble(nFro(iSym)+nIsh(iSym))
end do
Tot_El_Charge = Tot_El_Charge-dble(nActEl)
Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge

! Load bare nuclei Hamiltonian
call mma_allocate(Tmp1,nTot1,Label='Tmp1')
iComp = 1
iSyLbl = 1
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'OneHam  '
call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
if (iRc /= 0) then
  write(LF,*) 'SGFCIN: iRc from Call RdOne not 0'
  write(LF,*) 'Label = ',Label
  write(LF,*) 'iRc = ',iRc
  call Abend()
end if
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' CMO in SGFCIN'
  write(LF,*) ' ---------------------'
  write(LF,*)
  ioff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas /= 0) then
      write(6,*) 'Sym =',iSym
      do i=1,iBas
        write(6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
      end do
      iOff = iOff+(iBas*iBas)
    end if
  end do

  write(LF,*)
  write(LF,*) ' D1I in AO basis in SGFCIN'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1I(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do

  write(LF,*)
  write(LF,*) ' D1A in AO basis in SGFCIN'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do
end if
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' OneHam in AO basis in SGFCIN'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if

! Load the nuclear repulsion energy

!call Get_PotNuc(PotNuc)
call Get_dScalar('PotNuc',potNuc)

! modify the one electron Hamiltonian for reaction field calculations
ERFX = Zero
iCharge = int(Tot_Charge)

call mma_allocate(Tmp3,nTot1,Label='Tmp3')
call mma_allocate(Tmp4,nTot1,Label='Tmp4')
call mma_allocate(Tmp7,nTot1,Label='Tmp7')

call DecideOnESPF(Do_ESPF)

! Generate total density

call Fold(nSym,nBas,D1I,Tmp3)
call Fold(nSym,nBas,D1A,Tmp4)
call Daxpy_(nTot1,1.0d0,Tmp4,1,Tmp3,1)
call Put_dArray('D1ao',Tmp3,nTot1)
!write(LF,*)
!write(LF,*) ' D1ao in AO basis in SGFCIN'
!write(LF,*) ' ---------------------'
!write(LF,*)
!iOff = 1
!do iSym=1,nSym
!  iBas = nBas(iSym)
!  call TriPrt(' ','(5G17.11)',Tmp3(iOff),iBas)
!  iOff = iOff+(iBas*iBas+iBas)/2
!end do

! Generate spin-density

call Fold(nSym,nBas,D1S,Tmp7)
call Put_dArray('D1sao',Tmp7,nTot1)

if ((KSDFT(1:3) /= 'SCF') .or. Do_OFemb) then
  call Put_iArray('nFro',nFro,nSym)
  call Put_iArray('nAsh',nAsh,nSym)
  call Put_iArray('nIsh',nIsh,nSym)
end if

if (Do_ESPF .or. lRF .or. (KSDFT /= 'SCF') .or. Do_OFemb) then

  ! Scratch for one- and two-electron type contributions

  call mma_allocate(Tmp5,nTot1,Label='Tmp5')
  Tmp5(:) = 0.0d0
  call mma_allocate(Tmp6,nTot1,Label='Tmp6')
  Tmp6(:) = 0.0d0

  First = .true.
  Dff = .false.
  Do_DFT = .true.

  call Timing(Time(1),dum1,dum2,dum3)

  call DrvXV(Tmp5,Tmp6,Tmp3,PotNuc,nTot1,First,Dff,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,DFTFOCK,Do_DFT)
  if (IPRLEV >= DEBUG) then
    write(LF,*)
    write(LF,*) ' Tmp5, h1 (DFT), in AO basis in SGFCIN'
    write(LF,*) ' ---------------------'
    write(LF,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call TriPrt(' ','(5G17.11)',Tmp5(iOff),iBas)
      iOff = iOff+(iBas*iBas+iBas)/2
    end do
  end if
  call Timing(Time(2),dum1,dum2,dum3)
  TimeDens = TimeDens+Time(2)-Time(1)

  ERF1 = Zero
  ERF2 = dDot_(nTot1,Tmp6,1,Tmp4,1)
  ERFX = ERF1-0.5d0*ERF2
  call Daxpy_(nTot1,1.0d0,Tmp5,1,Tmp1,1)

  call Daxpy_(nTot1,1.0d0,Tmp6,1,FI,1)

  call mma_deallocate(Tmp6)
  call mma_deallocate(Tmp5)
end if

call mma_deallocate(Tmp7)
call mma_deallocate(Tmp4)
if (.not. Do_OFemb) call mma_deallocate(Tmp3)

if (RFpert) then

  ! Read the reaction field from RunFile or RunOld

  call f_Inquire('RUNOLD',Found)
  if (Found) call NameRun('RUNOLD')
  call mma_allocate(TmpZ,nTot1,Label='TmpZ')
  call Get_dScalar('RF Self Energy',ERFX)
  call Get_dArray('Reaction field',TmpZ,nTot1)
  call Daxpy_(nTot1,1.0d0,TmpZ,1,Tmp1,1)
  call mma_deallocate(TmpZ)
  if (Found) call NameRun('#Pop')
end if
call mma_allocate(Tmp2,nTot1,Label='Tmp2')
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' D1I in AO basis in SGFCIN'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',D1I(ioff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if
call Fold(nSym,nBas,D1I,Tmp2)

if (Do_OFemb) then
  if (OFE_first) then
    call mma_allocate(FMaux,nTot1,Label='FMAux')
    call Coul_DMB(.true.,1,Rep_EN,FMaux,Tmp3,Dumm,nTot1)
    OFE_first = .false.
  else
    call Coul_DMB(.false.,1,Rep_EN,FMaux,Tmp3,Dumm,nTot1)
  end if
  call DaXpY_(nTot1,One,FMaux,1,Tmp1,1)

  call NameRun('AUXRFIL') ! switch the RUNFILE name
  call Get_dExcdRa(Tmpx,nVxc)
  call DaXpY_(nTot1,One,Tmpx,1,Tmp1,1)
  if (nVxc == 2*nTot1) then ! Nuc Attr added twice
    call DaXpY_(nTot1,One,Tmpx(1+nTot1),1,Tmp1,1)
    call Get_dArray('Nuc Potential',Tmpx,nTot1)
    call DaXpY_(nTot1,-One,Tmpx,1,Tmp1,1)
  end if
  call mma_deallocate(Tmpx)
  call mma_deallocate(Tmp3)
  call NameRun('#Pop')    ! switch back to old RUNFILE
end if

! Compute energy contributions
Eone = dDot_(nTot1,Tmp2,1,Tmp1,1)
!call Get_PotNuc(PotNuc_Ref)
call Get_dScalar('PotNuc',PotNuc_Ref)
Eone = Eone+(PotNuc-PotNuc_Ref)
Etwo = dDot_(nTot1,Tmp2,1,FI,1)
call mma_deallocate(Tmp2)
EMY = PotNuc_Ref+Eone+0.5d0*Etwo+ERFX
CASDFT_Funct = Zero
if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) call Get_dScalar('CASDFT energy',CASDFT_Funct)
if (IPRLEV >= DEBUG) then
  write(LF,*) ' Nuclear repulsion 1      :',PotNuc
  write(LF,*) ' Nuclear repulsion 2 Ref  :',PotNuc_Ref
  write(LF,*) ' One-electron core energy :',Eone
  write(LF,*) ' Two-electron core energy :',Etwo
  write(LF,*) ' Total core energy        :',EMY
  if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT /= 'PAM')) write(LF,*) ' CASDFT Energy            :',CASDFT_Funct
end if

if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' FI matrix in AO in SGFCIN only 2-electron terms'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' Tmp1 matrix in SGFCIN, one-electron term'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if

! Assemble the one-electron Tmp1 and two-electron contribution to AO Fock matrix
call DaXpY_(nTot1,One,Tmp1,1,FI,1)
call mma_deallocate(Tmp1)

if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' Inactive Fock matrix in AO basis in SGFCIN'
  write(LF,*) '(it already contains OneHam and TwoEl contrib.)'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if

! Transform FI to active orbital basis and move it over to F.
! Remove also the symmetry blocking.
MXNB = 0
MXNA = 0
do ISYM=1,NSYM
  MXNB = max(MXNB,NBAS(ISYM))
  MXNA = max(MXNA,NASH(ISYM))
end do
call mma_allocate(X0,NTOT1,Label='X0')
call mma_allocate(X1,NTOT1,Label='X1')
call mma_allocate(X2,MXNB*MXNB,Label='X2')
call mma_allocate(X3,MXNB*MXNA,Label='X3')
call DCOPY_(NTOT1,FI,1,X1,1)
if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) then
  call Get_dExcdRa(TmpFckI,nTmpFck)
  call DaXpY_(NTOT1,1.0d0,TmpFckI,1,X1,1)
  if (IPRLEV >= DEBUG) then
    write(LF,*)
    write(LF,*) ' Exchange correlation in AO basis in SGFCIN'
    write(LF,*) ' ---------------------'
    write(LF,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call TriPrt(' ','(5G17.11)',TmpFckI(ioff),iBas)
      iOff = iOff+(iBas*iBas+iBas)/2
    end do
  end if
  call mma_deallocate(TmpFckI)
end if
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' Modified FI in AO basis in SGFCIN'
  write(LF,*) ' ---------------------'
  write(LF,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',X1(ioff),iBas)
    iOff = iOff+(iBas*iBas+iBas)/2
  end do
end if
call MOTRAC(CMO,X1,X2,X3)
call mma_deallocate(X3)
call mma_deallocate(X2)
call DCOPY_(NACPAR,[ZERO],0,F,1)
NTU = 0
ITU = 0
IADD = 0
!bjp
if (NACTEL /= 0) then
  EMYN = EMY/dble(NACTEL)
else
  EMYN = Zero
end if
do NST=1,NSYM
  NAT = NASH(NST)
  if (NAT /= 0) then
    do NT=1,NAT
      NTU = NTU+IADD
      do NU=1,NT
        NTU = NTU+1
        ITU = ITU+1
        F(NTU) = X1(ITU)
        if (NT == NU) F(NTU) = F(NTU)+EMYN
        X0(ITU) = F(NTU)
      end do
    end do
    IADD = IADD+NAT
  end if
end do

!Quan: Fix bug, skip Lucia stuff with DMRG
! and other external CI solvers.
if (.not. any([DoNECI,Do_CC_CI,DumpOnly,doDMRG,doBlockDMRG])) then
  INT1(1:ITU) = X0(1:ITU)
  INT1(ITU+1:) = Zero
  INT1O(1:ITU) = X0(1:ITU)
  INT1O(ITU+1:) = Zero
end if

call mma_deallocate(X1)
call mma_deallocate(X0)

! print h0
if (IPRLEV >= DEBUG) then
  write(LF,*)
  write(LF,*) ' Inactive Fock mat in act MO basis, h0, in SGFCIN'
  write(LF,*) ' ------------'
  write(LF,*)
  call TriPrt(' ',' ',F,NAC)
end if

return

end subroutine SGFCIN
