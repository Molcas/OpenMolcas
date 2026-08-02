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

use Index_Functions, only: nTri_Elem
use RunFile_procedures, only: Get_dExcdRa
use fcidump, only: DumpOnly
use fciqmc, only: DoNECI
use CC_CI_mod, only: Do_CC_CI
use timers, only: TimeDens
use lucia_data, only: INT1, INT1O
use rasscf_global, only: dftfock, doBlockDMRG, doDMRG, EMY, exfac, KSDFT, nac, nacpar, noneq, potnuc, rfpert, tot_charge, &
                         tot_el_charge, tot_nuc_charge
use OneDat, only: sNoNuc, sNoOri
use general_data, only: iSpin, nActEl, nAsh, nBas, nFro, nIsh, nSym, nTot1
use OFEmbed, only: Do_OFemb, FMaux, OFE_first, Rep_EN
use rctfld_module, only: lRF
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO(*), D1I(*), D1A(*), D1S(*)
real(kind=wp), intent(out) :: F(NACPAR)
real(kind=wp), intent(inout) :: FI(nTot1)
integer(kind=iwp) :: i, iadd, ibas, icharge, iComp, ioff, iopt, iprlev, irc, iSyLbl, iSym, iTu, j, mxna, mxnb, nAt, nst, nt, &
                     ntmpfck, ntu, nu, nvxc
real(kind=wp) :: CASDFT_Funct, dum1, dum2, dum3, dumm(1), Emyn, Eone, Erf1, Erf2, Erfx, Etwo, potnuc_ref, Time(2)
character(len=8) :: Label
logical(kind=iwp) :: Dff, Do_DFT, Do_ESPF, First, Found
real(kind=wp), allocatable :: Tmp0(:), Tmp1(:), Tmp2(:), Tmp3(:), Tmp4(:), Tmp5(:), Tmp6(:), Tmp7(:), TmpFckI(:), Tmpx(:), &
                              Tmpz(:), X0(:), X1(:), X2(:), X3(:)
real(kind=wp), external :: dDot_

! Local print level (if any)
IPRLEV = IPRLOC(3)
IPRLEV = 0000
if (IPRLEV >= DEBUG) write(u6,*) ' Entering SGFCIN'

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
  write(u6,*) 'SGFCIN: iRc from Call RdOne not 0'
  write(u6,*) 'Label = ',Label
  write(u6,*) 'iRc = ',iRc
  call Abend()
end if
call mma_deallocate(Tmp0)
Tot_El_Charge = -Two*sum(nFro(1:nSym)+nIsh(1:nSym))
Tot_El_Charge = Tot_El_Charge-real(nActEl,kind=wp)
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
  write(u6,*) 'SGFCIN: iRc from Call RdOne not 0'
  write(u6,*) 'Label = ',Label
  write(u6,*) 'iRc = ',iRc
  call Abend()
end if
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' CMO in SGFCIN'
  write(u6,*) ' ---------------------'
  write(u6,*)
  ioff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    if (iBas /= 0) then
      write(u6,*) 'Sym =',iSym
      do i=1,iBas
        write(u6,*) (CMO(ioff+iBas*(i-1)+j),j=0,iBas-1)
      end do
      iOff = iOff+(iBas*iBas)
    end if
  end do

  write(u6,*)
  write(u6,*) ' D1I in AO basis in SGFCIN'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1I(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do

  write(u6,*)
  write(u6,*) ' D1A in AO basis in SGFCIN'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do
end if
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' OneHam in AO basis in SGFCIN'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
    iOff = iOff+nTri_Elem(iBas)
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
Tmp3(:) = Tmp3(:)+Tmp4(:)
call Put_dArray('D1ao',Tmp3,nTot1)
!write(u6,*)
!write(u6,*) ' D1ao in AO basis in SGFCIN'
!write(u6,*) ' ---------------------'
!write(u6,*)
!iOff = 1
!do iSym=1,nSym
!  iBas = nBas(iSym)
!  call TriPrt(' ','(5G17.11)',Tmp3(iOff),iBas)
!  iOff = iOff+nTri_Elem(iBas)
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
  Tmp5(:) = Zero
  call mma_allocate(Tmp6,nTot1,Label='Tmp6')
  Tmp6(:) = Zero

  First = .true.
  Dff = .false.
  Do_DFT = .true.

  call Timing(Time(1),dum1,dum2,dum3)

  call DrvXV(Tmp5,Tmp6,Tmp3,PotNuc,nTot1,First,Dff,NonEq,lRF,KSDFT,ExFac,iCharge,iSpin,DFTFOCK,Do_DFT)
  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' Tmp5, h1 (DFT), in AO basis in SGFCIN'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call TriPrt(' ','(5G17.11)',Tmp5(iOff),iBas)
      iOff = iOff+nTri_Elem(iBas)
    end do
  end if
  call Timing(Time(2),dum1,dum2,dum3)
  TimeDens = TimeDens+Time(2)-Time(1)

  ERF1 = Zero
  ERF2 = dDot_(nTot1,Tmp6,1,Tmp4,1)
  ERFX = ERF1-Half*ERF2
  Tmp1(:) = Tmp1(:)+Tmp5(:)

  FI(:) = FI(:)+Tmp6(:)

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
  Tmp1(:) = Tmp1(:)+TmpZ(:)
  call mma_deallocate(TmpZ)
  if (Found) call NameRun('#Pop')
end if
call mma_allocate(Tmp2,nTot1,Label='Tmp2')
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' D1I in AO basis in SGFCIN'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',D1I(ioff),iBas)
    iOff = iOff+nTri_Elem(iBas)
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
  Tmp1(:) = Tmp1(:)+FMaux(:)

  call NameRun('AUXRFIL') ! switch the RUNFILE name
  call Get_dExcdRa(Tmpx,nVxc)
  Tmp1(:) = Tmp1(:)+Tmpx(1:nTot1)
  if (nVxc == 2*nTot1) then ! Nuc Attr added twice
    Tmpx(:) = Tmp1(:)+Tmpx(nTot1+1:2*nTot1)
    call Get_dArray('Nuc Potential',Tmpx,nTot1)
    Tmp1(:) = Tmp1(:)-Tmpx(nTot1+1:2*nTot1)
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
EMY = PotNuc_Ref+Eone+Half*Etwo+ERFX
CASDFT_Funct = Zero
if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) call Get_dScalar('CASDFT energy',CASDFT_Funct)
if (IPRLEV >= DEBUG) then
  write(u6,*) ' Nuclear repulsion 1      :',PotNuc
  write(u6,*) ' Nuclear repulsion 2 Ref  :',PotNuc_Ref
  write(u6,*) ' One-electron core energy :',Eone
  write(u6,*) ' Two-electron core energy :',Etwo
  write(u6,*) ' Total core energy        :',EMY
  if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT /= 'PAM')) write(u6,*) ' CASDFT Energy            :',CASDFT_Funct
end if

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' FI matrix in AO in SGFCIN only 2-electron terms'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
    iOff = iOff+nTri_Elem(iBas)
  end do
end if
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' Tmp1 matrix in SGFCIN, one-electron term'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
    iOff = iOff+nTri_Elem(iBas)
  end do
end if

! Assemble the one-electron Tmp1 and two-electron contribution to AO Fock matrix
FI(:) = FI(:)+Tmp1(:)
call mma_deallocate(Tmp1)

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' Inactive Fock matrix in AO basis in SGFCIN'
  write(u6,*) '(it already contains OneHam and TwoEl contrib.)'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
    iOff = iOff+nTri_Elem(iBas)
  end do
end if

! Transform FI to active orbital basis and move it over to F.
! Remove also the symmetry blocking.
MXNB = maxval(NBAS(1:NSYM))
MXNA = maxval(NASH(1:NSYM))
call mma_allocate(X0,NTOT1,Label='X0')
call mma_allocate(X1,NTOT1,Label='X1')
call mma_allocate(X2,MXNB*MXNB,Label='X2')
call mma_allocate(X3,MXNB*MXNA,Label='X3')
X1(:) = FI(:)
if ((KSDFT(1:3) /= 'SCF') .and. (KSDFT(1:3) /= 'PAM')) then
  call Get_dExcdRa(TmpFckI,nTmpFck)
  X1(:) = X1(:)+TmpFckI(:)
  if (IPRLEV >= DEBUG) then
    write(u6,*)
    write(u6,*) ' Exchange correlation in AO basis in SGFCIN'
    write(u6,*) ' ---------------------'
    write(u6,*)
    iOff = 1
    do iSym=1,nSym
      iBas = nBas(iSym)
      call TriPrt(' ','(5G17.11)',TmpFckI(ioff),iBas)
      iOff = iOff+nTri_Elem(iBas)
    end do
  end if
  call mma_deallocate(TmpFckI)
end if
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' Modified FI in AO basis in SGFCIN'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',X1(ioff),iBas)
    iOff = iOff+nTri_Elem(iBas)
  end do
end if
call MOTRAC(CMO,X1,X2,X3)
call mma_deallocate(X3)
call mma_deallocate(X2)
F(:) = Zero
NTU = 0
ITU = 0
IADD = 0
!bjp
if (NACTEL /= 0) then
  EMYN = EMY/real(NACTEL,kind=wp)
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
  write(u6,*)
  write(u6,*) ' Inactive Fock mat in act MO basis, h0, in SGFCIN'
  write(u6,*) ' ------------'
  write(u6,*)
  call TriPrt(' ',' ',F,NAC)
end if

return

end subroutine SGFCIN
