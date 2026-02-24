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
!               2013, Giovanni Li Manni                                *
!***********************************************************************

subroutine CASDFT_terms(CMO,F,FI,D1I,D1A,D1S)
! This routine is a modification of SGFCIN, adapted to a CASDFT
! implementation in which the CI step of a CASDFT calculation is
!     not corrected by DFT. DFT will play a role only in the Orbital
!     optimization step.
! Purpose:
! Generate the Fock-matrix for the frozen and inactive orbitals.
! Compute also the core energy. Finally, transform the Fock-matrix
! into the basis of the active orbitals.
!
! M.P. Fuelscher, Lund, July 1990
! GLM, Minneapolis,   May 2013

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use rctfld_module, only: lRF
use PrintLevel, only: DEBUG
use output_ras, only: IPRLOC
use general_data, only: ISPIN, NACTEL, NASH, NBAS, NFRO, NISH, NSYM, NTOT1
use rasscf_global, only: DFTFOCK, Emy, ExFac, KSDFT_temp, NAC, NACPAR, NONEQ, PotNuc, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge
#ifdef _DMRG_
use lucia_data, only: INT1, INT1O
use rasscf_global, only: DoDMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO(*), D1I(*), D1A(*), D1S(*)
real(kind=wp), intent(out) :: F(NACPAR)
real(kind=wp), intent(inout) :: FI(*)
integer(kind=iwp) :: i, IADD, iBas, iCharge, iCOmp, iOff, iOpt, iPrLev, iRC, iSyLbl, iSym, ITU, j, MXNA, MXNB, NAT, NST, NT, NTU, NU
real(kind=wp) :: CASDFT_Funct, Emyn, Eone, ETwo, PotNuc_Ref
logical(kind=iwp) :: Dff, Do_DFT, First
character(len=8) :: Label
real(kind=wp), allocatable :: Tmp0(:), Tmp1(:), Tmp2(:), Tmp3(:), Tmp4(:), Tmp5(:), Tmp6(:), Tmp7(:), X0(:), X1(:), X2(:), X3(:)
real(kind=wp), external :: DDot_

!**********************************************************
! Local print level (if any)
!**********************************************************
IPRLEV = IPRLOC(3)
!IPRLEV = 100
if (IPRLEV >= DEBUG) then
  write(u6,*) 'Printing matrices in CASDFT_Terms'
  write(u6,*)
  write(u6,*) ' CMO in CASDFT_terms'
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
  write(u6,*) ' D1I in AO basis in CASDFT_Terms'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1I(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do

  write(u6,*)
  write(u6,*) ' D1S in AO basis in CASDFT_Terms'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1S(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do

  write(u6,*)
  write(u6,*) ' D1A in AO basis in CASDFT_Terms'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call wrtmat(D1A(ioff),iBas,iBas,iBas,iBas)
    iOff = iOff+iBas*iBas
  end do

end if

!**********************************************************
! Generate molecular charges
!**********************************************************
call mma_allocate(Tmp0,nTot1+4,Label='Tmp0')
iRc = -1
iOpt = ibset(0,sNoOri)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(iRc,iOpt,Label,iComp,Tmp0,iSyLbl)
Tot_Nuc_Charge = Tmp0(nTot1+4)
if (iRc /= 0) then
  write(u6,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
  write(u6,*) 'Label = ',Label
  write(u6,*) 'iRc = ',iRc
  call Abend()
end if
call mma_deallocate(Tmp0)

Tot_El_Charge = -Two*sum(nFro(1:nSym)+nIsh(1:nSym))
Tot_El_Charge = Tot_El_Charge-real(nActEl,kind=wp)
Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
!if (IPRLEV >= DEBUG) then
!  write(u6,*)
!  write(u6,*) 'Total Charge :',Tot_Charge
!end if

!**********************************************************
! Load bare nuclei Hamiltonian
!**********************************************************
call mma_Allocate(Tmp1,nTot1,Label='Tmp1')
iComp = 1
iSyLbl = 1
iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
Label = 'OneHam  '
call RdOne(iRc,iOpt,Label,iComp,Tmp1,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'CASDFT_Terms: iRc from Call RdOne not 0'
  write(u6,*) 'Label = ',Label
  write(u6,*) 'iRc = ',iRc
  call Abend()
end if
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' OneHam in AO basis in CASDFT_Terms'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',Tmp1(iOff),iBas)
    iOff = iOff+nTri_Elem(iBas)
  end do
end if

!**********************************************************
! Load the nuclear repulsion energy
!**********************************************************
call Get_dScalar('PotNuc',potNuc)
!write(u6,*)
!write(u6,*) 'PotNuc in casdft_terms.f:',PotNuc

!**********************************************************
! Generate total density
!**********************************************************

!if (IPRLEV >= DEBUG) then
!  call mma_allocate(Tmp31,nBas*nBas,Label='Tmp31')
!  Tmp31(:) = D1I(1:nBas**2)+D1A(1:nBas**2)
!  write(u6,*)
!  write(u6,*) ' DMAT not folded in AO basis in CASDFT_Terms'
!  write(u6,*) ' ---------------------'
!  write(u6,*)
!  call wrtmat(Tmp31,nBas,nBas,nBas,nBas)
!  call mma_deallocate(Tmp3)
!end if

call mma_allocate(Tmp3,nTot1,Label='Tmp3')
call mma_allocate(Tmp4,nTot1,Label='Tmp4')
call Fold(nSym,nBas,D1I,Tmp3)
call Fold(nSym,nBas,D1A,Tmp4)
Tmp3(:) = Tmp3(:)+Tmp4(:)
call Put_dArray('D1ao',Tmp3,nTot1)
!**********************************************************
! Generate spin-density
!**********************************************************
call mma_allocate(Tmp7,nTot1,Label='Tmp7')
call Fold(nSym,nBas,D1S,Tmp7)
call Put_dArray('D1sao',Tmp7,nTot1)
call mma_deallocate(Tmp7)

!**********************************************************
! One- and two-electron type contributions
!**********************************************************

call mma_allocate(Tmp5,nTot1,Label='Tmp5')
Tmp5(:) = Zero
call mma_allocate(Tmp6,nTot1,Label='Tmp6')
Tmp6(:) = Zero

First = .true.
Dff = .false.
Do_DFT = .true.

call Put_iArray('nFro',nFro,nSym)
call Put_iArray('nAsh',nAsh,nSym)
call Put_iArray('nIsh',nIsh,nSym)

iCharge = int(Tot_Charge)
! Tmp5 and Tmp6 are not updated in DrvXV...
call DrvXV(Tmp5,Tmp6,Tmp3,PotNuc,nTot1,First,Dff,NonEq,lRF,KSDFT_TEMP,ExFac,iCharge,iSpin,DFTFOCK,Do_DFT)

Tmp1(:) = Tmp1(:)+Tmp5(:)
FI(1:nTot1) = FI(1:nTot1)+Tmp6(:)

call mma_deallocate(Tmp6)
call mma_deallocate(Tmp5)
call mma_deallocate(Tmp4)
call mma_deallocate(Tmp3)

!**********************************************************
!     Compute energy contributions
!**********************************************************
call mma_allocate(Tmp2,nTot1,Label='Tmp2')

call Fold(nSym,nBas,D1I,Tmp2)

Eone = dDot_(nTot1,Tmp2,1,Tmp1,1)
call Get_dScalar('PotNuc',PotNuc_Ref)
Eone = Eone+(PotNuc-PotNuc_Ref)
Etwo = dDot_(nTot1,Tmp2,1,FI,1)
call mma_deallocate(Tmp2)
EMY = PotNuc_Ref+Eone+Half*Etwo

CASDFT_Funct = Zero
call Get_dScalar('CASDFT energy',CASDFT_Funct)
if (IPRLEV >= DEBUG) then
  write(u6,*) ' Nuclear repulsion energy :',PotNuc
  !write(u6,*) ' Nuclear repulsion energy Ref :',PotNuc_Ref
  write(u6,*) ' One-electron core energy :',Eone
  write(u6,*) ' Two-electron core energy :',Etwo
  write(u6,*) ' Total core energy        :',EMY
  write(u6,*) ' CASDFT Energy            :',CASDFT_Funct
end if

!**********************************************************
! Printing matrices
!**********************************************************
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' FI matrix in CASDFT_Terms only 2-electron terms'
  write(u6,*) ' ---------------------'
  write(u6,*)
  iOff = 1
  do iSym=1,nSym
    iBas = nBas(iSym)
    call TriPrt(' ','(5G17.11)',FI(iOff),iBas)
    iOff = iOff+nTri_Elem(iBas)
  end do
end if

FI(1:nTot1) = FI(1:nTot1)+Tmp1(:)

if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*) ' Inactive Fock matrix in AO basis in CASDFT_terms'
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

call mma_deallocate(Tmp1)

!**********************************************************
!     Transform FI to active orbital basis and move it over to F.
!     Remove also the symmetry blocking.
!**********************************************************
!**********************************************************
! Shall I add here the DFT contribution? Maybe not yet!
! I am commenting off... if needed we can always re-activate.
!**********************************************************
MXNB = maxval(NBAS(1:NSYM))
MXNA = maxval(NASH(1:NSYM))
call mma_allocate(X0,NTOT1,Label='X0')
call mma_allocate(X1,NTOT1,Label='X1')
call mma_allocate(X2,MXNB*MXNB,Label='X2')
call mma_allocate(X3,MXNB*MXNA,Label='X3')
X1(:) = FI(1:NTOT1)
!call Get_dExcdRa(TmpFckI,nTmpFck)
!X1(:) = X1(:)+TmpFckI(1:NTOT1)
!if (IPRLEV >= DEBUG) then
!  write(u6,*)
!  write(u6,*) ' Exchange Corr. in AO basis in CASDFT_Terms'
!  write(u6,*) ' ---------------------'
!  write(u6,*)
!  iOff = 1
!  do iSym=1,nSym
!    iBas = nBas(iSym)
!    call TriPrt(' ','(5G17.11)',TmpFckI(ioff),iBas)
!    iOff = iOff+nTri_Elem(iBas)
!  end do
!end if
!Call mma_deallocate(TmpFckI)
!if (IPRLEV >= DEBUG) then
!  write(u6,*)
!  write(u6,*) ' Modified FI in AO basis in CASDFT_Terms'
!  write(u6,*) ' ---------------------'
!  write(u6,*)
!  iOff = 1
!  do iSym=1,nSym
!    iBas = nBas(iSym)
!    call TriPrt(' ','(5G17.11)',X1(ioff),iBas)
!    iOff = iOff+nTri_Elem(iBas)
!  end do
!end if

call MOTRAC(CMO,X1,X2,X3)
call mma_deallocate(X3)
call mma_deallocate(X2)
F(:) = Zero
NTU = 0
ITU = 0
IADD = 0

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
#ifdef _DMRG_
if (.not. doDMRG) then
  INT1(1:ITU) = X0(1:ITU)
  INT1(ITU+1:) = Zero
  INT1O(1:ITU) = X0(1:ITU)
  INT1O(ITU+1:) = Zero
end if
#endif
call mma_deallocate(X1)
call mma_deallocate(X0)

! print h0
if (IPRLEV >= DEBUG) then
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,*) ' Fock matrix in MO basis, h0, in CASDFT_TERMS'
  write(u6,*) ' ------------'
  call TriPrt(' ',' ',F,NAC)
end if

end subroutine CASDFT_terms
