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

subroutine InpOne()

use Index_Functions, only: nTri_Elem
use OneDat, only: sOpSiz
use rctfld_module, only: lRF
use MCLR_Data, only: CMO, Int1, KAIN1, nDens
use input_mclr, only: iSpin, nActEl, nAtoms, nBas, nFro, nIsh, nOrb, nSym, PotNuc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iCharge, iComp, idum(1), iiSym, iOpt, ip, ip2, iRC, iS, Leng
real(kind=wp) :: ExFac, Tot_Charge, Tot_El_Charge, Tot_Nuc_Charge
logical(kind=iwp) :: Dff, Do_DFT, Do_ESPF, First, NonEq
character(len=8) :: Label
real(kind=wp), allocatable :: D1ao(:), GTmp(:), HTmp(:), Nuc(:), Temp1(:), Temp2(:), Temp3(:)

iRc = -1
iOpt = ibset(0,sOpSiz)
nDens = sum(nBas(1:nSym)**2)
Label = 'ONEHAM'
iComp = 1
iisym = ibset(0,0)
call iRdOne(iRc,iOpt,Label,iComp,idum,iisym)
leng = idum(1)
if (iRC /= 0) then
  write(u6,*) 'InpOne: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
iRc = -1
iOpt = 0
call mma_allocate(Int1,nDens,Label='Int1')
kain1 => Int1

call mma_allocate(Temp1,leng+10,Label='Temp1')
call mma_allocate(Temp2,nDens,Label='Temp2')
call mma_allocate(Temp3,nDens,Label='Temp3')

call RdOne(iRc,iOpt,Label,iComp,Temp1,iisym)
if (iRC /= 0) then
  write(u6,*) 'InpOne: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
!nf

! Modify the one electron Hamiltonian for reaction
! field and ESPF calculations

call mma_allocate(Nuc,nAtoms,Label='Nuc')
call Get_dArray('Effective nuclear Charge',Nuc,nAtoms)
Tot_Nuc_Charge = sum(Nuc(1:nAtoms))
call mma_deallocate(Nuc)
Tot_El_Charge = -Two*sum(real(nFro(1:nSym)+nIsh(1:nSym),kind=wp))-real(nActEl,kind=wp)
Tot_Charge = Tot_Nuc_Charge+Tot_El_Charge
iCharge = int(Tot_Charge)
call DecideOnESPF(Do_ESPF)
if (Do_ESPF .or. lRF) then
  if (lRF) then
    write(u6,*) 'Sorry, MCLR+RF NYI'
    call Quit_OnUserError()
  end if

  ! Scratch for one- and two-electron type contributions
  ! + variational density-matrix

  call mma_allocate(Htmp,leng,Label='Htmp')
  call mma_allocate(Gtmp,leng,Label='Gtmp')
  Htmp(:) = Zero
  Gtmp(:) = Zero
  call mma_allocate(D1ao,leng,Label='D1ao')
  call Get_dArray_chk('D1ao',D1ao,leng)

  NonEq = .false.
  First = .true.
  Dff = .false.
  Do_DFT = .true.
  ExFac = Zero
  call Get_dScalar('PotNuc',PotNuc)

  ! Don't care about the last arguments: no (CAS-)DFT here I guess)

  call DrvXV(Htmp,Gtmp,D1ao,PotNuc,leng,First,Dff,NonEq,lRF,'SCF',ExFac,iCharge,iSpin,'1234',Do_DFT)
  Temp1(1:leng) = Temp1(1:leng)+Htmp(:)

  ! Hum, where the hell is FI (Fock Inactive) ???

  !FI(:) = FI(:)+Gtmp(:)
  call mma_deallocate(Gtmp)
  call mma_deallocate(Htmp)
  call mma_deallocate(D1ao)
end if
!nf
ip = 1
ip2 = 1
do iS=1,nSym
  if ((nBas(is) /= 0) .and. (nOrb(iS) /= 0)) then
    call Square(Temp1(ip),Temp2,1,nBas(is),nBas(is))
    ip = ip+nTri_Elem(nBas(is))
    call DGEMM_('T','N',nOrb(iS),nBas(iS),nBas(iS),One,CMO(ip2),nBas(iS),Temp2,nBas(iS),Zero,Temp3,nOrb(iS))
    call DGEMM_('N','N',nOrb(is),nOrb(iS),nBas(iS),One,Temp3,nOrb(iS),CMO(ip2),nBas(iS),Zero,Int1(ip2),nOrb(iS))
    ip2 = ip2+nBas(is)**2
  end if
end do
call mma_deallocate(Temp1)
call mma_deallocate(Temp2)
call mma_deallocate(Temp3)

end subroutine InpOne
