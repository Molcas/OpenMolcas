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

! NOT TESTED (used for OFEmbed below)
!#define _OFEmbed_
subroutine ADD1HAM(H1EFF,nH1Eff)
! ----------------------------------------------------------------
! Purpose: Reads and adds one-electron naked Hamiltonian into H1EFF.
! Dress it with reaction field (if any).
! Also get POTNUC at the same time.

#ifdef _OFEmbed_
use RunFile_procedures, only: Get_dExcdRa
use OFembed, only: Do_OFemb, FMAux, OFE_First
#endif
use OneDat, only: sNoNuc, sNoOri
use caspt2_module, only: ERFSelf, nBas, NBTRI, nSym, PotNuc, RFpert
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nH1EFF
real(kind=wp), intent(inout) :: H1EFF(nH1Eff)
integer(kind=iwp) :: ICOMP, IOPT, IRC, ISYLBL, iSym, nTemp
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: ISTLT
#endif
logical(kind=iwp) :: Found
character(len=8) :: Label
real(kind=wp), allocatable :: ONEHAM(:), Temp(:)
#ifdef _OFEmbed_
real(kind=wp), allocatable :: Coul(:)
#endif

! Add naked one-el Hamiltonian in AO basis to H1EFF:
call mma_allocate(ONEHAM,NBTRI,Label='OneHam')
IRC = -1
IOPT = ibset(ibset(0,sNoOri),sNoNuc)
ICOMP = 1
ISYLBL = 1
Label = 'OneHam'
call RDONE(IRC,IOPT,Label,ICOMP,ONEHAM,ISYLBL)
H1EFF(1:NBTRI) = H1EFF(1:NBTRI)+ONEHAM(:)
call mma_deallocate(ONEHAM)

! Read nuclear repulsion energy:
IRC = -1
IOPT = 0
ICOMP = 0
ISYLBL = 1
call Get_dScalar('PotNuc',PotNuc)

! If this is a perturbative reaction field calculation then
! modify the one-electron Hamiltonian by the reaction field and
! the nuclear attraction by the cavity self-energy
if (RFpert) then
  nTemp = 0
  do iSym=1,nSym
    nTemp = nTemp+nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call f_Inquire('RUNOLD',Found)
  if (Found) call NameRun('RUNOLD')
  call mma_allocate(Temp,nTemp,Label='Temp')
  call Get_dScalar('RF Self Energy',ERFSelf)
  call Get_dArray('Reaction field',Temp,nTemp)
  if (Found) call NameRun('#Pop')
  PotNuc = PotNuc+ERFself
  H1EFF(1:nTemp) = H1EFF(1:nTemp)+Temp(:)
  call mma_deallocate(Temp)
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' 1-EL HAMILTONIAN (MAY INCLUDE REACTION FIELD)'
ISTLT = 0
do ISYM=1,NSYM
  if (NBAS(ISYM) > 0) then
    write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
    call TRIPRT(' ',' ',H1EFF(1+ISTLT),NBAS(ISYM))
    ISTLT = ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
  end if
end do
#endif

#ifdef _OFEmbed_
! If this is a perturbative Orbital-Free Embedding (OFE) calculation
! then modify the one-electron Hamiltonian by the OFE potential and
! the nuclear attraction by the Rep_EN
if (Do_OFEmb) then
  nTemp = 0
  do iSym=1,nSym
    nTemp = nTemp+nBas(iSym)*(nBas(iSym)+1)/2
  end do
  call mma_allocate(Coul,nTemp,Label='Coul')
  Coul(:) = Zero
  if (OFE_First) then
    call mma_allocate(FMaux,nTemp,Label='FMaux')
    call Coul_DMB(.true.,1,Rep_EN,FMaux,Coul,Coul,nTemp)
  end if
  H1EFF(1:nTemp) = H1EFF(1:nTemp)+FMaux(:)
  call mma_deallocate(Coul)
  OFE_First = .false.

  call NameRun('AUXRFIL') ! switch the RUNFILE name
  call Get_dExcdRa(Vxc,nVxc)
  H1EFF(1:nTemp) = H1EFF(1:nTemp)+Vxc(1:nTemp)
  if (nVxc == 2*nTemp) then ! but fix for Nuc Attr added twice
    H1EFF(1:nTemp) = H1EFF(1:nTemp)+Vxc(1+nTemp:2*nTemp)
    call Get_dArray('Nuc Potential',Vxc,nTemp)
    H1EFF(1:nTemp) = H1EFF(1:nTemp)-Vxc(1:nTemp)
  end if
  call mma_deallocate(Vxc)
  call NameRun('#Pop')    ! switch back to old RUNFILE
# ifdef _DEBUGPRINT_
  write(u6,*) ' 1-EL HAMILTONIAN INCLUDING OFE POTENTIAL'
  ISTLT = 0
  do ISYM=1,NSYM
    if (NBAS(ISYM) > 0) then
      write(u6,'(6X,A,I2)') ' SYMMETRY SPECIES:',ISYM
      call TRIPRT(' ',' ',H1EFF(1+ISTLT),NBAS(ISYM))
      ISTLT = ISTLT+NBAS(ISYM)*(NBAS(ISYM)+1)/2
    end if
  end do
# endif
end if
#endif

end subroutine ADD1HAM
