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
! Copyright (C) 1991, Roland Lindh                                     *
!               2018, Sijia S. Dong                                    *
!***********************************************************************

subroutine Prpt()
!***********************************************************************
!                                                                      *
! Purpose: To set up all calling arguments for the subroutine          *
!          Prpt_ . For RASSCF to work MO coefficents and Occupation    *
!          numbers must be available on TMPORB.                        *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: iDummy(1), iError, iSym, iUHF, iWFType, Lu, n2Tot, nBas(8), nDim, nIrrep, nTriDim
real(kind=wp) :: Dummy(1)
logical(kind=iwp) :: ifallorb, Short, var
character(len=81) :: note
character(len=8) :: Method
character(len=4) :: PrpLst
character(len=2) :: lbl
real(kind=wp), allocatable :: Occ(:,:), Vec(:,:)
integer(kind=iwp), external :: isFreeUnit

call GetEnvf('MOLCAS_PROPERTIES',PrpLst)
call UpCase(PrpLst)
if (PrPlst(1:3) == 'LON') then
  Short = .false.
  !ifallorb = .True.
else
  Short = .true.
  ifallorb = .false.
end if

! This variable is used so we know if the density we search for is labeled
! variational or not.
var = .false.

call Get_cArray('Relax Method',Method,8)

call Get_iScalar('nSym',nIrrep)

call Get_iArray('nBas',nBas,nIrrep)

nDim = 0
nTriDim = 0
n2Tot = 0
do iSym=1,nIrrep
  nDim = nDim+nBas(iSym)
  nTriDim = nTriDim+nBas(iSym)*(nBas(iSym)+1)/2
  n2Tot = n2Tot+nBas(iSym)**2
end do

if ((Method == 'RHF-SCF ') .or. &
    (Method == 'IVO-SCF ') .or. &
    (Method == 'KS-DFT  ') .or. &
    (Method == 'UHF-SCF ')) then
  call Get_iScalar('SCF mode',iUHF)
else
  iUHF = 0
end if

if ((iUHF == 1) .or. (Method == 'RASSCFSA')) then
  call mma_allocate(Occ,nDim,2,label='Occ')
else
  call mma_allocate(Occ,nDim,1,label='Occ')
end if
if (Short) then
  call mma_allocate(Vec,0,2,label='Vec')
  lbl = 'O '
  n2Tot = 0
else
  if ((iUHF /= 1) .or. (Method == 'RASSCFSA')) then
    call mma_allocate(Vec,n2Tot,1,label='Vec')
  else
    call mma_allocate(Vec,n2Tot,2,label='Vec')
  end if
  lbl = 'CO'
end if

Lu = 10
Lu = IsFreeUnit(Lu)
if ((Method == 'RHF-SCF ') .or. &
    (Method == 'IVO-SCF ') .or. &
    (Method == 'KS-DFT  ') .or. &
    (Method == 'UHF-SCF ')) then
  if (iUHF /= 1) then
    call RdVec('SCFORB',Lu,Lbl,nIrrep,nBas,nBas,Vec(:,1),Occ(:,1),Dummy,iDummy,note,0,iError)
  else
    call RdVec_('UHFORB',Lu,Lbl,iUHF,nIrrep,nBas,nBas,Vec(:,1),Vec(:,2),Occ(:,1),Occ(:,2),Dummy,Dummy,iDummy,note,1,iError,iWFtype)
    if (Short) then
      Occ(:,1) = Occ(:,1)+Occ(:,2)
    end if
  end if
else if ((Method == 'RASSCF  ') .or. &
         (Method == 'CASSCF  ') .or. &
         (Method == 'CASDFT  ') .or. &
         (Method == 'CASSCFSA') .or. &
         (Method == 'CASPT2  ') .or. &
         (Method == 'RASSCFSA')) then
  call RdVec('TMPORB',Lu,Lbl,nIrrep,nBas,nBas,Vec(:,1),Occ(:,1),Dummy,iDummy,note,0,iError)
  if (Note(2:4) == 'var') var = .true.
else if (Method == 'MBPT2   ') then
  ! MBPT2 has no occupation-numbers at the moment.
  Occ(:,:) = Zero
  var = .true.
else
  write(u6,'(6X,2A)') 'Properties not supported for ',Method
end if

call Prpt_(nIrrep,nBas,nDim,Occ,n2Tot,Vec,var,Short,iUHF,ifallorb)

call mma_deallocate(Occ)
call mma_deallocate(Vec)

return

end subroutine Prpt
