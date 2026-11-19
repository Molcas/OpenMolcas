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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

subroutine UMP2_Setup_Spaces(ierr,Message)

use Definitions, only: iwp
use stdalloc, only: mma_allocate, mma_deallocate
use UMP2_Global, only: nSym, nBas, nOrb, nOccA, nOccB, &
                       nFro, ExplicitFrozen, nOA, nOB, nVA, nVB, &
                       CAlpha, CBeta, EOrbA, EOrbB, &
                       EOccA, EOccB, EVirA, EVirB, EnergyReady

implicit none

integer(kind=iwp), intent(out) :: ierr
character(len=*), intent(out) :: Message

ierr = 0
Message = ''

! Any previous energy becomes invalid when spaces are rebuilt.

EnergyReady = .false.

call mma_deallocate(EOccA,safe='*')
call mma_deallocate(EOccB,safe='*')
call mma_deallocate(EVirA,safe='*')
call mma_deallocate(EVirB,safe='*')

nOA = 0
nOB = 0
nVA = 0
nVB = 0

if (nSym /= 1) then
  call Fail('UMP2 orbital-space setup requires C1 symmetry.')
  return
end if

if ((nBas <= 0) .or. (nOrb <= 0) .or. (nOrb > nBas)) then
  call Fail('Invalid basis or retained-orbital dimensions.')
  return
end if

if (.not. ExplicitFrozen) then
  call Fail('Specify Frozen explicitly, including Frozen = 0.')
  return
end if

if ((nOccA < 0) .or. (nOccA > nOrb)) then
  call Fail('Invalid alpha occupied-orbital count.')
  return
end if

if ((nOccB < 0) .or. (nOccB > nOrb)) then
  call Fail('Invalid beta occupied-orbital count.')
  return
end if

if ((nFro < 0) .or. (nFro > min(nOccA,nOccB))) then
  call Fail('Frozen count must lie between zero and min(Nalpha,Nbeta).')
  return
end if

! Check allocation before querying array dimensions.
! Fortran does not guarantee short-circuit logical evaluation.

if (.not. allocated(CAlpha)) then
  call Fail('Alpha coefficient matrix has not been loaded.')
  return
end if

if (.not. allocated(CBeta)) then
  call Fail('Beta coefficient matrix has not been loaded.')
  return
end if

if (.not. allocated(EOrbA)) then
  call Fail('Alpha orbital energies have not been loaded.')
  return
end if

if (.not. allocated(EOrbB)) then
  call Fail('Beta orbital energies have not been loaded.')
  return
end if

if ((size(CAlpha,1) /= nBas) .or. &
    (size(CAlpha,2) /= nOrb)) then
  call Fail('Alpha coefficient-matrix dimensions are inconsistent.')
  return
end if

if ((size(CBeta,1) /= nBas) .or. &
    (size(CBeta,2) /= nOrb)) then
  call Fail('Beta coefficient-matrix dimensions are inconsistent.')
  return
end if

if ((size(EOrbA) /= nOrb) .or. &
    (size(EOrbB) /= nOrb)) then
  call Fail('Orbital-energy dimensions are inconsistent.')
  return
end if

! Freeze the lowest nFro occupied orbitals in each spin set.
! Frozen occupied orbitals never become virtual orbitals.

nOA = nOccA-nFro
nOB = nOccB-nFro

nVA = nOrb-nOccA
nVB = nOrb-nOccB

call mma_allocate(EOccA,nOA,label='UMP2 EOccA')
call mma_allocate(EOccB,nOB,label='UMP2 EOccB')
call mma_allocate(EVirA,nVA,label='UMP2 EVirA')
call mma_allocate(EVirB,nVB,label='UMP2 EVirB')

! Empty spaces are valid, for example the beta occupied space
! for a one-electron alpha-spin reference.

if (nOA > 0) EOccA(:) = EOrbA(nFro+1:nOccA)
if (nOB > 0) EOccB(:) = EOrbB(nFro+1:nOccB)

if (nVA > 0) EVirA(:) = EOrbA(nOccA+1:nOrb)
if (nVB > 0) EVirB(:) = EOrbB(nOccB+1:nOrb)

contains

subroutine Fail(Text)

  character(len=*), intent(in) :: Text

  ierr = 1
  Message = Text

end subroutine Fail

end subroutine UMP2_Setup_Spaces