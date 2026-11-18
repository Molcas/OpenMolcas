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

subroutine UMP2_Read_Reference(ierr,Message)

use Definitions, only: wp, iwp
use MBPT2_Global, only: MBPT2_All_Finite
use stdalloc, only: mma_allocate, mma_deallocate
use UMP2_Global, only: UMP2_Clean, ReferenceMethod, &
                       nSym, nBas, nOrb, nOccA, nOccB, &
                       CAlpha, CBeta, EOrbA, EOrbB, ESCF

implicit none

integer(kind=iwp), intent(out) :: ierr
character(len=*), intent(out) :: Message

integer(kind=iwp) :: nFrozenSCF, nCoeff, nFull
logical(kind=iwp) :: Found
real(kind=wp), allocatable :: Buffer(:)

ierr = 0
Message = ''

! The reader starts a fresh UMP2 calculation.
! Input options must be parsed after this call.

call UMP2_Clean()

! Relax Method is a required field in the existing MBPT2 driver.
! Require a direct ordinary UHF-SCF reference for this first version.

call Get_cArray('Relax Method',ReferenceMethod,8)

if (ReferenceMethod(1:7) /= 'UHF-SCF') then
  call Fail('UMP2 requires a direct UHF-SCF reference.')
  return
end if

call Qpg_iScalar('nSym',Found)
if (.not. Found) then
  call Fail('Missing RUNFILE integer scalar: nSym')
  return
end if

call Get_iScalar('nSym',nSym)

if (nSym /= 1) then
  call Fail('This UMP2 implementation requires C1 symmetry.')
  return
end if

call ReadCount('nBas',nBas)
if (ierr /= 0) return

call ReadCount('nOrb',nOrb)
if (ierr /= 0) return

call ReadCount('nIsh',nOccA)
if (ierr /= 0) return

! FinalSCF writes the final beta occupation counts under nIsh_ab.

call ReadCount('nIsh_ab',nOccB)
if (ierr /= 0) return

call ReadCount('nFro',nFrozenSCF)
if (ierr /= 0) return

! As in the existing RdMBPT reader, exclude SCF-frozen orbitals.
! Post-SCF correlation freezing is handled separately.

if (nFrozenSCF /= 0) then
  call Fail('Orbitals frozen during SCF are not supported.')
  return
end if

if ((nBas <= 0) .or. (nOrb <= 0)) then
  call Fail('Invalid basis or retained-orbital dimension.')
  return
end if

if (nOrb > nBas) then
  call Fail('The retained-orbital count exceeds the basis size.')
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

! Prevent overflow when forming the full coefficient-record length.

if (nBas > huge(nBas)/nBas) then
  call Fail('Coefficient-record length exceeds integer capacity.')
  return
end if

nCoeff = nBas*nOrb
nFull = nBas*nBas

call Qpg_dScalar('SCF energy',Found)
if (.not. Found) then
  call Fail('Missing RUNFILE scalar: SCF energy')
  return
end if

call Get_dScalar('SCF energy',ESCF)

if (.not. MBPT2_All_Finite(ESCF)) then
  call Fail('The SCF energy is not finite.')
  return
end if

! In C1, retained coefficient columns occupy the leading
! nBas*nOrb entries. FinalSCF can write an nBas*nBas record.

call ReadRealRecord('SCF orbitals',nCoeff,nFull,Buffer)
if (ierr /= 0) return

call mma_allocate(CAlpha,nBas,nOrb,label='UMP2 CAlpha')
CAlpha(:,:) = reshape(Buffer(1:nCoeff),shape(CAlpha))
call mma_deallocate(Buffer)

call ReadRealRecord('SCF orbitals_ab',nCoeff,nFull,Buffer)
if (ierr /= 0) return

call mma_allocate(CBeta,nBas,nOrb,label='UMP2 CBeta')
CBeta(:,:) = reshape(Buffer(1:nCoeff),shape(CBeta))
call mma_deallocate(Buffer)

! Orbital-energy records can similarly contain trailing entries.
! Only the retained nOrb energies belong to the MO spaces.

call ReadRealRecord('OrbE',nOrb,nBas,Buffer)
if (ierr /= 0) return

call mma_allocate(EOrbA,nOrb,label='UMP2 EOrbA')
EOrbA(:) = Buffer(1:nOrb)
call mma_deallocate(Buffer)

call ReadRealRecord('OrbE_ab',nOrb,nBas,Buffer)
if (ierr /= 0) return

call mma_allocate(EOrbB,nOrb,label='UMP2 EOrbB')
EOrbB(:) = Buffer(1:nOrb)
call mma_deallocate(Buffer)

contains

subroutine ReadCount(Label,Value)

  character(len=*), intent(in) :: Label
  integer(kind=iwp), intent(out) :: Value

  integer(kind=iwp) :: nData
  integer(kind=iwp) :: Values(8)
  logical(kind=iwp) :: Exists

  Value = 0
  nData = 0
  Values(:) = 0

  call Qpg_iArray(Label,Exists,nData)

  if (.not. Exists) then
    call Fail('Missing RUNFILE integer array: '//trim(Label))
    return
  end if

  ! Accept either a C1 record or a padded symmetry array.
  ! Only the first element applies after the nSym == 1 check.

  if ((nData < 1) .or. (nData > size(Values))) then
    call Fail('Unexpected RUNFILE array length: '//trim(Label))
    return
  end if

  call Get_iArray(Label,Values,nData)
  Value = Values(1)

end subroutine ReadCount

subroutine ReadRealRecord(Label,nKeep,nPadded,Values)

  character(len=*), intent(in) :: Label
  integer(kind=iwp), intent(in) :: nKeep,nPadded
  real(kind=wp), allocatable, intent(out) :: Values(:)

  integer(kind=iwp) :: nData
  logical(kind=iwp) :: Exists

  nData = 0
  call Qpg_dArray(Label,Exists,nData)

  if (.not. Exists) then
    call Fail('Missing RUNFILE real array: '//trim(Label))
    return
  end if

  if ((nData /= nKeep) .and. (nData /= nPadded)) then
    call Fail('Unexpected RUNFILE array length: '//trim(Label))
    return
  end if

  call mma_allocate(Values,nData,label='UMP2 read buffer')
  call Get_dArray(Label,Values,nData)

  ! Trailing padding is not part of the retained reference.
  ! Check only entries that will actually be used.

  if (.not. MBPT2_All_Finite(Values(1:nKeep))) then
    call mma_deallocate(Values)
    call Fail('Nonfinite RUNFILE data: '//trim(Label))
    return
  end if

end subroutine ReadRealRecord

subroutine Fail(Text)

  character(len=*), intent(in) :: Text

  ierr = 1
  Message = Text

  ! Discard any partially loaded reference.
  call UMP2_Clean()

end subroutine Fail

end subroutine UMP2_Read_Reference