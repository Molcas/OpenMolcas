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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine RPA_RdRun()

! Thomas Bondo Pedersen
!
! Read data from Runfile.

use RPA_globals, only: DFTFunctional, nBas, nDel, nFro, nOcc, nOrb, nSym, NuclearRepulsionEnergy, nVir, Reference
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iUHF, iSym, i
logical(kind=iwp) :: Warn
character(len=8) :: Model
character(len=*), parameter :: SecNam = 'RPA_RdRun'
integer(kind=iwp), external :: RPA_iUHF

! Set type of SCF reference wave function
! Note: in RPA, iUHF=1 means restricted, 2 means unrestricted.
call Get_cArray('Relax Method',Model,8)
call Get_iScalar('SCF mode',iUHF)
iUHF = iUHF+1 ! correct for convention in SCF
if (Model(1:7) == 'RHF-SCF') then
  if (iUHF == 1) then
    Reference = 'RHF'
    Warn = .false.
  else
    Reference = 'UHF'
    warn = .true.
  end if
else if (Model(1:7) == 'UHF-SCF') then
  if (iUHF == 1) then
    Reference = 'RHF'
    Warn = .true.
  else
    Reference = 'UHF'
    Warn = .false.
  end if
else if (Model(1:6) == 'KS-DFT') then
  Warn = .false.
  if (iUHF == 1) then
    Reference = 'RKS'
  else
    Reference = 'UKS'
  end if
else if ((Model(1:8) == 'dRPA@RHF') .or. (Model(1:8) == 'SOSX@RHF')) then
  if (iUHF == 1) then
    Reference = 'RHF'
    Warn = .false.
  else
    Reference = 'UHF'
    Warn = .true.
  end if
else if ((Model(1:8) == 'dRPA@UHF') .or. (Model(1:8) == 'SOSX@UHF')) then
  if (iUHF == 1) then
    Reference = 'RHF'
    Warn = .true.
  else
    Reference = 'UHF'
    Warn = .false.
  end if
else if ((Model(1:8) == 'dRPA@RKS') .or. (Model(1:8) == 'SOSX@RKS')) then
  if (iUHF == 1) then
    Reference = 'RKS'
    Warn = .false.
  else
    Reference = 'UKS'
    Warn = .true.
  end if
else if ((Model(1:8) == 'dRPA@UKS') .or. (Model(1:8) == 'SOSX@UKS')) then
  if (iUHF == 1) then
    Reference = 'RKS'
    Warn = .true.
  else
    Reference = 'UKS'
    Warn = .false.
  end if
else
  write(u6,'(A,A)') 'Reference model from Runfile: ',Model
  write(u6,'(A,I8)') 'iUHF from Runfile:            ',iUHF-1
  call RPA_Warn(2,'Illegal reference wave function in RPA')
  Reference = 'Non'
  Warn = .false.
end if
if (Warn) then
  call RPA_Warn(1,'Runfile restricted/unrestricted conflict in RPA')
  write(u6,'(A,A)') 'Reference model from Runfile: ',Model
  write(u6,'(A,I8)') 'iUHF from Runfile:            ',iUHF-1
  write(u6,'(A,A,A)') 'Assuming ',Reference,' reference!'
  call xFlush(u6)
end if

! Get nuclear potential energy
call Get_dScalar('PotNuc',NuclearRepulsionEnergy(1))

! Get DFT functional
if (Reference(2:3) == 'KS') then
  call Get_cArray('DFT functional',DFTFunctional,80)
else
  DFTFunctional = 'Hartree-Fock'
end if

! Get number of irreps
call Get_iScalar('nSym',nSym)
if ((nSym < 1) .or. (nSym > 8)) then
  call RPA_Warn(3,'nSym out of bonds in RPA')
end if

! Get iUHF (RPA convention: 1->RHF, 2->UHF)
!          (as opposed to SCF: 0->RHF, 1->UHF)
iUHF = RPA_iUHF()

! Get number of basis functions, orbitals, occupied,
! frozen (in SCF), deleted
call Get_iArray('nBas',nBas,nSym)
call Get_iArray('nOrb',nOrb,nSym)
call Get_iArray('nDel',nDel(1,1),nSym)
call Get_iArray('nFro',nFro(1,1),nSym)
call Get_iArray('nIsh',nOcc(1,1),nSym)
if (iUHF == 2) then
  ! unrestricted: read data for beta spin
  call Get_iArray('nIsh_ab',nOcc(1,2),nSym)
end if

! Check for orbitals frozen in SCF and check consistency.
do iSym=1,nSym
  if (nFro(iSym,1) /= 0) then
    write(u6,'(A,8I8)') 'nFro=',(nFro(i,1),i=1,nSym)
    call RPA_Warn(4,SecNam//': Some orbitals were frozen in SCF!')
  end if
  if (nDel(iSym,1) /= (nBas(iSym)-nOrb(iSym))) then
    write(u6,'(A,8I8)') 'nBas=     ',(nBas(i),i=1,nSym)
    write(u6,'(A,8I8)') 'nOrb=     ',(nOrb(i),i=1,nSym)
    write(u6,'(A,8I8)') 'nBas-nOrb=',((nBas(i)-nOrb(i)),i=1,nSym)
    write(u6,'(A,8I8)') 'nDel=     ',(nDel(i,1),i=1,nSym)
    call RPA_Warn(4,SecNam//': nDel != nBas-nOrb')
  end if
end do

! Set default frozen (core) orbitals
call Get_iArray('Non valence orbitals',nFro(1,1),nSym)
do iSym=1,nSym
  if (nFro(iSym,1) > nOcc(iSym,1)) then
    if (iUHF == 1) then
      write(u6,'(A,8I8)') 'nOcc=',(nOcc(i,1),i=1,nSym)
      write(u6,'(A,8I8)') 'nFro=',(nFro(i,1),i=1,nSym)
      call RPA_Warn(4,SecNam//': nFro > nOrb')
    else
      write(u6,'(A,8I8)') 'nOcc(alpha)=',(nOcc(i,1),i=1,nSym)
      write(u6,'(A,8I8)') 'nFro(alpha)=',(nFro(i,1),i=1,nSym)
      call RPA_Warn(4,SecNam//': nFro > nOrb [alpha]')
    end if
  end if
end do
if (iUHF == 2) then
  do iSym=1,nSym
    nFro(iSym,2) = nFro(iSym,1)
    if (nFro(iSym,2) > nOcc(iSym,2)) then
      write(u6,'(A,8I8)') 'nOcc(beta)=',(nOcc(i,2),i=1,nSym)
      write(u6,'(A,8I8)') 'nFro(beta)=',(nFro(i,2),i=1,nSym)
      call RPA_Warn(4,SecNam//': nFro > nOrb [beta]')
    end if
  end do
else
  do iSym=1,nSym
    nFro(iSym,2) = 0
  end do
end if

! Compute number of virtual orbitals
do i=1,iUHF
  do iSym=1,nSym
    nVir(iSym,i) = nOrb(iSym)-nOcc(iSym,i)
  end do
end do

! Set default deleted (virtual) orbitals
nDel(:,:) = 0

end subroutine RPA_RdRun
