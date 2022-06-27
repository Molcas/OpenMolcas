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

subroutine GetInfo_Localisation_0()
! Author: T.B. Pedersen
!
! Purpose: read basic info from runfile and INPORB.

use Localisation_globals, only: BName, CMO, EOrb, Ind, LC_FileOrb, nAtoms, nBas, nCMO, nOccInp, nOrb, nSym, nVirInp, Occ
use stdalloc, only: mma_allocate
use Constants, only: Zero
use Definitions, only: iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: i, iSym, kOff, n2Bas, nBasT, nOrbT
character(len=80) :: Txt
character(len=512) :: FName
character(len=*), parameter :: SecNam = 'GetInfo_Localisation_0'

! Read number of irreps.
! ----------------------

call Get_iScalar('nSym',nSym)
if ((nSym < 1) .or. (nSym > MxSym)) then
  write(Txt,'(A,I9)') 'nSym =',nSym
  call SysAbendMsg(SecNam,'Number of irreps out of bounds!',Txt)
end if

! Read number of basis functions.
! -------------------------------

call Get_iArray('nBas',nBas,nSym)
nBasT = nBas(1)
do iSym=2,nSym
  nBasT = nBasT+nBas(iSym)
end do
if ((nBasT < 1) .or. (nBasT > MxBas)) then
  write(Txt,'(A,I9)') 'nBasT =',nBasT
  call SysAbendMsg(SecNam,'Basis set limits exceeded!',Txt)
end if

! Set number of orbitals equal to nBas.
! -------------------------------------

nOrb(:) = nBas(:)
nOrbT = nOrb(1)
do iSym=2,nSym
  nOrbT = nOrbT+nOrb(iSym)
end do
if ((nOrbT < 1) .or. (nOrbT > MxBas)) then
  write(Txt,'(A,I9)') 'nOrbT =',nOrbT
  call SysAbendMsg(SecNam,'Orbital limits exceeded!',Txt)
end if
do iSym=1,nSym
  if (nOrb(iSym) > nBas(iSym)) then
    write(Txt,'(A,I2,2(1X,I9))') 'iSym,nOrb,nBas:',iSym,nOrb(iSym),nBas(iSym)
    call SysAbendMsg(SecNam,'#orb > #bas:',Txt)
  end if
end do

! Read MO coefficients, orbital occupations, and orbital energies
! from INPORB.
! ---------------------------------------------------------------

n2Bas = nBas(1)**2
do iSym=2,nSym
  n2Bas = n2Bas+nBas(iSym)**2
end do

nCMO = n2Bas
call mma_allocate(CMO,nCMO,label='CMO')
call mma_allocate(Occ,nBasT,label='Occup')
call mma_allocate(Eorb,nBasT,label='OrbEn')
call mma_allocate(Ind,nBasT,label='IndT')
FName = LC_FileOrb
if (FName == '') FName = 'INPORB' ! file name
call RdVec_Localisation(nSym,nBas,nOrb,Ind,CMO,Occ,EOrb,trim(FName))

! Set number of occupied and virtual orbitals according to the
! occupation numbers from INPORB (assuming that occupied orbitals
! precede virtual ones on file).
! ---------------------------------------------------------------

kOff = 0
do iSym=1,nSym
  nOccInp(iSym) = 0
  i = 0
  do while (i < nOrb(iSym))
    i = i+1
    if (Occ(kOff+i) > Zero) then
      nOccInp(iSym) = nOccInp(iSym)+1
    else
      i = nOrb(iSym) ! break while loop
    end if
  end do
  nVirInp(iSym) = nOrb(iSym)-nOccInp(iSym)
  if (nVirInp(iSym) < 0) then
    write(Txt,'(3(A,I9))') 'No. of occupied: ',nOccInp(iSym),' No. of orbitals: ',nOrb(iSym),' Symmetry: ',iSym
    call SysAbendMsg(SecNam,'#occ > #orb:',Txt)
  end if
  kOff = kOff+nBas(iSym)
end do

! Read number of atoms, atomic labels, and basis function labels
! from runfile.
! --------------------------------------------------------------

call Get_nAtoms_All(nAtoms)
if ((nAtoms < 1) .or. (nAtoms > MxAtom)) then
  write(Txt,'(A,I9)') 'nAtoms =',nAtoms
  call SysAbendMsg(SecNam,'Atom limit exceeded!',Txt)
end if
call mma_allocate(BName,nBasT,label='BName')
!call Get_cArray('Unique Atom Names',AtomLbl,4*nAtoms)
call Get_cArray('Unique Basis Names',BName,len(BName)*nBasT)

end subroutine GetInfo_Localisation_0
