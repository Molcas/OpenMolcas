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

subroutine Print_MPPROP(rMP,xrMP,xnrMP,nij,nElem,lMax,EC,Polar,Lbl,nAtoms,iANr,NoField,CoC,Coor,nOcOb,Energy_Without_FFPT,Ene_Occ, &
                        MpProp_Level,Bond_Threshold,nReal_Centers)

!  rMP : Multipole moments moved to center of charge
! xrMP : Multipole moments kept on the expansion center

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nij, nElem, lMax, nAtoms, iANr(nAtoms), nOcOb, MpProp_Level, nReal_Centers
real(kind=wp), intent(in) :: rMP(nij,nElem), xrMP(nij,nElem), EC(3,nij), Polar(6,nij), CoC(3), Coor(3,nAtoms), &
                             Energy_Without_FFPT, Ene_Occ(nOcOb), Bond_Threshold
real(kind=wp), intent(inout) :: xnrMP(nij,nElem)
character(len=LenIn), intent(in) :: Lbl(nAtoms)
logical(kind=iwp), intent(in) :: NoField
integer(kind=iwp) :: i, iAtom, iElem, iEnd, ii, ij, iSize, iStrt, iUnit, j, jAtom, l, Last_NonBlank, lMax_for_size, MxMP, nBas(8), &
                     nCenters, nSym
real(kind=wp) :: Polar_M(6)
integer(kind=iwp), parameter :: lHeader = 144
character(len=lHeader) :: Header
character(len=8) :: Method, Label
character(len=6) :: fName
logical(kind=iwp) :: Exists, Text, Bond_OK, Check_Bond, Found
logical(kind=iwp), parameter :: Use_Two_Centers_For_Atoms = .false.
real(kind=wp), allocatable :: Scratch_ele(:), Scratch_nuc(:)
real(kind=wp), external :: DDot_

MxMP = lMax
if (lMax > MpProp_Level) MxMP = MpProp_Level

call Qpg_carray('Relax Method',Found,iSize)
if (Found) then
  call Get_cArray('Relax Method',Method,8)
else
  Method = 'UNDEF   '
end if
call Get_cArray('Seward Title',Header,lHeader)
j = 0
Text = .false.
Last_NonBlank = 0
do i=1,72
  if (Text .or. (Header(i:i) /= ' ')) then
    Text = .true.
    j = j+1
    Header(j:j) = Header(i:i)
    if (Header(i:i) /= ' ') Last_NonBlank = j
  end if
end do

fName = 'MPPROP'
iUnit = 11
call OpnFl(fName,iUnit,Exists)

nCenters = nReal_Centers
if (Use_Two_Centers_For_Atoms) then
  nCenters = nCenters+nAtoms
end if

write(iUnit,'(A)') '**************************************************'
write(iUnit,'(A)') '* Molecule'
write(iUnit,'(A)') Header(1:Last_NonBlank)
write(iUnit,'(A)') '* Method'
write(iUnit,'(A)') Method
write(iUnit,'(A)') '* Level of Multipoles and Polarizabilities'
write(iUnit,'(2I5)') MxMP,1
write(iUnit,'(A)') '* All centers'
write(iUnit,'(I5)') nCenters

! Insert atoms
do iAtom=1,nAtoms
  ii = iAtom*(iAtom+1)/2
  call ReExpand(xnrMP,nij,nElem,CoC,EC(1,ii),ii,lMax)

  ! Informations on the atom
  write(iUnit,'(2I5,4X,A)') 2,1,Lbl(iAtom)
  ! The expansion center
  write(iUnit,'(3F20.10)') (EC(iElem,ii),iElem=1,3)
  ! The multipole moments
  iStrt = 1
  do l=0,MxMP
    iEnd = iStrt+(l+1)*(l+2)/2-1
    write(iUnit,'(3F20.10)') (xrMP(ii,iElem)+xnrMP(ii,iElem),iElem=iStrt,iEnd)
    iStrt = iEnd+1
  end do
  ! The size parameters (never higher than quadrupole)
  lMax_for_size = min(MxMP,2)
  iStrt = 1
  do l=0,lMax_for_size
    iEnd = iStrt+(l+1)*(l+2)/2-1
    write(iUnit,'(3F20.10)') (Zero,iElem=iStrt,iEnd)
    iStrt = iEnd+1
  end do
  ! Polarizabilities
  if (NoField) then
    write(iUnit,'(3F20.10)') (Zero,iElem=1,6)
  else
    Polar_M(1) = Polar(1,ii)
    Polar_M(2) = Polar(2,ii)
    Polar_M(3) = Polar(4,ii)
    Polar_M(4) = Polar(3,ii)
    Polar_M(5) = Polar(5,ii)
    Polar_M(6) = Polar(6,ii)
    write(iUnit,'(3F20.10)') (Polar_M(iElem),iElem=1,6)
  end if

  call ReExpand(xnrMP,nij,nElem,EC(1,ii),CoC,ii,lMax)
end do
if (Use_Two_Centers_For_Atoms) then
  ! Insert real atoms and fake values
  do iAtom=1,nAtoms

    ! Informations on the center
    write(iUnit,'(2I5,4X,A)') iANr(iAtom),1,Lbl(iAtom)
    ! The expansion center
    write(iUnit,'(3F20.10)') (Coor(iElem,iAtom),iElem=1,3)
    ! The multipole moments
    iStrt = 1
    do l=0,MxMP
      iEnd = iStrt+(l+1)*(l+2)/2-1
      write(iUnit,'(3F20.10)') (Zero,iElem=iStrt,iEnd)
      iStrt = iEnd+1
    end do
    ! The size parameters (never higher than quadrupole)
    lMax_for_size = min(MxMP,2)
    iStrt = 1
    do l=0,lMax_for_size
      iEnd = iStrt+(l+1)*(l+2)/2-1
      write(iUnit,'(3F20.10)') (One,iElem=iStrt,iEnd)
      iStrt = iEnd+1
    end do
    ! Polarizabilities
    write(iUnit,'(3F20.10)') (Zero,iElem=1,6)

  end do
end if
! Insert bonds
do jAtom=1,nAtoms
  do iAtom=jAtom+1,nAtoms
    ij = iAtom*(iAtom-1)/2+jAtom

    Bond_OK = Check_Bond(Coor(1,iAtom),Coor(1,jAtom),iANr(iAtom),iANr(jAtom),Bond_Threshold)
    if (Bond_OK) then
      Label = '       '
      write(Label,'(I3,A,I3)') jAtom,'-',iAtom
      j = 0
      do i=1,8
        if (Label(i:i) /= ' ') then
          j = j+1
          Label(j:j) = Label(i:i)
        end if
      end do
      ! Informations on the center
      write(iUnit,'(2I5,4X,A)') 2,1,Label(1:j)
      ! The expansion center
      write(iUnit,'(3F20.10)') (EC(iElem,ij),iElem=1,3)
      ! The multipole moments
      iStrt = 1
      do l=0,MxMP
        iEnd = iStrt+(l+1)*(l+2)/2-1
        write(iUnit,'(3F20.10)') (xrMP(ij,iElem),iElem=iStrt,iEnd)
        iStrt = iEnd+1
      end do
      ! The size parameters (never higher than quadrupole)
      lMax_for_size = min(MxMP,2)
      iStrt = 1
      do l=0,lMax_for_size
        iEnd = iStrt+(l+1)*(l+2)/2-1
        write(iUnit,'(3F20.10)') (One,iElem=iStrt,iEnd)
        iStrt = iEnd+1
      end do
      ! Polarizabilities
      if (NoField) then
        write(iUnit,'(3F20.10)') (Zero,iElem=1,6)
      else
        Polar_M(1) = Polar(1,ij)
        Polar_M(2) = Polar(2,ij)
        Polar_M(3) = Polar(4,ij)
        Polar_M(4) = Polar(3,ij)
        Polar_M(5) = Polar(5,ij)
        Polar_M(6) = Polar(6,ij)
        write(iUnit,'(3F20.10)') (Polar_M(iElem),iElem=1,6)
      end if
    end if

  end do
end do

!*** Molecular properties

write(iUnit,'(A)') '* Molecule properties'
call mma_allocate(Scratch_ele,nElem,label='Scratch_ele')
call mma_allocate(Scratch_nuc,nElem,label='Scratch_nuc')
do iElem=1,nElem
  Scratch_ele(iElem) = sum(rMP(:,iElem))
  Scratch_nuc(iElem) = sum(xnrMP(:,iElem))
end do
iStrt = 1
do l=0,MxMP
  iEnd = iStrt+(l+1)*(l+2)/2-1
  write(iUnit,'(3F20.10)') (Scratch_ele(iElem)+Scratch_nuc(iElem),iElem=iStrt,iEnd)
  iStrt = iEnd+1
end do
call mma_deallocate(Scratch_ele)
call mma_deallocate(Scratch_nuc)

if (NoField) then
  write(iUnit,'(3F20.10)') (Zero,iElem=1,6)
else
  Polar_M(1) = DDot_(nij,[One],0,Polar(1,1),6)
  Polar_M(2) = DDot_(nij,[One],0,Polar(2,1),6)
  Polar_M(3) = DDot_(nij,[One],0,Polar(4,1),6)
  Polar_M(4) = DDot_(nij,[One],0,Polar(3,1),6)
  Polar_M(5) = DDot_(nij,[One],0,Polar(5,1),6)
  Polar_M(6) = DDot_(nij,[One],0,Polar(6,1),6)
  write(iUnit,'(3F20.10)') (Polar_M(iElem),iElem=1,6)
end if

!*** Orbital information

call Get_iScalar('nSym',nSym)
call Get_iArray('nBas',nBas,nSym)

write(iUnit,'(A)') '* Orbital information'
write(iUnit,'(2I5)') nBas(1),nOcOb
write(iUnit,'(F20.10)') Energy_Without_FFPT
write(iUnit,'(3F20.10)') (Ene_Occ(i),i=1,nOcOb)

close(iUnit)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Print_MPPROP
