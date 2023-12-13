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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_PrInp()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Print RPA configuration after input processing.

use RPA_globals, only: CMO, DFTFunctional, EMO, iPrint, nBas, nDel, nFro, nOcc, nOrb, nSym, nTitle, nVir, OccEn, Reference, Title, &
                       VirEn
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iUHF, lLine, nLine, l_orbitals, i, j, k, left, iSym, nB, l_Name, iCount
real(kind=wp) :: Dummy(1)
character(len=3) :: lIrrep(8)
character(len=7) :: spin(2)
character(len=8) :: Fmt1, Fmt2
character(len=13) :: orbitals
character(len=120) :: Line, BlLine, StLine
character(len=LenIn8), allocatable :: UBName(:)
integer(kind=iwp), parameter :: lPaper = 132
character(len=*), parameter :: SecNam = 'RPA_PrInp'
integer(kind=iwp), external :: RPA_iUHF

! set restricted(1)/unrestricted(2)
iUHF = RPA_iUHF()

! set labels
if (iUHF == 1) then
  orbitals = 'orbitals'
  l_orbitals = 8
  spin(1) = ' '
  spin(2) = ' '
else if (iUHF == 2) then
  orbitals = 'spin-orbitals'
  l_orbitals = 13
  spin(1) = '(alpha)'
  spin(2) = '(beta)'
else
  write(u6,'(A,I6)') 'iUHF=',iUHF
  call RPA_Warn(3,SecNam//': iUHF error')
  orbitals = ' '
  l_orbitals = 1
  spin(1) = ' '
  spin(2) = ' '
end if

! set irrep names
call Get_cArray('Irreps',lIrrep,24)
do iSym=1,nSym
  lIrrep(iSym) = adjustr(lIrrep(iSym))
end do

! init blank line and "star" line
lLine = len(Line)
do i=1,lLine
  BlLine(i:i) = ' '
  StLine(i:i) = '*'
end do
left = (lPaper-lLine)/2
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'

! print title from input
if (nTitle > 0) then
  write(u6,*)
  nLine = nTitle+5
  do i=1,nLine
    Line = BlLine
    if ((i == 1) .or. (i == nLine)) Line = StLine
    if (i == 3) Line = 'Title:'
    if ((i >= 4) .and. (i <= nLine-2)) Line = Title(i-3)
    call Center_Text(Line)
    write(u6,Fmt1) '*'//Line//'*'
  end do
  write(u6,*)
end if

! print coordinates of the molecule
if (iPrint >= 2) then
  call PrCoor()
end if

! print orbital info
if (iPrint >= 2) then
  write(u6,*)
  write(u6,Fmt2//'A,2(1X,A))') Reference,'reference',orbitals
  j = len(Reference)+11+l_orbitals
  write(u6,Fmt2//'A)') repeat('-',j)
  if (Reference(2:3) == 'KS') then
    write(u6,Fmt2//'A,A)') 'DFT functional: ',DFTFunctional
  end if
  write(u6,*)
  write(u6,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Number of orbitals',(nOrb(iSym),iSym=1,nSym)
  do k=1,iUHF
    write(u6,Fmt2//'A,2(1X,A),T47,8I4)') 'Frozen occupied',orbitals,spin(k),(nFro(iSym,k),iSym=1,nSym)
  end do
  do k=1,iUHF
    write(u6,Fmt2//'A,2(1X,A),T47,8I4)') 'Active occupied',orbitals,spin(k),(nOcc(iSym,k),iSym=1,nSym)
  end do
  do k=1,iUHF
    write(u6,Fmt2//'A,2(1X,A),T47,8I4)') 'Active virtual',orbitals,spin(k),(nVir(iSym,k),iSym=1,nSym)
  end do
  do k=1,iUHF
    write(u6,Fmt2//'A,2(1X,A),T47,8I4)') 'Deleted virtual',orbitals,spin(k),(nDel(iSym,k),iSym=1,nSym)
  end do
end if

! print orbital energies
if (iPrint >= 2) then
  iCount = 0
  do k=1,iUHF
    do iSym=1,nSym
      iCount = iCount+nFro(iSym,k)
    end do
  end do
  if (iCount > 0) then
    write(u6,*)
    write(u6,*)
    write(u6,Fmt2//'A,1X,A,T47)') 'Energies of the frozen occupied',orbitals
    do k=1,iUHF
      i = 0
      do iSym=1,nSym
        if (nFro(iSym,k) > 0) then
          write(u6,*)
          write(u6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))') 'symmetry species',iSym,lIrrep(iSym),spin(k),(OccEn(i+j,k),j=1,nFro(iSym,k))
          i = i+nFro(iSym,k)+nOcc(iSym,k)
        end if
      end do
    end do
  end if
  write(u6,*)
  write(u6,*)
  write(u6,Fmt2//'A,1X,A,T47)') 'Energies of the active occupied',orbitals
  do k=1,iUHF
    i = 0
    do iSym=1,nSym
      if (nOcc(iSym,k) > 0) then
        write(u6,*)
        write(u6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))') &
          'symmetry species',iSym,lIrrep(iSym),spin(k),(OccEn(i+nFro(iSym,k)+j,k),j=1,nOcc(iSym,k))
        i = i+nFro(iSym,k)+nOcc(iSym,k)
      end if
    end do
  end do
  write(u6,*)
  write(u6,*)
  write(u6,Fmt2//'A,1X,A,T47)') 'Energies of the active virtual',orbitals
  do k=1,iUHF
    i = 0
    do iSym=1,nSym
      if (nVir(iSym,k) > 0) then
        write(u6,*)
        write(u6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))') 'symmetry species',iSym,lIrrep(iSym),spin(k),(VirEn(i+j,k),j=1,nVir(iSym,k))
        i = i+nVir(iSym,k)+nDel(iSym,k)
      end if
    end do
  end do
  iCount = 0
  do k=1,iUHF
    do iSym=1,nSym
      iCount = iCount+nDel(iSym,k)
    end do
  end do
  if (iCount > 0) then
    write(u6,*)
    write(u6,*)
    write(u6,Fmt2//'A,1X,A,T47)') 'Energies of the deleted virtual',orbitals
    do k=1,iUHF
      i = 0
      do iSym=1,nSym
        if (nDel(iSym,k) > 0) then
          write(u6,*)
          write(u6,Fmt2//'A,I2,2(1X,A),(T40,5F14.6))') &
            'symmetry species',iSym,lIrrep(iSym),spin(k),(VirEn(i+nVir(iSym,k)+j,k),j=1,nDel(iSym,k))
          i = i+nVir(iSym,k)+nDel(iSym,k)
        end if
      end do
    end do
  end if
end if

! print orbitals
if (iPrint >= 2) then
  nB = nBas(1)
  do iSym=2,nSym
    nB = nB+nBas(iSym)
  end do
  l_Name = nB
  call mma_allocate(UBName,l_Name,label='Name')
  call Get_cArray('Unique Basis Names',UBName,LenIn8*nB)
  do k=1,iUHF
    call PriMO(Reference//' '//orbitals//' '//spin(k),.false.,.true.,-9.9e9_wp,9.9e9_wp,nSym,nBas,nOrb,UBName,EMO(:,k),Dummy, &
               CMO(:,k),-1)
  end do
  call mma_deallocate(UBName)
end if

! flush output buffer
call xFlush(u6)

end subroutine RPA_PrInp
