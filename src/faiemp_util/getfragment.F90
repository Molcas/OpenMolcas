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
! Copyright (C) Ben Swerts                                             *
!               2016, Liviu Ungur                                      *
!***********************************************************************

subroutine GetFragment(lUnit,iCnttp)
!***********************************************************************
!                                                                      *
!    Objective: To read frozen fragment information.                   *
!               This means that we read (from input stream) the unique *
!               centers with their basis sets, the coordinates of all  *
!               atoms of the fragment, the orbital energies, the       *
!               orbital coefficients and the Mulliken charges.         *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
!     Author: Ben Swerts                                               *
!   Modified: Liviu Ungur                                              *
!***********************************************************************

use Basis_Info

implicit none
#include "stdalloc.fh"
#include "angstr.fh"
integer lUnit, iCnttp
character*180 Line, Get_Ln
integer storageSize, LineWords
parameter(storageSize=200,LineWords=storageSize/8)
integer nFragType, nFragCoor, nFragEner, nFragDens
character*(storageSize) sBasis
real*8 eqBasis(LineWords)
integer iPrint, i, j, iBasis, ierr
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Read_v(lunit,Work,istrt,iend,inc,ierr)
    integer lUnit, iStrt, iEnd, Inc, iErr
    real*8 Work(iend)
  end subroutine Read_v
end interface
!LineWords = 25
!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
iPrint = 99
#else
iPrint = 5
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
nFrag_LineWords = LineWords ! Needed in Module Basis_Info
!                                                                      *
!***********************************************************************
!                                                                      *
! Keyword: LBASIS
!                                                                      *
!***********************************************************************
!                                                                      *
! Local basis sets for unique centers
!
if (iPrint >= 99) write(6,*) 'Reading LBASIS'
Line = Get_Ln(lUnit)
if (index(Line,'LBASIS') == 0) then
  write(6,*) 'ERROR: Keyword LBASIS expected, offending line:'
  write(6,*) Line
  call Quit_OnUserError()
end if

! Fragment types

Line = Get_Ln(lUnit)
call Get_i1(1,nFragType)
dbsc(iCnttp)%nFragType = nFragType
call mma_allocate(dbsc(iCnttp)%FragType,LineWords,nFragType,Label='FragType')
if (iPrint >= 99) write(6,*) 'number of LBASIS = ',nFragType
!
!     read the basis sets labels
!
do i=1,nFragType
  sBasis = Get_Ln(lUnit)
  eqBasis = transfer(sBasis,eqBasis) ! ???
  do j=1,LineWords
    dbsc(iCnttp)%FragType(j,i) = eqBasis(j)
  end do
  if (iPrint >= 49) write(6,*) 'GetFragment: basis set ',sBasis
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Keyword: RELCOORDS
!                                                                      *
!***********************************************************************
!                                                                      *
! All atoms: index of the associated basis set and coordinates

if (iPrint >= 99) write(6,*) 'Reading RELCOORDS'

Line = Get_Ln(lUnit)
if (index(Line,'RELCOORDS') == 0) then
  write(6,*) 'ERROR: Keyword RELCOORDS expected, offending line:'
  write(6,*) Line
  call Quit_OnUserError()
end if

Line = Get_Ln(lUnit)
call Get_i1(1,nFragCoor)
dbsc(iCnttp)%nFragCoor = nFragCoor
if (iPrint >= 99) write(6,*) 'number of RELCOORDS = ',nFragCoor

! read all centers, but reserve space for the Mulliken charges
!
! FragCoor(1,i): Index of the FragType
! FragCoor(2,i): x coordinate
! FragCoor(3,i): y coordinate
! FragCoor(4,i): z coordinate
! FragCoor(5,i): Mulliken charge

call mma_allocate(dbsc(iCnttp)%FragCoor,5,nFragCoor,Label='FragCoor')
do i=1,nFragCoor
  Line = Get_Ln(lUnit)
  call Get_i1(1,iBasis)
  dbsc(iCnttp)%FragCoor(1,i) = dble(iBasis)
  call Get_f(2,dbsc(iCnttp)%FragCoor(2,i),3)
  if (index(Line,'ANGSTROM') /= 0) then
    if (iPrint >= 49) write(6,*) 'Reading the relcoords in Angstrom'
    dbsc(iCnttp)%FragCoor(2:4,i) = dbsc(iCnttp)%FragCoor(2:4,i)/Angstr
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! keyword: ENERGIES
!                                                                      *
!***********************************************************************
!                                                                      *
! Orbital energies (taken from the ONE ELECTRON ENERGIES in ScfOrb)

if (iPrint >= 99) write(6,*) 'Reading ENERGIES'
Line = Get_Ln(lUnit)
if (index(Line,'ENERGIES') == 0) then
  write(6,*) 'ERROR: Keyword ENERGIES expected, offending line:'
  write(6,*) Line
  call Quit_OnUserError()
end if

Line = Get_Ln(lUnit)
call Get_i1(1,nFragEner)
dbsc(iCnttp)%nFragEner = nFragEner
call mma_Allocate(dbsc(iCnttp)%FragEner,nFragEner,Label='FragEner')

call Read_v(lUnit,dbsc(iCnttp)%FragEner,1,nFragEner,1,ierr)
if (iPrint >= 99) call RecPrt('Fragment MO energies',' ',dbsc(iCnttp)%FragEner,1,nFragEner)
if (ierr /= 0) then
  write(6,*) 'ERROR: number of energy values is not correct'
  write(6,*) ierr
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! keyword: MOCOEFF
!                                                                      *
!***********************************************************************
!                                                                      *
! MO coefficients (taken from the ORBITALs in ScfOrb)

if (iPrint >= 99) write(6,*) 'Reading MOCOEFF'
Line = Get_Ln(lUnit)
if (index(Line,'MOCOEFF') == 0) then
  write(6,*) 'ERROR: Keyword MOCOEFF expected, offending line:'
  write(6,*) Line
  call Quit_OnUserError()
end if

Line = Get_Ln(lUnit)
call Get_i1(1,nFragDens)
dbsc(iCnttp)%nFragDens = nFragDens
call mma_allocate(dbsc(iCnttp)%FragCoef,nFragDens,nFragEner,Label='FragCoef')

call Read_v(lUnit,dbsc(iCnttp)%FragCoef,1,nFragDens*nFragEner,1,iErr)
if (iPrint >= 99) call RecPrt('Fragment MO coefficients',' ',dbsc(iCnttp)%FragCoef,nFragDens,nFragEner)
if (ierr /= 0) then
  write(6,*) 'ERROR: number of coefficients is not correct'
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! keyword: MULLIKEN
!                                                                      *
!***********************************************************************
!                                                                      *
! Mulliken charges

if (iPrint >= 99) write(6,*) 'Reading MULLIKEN'
Line = Get_Ln(lUnit)
if (index(Line,'MULLIKEN') == 0) then
  write(6,*) 'ERROR: Keyword MULLIKEN expected, offending line:'
  write(6,*) Line
  call Quit_OnUserError()
end if
call Read_v(lUnit,dbsc(iCnttp)%FragCoor,5,5*nFragCoor,5,ierr)
if (iPrint >= 99) call RecPrt('Fragment Mulliken charges',' ',dbsc(iCnttp)%FragCoor,5,nFragCoor)
if (ierr /= 0) then
  write(6,*) 'ERROR: number of Mulliken charges is not correct'
  call Quit_OnUserError()
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine GetFragment
