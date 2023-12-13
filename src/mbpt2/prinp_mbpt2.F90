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
! Copyright (C) 1992, Markus P. Fuelscher                              *
!***********************************************************************

subroutine PrInp_MBPT2(Eocc,Eext,iTst)
!***********************************************************************
!                                                                      *
!     Print the program banner, date and time of execution             *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1992                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use MBPT2_Global, only: iDel, iFro, iPL, nBas, nDel1, nDel2, nDsto, nFro1, nFro2, nTit, Title
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: Eocc(*), Eext(*)
integer(kind=iwp), intent(in) :: iTst
integer(kind=iwp) :: i, ii, iOrb, iSym, k, left, lLine, lPaper, nLine
logical(kind=iwp) :: lFro, lDel
character(len=3) :: lIrrep(8)
character(len=8) :: Fmt1, Fmt2
character(len=120) :: Line, StLine
character(len=102) :: XLine
#include "corbinf.fh"

!----------------------------------------------------------------------*
!     Start and define the paper width                                 *
!----------------------------------------------------------------------*
lPaper = 132
!----------------------------------------------------------------------*
!     Initialize header lines                                          *
!----------------------------------------------------------------------*
lLine = len(Line)
do i=1,lLine
  StLine(i:i) = '*'
end do
lPaper = 132
left = (lPaper-lLine)/2
write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
!----------------------------------------------------------------------*
!     Print the project title                                          *
!----------------------------------------------------------------------*
if (nTit > 0) then
  write(u6,*)
  nLine = nTit+5
  do i=1,nLine
    Line = ''
    if ((i == 1) .or. (i == nLine)) Line = StLine
    if (i == 3) Line = 'Project:'
    if ((i >= 4) .and. (i <= nLine-2)) then
      Line = Title(i-3)
    end if
    call Center_Text(Line)
    write(u6,Fmt1) '*'//Line//'*'
  end do
  write(u6,*)
end if
!----------------------------------------------------------------------*
!     Print the coordinates of the system                              *
!----------------------------------------------------------------------*
if (iPL >= 2) call PrCoor()
!----------------------------------------------------------------------*
!     Print contents of the runfile RUNFILE                            *
!----------------------------------------------------------------------*

call Get_cArray('Irreps',lIrrep,24)
do iSym=1,nSym
  lIrrep(iSym) = adjustr(lIrrep(iSym))
end do

if (iPL >= 2) then
  write(u6,*)
  write(u6,Fmt2//'A)') 'Contents of RUNFILE file:'
  write(u6,Fmt2//'A)') '-------------------------'
  write(u6,*)
  write(u6,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Frozen occupied orbitals',(nFro(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Active occupied orbitals',(nOcc(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Active external orbitals',(nExt(iSym),iSym=1,nSym)
  write(u6,Fmt2//'A,T47,8I4)') 'Deleted external orbitals',(nDel(iSym),iSym=1,nSym)
end if
!----------------------------------------------------------------------*
!     ordering and orbital energies of the frozen occupied orbitals    *
!----------------------------------------------------------------------*
lFro = .false.
do iSym=1,nSym
  if (nFro1(iSym)+nFro2(iSym) /= 0) lFro = .true.
end do
if (iPL >= 2) then
  if (lFro) then
    write(u6,*)
    write(u6,*)
    write(u6,Fmt2//'A,T47)') 'Reference numbers of frozen occupied orbitals according to the original input sequence'
    do iSym=1,nSym
      if (nFro1(iSym) /= 0) then
        write(u6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym,(iOrb,iOrb=1,nFro1(iSym))
      end if
      if (nFro2(iSym) /= 0) then
        write(u6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym,(iFro(iSym,iOrb),iOrb=1,nFro2(iSym))
      end if
    end do
  end if
  write(u6,*)
  write(u6,*)
  write(u6,Fmt2//'A,T47)') 'Energies of the active occupied orbitals'
  ii = 0
  do iSym=1,nSym
    if (nOcc(iSym) /= 0) then
      write(u6,*)
      write(u6,Fmt2//'A,I2,(T40,5F14.6))') 'symmetry species',iSym,(Eocc(ii+k),k=1,nOcc(iSym))
      ii = ii+nOcc(iSym)
    end if
  end do
end if
!----------------------------------------------------------------------*
!     ordering and orbital energies of the deleted external orbitals   *
!----------------------------------------------------------------------*
lDel = .false.
do iSym=1,nSym
  if (nDel(iSym)+nDel1(iSym)+nDel2(iSym) /= 0) lDel = .true.
end do
if (iPL >= 2) then
  if (lDel) then
    write(u6,*)
    write(u6,*)
    write(u6,Fmt2//'A,T47)') 'Reference numbers of deleted external orbitals according to the original input sequence'
    do iSym=1,nSym
      if ((nDel1(iSym) /= 0) .or. (nDsto(iSym) /= 0)) then
        write(u6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym, &
                                                     (nBas(iSym)-nFro(iSym)-nOcc(iSym)-iOrb+1,iOrb=nDsto(iSym)+nDel1(iSym),1,-1)
      end if
      if (nDel2(iSym) /= 0) then
        write(u6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym,(iDel(iSym,iOrb),iOrb=1,nDel2(iSym))
      end if
    end do
  end if
  write(u6,*)
  write(u6,*)
  write(u6,Fmt2//'A,T47)') 'Energies of the active external orbitals'
  ii = 0
  do iSym=1,nSym
    if (nExt(iSym) /= 0) then
      write(u6,*)
      write(u6,Fmt2//'A,I2,(T40,5F14.6))') 'symmetry species',iSym,(Eext(ii+k),k=1,nExt(iSym))
      ii = ii+nExt(iSym)
    end if
  end do
end if
!----------------------------------------------------------------------*
!     print header for final results                                   *
!----------------------------------------------------------------------*
if ((iTst == 0) .and. (iPL >= 2)) then
  write(u6,*)
  write(u6,*)
  nLine = 3
  do i=1,nLine
    XLine = ''
    if ((i == 1) .or. (i == nLine)) XLine = trim(StLine)
    if (i == 2) XLine = 'Results'
    call Center_Text(XLine)
    write(u6,Fmt1) '*'//XLine//'*'
  end do
  write(u6,*)
end if
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*
return

end subroutine PrInp_MBPT2
