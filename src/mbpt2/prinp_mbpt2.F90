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

subroutine PrInp_MBPT2(Eocc,Eext)
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

implicit real*8(A-H,O-Z)
! declare local variables...
dimension Eocc(*), Eext(*)
character*8 Fmt1, Fmt2
character*120 Line, BlLine, StLine
character*102 XLine
character*3 lIrrep(8)
logical lFro, lDel
#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "cdtfaux.fh"
#include "print_mbpt2.fh"

!----------------------------------------------------------------------*
!     Start and define the paper width                                 *
!----------------------------------------------------------------------*
lPaper = 132
!----------------------------------------------------------------------*
!     Initialize blank and header lines                                *
!----------------------------------------------------------------------*
lLine = len(Line)
do i=1,lLine
  BlLine(i:i) = ' '
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
  write(6,*)
  nLine = nTit+5
  do i=1,nLine
    Line = BlLine
    if ((i == 1) .or. (i == nLine)) Line = StLine
    if (i == 3) Line = 'Project:'
    if ((i >= 4) .and. (i <= nLine-2)) then
      Line = Title(i-3)
    end if
    call Center(Line)
    write(6,Fmt1) '*'//Line//'*'
  end do
  write(6,*)
end if
!----------------------------------------------------------------------*
!     Print the coordinates of the system                              *
!----------------------------------------------------------------------*
if (iPL >= 2) call PrCoor()
!----------------------------------------------------------------------*
!     Print contents of the runfile RUNFILE                            *
!----------------------------------------------------------------------*
!
call Get_cArray('Irreps',lIrrep,24)
do iSym=1,nSym
  call RightAd(lIrrep(iSym))
end do
!
if (iPL >= 2) then
  write(6,*)
  write(6,Fmt2//'A)') 'Contents of RUNFILE file:'
  write(6,Fmt2//'A)') '-------------------------'
  write(6,*)
  write(6,Fmt2//'A,T47,8I4)') 'Symmetry species',(iSym,iSym=1,nSym)
  write(6,Fmt2//'A,T47,8(1X,A))') '                ',(lIrrep(iSym),iSym=1,nSym)
  write(6,Fmt2//'A,T47,8I4)') 'Number of basis functions',(nBas(iSym),iSym=1,nSym)
  write(6,Fmt2//'A,T47,8I4)') 'Frozen occupied orbitals',(nFro(iSym),iSym=1,nSym)
  write(6,Fmt2//'A,T47,8I4)') 'Active occupied orbitals',(nOcc(iSym),iSym=1,nSym)
  write(6,Fmt2//'A,T47,8I4)') 'Active external orbitals',(nExt(iSym),iSym=1,nSym)
  write(6,Fmt2//'A,T47,8I4)') 'Deleted external orbitals',(nDel(iSym),iSym=1,nSym)
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
    write(6,*)
    write(6,*)
    write(6,Fmt2//'A,T47)') 'Reference numbers of frozen occupied orbitals according to the original input sequence'
    do iSym=1,nSym
      if (nFro1(iSym) /= 0) then
        write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym,(iOrb,iOrb=1,nFro1(iSym))
      end if
      if (nFro2(iSym) /= 0) then
        write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym,(iFro(iSym,iOrb),iOrb=1,nFro2(iSym))
      end if
    end do
  end if
  write(6,*)
  write(6,*)
  write(6,Fmt2//'A,T47)') 'Energies of the active occupied orbitals'
  ii = 0
  do iSym=1,nSym
    if (nOcc(iSym) /= 0) then
      write(6,*)
      write(6,Fmt2//'A,I2,(T40,5F14.6))') 'symmetry species',iSym,(Eocc(ii+k),k=1,nOcc(iSym))
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
    write(6,*)
    write(6,*)
    write(6,Fmt2//'A,T47)') 'Reference numbers of deleted external orbitals according to the original input sequence'
    do iSym=1,nSym
      if ((nDel1(iSym) /= 0) .or. (nDsto(iSym) /= 0)) then
        write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym, &
                                                    (nBas(iSym)-nFro(iSym)-nOcc(iSym)-iOrb+1,iOrb=nDsto(iSym)+nDel1(iSym),1,-1)
      end if
      if (nDel2(iSym) /= 0) then
        write(6,Fmt2//'A,I2,T47,20I3,(/T49,20I3))') 'symmetry species',iSym,(iDel(iSym,iOrb),iOrb=1,nDel2(iSym))
      end if
    end do
  end if
  write(6,*)
  write(6,*)
  write(6,Fmt2//'A,T47)') 'Energies of the active external orbitals'
  ii = 0
  do iSym=1,nSym
    if (nExt(iSym) /= 0) then
      write(6,*)
      write(6,Fmt2//'A,I2,(T40,5F14.6))') 'symmetry species',iSym,(Eext(ii+k),k=1,nExt(iSym))
      ii = ii+nExt(iSym)
    end if
  end do
end if
!----------------------------------------------------------------------*
!     # of orbitals in 1st index to be skipped (restart)               *
!----------------------------------------------------------------------*
if ((iRest /= 0) .and. (iPL >= 2)) then
  write(6,*)
  write(6,Fmt2//'A)') 'restart information...'
  write(6,Fmt2//'A,T47,I8)') 'max # integral passes',nPass
  write(6,Fmt2//'A,T47,8I4)') 'symmetry species',(iSym,iSym=1,nSym)
  write(6,Fmt2//'A,T47,8I4)') '# of orbitals in 1st index skipped',(MOSkip(iSym-1),iSym=1,nSym)
end if
!----------------------------------------------------------------------*
!     print header for final results                                   *
!----------------------------------------------------------------------*
if ((iTst == 0) .and. (iPL >= 2)) then
  write(6,*)
  write(6,*)
  nLine = 3
  do i=1,nLine
    XLine = trim(BlLine)
    if ((i == 1) .or. (i == nLine)) XLine = trim(StLine)
    if (i == 2) XLine = 'Results'
    call Center(XLine)
    write(6,Fmt1) '*'//XLine//'*'
  end do
  write(6,*)
end if
!----------------------------------------------------------------------*
!     Normal termination                                               *
!----------------------------------------------------------------------*
return

end subroutine PrInp_MBPT2
