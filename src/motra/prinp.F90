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
! Copyright (C) 1991, Markus P. Fuelscher                              *
!***********************************************************************

subroutine PrInp(CMO)
!***********************************************************************
!                                                                      *
! Purpose:                                                             *
! Echo all input                                                       *
!                                                                      *
!**** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************

implicit real*8(A-H,O-Z)
#include "motra_global.fh"
real*8 CMO(*)
character*120 Line, BlLine, StLine
character*8 Fmt
logical PrOcc, PrEne

!----------------------------------------------------------------------*
! Start and define the paper width                                     *
!----------------------------------------------------------------------*
lPaper = 132
!----------------------------------------------------------------------*
! Initialize blank and header lines                                    *
!----------------------------------------------------------------------*
lLine = len(Line)
do i=1,lLine
  BlLine(i:i) = ' '
  StLine(i:i) = '*'
end do
left = (lPaper-lLine)/2
write(Fmt,'(A,I3.3,A)') '(',left,'X,A)'
!----------------------------------------------------------------------*
! Print the project title                                              *
!----------------------------------------------------------------------*
if (nTit > 0) then
  write(6,*)
  nLine = nTit+5
  do i=1,nLine
    Line = BlLine
    if (i == 1 .or. i == nLine) Line = StLine
    if (i == 3) Line = 'Project:'
    if (i >= 4 .and. i <= nLine-2) then
      Line = Title(i-3)
    end if
    call Center(Line)
    write(6,Fmt) '*'//Line//'*'
  end do
  write(6,*)
end if
!----------------------------------------------------------------------*
! Print the integral file header                                       *
!----------------------------------------------------------------------*
write(6,*)
write(6,'(6X,A)') 'Header of the integral files:'
write(Line,'(A)') Header(1:72)
write(6,'(6X,A)') Line(:mylen(Line))
write(Line,'(A)') Header(73:144)
write(6,'(6X,A)') Line(:mylen(Line))
write(6,*)
!----------------------------------------------------------------------*
! Print the header of the source file of MO coefficients               *
!----------------------------------------------------------------------*
write(6,*)
write(6,'(6X,A)') 'Header of MO coefficients source file:'
write(6,'(6X,A)') VecTit
write(6,*)
!----------------------------------------------------------------------*
! Print coordinates of the system                                      *
!----------------------------------------------------------------------*
call PrCoor()
!----------------------------------------------------------------------*
! Print the orbital specifications                                     *
!----------------------------------------------------------------------*
write(6,*)
write(6,'(6X,A)') 'Orbital specifications:'
write(6,'(6X,A)') '-----------------------'
write(6,*)
write(6,'(6X,A,T35,8i4)') 'Symmetry species:',(i,i=1,nSym)
write(6,'(6X,A,T35,8i4)') 'Number of basis functions:',(nBas(i),i=1,nSym)
write(6,'(6X,A,T35,8i4)') 'Frozen orbitals:',(nFro(i),i=1,nSym)
write(6,'(6X,A,T35,8i4)') 'Deleted orbitals:',(nDel(i),i=1,nSym)
write(6,'(6X,A,T35,8i4)') 'Number of orbitals used:',(nOrb(i),i=1,nSym)
if (iAutoCut == 1) then
  write(6,'(6X,A)') 'Automatic orbital deletion is turned on'
  if (iAutoCut == 1) then
    write(6,'(6X,A,T35,8F10.8)') 'Cutting thresholds:',(CutThrs(i),i=1,nSym)
  end if
end if
if (iRFpert /= 0) then
  write(6,*)
  write(6,*)
  write(6,'(6X,A)') 'Reaction field specifications:'
  write(6,'(6X,A)') '------------------------------'
  write(6,*)
  write(6,'(6X,A)') 'The Reaction field is added as a perturbation and has been determined in a previos calculation'
  write(6,*)
end if

!----------------------------------------------------------------------*
! Print MO coefficients                                                *
!----------------------------------------------------------------------*
if ((iPrint >= 2) .or. (Debug == 1)) then
  PrEne = .false.
  PrOcc = .false.
  if (iAutoCut == 1) PrOcc = .true.
  ThrOcc = 0.0d0
  ThrEne = 0.0d0
  Ene = 0.0d0
  Line = 'Input orbitals after orthogonalization'
  call PRIMO(Line,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nBas,BsLbl,[Ene],Occ,Cmo,-1)
end if
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*

return

end subroutine PrInp
