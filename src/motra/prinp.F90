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

use motra_global, only: BsLbl, CutThrs, Debug, Header, iAutoCut, iPrint, iRFpert, nBas, nDel, nFro, nOrb, nSym, nTit, Occ, Title, &
                        VecTit
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: CMO(*)
integer(kind=iwp) :: i, left, lLine, lPaper, nLine
real(kind=wp) :: Ene, ThrEne, ThrOcc
character(len=120) :: Line, StLine
character(len=8) :: Frmt
logical(kind=iwp) :: PrOcc, PrEne

!----------------------------------------------------------------------*
! Start and define the paper width                                     *
!----------------------------------------------------------------------*
lPaper = 132
!----------------------------------------------------------------------*
! Initialize blank and header lines                                    *
!----------------------------------------------------------------------*
lLine = len(Line)
do i=1,lLine
  StLine(i:i) = '*'
end do
left = (lPaper-lLine)/2
write(Frmt,'(A,I3.3,A)') '(',left,'X,A)'
!----------------------------------------------------------------------*
! Print the project title                                              *
!----------------------------------------------------------------------*
if (nTit > 0) then
  write(u6,*)
  nLine = nTit+5
  do i=1,nLine
    Line = ''
    if (i == 1 .or. i == nLine) Line = StLine
    if (i == 3) Line = 'Project:'
    if (i >= 4 .and. i <= nLine-2) then
      Line = Title(i-3)
    end if
    call Center_Text(Line)
    write(u6,Frmt) '*'//Line//'*'
  end do
  write(u6,*)
end if
!----------------------------------------------------------------------*
! Print the integral file header                                       *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'Header of the integral files:'
write(Line,'(A)') Header(1:72)
write(u6,'(6X,A)') trim(Line)
write(Line,'(A)') Header(73:144)
write(u6,'(6X,A)') trim(Line)
write(u6,*)
!----------------------------------------------------------------------*
! Print the header of the source file of MO coefficients               *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'Header of MO coefficients source file:'
write(u6,'(6X,A)') VecTit
write(u6,*)
!----------------------------------------------------------------------*
! Print coordinates of the system                                      *
!----------------------------------------------------------------------*
call PrCoor()
!----------------------------------------------------------------------*
! Print the orbital specifications                                     *
!----------------------------------------------------------------------*
write(u6,*)
write(u6,'(6X,A)') 'Orbital specifications:'
write(u6,'(6X,A)') '-----------------------'
write(u6,*)
write(u6,'(6X,A,T35,8i4)') 'Symmetry species:',(i,i=1,nSym)
write(u6,'(6X,A,T35,8i4)') 'Number of basis functions:',(nBas(i),i=1,nSym)
write(u6,'(6X,A,T35,8i4)') 'Frozen orbitals:',(nFro(i),i=1,nSym)
write(u6,'(6X,A,T35,8i4)') 'Deleted orbitals:',(nDel(i),i=1,nSym)
write(u6,'(6X,A,T35,8i4)') 'Number of orbitals used:',(nOrb(i),i=1,nSym)
if (iAutoCut == 1) then
  write(u6,'(6X,A)') 'Automatic orbital deletion is turned on'
  if (iAutoCut == 1) then
    write(u6,'(6X,A,T35,8F10.8)') 'Cutting thresholds:',(CutThrs(i),i=1,nSym)
  end if
end if
if (iRFpert /= 0) then
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') 'Reaction field specifications:'
  write(u6,'(6X,A)') '------------------------------'
  write(u6,*)
  write(u6,'(6X,A)') 'The Reaction field is added as a perturbation and has been determined in a previos calculation'
  write(u6,*)
end if

!----------------------------------------------------------------------*
! Print MO coefficients                                                *
!----------------------------------------------------------------------*
if ((iPrint >= 2) .or. (Debug == 1)) then
  PrEne = .false.
  PrOcc = .false.
  if (iAutoCut == 1) PrOcc = .true.
  ThrOcc = Zero
  ThrEne = Zero
  Ene = Zero
  Line = 'Input orbitals after orthogonalization'
  call PRIMO(Line,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nBas,BsLbl,[Ene],Occ,Cmo,-1)
end if
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*

return

end subroutine PrInp
