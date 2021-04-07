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

subroutine RdInp_Motra()
!***********************************************************************
!                                                                      *
! Subroutine RdInp                                                     *
!                                                                      *
! Purpose: Read input from inputstream                                 *
!                                                                      *
!**** M.P. Fuelscher, University of Lund, Sweden, 1991 *****************

implicit real*8(A-H,O-Z)
parameter(nCmd=16)
parameter(lCmd=4)
character*4 CmdTab(nCmd)
#include "motra_global.fh"
#include "trafo_motra.fh"
#include "files_motra.fh"
character*180 Line, Blank
integer nDel2(8)
character*(180) Get_Ln
external Get_Ln
data CmdTab/'TITL','FROZ','DELE','PRIN','MOLO','LUMO','JOBI','ONEL','FILE','AUTO','EXTR','RFPE','CTON','DIAG','HDF5','END '/
#include "cho_minp.fh"
#include "chotraw.fh"

iCTonly = 0
iDoInt = 0
ihdf5 = 0
tv2disk = 'PQK'
!----------------------------------------------------------------------*
! Initialize some arrays                                           *
!----------------------------------------------------------------------*
do iSym=1,mxSym
  nDel(iSym) = 0
  nOrb(iSym) = 0
  CutThrs(iSym) = 0.0d0
end do
call Get_iArray('Non valence orbitals',nFro,nSym)
do i=1,72
  Blank(i:i) = ' '
end do
!----------------------------------------------------------------------*
! Locate "start of input"                                          *
!----------------------------------------------------------------------*
LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
call RdNLst(LuSpool,'MOTRA')
!----------------------------------------------------------------------*
! Read the input stream line by line and identify key command      *
!----------------------------------------------------------------------*
100 read(LuSpool,'(A)') Line
if (Line(1:72) == Blank(1:72) .or. Line(1:1) == '*') goto 100
call UpCase(Line)
110 jCmd = 0
do iCmd=1,nCmd
  if (Line(1:lCmd) == CmdTab(iCmd)(1:lCmd)) jCmd = iCmd
end do
if (jCmd == 0) then
  write(6,*) 'RdInp: Unknown command at line: ',trim(Line)
  call Abend()
end if
!----------------------------------------------------------------------*
! Branch to the processing of the command sections                 *
!----------------------------------------------------------------------*
goto(1010,1020,1030,1040,1050,1060,1070,1080,1090,1100,1110,1120,1130,1140,1150,2000),jCmd
!---  Process the "TITLe" command -------------------------------------*
1010 nTit = 0
15 read(LuSpool,'(A)') Line
if (Line(1:72) == Blank(1:72) .or. Line(1:1) == '*') goto 15
call UpCase(Line)
if (nTit > 0) then
  do iCmd=1,nCmd
    if (Line(1:lCmd) == CmdTab(iCmd)(1:lCmd)) goto 110
  end do
end if
nTit = nTit+1
if (nTit <= mxTit) Title(nTit) = trim(Line)
goto 15
!---  Process the "FROZen orbitals" command ---------------------------*
1020 continue

if (iPrint >= 0) then
  write(6,*)
  write(6,'(6X,A)') '*** WARNING: Default frozen orbitals is overwritten by user input.'
  write(6,'(6X,A,8I4)') '*** Default values:',(nFro(iSym),iSym=1,nSym)
end if

25 read(LuSpool,'(A)') Line
if (Line(1:72) == Blank(1:72) .or. Line(1:1) == '*') goto 25
read(Line,*,Err=994) (nFro(iSym),iSym=1,nSym)
goto 100
!---  Process the "DELEted orbitals" command --------------------------*
1030 continue
35 read(LuSpool,'(A)') Line
if (Line(1:72) == Blank(1:72) .or. Line(1:1) == '*') goto 35
read(Line,*,Err=994) (nDel(iSym),iSym=1,nSym)
goto 100
!---  Process the "PRINt level" command -------------------------------*
1040 continue
45 read(LuSpool,'(A)') Line
if (Line(1:72) == Blank(1:72) .or. Line(1:1) == '*') goto 45
read(Line,*,Err=994) iPrint
goto 100
!---  Process the "MOLOrb" command ------------------------------------*
1050 continue
iVecTyp = 1
goto 100
!---  Process the "LUMOrb" command ------------------------------------*
1060 continue
iVecTyp = 2
goto 100
!---  Process the "JOBIph" command ------------------------------------*
1070 continue
iVecTyp = 3
goto 100
!---  Process the "ONEL only" command ---------------------------------*
1080 continue
iOneOnly = 1
goto 100
!---  Process the "FILEORB" command------------------------------------*
1090 continue
iVecTyp = 2
Line = Get_Ln(LuSpool)
write(6,*) ' RdInp_Motra before calling fileorb.'
write(6,*) ' Line:'//line(1:60)
write(6,*) '   Calling fileorb now...'
call fileorb(Line,FnInpOrb)
write(6,*) '   Back from fileorb.'
write(6,*) ' Line:'//line(1:60)
write(6,*) FnInpOrb
goto 100
!---  Process the "AUTO delete" command--------------------------------*
1100 continue
iAutoCut = 1
105 read(LuSpool,'(A)') Line
if ((Line(1:72) == Blank(1:72)) .or. (Line(1:1) == '*')) goto 105
read(Line,*,Err=994) (CutThrs(iSym),iSym=1,nSym)
goto 100
!---  Process the "EXTRact" command------------------------------------*
1110 write(6,*) 'The EXTRACT option is redundant and is ignored!'
goto 100
!---  Process the "RFperturbation" command ----------------------------*
1120 continue
iRFpert = 1
goto 100
!---  Process the "CTonly" to perform exclusively CD vectors transform-*
1130 continue
Line = Get_Ln(LuSpool)
call UpCase(Line)
call LeftAd(Line)
tv2disk = Line(1:3)
if ((tv2disk /= 'PQK') .and. (tv2disk /= 'KPQ')) tv2disk = 'PQK' !def
iCTonly = 1
goto 100
!---  Process the "DIAGonal ERI evaluation" command -------------------*
1140 continue
iDoInt = 1
goto 100
!---  Process the "HDF5 output file" command --------------------------*
1150 continue
ihdf5 = 1
goto 100
!---  Process the "END of input" command ------------------------------*
2000 continue

! New rules for title lines...warning needed?
if (nTit > mxTit) then
  write(6,*) ' NOTE: New input specifications says TITLE keyword'
  write(6,*) ' must be followed by exactly one title line.'
  nTit = mxTit
end if

! Check for deleted orbitals

call Get_iArray('nDel',nDel2,nSym)
do iSym=1,nSym
  if (nDel2(iSym) > nDel(iSym)) then
    write(6,*)
    write(6,*) 'Orbitals deleted at an earlier stage.'
    write(6,*) 'iSym=',iSym
    write(6,*) 'Input value  :',nDel(iSym)
    write(6,*) 'Earlier value:',nDel2(iSym)
    write(6,*)
    write(6,*) 'Input value is now updated to earlier value!'
    write(6,*)
    nDel(iSym) = nDel2(iSym)
  end if
end do
!----------------------------------------------------------------------*
! Normal termination                                                   *
!----------------------------------------------------------------------*
nOrbt = 0
nOrbtt = 0
do iSym=1,nSym
  nOrb(iSym) = nBas(iSym)-nFro(iSym)-nDel(iSym)
  nOrbt = nOrbt+nOrb(iSym)
  nOrbtt = nOrbtt+nOrb(iSym)*(nOrb(iSym)+1)/2
end do
call Put_iArray('nFro',nFro,nSym)
close(LuSpool)
return
!----------------------------------------------------------------------*
! Error Exit                                                           *
!----------------------------------------------------------------------*
994 write(6,*) 'RdInp: error readin input file!'
write(6,*) 'Command=',CmdTab(jCmd)
call Abend()

end subroutine RdInp_Motra
