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

use motra_global, only: CutThrs, FnInpOrb, iAutoCut, iCTonly, iDoInt, ihdf5, iOneOnly, iortho, iPrint, iRFpert, iVecTyp, nBas, &
                        nDel, nFro, nOrb, nOrbt, nOrbtt, nSym, nTit, Title
use Cholesky, only: tv2disk
use Constants, only: Zero
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iCmd, istatus, iSym, jCmd, LuSpool, mxTit, nDel2(nSym)
character(len=180) :: Line
logical(kind=iwp) :: Skip
integer(kind=iwp), parameter :: nCmd = 17, lCmd = 4
character(len=lCmd), parameter :: CmdTab(nCmd) = ['TITL','FROZ','DELE','PRIN','MOLO','LUMO','JOBI','ONEL','FILE','AUTO', &
                                                  'EXTR','RFPE','CTON','DIAG','HDF5','NOOR','END ']
character(len=180), external :: Get_Ln

iortho = 0
iCTonly = 0
iDoInt = 0
ihdf5 = 0
tv2disk = 'PQK'
!----------------------------------------------------------------------*
! Initialize some arrays                                               *
!----------------------------------------------------------------------*
nDel(:) = 0
nOrb(:) = 0
CutThrs(:) = Zero
mxTit = size(Title)
call Get_iArray('Non valence orbitals',nFro,nSym)
!----------------------------------------------------------------------*
! Locate "start of input"                                              *
!----------------------------------------------------------------------*
LuSpool = 17
call SpoolInp(LuSpool)
rewind(LuSpool)
call RdNLst(LuSpool,'MOTRA')
!----------------------------------------------------------------------*
! Read the input stream line by line and identify key command          *
!----------------------------------------------------------------------*
Skip = .false.
input: do
  if (.not. Skip) then
    read(LuSpool,'(A)') Line
    if ((Line(1:72) == ' ') .or. (Line(1:1) == '*')) cycle
    call UpCase(Line)
  end if
  Skip = .false.
  jCmd = 0
  do iCmd=1,nCmd
    if (Line(1:lCmd) == CmdTab(iCmd)(1:lCmd)) jCmd = iCmd
  end do
  if (jCmd == 0) then
    write(u6,*) 'RdInp: Unknown command at line: ',trim(Line)
    call Abend()
  end if
  !----------------------------------------------------------------------*
  ! Branch to the processing of the command sections                     *
  !----------------------------------------------------------------------*
  select case (jCmd)
    case (1)
      !---  Process the "TITLe" command -------------------------------*
      nTit = 0
      do
        do
          read(LuSpool,'(A)') Line
          if ((Line(1:72) /= ' ') .and. (Line(1:1) /= '*')) exit
        end do
        call UpCase(Line)
        if (nTit > 0) then
          do iCmd=1,nCmd
            if (Line(1:lCmd) == CmdTab(iCmd)(1:lCmd)) then
              Skip = .true.
              cycle input
            end if
          end do
        end if
        nTit = nTit+1
        if (nTit <= mxTit) Title(nTit) = trim(Line)
      end do
    case (2)
      !---  Process the "FROZen orbitals" command ---------------------*

      if (iPrint >= 0) then
        write(u6,*)
        write(u6,'(6X,A)') '*** WARNING: Default frozen orbitals is overwritten by user input.'
        write(u6,'(6X,A,8I4)') '*** Default values:',(nFro(iSym),iSym=1,nSym)
      end if

      do
        read(LuSpool,'(A)') Line
        if ((Line(1:72) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) nFro(1:nSym)
      if (istatus /= 0) call Error()
    case (3)
      !---  Process the "DELEted orbitals" command --------------------*
      do
        read(LuSpool,'(A)') Line
        if ((Line(1:72) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) nDel(1:nSym)
      if (istatus /= 0) call Error()
    case (4)
      !---  Process the "PRINt level" command -------------------------*
      do
        read(LuSpool,'(A)') Line
        if ((Line(1:72) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) iPrint
      if (istatus /= 0) call Error()
    case (5)
      !---  Process the "MOLOrb" command ------------------------------*
      iVecTyp = 1
    case (6)
      !---  Process the "LUMOrb" command ------------------------------*
      iVecTyp = 2
    case (7)
      !---  Process the "JOBIph" command ------------------------------*
      iVecTyp = 3
    case (8)
      !---  Process the "ONEL only" command ---------------------------*
      iOneOnly = 1
    case (9)
      !---  Process the "FILEORB" command------------------------------*
      iVecTyp = 2
      Line = Get_Ln(LuSpool)
      write(u6,*) ' RdInp_Motra before calling fileorb.'
      write(u6,*) ' Line:'//line(1:60)
      write(u6,*) '   Calling fileorb now...'
      call fileorb(Line,FnInpOrb)
      write(u6,*) '   Back from fileorb.'
      write(u6,*) ' Line:'//line(1:60)
      write(u6,*) FnInpOrb
    case (10)
      !---  Process the "AUTO delete" command--------------------------*
      iAutoCut = 1
      do
        read(LuSpool,'(A)') Line
        if ((Line(1:72) /= ' ') .and. (Line(1:1) /= '*')) exit
      end do
      read(Line,*,iostat=istatus) CutThrs(1:nSym)
      if (istatus /= 0) call Error()
    case (11)
      !---  Process the "EXTRact" command------------------------------*
      write(u6,*) 'The EXTRACT option is redundant and is ignored!'
    case (12)
      !---  Process the "RFperturbation" command ----------------------*
      iRFpert = 1
    case (13)
      !---  Process the "CTonly" to perform exclusively CD vectors transform-*
      Line = Get_Ln(LuSpool)
      call UpCase(Line)
      Line = adjustl(Line)
      tv2disk = Line(1:3)
      if ((tv2disk /= 'PQK') .and. (tv2disk /= 'KPQ')) tv2disk = 'PQK' !def
      iCTonly = 1
    case (14)
      !---  Process the "DIAGonal ERI evaluation" command -------------*
      iDoInt = 1
    case (15)
      !---  Process the "HDF5 output file" command --------------------*
      ihdf5 = 1
    case (16)
      !---  Process the "NOORthogonalization" command -----------------*
      iortho = 1
    case (17)
      !---  Process the "END of input" command ------------------------*
      exit
    case default
  end select
end do input

! New rules for title lines... warning needed?
if (nTit > mxTit) then
  write(u6,*) ' NOTE: New input specifications says TITLE keyword'
  write(u6,*) ' must be followed by exactly one title line.'
  nTit = mxTit
end if

! Check for deleted orbitals

call Get_iArray('nDel',nDel2,nSym)
do iSym=1,nSym
  if (nDel2(iSym) > nDel(iSym)) then
    write(u6,*)
    write(u6,*) 'Orbitals deleted at an earlier stage.'
    write(u6,*) 'iSym=',iSym
    write(u6,*) 'Input value  :',nDel(iSym)
    write(u6,*) 'Earlier value:',nDel2(iSym)
    write(u6,*)
    write(u6,*) 'Input value is now updated to earlier value!'
    write(u6,*)
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
! Bug in original code?? (tps/cdg 20210430)
call Put_iArray('nDel',nDel,nSym)
close(LuSpool)

return

contains

subroutine Error()
  write(u6,*) 'RdInp: error readin input file!'
  write(u6,*) 'Command=',CmdTab(jCmd)
  call Abend()
end subroutine Error

end subroutine RdInp_Motra
