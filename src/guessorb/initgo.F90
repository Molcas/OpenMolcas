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
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine initialize guessorb.                                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************

subroutine InitGO(StandAlone)

use GuessOrb_Global, only: GapThr, iPrFmt, Label, LenIn, LenIn1, LenIn8, MxAtom,  MxSym, Name, nBas, nDel, nNuc, nOcc, nSym, nVir, &
                           PrintEor, PrintMOs, PrintPop, PrThr, SThr, TThr, xCharge

implicit none
logical StandAlone
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
logical Debug
logical Trace
logical Found
integer nBasTot
integer iSym
integer iBas
integer iPrt
integer i
!----------------------------------------------------------------------*
! External entry points                                                *
!----------------------------------------------------------------------*
integer iPrintLevel
external iPrintLevel
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
Debug = .false.
Trace = .false.
if (Trace) write(6,*) '>>> Entering initgo'
!----------------------------------------------------------------------*
! Set default for MO printing.                                         *
!----------------------------------------------------------------------*
iPrt = iPrintLevel(-1)
if (iPrt >= 4) then
  PrintMOs = .true.
  PrintEor = .true.
  PrintPop = .true.
  iPrFmt = 3
  PrThr = 1.0d6
else if (iPrt >= 3) then
  PrintMOs = .false.
  PrintEor = .false.
  PrintPop = .false.
  iPrFmt = 1
  PrThr = 5.0d0
else
  PrintMOs = .false.
  PrintEor = .false.
  PrintPop = .false.
  PrThr = 5.0d0
end if
!----------------------------------------------------------------------*
! Set other defaults.                                                  *
!----------------------------------------------------------------------*
call Qpg_dScalar('S delete thr',Found)
if (Found) then
  call Get_dScalar('S delete thr',SThr)
else
  SThr = 1.0d-9
  call Put_dScalar('S delete thr',SThr)
end if
call Qpg_dScalar('T delete thr',Found)
if (Found) then
  call Get_dScalar('T delete thr',TThr)
else
  TThr = 1.0d+6
  call Put_dScalar('T delete thr',TThr)
end if
GapThr = 0.01d0
!----------------------------------------------------------------------*
! Get basic data from runfile.                                         *
!----------------------------------------------------------------------*
call get_iscalar('nSym',nSym)
call get_iarray('nBas',nBas,nSym)
do iSym=1,MxSym
  nOcc(iSym) = 0
  nVir(iSym) = 0
  nDel(iSym) = 0
end do
nBasTot = 0
do iSym=1,nSym
  nBasTot = nBasTot+nBas(iSym)
end do
if (Debug) then
  write(6,'(a,8i5)') 'initgo: nSym',nSym
  write(6,'(a,8i5)') 'initgo: nBas',nBas
  write(6,'(a,8i5)') 'initgo: nOcc',nOcc
  write(6,'(a,8i5)') 'initgo: nVir',nVir
  write(6,'(a,8i5)') 'initgo: nBasTot',nBasTot
end if
call get_iscalar('Unique Atoms',nNuc)
if (nNuc > MxAtom) then
  call SysAbendMsg('initgo','Fatal:','Too many atoms, increase MxAtom')
end if
call get_carray('Unique Atom Names',Name,LENIN*nNuc)
call get_carray('Unique Basis Names',Label,(LENIN8)*nBasTot)
call get_darray('Nuclear Charge',xCharge,nNuc)
if (Debug) then
  write(6,'(a,8i5)') 'initgo: nNuc',nNuc
  write(6,'(a,8i5)') 'initgo: nBasTot',nBasTot
  write(6,'(a,8f12.6)') 'initgo: Charge',(xCharge(i),i=1,nNuc)
  write(6,'(a,8a4)') 'initgo: Name ',(Name(i),i=1,nNuc)
  write(6,'(a)') 'initgo: Basis functions'
  do iBas=1,nBasTot
    write(6,'(2a)') Label(iBas)(1:LENIN),Label(iBas)(LENIN1:LENIN8)
  end do
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
if (Trace) write(6,*) '<<< Exiting initgo'

return

! Avoid unused argument warnings
if (.false.) call Unused_logical(StandAlone)

end subroutine InitGO
