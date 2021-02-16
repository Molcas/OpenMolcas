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

subroutine InitGO()

use GuessOrb_Global, only: GapThr, iPrFmt, Label, MxAtom, MxSym, AtName, nBas, nDel, nNuc, nOcc, nSym, nVir, PrintEor, PrintMOs, &
                           PrintPop, PrThr, SThr, TThr, xCharge
use Constants, only: Five
use Definitions, only: wp, iwp, u6

implicit none
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
logical(kind=iwp) :: Debug, Trace, Found
integer(kind=iwp) :: nBasTot, iSym, iBas, iPrt, i, lenName
!----------------------------------------------------------------------*
! External entry points                                                *
!----------------------------------------------------------------------*
integer(kind=iwp), external :: iPrintLevel
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
Debug = .false.
Trace = .false.
lenName = len(AtName)
if (Trace) write(u6,*) '>>> Entering initgo'
!----------------------------------------------------------------------*
! Set default for MO printing.                                         *
!----------------------------------------------------------------------*
iPrt = iPrintLevel(-1)
if (iPrt >= 4) then
  PrintMOs = .true.
  PrintEor = .true.
  PrintPop = .true.
  iPrFmt = 3
  PrThr = huge(PrThr)
else if (iPrt >= 3) then
  PrintMOs = .false.
  PrintEor = .false.
  PrintPop = .false.
  iPrFmt = 1
  PrThr = Five
else
  PrintMOs = .false.
  PrintEor = .false.
  PrintPop = .false.
  PrThr = Five
end if
!----------------------------------------------------------------------*
! Set other defaults.                                                  *
!----------------------------------------------------------------------*
call Qpg_dScalar('S delete thr',Found)
if (Found) then
  call Get_dScalar('S delete thr',SThr)
else
  SThr = 1.0e-9_wp
  call Put_dScalar('S delete thr',SThr)
end if
call Qpg_dScalar('T delete thr',Found)
if (Found) then
  call Get_dScalar('T delete thr',TThr)
else
  TThr = huge(TThr)
  call Put_dScalar('T delete thr',TThr)
end if
GapThr = 0.01_wp
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
  write(u6,'(a,8i5)') 'initgo: nSym',nSym
  write(u6,'(a,8i5)') 'initgo: nBas',nBas
  write(u6,'(a,8i5)') 'initgo: nOcc',nOcc
  write(u6,'(a,8i5)') 'initgo: nVir',nVir
  write(u6,'(a,8i5)') 'initgo: nBasTot',nBasTot
end if
call get_iscalar('Unique Atoms',nNuc)
if (nNuc > MxAtom) then
  call SysAbendMsg('initgo','Fatal:','Too many atoms, increase MxAtom')
end if
call get_carray('Unique Atom Names',AtName,lenName*nNuc)
call get_carray('Unique Basis Names',Label,len(Label)*nBasTot)
call get_darray('Nuclear Charge',xCharge,nNuc)
if (Debug) then
  write(u6,'(a,8i5)') 'initgo: nNuc',nNuc
  write(u6,'(a,8i5)') 'initgo: nBasTot',nBasTot
  write(u6,'(a,8f12.6)') 'initgo: Charge',(xCharge(i),i=1,nNuc)
  write(u6,'(a,8a4)') 'initgo: Name ',(AtName(i),i=1,nNuc)
  write(u6,'(a)') 'initgo: Basis functions'
  do iBas=1,nBasTot
    write(u6,'(2a)') Label(iBas)(1:lenName),Label(iBas)(lenName+1:)
  end do
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
if (Trace) write(u6,*) '<<< Exiting initgo'

return

end subroutine InitGO
