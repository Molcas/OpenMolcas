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

!#define _DEBUGPRINT_
subroutine InitGO()

use GuessOrb_Global, only: AtName, GapThr, iPrFmt, Label, nBas, nDel, nNuc, nOcc, nSym, nVir, PrintEor, PrintMOs, PrintPop, PrThr, &
                           SThr, TThr
#ifdef _OLD_
use GuessOrb_Global, only: xCharge
#endif
use Molcas, only: MxAtom
use Constants, only: Five
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: iPrt, lenName, nBasTot
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i, iBas
#endif
logical(kind=iwp) :: Found
integer(kind=iwp), external :: iPrintLevel

!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
lenName = len(AtName)
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
nOcc(:) = 0
nVir(:) = 0
nDel(:) = 0
nBasTot = sum(nBas(1:nSym))
#ifdef _DEBUGPRINT_
write(u6,'(a,8i5)') 'initgo: nSym',nSym
write(u6,'(a,8i5)') 'initgo: nBas',nBas
write(u6,'(a,8i5)') 'initgo: nOcc',nOcc
write(u6,'(a,8i5)') 'initgo: nVir',nVir
write(u6,'(a,8i5)') 'initgo: nBasTot',nBasTot
#endif
call get_iscalar('Unique Atoms',nNuc)
if (nNuc > MxAtom) call SysAbendMsg('initgo','Fatal:','Too many atoms, increase MxAtom')
call get_carray('Unique Atom Names',AtName,lenName*nNuc)
call get_carray('Unique Basis Names',Label,len(Label)*nBasTot)
#ifdef _DEBUGPRINT_
write(u6,'(a,8i5)') 'initgo: nNuc',nNuc
write(u6,'(a,8i5)') 'initgo: nBasTot',nBasTot
#ifdef _OLD_
call get_darray('Nuclear Charge',xCharge,nNuc)
write(u6,'(a,8f12.6)') 'initgo: Charge',(xCharge(i),i=1,nNuc)
#endif
write(u6,'(a,8a4)') 'initgo: Name ',(AtName(i),i=1,nNuc)
write(u6,'(a)') 'initgo: Basis functions'
do iBas=1,nBasTot
  write(u6,'(2a)') Label(iBas)(1:lenName),Label(iBas)(lenName+1:)
end do
#endif
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*

end subroutine InitGO
