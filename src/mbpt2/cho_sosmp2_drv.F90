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
! Copyright (C) 2007, Francesco Aquilante                              *
!***********************************************************************

subroutine Cho_SOSmp2_Drv(irc,EMP2,CMO,EOcc,EVir)
! Francesco Aquilante, May 2007.
!
! Purpose: driver for computing the Scaled Opposite-Spin (SOS)
!          MP2 energy correction EMP2
!          using Cholesky (or RI) representation for the
!          two-electron integrals.
!          Input must have been processed and MO coefficients
!          and orbital energies must be passed as arguments.
!
! Notes:
!
!   - all MO Cholesky vector files generated here are deleted before
!     exit, except for error terminations (i.e. no cleanup actions
!     are taken!)

#include "implicit.fh"
dimension CMO(*), EOcc(*), EVir(*)
character*3 ThisNm
character*14 SecNam
parameter(SecNam='Cho_SOSmp2_Drv',ThisNm='Drv')
parameter(Chk_Mem_ChoMP2=0.123456789d0,Tol=1.0D-15)
parameter(iFmt=0)
logical Delete, Delete_def
parameter(Delete_def=.true.)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"

#if defined (_DEBUGPRINT_)
Verbose = .true.
#endif
if (Verbose) then
  call CWTime(CPUTot1,WallTot1)
end if

! Initializations.
! ----------------

irc = 0

EMP2 = 0.0d0

if (Verbose) then
  call CWTime(CPUIni1,WallIni1)
end if

l_Dum = 1
call GetMem('Dummy','Allo','Real',ip_Dum,l_Dum)
Work(ip_Dum) = Chk_Mem_ChoMP2

FracMem = 0.0d0 ! no buffer allocated
call Cho_X_Init(irc,FracMem)
if (irc /= 0) then
  write(6,*) SecNam,': Cho_X_Init returned ',irc
  call ChoMP2_Quit(SecNam,'Cholesky initialization error',' ')
end if

call Cho_SOSmp2_Setup(irc)
if (irc /= 0) then
  write(6,*) SecNam,': Cho_SOSmp2_Setup returned ',irc
  Go To 1  ! exit
end if

if (Verbose) then
  call Cho_SOSmp2_Setup_Prt(irc)
  if (irc /= 0) then
    write(6,*) SecNam,': Cho_SOSmp2_Setup_Prt returned ',irc
    Go To 1  ! exit
  end if
  call CWTime(CPUIni2,WallIni2)
  call Cho_PrtTim('Cholesky SOS-MP2 initialization',CPUIni2,CPUIni1,WallIni2,WallIni1,iFmt)
end if

! Transform Cholesky vectors directly from reduced set to MO
! representation. Result vectors are stored on disk.
! Compute also the (ai|ai)^2 diagonal here.
! ----------------------------------------------------------

if (Verbose) then
  call CWTime(CPUTra1,WallTra1)
end if
lDiag = nT1am(1)
do iSym=2,nSym
  lDiag = lDiag+nT1am(iSym)
end do
call GetMem('Diag','Allo','Real',ipDiag,lDiag)

call ChoMP2_TraDrv(irc,CMO,Work(ipDiag),.true.)
if (irc /= 0) then
  write(6,*) SecNam,': ChoMP2_TraDrv returned ',irc
  Go To 1  ! exit
end if

! Squaring each diagonal element
! ------------------------------
do ia=0,lDiag-1
  Work(ipDiag+ia) = Work(ipDiag+ia)**2
end do
if (set_cd_thr) ThrMP2 = ddot_(lDiag,[1.0d0],0,Work(ipDiag),1)/(5.0d0*lDiag)

if (Verbose) then
  call CWTime(CPUTra2,WallTra2)
  call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,WallTra2,WallTra1,iFmt)
end if

! Finalize Cholesky info (to release memory).
! Retain essential info: LuPri, nSym, and NumCho(*).
! --------------------------------------------------

nSym_Sav = nSym
call iCopy(nSym,NumCho,1,nMP2Vec,1)

call Cho_X_Final(irc)
if (irc /= 0) then
  write(6,*) SecNam,': Cho_X_Final returned ',irc
  Go To 1 ! exit
end if

LuPri = 6
nSym = nSym_Sav
call iCopy(nSym,nMP2Vec,1,NumCho,1)

! Decompose M(ai,bj) = (ai|bj)^2 .
! Set number of vectors to be used in energy calculation.
! -------------------------------------------------------

if (Verbose) then
  call CWTime(CPUDec1,WallDec1)
end if
Delete = Delete_def ! delete transf. vector files after dec.
call Cho_SOSmp2_DecDrv(irc,Delete,Work(ipDiag))
if (irc /= 0) then
  write(6,*) SecNam,': Cho_SOSmp2_DecDrv returned ',irc
  call ChoMP2_Quit(SecNam,'SOS-MP2 decomposition failed!',' ')
end if
if (Verbose) then
  call CWTime(CPUDec2,WallDec2)
  call Cho_PrtTim('Cholesky SOS-MP2 decomposition',CPUDec2,CPUDec1,WallDec2,WallDec1,iFmt)
end if
call GetMem('Diag','Free','Real',ipDiag,lDiag)

! Compute SOS-MP2 energy correction.
! ----------------------------------

if (Verbose) then
  call CWTime(CPUEnr1,WallEnr1)
end if
Delete = Delete_def
call Cho_SOSmp2_Energy(irc,EMP2,EOcc,EVir,Delete)
if (irc /= 0) then
  write(6,*) SecNam,': Cho_SOSmp2_Energy returned ',irc
  Go To 1 ! exit
end if
if (Verbose) then
  call CWTime(CPUEnr2,WallEnr2)
  call Cho_PrtTim('Cholesky SOS-MP2 energy',CPUEnr2,CPUEnr1,WallEnr2,WallEnr1,iFmt)
end if

! Exit.
! -----

1 continue
Diff = abs(Work(ip_Dum)-Chk_Mem_ChoMP2)
if (Diff > Tol) then
  write(6,*) SecNam,': Memory Boundary Error!'
  if (irc == 0) irc = -9999
end if
if (Verbose) then
  call CWTime(CPUTot2,WallTot2)
  call Cho_PrtTim('Cholesky SOS-MP2',CPUTot2,CPUTot1,WallTot2,WallTot1,iFmt)
end if
call GetMem('Flush','Flush','Real',ip_Dum,l_Dum)
call GetMem('Dummy','Free','Real',ip_Dum,l_Dum)

return

end subroutine Cho_SOSmp2_Drv
