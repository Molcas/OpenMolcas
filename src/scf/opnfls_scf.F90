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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!***********************************************************************

subroutine OpnFls_SCF()
!***********************************************************************
!                                                                      *
!     purpose: Open files needed by SCF                                *
!                                                                      *
!***********************************************************************

use InfSCF, only: DoCholesky, DSCF
use SCFFiles, only: FnDel, FnDGd, FnDst, FnGrd, FnOrd, FnOSt, FnTSt, Fnx, Fny, LuDel, LuDGd, LuDst, LuGrd, LuOrd, LuOSt, LuTSt, &
                    Lux, Luy
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iOpt, iRC
logical(kind=iwp) :: test

!---  open two-electron integral file ---------------------------------*
call f_Inquire(FnOrd,test)
call DecideOnDirect(.true.,test,DSCF,DoCholesky)
if ((.not. DSCF) .and. (.not. DoCholesky)) then
  !InVec = 0
  iRc = -1
  iOpt = 0
  call OpnOrd(iRC,iOpt,FnOrd,LuOrd)
  if (iRc /= 0) then
    write(u6,*) 'OpnFls: Error opening ORDINT'
    call Abend()
  end if
end if

!---  open DNSMAT, dVxcdR, TWOHAM and GRADIENT ------------------------*
call DAName(LuDSt,FnDSt)
call DAName(LuOSt,FnOSt)
call DAName(LuTSt,FnTSt)
call DAName(LuGrd,FnGrd)

!---  open 2nd order update files      --------------------------------*
call DAName(LuDGd,FnDGd)
call DAName(Lux,Fnx)
call DAName(LuDel,FnDel)
call DAName(Luy,Fny)

end subroutine OpnFls_SCF
