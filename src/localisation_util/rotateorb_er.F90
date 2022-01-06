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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RotateOrb_ER(R,CMO,nBasis,nOrb2Loc,Debug)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: rotate ER orbitals,
!          CMO -> CMO * U
!          U = R*[R^T*R]^(-1/2)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc
real(kind=wp), intent(in) :: R(nOrb2Loc,nOrb2Loc)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
logical(kind=iwp), intent(in) :: Debug
#include "WrkSpc.fh"
integer(kind=iwp) :: ipCMO, ipU, irc, lCMO, lU
real(kind=wp) :: ThrU
character(len=*), parameter :: SecNam = 'RotateOrb_ER'

if ((nOrb2Loc < 1) .or. (nBasis < 1)) return

! Allocate transformation matrix U.
! ---------------------------------

lU = nOrb2Loc**2
call GetMem('Umat','Allo','Real',ipU,lU)

! Compute U.
! ----------

call GetU_ER(Work(ipU),R,nOrb2Loc)

! Debug: check that U is unitary.
! -------------------------------

if (Debug) then
  ThrU = 1.0e-10_wp
  irc = -1
  call Chk_Unitary(irc,Work(ipU),nOrb2Loc,ThrU)
  if (irc /= 0) then
    call SysAbendMsg(SecNam,'U matrix is not unitary!',' ')
  end if
end if

! Update C.
! ---------

lCMO = nBasis*nOrb2Loc
call GetMem('CMOscr','Allo','Real',ipCMO,lCMO)
call dCopy_(lCMO,CMO,1,Work(ipCMO),1)
call DGEMM_('N','N',nBasis,nOrb2Loc,nOrb2Loc,One,Work(ipCMO),nBasis,Work(ipU),nOrb2Loc,Zero,CMO,nBasis)
call GetMem('CMOscr','Free','Real',ipCMO,lCMO)

! De-allocate U.
! --------------

call GetMem('Umat','Free','Real',ipU,lU)

end subroutine RotateOrb_ER
