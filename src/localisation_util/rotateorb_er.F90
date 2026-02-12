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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBasis, nOrb2Loc
real(kind=wp), intent(in) :: R(nOrb2Loc,nOrb2Loc)
real(kind=wp), intent(inout) :: CMO(nBasis,nOrb2Loc)
logical(kind=iwp), intent(in) :: Debug
integer(kind=iwp) :: irc
real(kind=wp) :: ThrU
real(kind=wp), allocatable :: CMOscr(:,:), U(:,:)
character(len=*), parameter :: SecNam = 'RotateOrb_ER'

if ((nOrb2Loc < 1) .or. (nBasis < 1)) return

! Allocate transformation matrix U.
! ---------------------------------

call mma_allocate(U,nOrb2Loc,nOrb2Loc,label='Umat')

! Compute U.
! ----------

call GetU_ER(U,R,nOrb2Loc)

! Debug: check that U is unitary.
! -------------------------------

if (Debug) then
  ThrU = 1.0e-10_wp
  irc = -1
  call Chk_Unitary(irc,U,nOrb2Loc,ThrU)
  if (irc /= 0) then
    call SysAbendMsg(SecNam,'U matrix is not unitary!',' ')
  end if
end if

! Update C.
! ---------

call mma_allocate(CMOscr,nBasis,nOrb2Loc,label='CMOscr')
CMOscr(:,:) = CMO
call DGEMM_('N','N',nBasis,nOrb2Loc,nOrb2Loc,One,CMOscr,nBasis,U,nOrb2Loc,Zero,CMO,nBasis)
call mma_deallocate(CMOscr)

! De-allocate U.
! --------------

call mma_deallocate(U)

end subroutine RotateOrb_ER
