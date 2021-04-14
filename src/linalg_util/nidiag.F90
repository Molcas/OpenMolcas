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
! Copyright (C) 2005, Per-Olof Widmark                                 *
!***********************************************************************

subroutine NIdiag(H,U,n,nv,iOpt)
!***********************************************************************
!                                                                      *
! This routine is a wrapper that calls appropriate routines to         *
! perform diagonalization of of symmetric matrices stored in lower     *
! triangular form.                                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund university, Sweden                                     *
! Written November 2005                                                *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp, r8

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
! n    - Dimension of matrix                                           *
! nv   - Length of eigenvectors nv>=n                                  *
! H    - Matrix to be diagonalized                                     *
! U    - Eigenvectors                                                  *
! iOpt - Option flag, for future improvements.                         *
!----------------------------------------------------------------------*
integer(kind=iwp), intent(in) :: n, nv, iOpt
real(kind=wp), intent(inout) :: H(*), U(nv,n)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp) :: i, ierr
real(kind=wp) :: Tmp
real(kind=r8), external :: OrbPhase
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
if (n == 0) return
call Givens(H,U,n,nv)
call QLdiag(H,U,n,nv,ierr)
if (ierr == 1) then
  !write(u6,*) 'NIdiag: backup call to Jacob!'
  call Jacob(H,U,n,nv)
end if
do i=1,n
  Tmp = OrbPhase(U(1,i),nV)
end do
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
return
#ifdef _WARNING_WORKAROUND_
if (.false.) then
  call Unused_integer(iOpt)
  call Unused_real(Tmp)
end if
#endif

end subroutine NIdiag
