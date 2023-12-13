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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Tra2B(iSym,jSym,iBas,jBas,iAsh,jAsh,iOrb,jOrb,ikl,nkl,CMO_ip,CMO_jp,CMO_ia,CMO_ja,VXIJ,X2,X3_ip,X3_jp,PUVX_jp,PUVX_ip)
!***********************************************************************
!                                                                      *
!     run the second half of the AO --> MO transformation with         *
!     one active and one general index, i.e.                           *
!     (vx!ij) --> (vx!up)                                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSym, jSym, iBas, jBas, iAsh, jAsh, iOrb, jOrb, ikl, nkl
real(kind=wp), intent(in) :: CMO_ip(iBas*iOrb), CMO_jp(jBas*jOrb), CMO_ia(iBas*iAsh), CMO_ja(jBas*jAsh), VXIJ(*)
real(kind=wp), intent(out) :: X2(*), X3_ip(iAsh,jOrb), X3_jp(jAsh,iOrb)
real(kind=wp), intent(inout) :: PUVX_jp(iOrb,jAsh,nkl), PUVX_ip(jOrb,iAsh,nkl)
integer(kind=iwp) :: i, j

if (jAsh*iOrb /= 0) then
  call DGEMM_('T','N',iBas,jAsh,jBas,One,VXIJ,jBas,CMO_ja,jBas,Zero,X2,iBas)

  call DGEMM_('T','N',jAsh,iOrb,iBas,One,X2,iBas,CMO_ip,iBas,Zero,X3_jp,jAsh)

  do i=1,jAsh
    do j=1,iOrb
      PUVX_jp(j,i,ikl) = X3_jp(i,j)
    end do
  end do
end if

if ((iSym /= jSym) .and. (iAsh*jOrb /= 0)) then
  call DGEMM_('N','N',jBas,iAsh,iBas,One,VXIJ,jBas,CMO_ia,iBas,Zero,X2,jBas)

  call DGEMM_('T','N',iAsh,jOrb,jBas,One,X2,jBas,CMO_jp,jBas,Zero,X3_ip,iAsh)

  do i=1,iAsh
    do j=1,jOrb
      PUVX_ip(j,i,ikl) = X3_ip(i,j)
    end do
  end do
end if

return

end subroutine Tra2B
