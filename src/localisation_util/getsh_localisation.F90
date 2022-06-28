!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine GetSh_Localisation(X,nBas,m,XSh,nShell,iSO2Sh,iOpt,Norm)

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nBas, m, nShell, iSO2Sh(*), iOpt
real(kind=wp), intent(in) :: X(nBas,m)
real(kind=wp), intent(_OUT_) :: XSh(nShell,*)
character(len=3), intent(in) :: Norm  ! 'MAX' or 'FRO'
integer(kind=iwp) :: i, iShell, j, jShell
character(len=3) :: myNorm

if ((nBas < 1) .or. (nShell < 1)) return

myNorm = Norm
call UpCase(myNorm)

if (iOpt == 1) then
  XSh(:,1:m) = Zero
  if (myNorm == 'MAX') then
    do j=1,m
      do i=1,nBas
        iShell = iSO2Sh(i)
        XSh(iShell,j) = max(abs(X(i,j)),XSh(iShell,j))
      end do
    end do
  else if (myNorm == 'FRO') then
    do j=1,m
      do i=1,nBas
        iShell = iSO2Sh(i)
        XSh(iShell,j) = XSh(iShell,j)+X(i,j)**2
      end do
      do iShell=1,nShell
        XSh(iShell,j) = sqrt(XSh(iShell,j))
      end do
    end do
  end if
else
  if (m /= nBas) call SysAbendMsg('GetSh_Localisation','Fatal error','m != nBas')
  XSh(:,1:nShell) = Zero
  if (myNorm == 'MAX') then
    do j=1,nBas
      jShell = iSO2Sh(j)
      do i=1,nBas
        iShell = iSO2Sh(i)
        XSh(iShell,jShell) = max(abs(X(i,j)),XSh(iShell,jShell))
      end do
    end do
  else if (myNorm == 'FRO') then
    do j=1,nBas
      jShell = iSO2Sh(j)
      do i=1,nBas
        iShell = iSO2Sh(i)
        XSh(iShell,jShell) = XSh(iShell,jShell)+X(i,j)**2
      end do
      Xsh(:,jShell) = sqrt(Xsh(:,jShell))
    end do
  end if
end if

end subroutine GetSh_Localisation
