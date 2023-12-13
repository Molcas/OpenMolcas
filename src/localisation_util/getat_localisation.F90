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

subroutine GetAt_Localisation(X,nBas,m,XAt,nAtoms,iOpt,nBas_per_Atom,nBas_Start,Norm)

use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nBas, m, nAtoms, iOpt, nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
real(kind=wp), intent(in) :: X(nBas,m)
real(kind=wp), intent(_OUT_) :: XAt(nAtoms,*)
character(len=3), intent(in) :: Norm  ! 'MAX' or 'FRO'
integer(kind=iwp) :: i, i1, i2, iAt, j, j1, j2, jAt
character(len=3) :: myNorm

if ((nBas < 1) .or. (nAtoms < 1)) return

myNorm = Norm
call UpCase(myNorm)

if (iOpt == 1) then
  XAt(:,1:m) = Zero
  if (myNorm == 'MAX') then
    do j=1,m
      do iAt=1,nAtoms
        i1 = nBas_Start(iAt)
        i2 = i1+nBas_per_Atom(iAt)-1
        do i=i1,i2
          XAt(iAt,j) = max(abs(X(i,j)),XAt(iAt,j))
        end do
      end do
    end do
  else if (myNorm == 'FRO') then
    do j=1,m
      do iAt=1,nAtoms
        i1 = nBas_Start(iAt)
        i2 = i1+nBas_per_Atom(iAt)-1
        do i=i1,i2
          XAt(iAt,j) = XAt(iAt,j)+X(i,j)**2
        end do
        XAt(iAt,j) = sqrt(XAt(iAt,j))
      end do
    end do
  end if
else
  if (m /= nBas) then
    call SysAbendMsg('GetAt_Localisation','Fatal error','m != nBas')
  end if
  XAt(:,1:nAtoms) = Zero
  if (myNorm == 'MAX') then
    do jAt=1,nAtoms
      j1 = nBas_Start(jAt)
      j2 = j1+nBas_per_Atom(jAt)-1
      do j=j1,j2
        do iAt=1,nAtoms
          i1 = nBas_Start(iAt)
          i2 = i1+nBas_per_Atom(iAt)-1
          do i=i1,i2
            XAt(iAt,jAt) = max(abs(X(i,j)),XAt(iAt,jAt))
          end do
        end do
      end do
    end do
  else if (myNorm == 'FRO') then
    do jAt=1,nAtoms
      j1 = nBas_Start(jAt)
      j2 = j1+nBas_per_Atom(jAt)-1
      do j=j1,j2
        do iAt=1,nAtoms
          i1 = nBas_Start(iAt)
          i2 = i1+nBas_per_Atom(iAt)-1
          do i=i1,i2
            XAt(iAt,jAt) = XAt(iAt,jAt)+X(i,j)**2
          end do
        end do
      end do
      XAt(:,jAt) = sqrt(XAt(:,jAt))
    end do
  end if
end if

end subroutine GetAt_Localisation
