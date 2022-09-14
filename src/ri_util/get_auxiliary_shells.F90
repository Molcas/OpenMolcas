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

subroutine Get_Auxiliary_Shells(iSO,nSO,jOff,iSO2Shl,nSO2Shl,iPair,nPair)

use Index_Functions, only: iTri
use Definitions, only: iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nSO, iSO(2,nSO), jOff, nSO2Shl, iSO2Shl(nSO2Shl), nPair
integer(kind=iwp), intent(_OUT_) :: iPair(nPair)
integer(kind=iwp) :: i, k, kl, kSh, l, lSh

!write(u6,*) 'iSO'
!write(u6,*) '==='
!do i=1,nSO
!  write(u6,*) iSO(1,i),iSO(2,i)
!end do

!write(u6,*) 'iSO2Shl'
!write(u6,*) '======='
!do i=1,nSO2Shl
!  write(u6,*) i,iSO2Shl(i)
!end do
do i=1,nSO
  k = iSO(1,i)+jOff
  l = iSO(2,i)+jOff
  kSh = iSO2Shl(k)
  lSh = iSO2Shl(l)
  !write(u6,*) 'k,kSh=',k,kSh
  !write(u6,*) 'l,lSh=',k,lSh
  kl = iTri(kSh,lSh)
  iPair(kl) = 1
end do
!write(u6,*) 'iPairs'
!write(u6,*) '======'
!do i=1,nPair
!  write(u6,*) iPair(i)
!end do

return

end subroutine Get_Auxiliary_Shells
