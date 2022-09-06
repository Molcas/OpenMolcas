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

integer iSO(2,nSO), iSO2Shl(nSO2Shl), iPair(nPair)

!write(6,*) 'iSO'
!write(6,*) '==='
!do i=1,nSO
!  write(6,*) iSO(1,i),iSO(2,i)
!end do

!write(6,*) 'iSO2Shl'
!write(6,*) '======='
!do i=1,nSO2Shl
!  write(6,*) i,iSO2Shl(i)
!end do
do i=1,nSO
  k = iSO(1,i)+jOff
  l = iSO(2,i)+jOff
  kSh = iSO2Shl(k)
  lSh = iSO2Shl(l)
  !write(6,*) 'k,kSh=',k,kSh
  !write(6,*) 'l,lSh=',k,lSh
  kl = max(kSh,lSh)*(max(kSh,lSh)-1)/2+min(kSh,lSh)
  iPair(kl) = 1
end do
!write(6,*) 'iPairs'
!write(6,*) '======'
!do i=1,nPair
!  write(6,*) iPair(i)
!end do

return

end subroutine Get_Auxiliary_Shells
