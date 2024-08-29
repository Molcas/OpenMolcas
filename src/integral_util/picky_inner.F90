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
! Copyright (C) 1993, Roland Lindh                                     *
!***********************************************************************

subroutine Picky_inner(Din,n,m,nijPrm,nijCmp,nDCR,nSt,nEnd,mSt,mEnd,Dout)
!***********************************************************************
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, m, nijPrm, nijCmp, nDCR, nSt, nEnd, mSt, mEnd
real(kind=wp), intent(in) :: Din(((n*m+1)*nijCmp)+nijPrm+1,nDCR)
real(kind=wp), intent(out) :: Dout((((nEnd-nSt+1)*(mEnd-mSt+1)+1)*nijCmp)+nijPrm+1,nDCR)
integer(kind=iwp) :: i_n, iDCR, iIn, ijCmp, im, iOut, jIn, jm, jn, jOut
! Statement functions
integer(kind=iwp) :: i, j, k, ij1, ij2
ij1(i,j,k) = (k-1)*(n*m+1)+(j-1)*n+i
ij2(i,j,k) = (k-1)*((nEnd-nSt+1)*(mEnd-mSt+1)+1)+(j-1)*(nEnd-nSt+1)+i

if ((nSt == 1) .and. (nEnd == n) .and. (mSt == 1) .and. (mEnd == m)) then
  ! Copy the whole block
  call dcopy_((((n*m+1)*nijCmp)+nijPrm+1)*nDCR,Din,1,Dout,1)
else
  iIn = (n*m+1)*nijCmp
  iOut = (((nEnd-nSt+1)*(mEnd-mSt+1)+1)*nijCmp)
  ! Loop over desymmetrized density blocks
  do iDCR=1,nDCR
    ! Loop over angular combinations
    do ijCmp=1,nijCmp
      jm = 0
      ! Loop over subset of contracted basis
      do im=mSt,mEnd
        jm = jm+1
        jn = 0
        ! Loop over subset of contracted basis
        do i_n=nSt,nEnd
          jn = jn+1
          Dout(ij2(jn,jm,ijCmp),iDCR) = Din(ij1(i_n,im,ijCmp),iDCR)
        end do
      end do
      ! Move the largest density matrix element for
      ! this angular combination
      jIn = ij1(n,m,ijCmp)+1
      jOut = ij2(nEnd-nSt+1,mEnd-mSt+1,ijCmp)+1
      Dout(jOut,iDCR) = Din(jIn,iDCR)
    end do
    ! Move the largest density matrix element for
    ! each pair plus the overall largest element
    call dcopy_(nijPrm+1,Din(iIn+1,iDCR),1,Dout(iOut+1,iDCR),1)
  end do
end if

return

end subroutine Picky_inner
