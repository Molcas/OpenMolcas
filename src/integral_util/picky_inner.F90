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
integer(kind=iwp) :: i_n, iDCR, iIn, ijCmp, im, iOut, jIn, jm, jn, jOut, nl, ml

if ((nSt == 1) .and. (nEnd == n) .and. (mSt == 1) .and. (mEnd == m)) then
  ! Copy the whole block
  Dout(:,:) = Din(:,:)
else
  nl = nEnd-nSt+1
  ml = mEnd-mSt+1
  iIn = (n*m+1)*nijCmp
  iOut = (nl*ml+1)*nijCmp
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
          Dout(ij(jn,jm,ijCmp,nl,ml),iDCR) = Din(ij(i_n,im,ijCmp,n,m),iDCR)
        end do
      end do
      ! Move the largest density matrix element for
      ! this angular combination
      jIn = ij(n,m,ijCmp,n,m)+1
      jOut = ij(nl,ml,ijCmp,nl,ml)+1
      Dout(jOut,iDCR) = Din(jIn,iDCR)
    end do
    ! Move the largest density matrix element for
    ! each pair plus the overall largest element
    Dout(iOut+1:iOut+nijPrm+1,iDCR) = Din(iIn+1:iIn+nijPrm+1,iDCR)
  end do
end if

return

contains

pure function ij(i,j,k,nn,mm)

  integer(kind=iwp) :: ij
  integer(kind=iwp), intent(in) :: i, j, k, nn, mm

  ij = (k-1)*(nn*mm+1)+(j-1)*nn+i

end function ij

end subroutine Picky_inner
