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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine gaussj2_cvb(a,lrow,irows,ijs,oijs,n)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(inout) :: a(n,n)
integer(kind=iwp), intent(out) :: lrow(n), irows(n), ijs(2,n*n)
real(kind=wp), intent(out) :: oijs(n*n)
integer(kind=iwp) :: i, idum, ihad, ii, ii2, imain, imx, j, jmx, nij
real(kind=wp) :: amx, dum, oneovamx
logical(kind=iwp) :: done
integer(kind=iwp), allocatable :: lcol(:), ibook(:)
real(kind=wp), parameter :: thresh = 1.0e-10_wp

call mma_allocate(lcol,n,label='lcol')
call mma_allocate(ibook,n,label='ibook')

! initialize imx & jmx to suppress compiler warnings ...
imx = 0
jmx = 0

nij = n*n
do i=1,n
  irows(i) = i
end do
ibook(:) = 0
done = .false.
do imain=1,n
  amx = Zero
  do i=1,n
    if (ibook(i) /= 1) then
      do j=1,n
        if (ibook(j) == 0) then
          if (abs(a(i,j)) >= amx) then
            amx = abs(a(i,j))
            imx = i
            jmx = j
          end if
        else if (ibook(j) > 1) then
          write(u6,*) ' Singular matrix in GAUSSJ !'
          call abend_cvb()
        end if
      end do
    end if
  end do
  ibook(jmx) = ibook(jmx)+1
  if (imx /= jmx) then
    do ii=1,n
      dum = a(imx,ii)
      a(imx,ii) = a(jmx,ii)
      a(jmx,ii) = dum
    end do
    idum = irows(imx)
    irows(imx) = irows(jmx)
    irows(jmx) = idum
  end if
  lrow(imain) = imx
  lcol(imain) = jmx
  if (abs(a(jmx,jmx)) < thresh) then
    ihad = imain-1
    done = .true.
    exit
  end if
  ijs(1,nij) = irows(jmx)
  ijs(2,nij) = irows(jmx)
  oijs(nij) = a(jmx,jmx)
  nij = nij-1
  oneovamx = One/a(jmx,jmx)
  a(jmx,jmx) = One
  a(jmx,:) = oneovamx*a(jmx,:)
  do ii2=1,n
    if (ii2 /= jmx) then
      dum = a(ii2,jmx)
      a(ii2,jmx) = Zero
      a(ii2,:) = a(ii2,:)-dum*a(jmx,:)
      ijs(1,nij) = irows(ii2)
      ijs(2,nij) = irows(jmx)
      oijs(nij) = dum
      nij = nij-1
    end if
  end do
end do
if (done) then
  outer: do i=1,n
    do j=1,ihad
      if (lcol(j) == i) cycle outer
    end do
    ijs(1,nij) = irows(i)
    ijs(2,nij) = irows(i)
    oijs(nij) = Zero
    nij = nij-1
    do j=1,n
      if (j /= i) then
        ijs(1,nij) = irows(j)
        ijs(2,nij) = irows(i)
        oijs(nij) = a(j,i)
        nij = nij-1
      end if
    end do
  end do outer
else
  do ii=n,1,-1
    if (lrow(ii) /= lcol(ii)) then
      do ii2=1,n
        dum = a(ii2,lrow(ii))
        a(ii2,lrow(ii)) = a(ii2,lcol(ii))
        a(ii2,lcol(ii)) = dum
      end do
    end if
  end do
end if

call mma_deallocate(lcol)
call mma_deallocate(ibook)

return

end subroutine gaussj2_cvb
