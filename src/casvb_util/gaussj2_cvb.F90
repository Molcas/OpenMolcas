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

subroutine gaussj2_cvb(a,lrow,lcol,ibook,irows,ijs,oijs,n)

implicit real*8(a-h,o-z)
dimension a(n,n)
dimension lrow(n), lcol(n), ibook(n)
dimension irows(n), ijs(2,n*n), oijs(n*n)
save zero, one, thresh
data zero/0.d0/,one/1.d0/,thresh/1.d-10/

! initialize imx & jmx to suppress compiler warnings ...
imx = 0
jmx = 0

nij = n*n
do i=1,n
  irows(i) = i
  ibook(i) = 0
end do
do imain=1,n
  amx = zero
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
          write(6,*) ' Singular matrix in GAUSSJ !'
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
    goto 3000
  end if
  ijs(1,nij) = irows(jmx)
  ijs(2,nij) = irows(jmx)
  oijs(nij) = a(jmx,jmx)
  nij = nij-1
  oneovamx = one/a(jmx,jmx)
  a(jmx,jmx) = one
  do ii=1,n
    a(jmx,ii) = oneovamx*a(jmx,ii)
  end do
  do ii2=1,n
    if (ii2 /= jmx) then
      dum = a(ii2,jmx)
      a(ii2,jmx) = zero
      do ii=1,n
        a(ii2,ii) = a(ii2,ii)-dum*a(jmx,ii)
      end do
      ijs(1,nij) = irows(ii2)
      ijs(2,nij) = irows(jmx)
      oijs(nij) = dum
      nij = nij-1
    end if
  end do
end do
do ii=n,1,-1
  if (lrow(ii) /= lcol(ii)) then
    do ii2=1,n
      dum = a(ii2,lrow(ii))
      a(ii2,lrow(ii)) = a(ii2,lcol(ii))
      a(ii2,lcol(ii)) = dum
    end do
  end if
end do
return
3000 continue
do i=1,n
  do j=1,ihad
    if (lcol(j) == i) goto 3100
  end do
  ijs(1,nij) = irows(i)
  ijs(2,nij) = irows(i)
  oijs(nij) = zero
  nij = nij-1
  do j=1,n
    if (j /= i) then
      ijs(1,nij) = irows(j)
      ijs(2,nij) = irows(i)
      oijs(nij) = a(j,i)
      nij = nij-1
    end if
  end do
3100 continue
end do

return

end subroutine gaussj2_cvb
