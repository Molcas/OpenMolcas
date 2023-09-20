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

subroutine ddsol72_cvb(hp,eigval,eigvec,dum,itdav,maxdav,nfrdim1,solp,solp_res,eig,eig_res)
! Solve linear equation in Davidson subspace.

use casvb_global, only: ifollow, ipdd, iroot, jroot, nfrdim, nroot

implicit real*8(a-h,o-z)
dimension hp(maxdav,maxdav), eigval(itdav), eigvec(itdav,itdav)
dimension solp(maxdav), solp_res(maxdav)

if (ipdd >= 3) then
  write(6,*) ' HP matrix (b) :'
  call mxprint2_cvb(hp,itdav,maxdav,itdav,0)
end if

do it=1,itdav
  call fmove_cvb(hp(1,it),eigvec(1,it),itdav)
end do
call mxdiag_cvb(eigvec,eigval,itdav)

if (ifollow <= 2) then
  iroot = nroot
  jroot = mod(itdav,nroot)
  if (jroot == 0) jroot = nroot
  if ((itdav == maxdav) .or. (itdav == nfrdim)) jroot = nroot
  iroot = min(itdav,iroot)
  jroot = min(itdav,jroot)
  if (ifollow == 1) then
    iroot = itdav-iroot+1
    jroot = itdav-jroot+1
  end if
else if (ifollow == 3) then
  write(6,*) ' Overlap-based root following not yet implemented!'
  call abend_cvb()
else if (ifollow == 4) then
  ! Eigenvalue-based root following -- determine closest root:
  iroot = 1
  delmin = abs(eigval(1)-eig)
  do i=1,min(itdav,nroot)
    del = abs(eigval(i)-eig)
    if (del < delmin) then
      delmin = del
      iroot = i
    end if
  end do
  jroot = iroot
end if
eig = eigval(iroot)
call fmove_cvb(eigvec(1,iroot),solp,itdav)
eig_res = eigval(jroot)
call fmove_cvb(eigvec(1,jroot),solp_res,itdav)
if (ipdd >= 2) then
  write(6,'(a)') ' Eigenvalues :'
  call vecprint_cvb(eigval,itdav)
  write(6,'(a,i3,a)') ' Eigenvector number',iroot,' :'
  call vecprint_cvb(solp,itdav)
  if (jroot /= iroot) then
    write(6,'(a,i3,a)') ' Eigenvector number',jroot,' :'
    call vecprint_cvb(solp_res,itdav)
  end if
end if

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(dum)
  call Unused_integer(nfrdim1)
end if

end subroutine ddsol72_cvb
