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

subroutine applyt2_cvb(vec,gjorb,igjorb,i1alf,i1bet,iato,ibto,phato,phbto)
! Apply T(O) to the vector VEC. O is defined in terms of GJORB.

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: vec(nda,ndb), gjorb(norb*norb), phato(norb,nam1), phbto(norb,nbm1)
integer(kind=iwp) :: igjorb(2,norb*norb), i1alf(n1a,norb), i1bet(n1b,norb), iato(norb,0:nam1), ibto(norb,0:nbm1)
integer(kind=iwp) :: ia, iak, iax, iaxtmp, ib, ibk, ibx, ibxtmp, ij, iorb, jax, jbx, jorb
real(kind=wp) :: scl, tcof
real(kind=wp), parameter :: thresh = 1.0e-10_wp

do ij=1,norb*norb
  iorb = igjorb(2,ij)
  jorb = igjorb(1,ij)
  scl = gjorb(ij)
  if ((iorb /= jorb) .and. (abs(scl) > thresh)) then
    ! a) Alpha excitation
    if (absym(2)) then
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          tcof = scl*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          if (jax > iax) then
            call daxpy_(ndb-jax+1,tcof,vec(iax,jax),nda,vec(jax,jax),nda)
          else
            call daxpy_(iax-jax,tcof,vec(jax,iax),1,vec(jax,jax),nda)
            vec(jax,iax) = vec(jax,iax)+tcof*vec(iax,iax)
            if (ndb-iax > 0) call daxpy_(ndb-iax,tcof,vec(iax,iax+1),nda,vec(jax,iax+1),nda)
          end if
        end if
      end do
    else
      do ia=1,n1a
        iaxtmp = i1alf(ia,iorb)
        jax = iato(jorb,iaxtmp)
        if (jax /= 0) then
          iax = iato(iorb,iaxtmp)
          tcof = scl*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
          call daxpy_(ndb,tcof,vec(iax,1),nda,vec(jax,1),nda)
        end if
      end do
    end if
    ! b) Beta excitation
    if (absym(2)) then
      do ib=1,n1b
        ibxtmp = i1bet(ib,iorb)
        jbx = ibto(jorb,ibxtmp)
        if (jbx /= 0) then
          ibx = ibto(iorb,ibxtmp)
          tcof = scl*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
          if (jbx > ibx) then
            call daxpy_(ibx-1,tcof,vec(1,ibx),1,vec(1,jbx),1)
            vec(ibx,jbx) = vec(ibx,jbx)+tcof*vec(ibx,ibx)
            if (jbx-ibx > 0) call daxpy_(jbx-ibx,tcof,vec(ibx,ibx+1),nda,vec(ibx+1,jbx),1)
          else
            call daxpy_(jbx,tcof,vec(1,ibx),1,vec(1,jbx),1)
          end if
        end if
      end do
    else
      do ib=1,n1b
        ibxtmp = i1bet(ib,iorb)
        jbx = ibto(jorb,ibxtmp)
        if (jbx /= 0) then
          ibx = ibto(iorb,ibxtmp)
          tcof = scl*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
          call daxpy_(nda,tcof,vec(1,ibx),1,vec(1,jbx),1)
        end if
      end do
    end if
  else if ((iorb == jorb) .and. (abs(scl-One) > thresh)) then
    ! Alpha singly occupied
    if (absym(2)) then
      do ia=1,n1a
        iak = iato(iorb,i1alf(ia,iorb))
        call dscal_(ndb-iak+1,scl,vec(iak,iak),nda)
      end do
    else
      do ia=1,n1a
        iak = iato(iorb,i1alf(ia,iorb))
        call dscal_(ndb,scl,vec(iak,1),nda)
      end do
    end if
    ! Beta singly occupied
    if (absym(2)) then
      do ib=1,n1b
        ibk = ibto(iorb,i1bet(ib,iorb))
        call dscal_(ibk,scl,vec(1,ibk),1)
      end do
    else
      do ib=1,n1b
        ibk = ibto(iorb,i1bet(ib,iorb))
        call dscal_(nda,scl,vec(1,ibk),1)
      end do
    end if
  end if
end do
if (absym(2)) then
  do ia=1,nda
    do ib=ia+1,ndb
      vec(ib,ia) = vec(ia,ib)
    end do
  end do
end if

return

end subroutine applyt2_cvb
