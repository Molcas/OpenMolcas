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

use casvb_global, only: absym, n1a, n1b, nam1, nbm1, nda, ndb, norb
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vec(nda,ndb)
real(kind=wp), intent(in) :: gjorb(norb*norb), phato(norb,nam1), phbto(norb,nbm1)
integer(kind=iwp), intent(in) :: igjorb(2,norb*norb), i1alf(n1a,norb), i1bet(n1b,norb), iato(norb,0:nam1), ibto(norb,0:nbm1)
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
            vec(jax,jax:) = vec(jax,jax:)+tcof*vec(iax,jax:)
          else
            vec(jax,jax:iax) = vec(jax,jax:iax)+tcof*vec(jax:iax,iax)
            if (iax < ndb) vec(jax,iax+1:) = vec(jax,iax+1:)+tcof*vec(iax,iax+1:)
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
          vec(jax,:) = vec(jax,:)+tcof*vec(iax,:)
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
            vec(1:ibx,jbx) = vec(1:ibx,jbx)+tcof*vec(1:ibx,ibx)
            vec(ibx+1:jbx,jbx) = vec(ibx+1:jbx,jbx)+tcof*vec(ibx,ibx+1:jbx)
          else
            vec(1:jbx,jbx) = vec(1:jbx,jbx)+tcof*vec(1:jbx,ibx)
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
          vec(:,jbx) = vec(:,jbx)+tcof*vec(:,ibx)
        end if
      end do
    end if
  else if ((iorb == jorb) .and. (abs(scl-One) > thresh)) then
    ! Alpha singly occupied
    if (absym(2)) then
      do ia=1,n1a
        iak = iato(iorb,i1alf(ia,iorb))
        vec(iak,iak:) = scl*vec(iak,iak:)
      end do
    else
      do ia=1,n1a
        iak = iato(iorb,i1alf(ia,iorb))
        vec(iak,:) = scl*vec(iak,:)
      end do
    end if
    ! Beta singly occupied
    if (absym(2)) then
      do ib=1,n1b
        ibk = ibto(iorb,i1bet(ib,iorb))
        vec(1:ibk,ibk) = scl*vec(1:ibk,ibk)
      end do
    else
      do ib=1,n1b
        ibk = ibto(iorb,i1bet(ib,iorb))
        vec(:,ibk) = scl*vec(:,ibk)
      end do
    end if
  end if
end do
if (absym(2)) then
  do ia=1,nda
    vec(ia+1:ndb,ia) = vec(ia,ia+1:ndb)
  end do
end if

return

end subroutine applyt2_cvb
