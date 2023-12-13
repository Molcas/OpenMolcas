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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               1995,1996, Pavel Neogrady                              *
!***********************************************************************

subroutine mod1(nsym,nfro,nish,nssh,ndel,norb,nfror,ndelr,firas,fi,epsras,eps)
! this routine does:
! 1) reduce firas, epsras if nfror>nfro, ndelr>ndel
! 2) redefine nfro,nish,nssh,ndel,norb to proper ones

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nsym, ndel(8), nfror(8), ndelr(8)
integer(kind=iwp), intent(inout) :: nfro(8), nish(8), nssh(8), norb(8)
real(kind=wp), intent(in) :: firas(*), epsras(*)
real(kind=wp), intent(_OUT_) :: fi(*), eps(*)
integer(kind=iwp) :: isym, ndd, ndf, nlow, nup, p, pnew, pqnew, pqras, pras, q

!1 reduce fi

pqras = 0
pqnew = 0
do isym=1,nsym

  ndf = nfror(isym)-nfro(isym)
  ndd = ndelr(isym)-ndel(isym)
  nlow = ndf+1
  nup = norb(isym)-ndd

  do p=1,norb(isym)
    do q=1,p
      pqras = pqras+1

      if ((p >= nlow) .and. (p <= nup)) then
        if ((q >= nlow) .and. (q <= nup)) then
          pqnew = pqnew+1
          fi(pqnew) = firas(pqras)
        end if
      end if

    end do
  end do

end do

!2 reduce eps

pras = 0
pnew = 0
do isym=1,nsym

  ndf = nfror(isym)-nfro(isym)
  ndd = ndelr(isym)-ndel(isym)
  nlow = ndf+1
  nup = norb(isym)-ndd

  do p=1,norb(isym)
    pras = pras+1

    if ((p >= nlow) .and. (p <= nup)) then
      pnew = pnew+1
      eps(pnew) = epsras(pras)
    end if

  end do

end do

!3 define new nfro,nish,nssh,ndel,norb

nish(1:nsym) = nish(1:nsym)-nfror(1:nsym)+nfro(1:nsym)
nssh(1:nsym) = nssh(1:nsym)-ndelr(1:nsym)+ndel(1:nsym)
norb(1:nsym) = norb(1:nsym)-nfror(1:nsym)+nfro(1:nsym)-ndelr(1:nsym)+ndel(1:nsym)
nfro(1:nsym) = nfror(1:nsym)

return

end subroutine mod1
