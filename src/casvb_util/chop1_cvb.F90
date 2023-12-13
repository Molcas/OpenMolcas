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

subroutine chop1_cvb()

use casvb_global, only: absym, i1alf, i1bet, i1c, ia12ind, iafrm, iapr, iapr1, iato, ib12ind, ibfrm, ibpr, ibpr1, ibto, icfrm, &
                        iconfs, icto, idetvb, ixapr, ixapr1, ixbpr, ixbpr1, n1a, n1b, nalf, nalf_fr, nam1, naprodvb, nbet, &
                        nbet_fr, nbm1, nbprodvb, ncivb, nconf, nda, nda_fr, ndb, ndb_fr, ndet, ndetvb, ndetvb_fr, nel, nfrag, noe, &
                        norb, npvb, phato, phbto, phcto, release
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i, ifrag, ndavb, ndbvb

if (release(1)) then
  call mma_deallocate(i1alf)
  call mma_deallocate(iafrm)
  call mma_deallocate(iato)
  call mma_deallocate(phato)
  nullify(i1bet)
  nullify(ibfrm)
  nullify(ibto)
  nullify(phbto)
  if (allocated(i1c)) call mma_deallocate(i1c)
  if (allocated(icfrm)) call mma_deallocate(icfrm)
  if (allocated(icto)) call mma_deallocate(icto)
  if (allocated(phcto)) call mma_deallocate(phcto)
  call mma_deallocate(iapr)
  call mma_deallocate(ixapr)
  call mma_deallocate(ibpr)
  call mma_deallocate(ixbpr)
  call mma_deallocate(iconfs)
  call mma_deallocate(idetvb)
  call mma_deallocate(ia12ind)
  call mma_deallocate(ib12ind)
  call mma_deallocate(iapr1)
  call mma_deallocate(ixapr1)
  call mma_deallocate(ibpr1)
  call mma_deallocate(ixbpr1)
end if
release(1) = .true.
release(2) = .false.

! Dimensions
call icomb_cvb(norb,nalf,nda)
call icomb_cvb(norb,nbet,ndb)
!FIXME: There's some inconsistency in the *_fr indices,
!       the fragment index is supposed to be the second one
do i=1,nfrag
  call icomb_cvb(norb,nalf_fr(i,1),nda_fr(i,1))
  call icomb_cvb(norb,nbet_fr(i,1),ndb_fr(i,1))
end do
call icomb_cvb(norb-1,nalf-1,n1a)
call icomb_cvb(norb-1,nbet-1,n1b)
call icomb_cvb(norb,nalf-1,nam1)
call icomb_cvb(norb,nbet-1,nbm1)
ndet = nda*ndb
! Symmetry of determinant strings:
call getnci_cvb(ncivb,nel,nalf-nbet,0)
! Identical indexing arrays may share memory:
call mma_allocate(i1alf,n1a,norb,label='i1alf')
call mma_allocate(iafrm,norb,nda,label='iafrm')
call mma_allocate(iato,[1,norb],[0,nam1],label='iato')
call mma_allocate(phato,norb,nam1,label='phato')
if (absym(4)) then
  i1bet => i1alf
  ibfrm => iafrm
  ibto => iato
  phbto => phato
else
  call mma_allocate(i1c,n1b,norb,label='i1c')
  call mma_allocate(icfrm,norb,ndb,label='icfrm')
  call mma_allocate(icto,[1,norb],[0,nbm1],label='icto')
  call mma_allocate(phcto,norb,nbm1,label='phcto')
  i1bet => i1c
  ibfrm => icfrm
  ibto => icto
  phbto => phcto
end if

! Determinant dimensioning for VB wavefunction:
ndavb = 0
ndbvb = 0
naprodvb = 1
nbprodvb = 1
npvb = 1
do ifrag=1,nfrag
  ndavb = ndavb+nda_fr(1,ifrag)+1
  ndbvb = ndbvb+ndb_fr(1,ifrag)+1
  naprodvb = naprodvb*nda_fr(1,ifrag)
  nbprodvb = nbprodvb*ndb_fr(1,ifrag)
  npvb = npvb*ndetvb_fr(ifrag)
end do
if (nfrag <= 1) then
  naprodvb = 0
  nbprodvb = 0
end if

call mma_allocate(iapr,npvb,label='iapr')
call mma_allocate(ixapr,nda+1,label='ixapr')
call mma_allocate(ibpr,npvb,label='ibpr')
call mma_allocate(ixbpr,ndb+1,label='ixbpr')
call mma_allocate(iconfs,noe,nconf,label='iconfs')
call mma_allocate(idetvb,ndetvb,label='idetvb')
call mma_allocate(ia12ind,naprodvb,label='ia12ind')
call mma_allocate(ib12ind,nbprodvb,label='ib12ind')
call mma_allocate(iapr1,ndetvb,label='iapr1')
call mma_allocate(ixapr1,ndavb,label='ixapr1')
call mma_allocate(ibpr1,ndetvb,label='ibpr1')
call mma_allocate(ixbpr1,ndbvb,label='ixbpr1')

return

end subroutine chop1_cvb
