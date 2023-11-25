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

subroutine dpgendet_cvb()

use casvb_global, only: ia12ind, iapr1, ib12ind, ibpr1, iconfs, idetvb, ixapr1, ixbpr1, nalf, nalf_fr, nbet, nbet_fr, nconf_fr, &
                        nconfion_fr, nda_fr, ndb_fr, ndetvb_fr, nel_fr, noe, norb, nfrag
use Data_Structures, only: Alloc1DiArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iapr_add, ibpr_add, iconfs_add, idetvb_add, ifrag, ixapr_add, ixbpr_add, mxstack, nalf_l, nbet_l, nda_l, ndb_l
type(Alloc1DiArray_Type) :: astr_fr(nfrag), bstr_fr(nfrag)

iapr_add = 1
ixapr_add = 1
ibpr_add = 1
ixbpr_add = 1
iconfs_add = 1
idetvb_add = 1
do ifrag=1,nfrag
  nalf_l = nalf_fr(1,ifrag)
  nbet_l = nbet_fr(1,ifrag)
  call icomb_cvb(norb,nalf_l,nda_l)
  call icomb_cvb(norb,nbet_l,ndb_l)

  nda_fr(1,ifrag) = nda_l
  nda_fr(1,ifrag) = nda_l

  call vbgendet_cvb(iapr1(iapr_add),ixapr1(ixapr_add),ibpr1(ibpr_add),ixbpr1(ixbpr_add),iconfs(:,iconfs_add),idetvb(idetvb_add), &
                    nconf_fr(ifrag),nconfion_fr(0,ifrag),nda_l,ndb_l,ndetvb_fr(ifrag),nel_fr(ifrag),noe,nalf_l,nbet_l,norb)
  iapr_add = iapr_add+ndetvb_fr(ifrag)
  ibpr_add = ibpr_add+ndetvb_fr(ifrag)
  ixapr_add = ixapr_add+nda_fr(1,ifrag)+1
  ixbpr_add = ixbpr_add+ndb_fr(1,ifrag)+1
  idetvb_add = idetvb_add+ndetvb_fr(ifrag)
  iconfs_add = iconfs_add+nconf_fr(ifrag)
end do

do ifrag=1,nfrag
  call mma_allocate(astr_fr(ifrag)%A,nalf_fr(1,ifrag)*nda_fr(1,ifrag),label='astr_fr')
  call mma_allocate(bstr_fr(ifrag)%A,nbet_fr(1,ifrag)*ndb_fr(1,ifrag),label='bstr_fr')
  call stringen_cvb(nel_fr(ifrag),nalf_fr(1,ifrag),astr_fr(ifrag)%A,bstr_fr(ifrag)%A)
end do

ia12ind(:) = 0
ib12ind(:) = 0
mxstack = 100
call detsort2_cvb(norb,nalf,nfrag,nda_fr(1,1),nalf_fr(1,1),ia12ind,astr_fr,mxstack)
call detsort2_cvb(norb,nbet,nfrag,ndb_fr(1,1),nbet_fr(1,1),ib12ind,bstr_fr,mxstack)

do ifrag=1,nfrag
  call mma_deallocate(astr_fr(ifrag)%A)
  call mma_deallocate(bstr_fr(ifrag)%A)
end do
call setiaprtot_cvb()

return

end subroutine dpgendet_cvb
