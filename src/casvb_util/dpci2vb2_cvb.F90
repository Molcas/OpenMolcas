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

subroutine dpci2vb2_cvb(civec,cvbdet,dvbdet,evbdet,ic1,ret,ic)

use casvb_global, only: ia12ind, iapr1, ib12ind, ixapr1, nda, nda_fr, ndb, ndb_fr, ndetvb, ndetvb_fr, nfrag
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: civec(nda,ndb), cvbdet(ndetvb), evbdet(ndetvb)
real(kind=wp), intent(in) :: dvbdet(ndetvb)
integer(kind=iwp), intent(in) :: ic1, ic
real(kind=wp), intent(_OUT_) :: ret
integer(kind=iwp) :: i, ia, ia_ci, ib, ib_ci, idetvb, ifr, iter, ixa, jfr, kfr, mxiter, ndetvb_add, nestlevel, nloop
real(kind=wp) :: cf, cf2, cinrm, cnrm, fac, fac1
logical(kind=iwp) :: fail
integer(kind=iwp), allocatable :: idetind(:), ipr_off(:), istack(:), ixapr_off(:), ixbpr_off(:), nc_facalf(:), nc_facbet(:), &
                                  ncindalf(:), ncindbet(:)
real(kind=wp), allocatable :: coeff(:)
integer(kind=iwp), parameter :: mxstack = 100
real(kind=wp), external :: ddot_

if (ic == 0) then
  cvbdet(:) = Zero
else if ((ic == 1) .or. (ic == 4)) then
  civec(:,:) = Zero
else if (ic == 3) then
  ret = Zero
else if (ic == 5) then
  evbdet(:) = Zero
end if

call mma_allocate(nc_facalf,nfrag,label='nc_facalf')
call mma_allocate(nc_facbet,nfrag,label='nc_facbet')
call mma_allocate(ncindalf,[0,nfrag],label='ncindalf')
call mma_allocate(ncindbet,[0,nfrag],label='ncindbet')
call mma_allocate(coeff,[0,nfrag],label='coeff')
call mma_allocate(ipr_off,nfrag,label='ipr_off')
call mma_allocate(ixapr_off,nfrag,label='ixapr_off')
call mma_allocate(ixbpr_off,nfrag,label='ixbpr_off')

!FIXME: There's some inconsistency in the nd[ab]_fr indices,
!       the fragment index is supposed to be the second one
do ifr=1,nfrag
  if (ifr == 1) then
    ipr_off(ifr) = 0
    ixapr_off(ifr) = 0
    ixbpr_off(ifr) = 0
  else
    ipr_off(ifr) = ipr_off(ifr-1)+ndetvb_fr(ifr-1)
    ixapr_off(ifr) = ixapr_off(ifr-1)+nda_fr(ifr-1,1)+1
    ixbpr_off(ifr) = ixbpr_off(ifr-1)+ndb_fr(ifr-1,1)+1
  end if
end do

cinrm = Zero
coeff(0) = One
ncindalf(0) = 1
ncindbet(0) = 1
do i=1,nfrag
  if (i == 1) then
    nc_facalf(i) = 1
    nc_facbet(i) = 1
  else
    nc_facalf(i) = nc_facalf(i-1)*nda_fr(i-1,1)
    nc_facbet(i) = nc_facbet(i-1)*ndb_fr(i-1,1)
  end if
end do

nloop = nfrag
! MXITERS -> NDA_FR

! Following is code for a set of nested loops. To deal with the
! complication that the number of nested loops is not known at
! compile time, a simple integer stack is used.
! NESTLEVEL=1 signifies we are doing outermost loop and so on.

call mma_allocate(istack,mxstack,label='istack')
call mma_allocate(idetind,nfrag,label='idetind')

nestlevel = 0
call istkinit_cvb(istack,mxstack)

outer: do
  ! Here we go to the beginning of the next loop in the sequence:
  if (nestlevel < nloop) then
    nestlevel = nestlevel+1
    iter = 0
    mxiter = ndetvb_fr(nestlevel)
    call istkpush_cvb(istack,iter)
    call istkpush_cvb(istack,mxiter)
  end if

  ! Here we do the next loop iteration of the current loop:
  do
    if (nestlevel == 0) exit outer
    call istkpop_cvb(istack,mxiter)
    call istkpop_cvb(istack,iter)
    iter = iter+1
    if (iter > mxiter) then
      nestlevel = nestlevel-1
    else
      call istkpush_cvb(istack,iter)
      call istkpush_cvb(istack,mxiter)
      exit
    end if
  end do

  ! Here goes the code specific to this loop level.
  idetvb = 0
  ixa = 0 ! dummy initialize
  fail = .true.
  do ia=1,nda_fr(nestlevel,1)
    do ixa=ixapr1(ia+ixapr_off(nestlevel)),ixapr1(ia+1+ixapr_off(nestlevel))-1
      idetvb = idetvb+1
      if (idetvb == iter) then
        fail = .false.
        exit
      end if
    end do
  end do
  if (fail) then
    write(u6,*) ' Error in DPCI2VB '
    call abend_cvb()
  end if
  ib = iapr1(ixa+ipr_off(nestlevel))

  idetind(nestlevel) = iter+ipr_off(nestlevel)
  if (((ic == 1) .and. (ic1 == 0)) .or. (ic == 3)) then
    coeff(nestlevel) = coeff(nestlevel-1)*cvbdet(iter+ipr_off(nestlevel))
  else if (ic /= 4) then
    coeff(nestlevel) = dvbdet(iter+ipr_off(nestlevel))
  end if

  ncindalf(nestlevel) = ncindalf(nestlevel-1)+nc_facalf(nestlevel)*(ia-1)
  ncindbet(nestlevel) = ncindbet(nestlevel-1)+nc_facbet(nestlevel)*(ib-1)

  if (nestlevel == nfrag) then
    ia_ci = ia12ind(ncindalf(nestlevel))
    if (ia_ci /= 0) then
      ib_ci = ib12ind(ncindbet(nestlevel))
      if (ib_ci /= 0) then
        if (ic == 0) then
          if (ic1 == 0) then
            ! --  CI2VBC  --
            cinrm = cinrm+civec(abs(ia_ci),abs(ib_ci))*civec(abs(ia_ci),abs(ib_ci))
            if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
              do ifr=1,nfrag
                cvbdet(idetind(ifr)) = cvbdet(idetind(ifr))+civec(abs(ia_ci),abs(ib_ci))
              end do
            else
              do ifr=1,nfrag
                cvbdet(idetind(ifr)) = cvbdet(idetind(ifr))-civec(abs(ia_ci),abs(ib_ci))
              end do
            end if
          else if (ic1 == 2) then
            ! --  CI2VBG  --
            do ifr=1,nfrag
              cf = One
              do jfr=1,nfrag
                if (jfr /= ifr) cf = cf*coeff(jfr)
              end do
              if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
                cvbdet(idetind(ifr)) = cvbdet(idetind(ifr))+cf*civec(abs(ia_ci),abs(ib_ci))
              else
                cvbdet(idetind(ifr)) = cvbdet(idetind(ifr))-cf*civec(abs(ia_ci),abs(ib_ci))
              end if
            end do
          end if
        else if (ic == 1) then
          if (ic1 == 0) then
            ! --  VB2CIC  --
            if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
              civec(abs(ia_ci),abs(ib_ci)) = coeff(nestlevel)
            else
              civec(abs(ia_ci),abs(ib_ci)) = -coeff(nestlevel)
            end if
          else if (ic1 == 1) then
            ! --  VB2CIF  --
            do ifr=1,nfrag
              cf = One
              do jfr=1,nfrag
                if (jfr /= ifr) cf = cf*coeff(jfr)
              end do
              if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
                civec(abs(ia_ci),abs(ib_ci)) = civec(abs(ia_ci),abs(ib_ci))+cf*cvbdet(idetind(ifr))
              else
                civec(abs(ia_ci),abs(ib_ci)) = civec(abs(ia_ci),abs(ib_ci))-cf*cvbdet(idetind(ifr))
              end if
            end do
          end if
        else if (ic == 2) then
          ! --  VB2CIAF  --
          do ifr=1,nfrag
            cf = One
            do jfr=1,nfrag
              if (jfr /= ifr) cf = cf*coeff(jfr)
            end do
            if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
              civec(abs(ia_ci),abs(ib_ci)) = civec(abs(ia_ci),abs(ib_ci))+cf*cvbdet(idetind(ifr))
            else
              civec(abs(ia_ci),abs(ib_ci)) = civec(abs(ia_ci),abs(ib_ci))-cf*cvbdet(idetind(ifr))
            end if
          end do
        else if (ic == 3) then
          ! --   VB2CIDOT  --
          if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
            ret = ret+civec(abs(ia_ci),abs(ib_ci))*coeff(nestlevel)
          else
            ret = ret-civec(abs(ia_ci),abs(ib_ci))*coeff(nestlevel)
          end if
        else if (ic == 4) then
          ! --  SETIAPR  --
          civec(abs(ia_ci),abs(ib_ci)) = One
        else if (ic == 5) then
          ! --  CI2ORDR  --
          do ifr=1,nfrag
            cf = Zero
            do jfr=1,nfrag
              if (jfr /= ifr) then
                cf2 = cvbdet(idetind(jfr))
                do kfr=1,nfrag
                  if ((kfr /= ifr) .and. (kfr /= jfr)) cf2 = cf2*dvbdet(idetind(kfr))
                end do
                cf = cf+cf2
              end if
            end do
            if ((ia_ci > 0) .eqv. (ib_ci > 0)) then
              evbdet(idetind(ifr)) = evbdet(idetind(ifr))+cf*civec(abs(ia_ci),abs(ib_ci))
            else
              evbdet(idetind(ifr)) = evbdet(idetind(ifr))-cf*civec(abs(ia_ci),abs(ib_ci))
            end if
          end do
        end if
      end if
    end if
  end if

end do outer

call mma_deallocate(nc_facalf)
call mma_deallocate(nc_facbet)
call mma_deallocate(ncindalf)
call mma_deallocate(ncindbet)
call mma_deallocate(istack)
call mma_deallocate(coeff)
call mma_deallocate(idetind)
call mma_deallocate(ipr_off)
call mma_deallocate(ixapr_off)
call mma_deallocate(ixbpr_off)

! This is the end ...

if ((ic == 0) .and. (ic1 == 0)) then
  ! "Normalize" the coefficients for each fragment:
  fac = One/sqrt(cinrm**(One/real(nfrag,kind=wp)))
  ndetvb_add = 0
  do ifr=1,nfrag
    cnrm = ddot_(ndetvb_fr(ifr),cvbdet(ndetvb_add+1),1,cvbdet(ndetvb_add+1),1)
    fac1 = fac*cnrm
    cvbdet(ndetvb_add+1:ndetvb_add+ndetvb_fr(ifr)) = fac1*sqrt(cnrm)*cvbdet(ndetvb_add+1:ndetvb_add+ndetvb_fr(ifr))
    ndetvb_add = ndetvb_add+ndetvb_fr(ifr)
  end do
end if

return

end subroutine dpci2vb2_cvb
