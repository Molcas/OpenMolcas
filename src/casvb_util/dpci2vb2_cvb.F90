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

subroutine dpci2vb2_cvb(civec,cvbdet,dvbdet,evbdet,ic1,ret,ic,nda,ndb,ndetvb,nfrag,nda_fr,ndb_fr,ia12ind,ib12ind,nc_facalf, &
                        nc_facbet,ncindalf,ncindbet,istack,mxstack,coeff,idetind,iapr,ixapr,ipr_off,ixapr_off,ixbpr_off,ndetvb_fr, &
                        ndavb)

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: ic1, ic, nda, ndb, ndetvb, nfrag, nda_fr(nfrag), ndb_fr(nfrag), ia12ind(*), ib12ind(*), nc_facalf(nfrag), &
                     nc_facbet(nfrag), ncindalf(0:nfrag), ncindbet(0:nfrag), mxstack, istack(mxstack), idetind(nfrag), &
                     iapr(ndetvb), ndavb, ixapr(ndavb), ipr_off(nfrag), ixapr_off(nfrag), ixbpr_off(nfrag), ndetvb_fr(nfrag)
real(kind=wp) :: civec(nda,ndb), cvbdet(ndetvb), dvbdet(ndetvb), evbdet(ndetvb), ret, coeff(0:nfrag)
integer(kind=iwp) :: i, ia, ia_ci, ib, ib_ci, idetvb, ifr, iter, ixa, jfr, kfr, mxiter, ndetvb_add, nestlevel, nloop
real(kind=wp) :: cf, cf2, cinrm, cnrm, fac, fac1
logical(kind=iwp) :: fail
real(kind=wp), external :: ddot_

if (ic == 0) then
  call fzero(cvbdet,ndetvb)
else if ((ic == 1) .or. (ic == 4)) then
  call fzero(civec,nda*ndb)
else if (ic == 3) then
  ret = Zero
else if (ic == 5) then
  call fzero(evbdet,ndetvb)
end if

do ifr=1,nfrag
  if (ifr == 1) then
    ipr_off(ifr) = 0
    ixapr_off(ifr) = 0
    ixbpr_off(ifr) = 0
  else
    ipr_off(ifr) = ipr_off(ifr-1)+ndetvb_fr(ifr-1)
    ixapr_off(ifr) = ixapr_off(ifr-1)+nda_fr(ifr-1)+1
    ixbpr_off(ifr) = ixbpr_off(ifr-1)+ndb_fr(ifr-1)+1
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
    nc_facalf(i) = nc_facalf(i-1)*nda_fr(i-1)
    nc_facbet(i) = nc_facbet(i-1)*ndb_fr(i-1)
  end if
end do

nloop = nfrag
! MXITERS -> NDA_FR

! Following is code for a set of nested loops. To deal with the
! complication that the number of nested loops is not known at
! compile time, a simple integer stack is used.
! NESTLEVEL=1 signifies we are doing outermost loop and so on.

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
  do ia=1,nda_fr(nestlevel)
    do ixa=ixapr(ia+ixapr_off(nestlevel)),ixapr(ia+1+ixapr_off(nestlevel))-1
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
  ib = iapr(ixa+ipr_off(nestlevel))

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

! This is the end ...

if ((ic == 0) .and. (ic1 == 0)) then
  ! "Normalize" the coefficients for each fragment:
  fac = One/sqrt(cinrm**(One/real(nfrag,kind=wp)))
  ndetvb_add = 1
  do ifr=1,nfrag
    cnrm = ddot_(ndetvb_fr(ifr),cvbdet(ndetvb_add),1,cvbdet(ndetvb_add),1)
    fac1 = fac*cnrm
    call dscal_(ndetvb_fr(ifr),fac1*sqrt(cnrm),cvbdet(ndetvb_add),1)
    ndetvb_add = ndetvb_add+ndetvb_fr(ifr)
  end do
end if

return

end subroutine dpci2vb2_cvb
