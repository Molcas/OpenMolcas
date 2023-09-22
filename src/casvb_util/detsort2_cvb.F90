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

subroutine detsort2_cvb(xalf,norb,nalf,nfrag,nda_fr,nalf_fr,nalf_acc,ia12ind,iphase,nc_fac,ncombindex,iastr,iastr_off,iastr_acc, &
                        istack,mxstack)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: norb, nalf, xalf(0:norb,0:nalf), nfrag, nda_fr(nfrag), nalf_fr(nfrag), nalf_acc(nfrag), ia12ind(*), &
                     iphase(nfrag), nc_fac(nfrag), ncombindex(0:nfrag), iastr(*), iastr_off(nfrag), iastr_acc(norb,nfrag), &
                     mxstack, istack(mxstack)
integer(kind=iwp) :: i, iatotindx, iter, mxiter, nestlevel, nloop
integer(kind=iwp), external :: ioemrg2_cvb, minind_cvb

call weightfl_cvb(xalf,nalf,norb)

ncombindex(0) = 1
do i=1,nfrag
  if (i == 1) then
    nalf_acc(i) = nalf_fr(i)
    nc_fac(i) = 1
  else
    nalf_acc(i) = nalf_acc(i-1)+nalf_fr(i)
    nc_fac(i) = nc_fac(i-1)*nda_fr(i-1)
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
    mxiter = nda_fr(nestlevel)
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
  if (nestlevel == 1) then
    call imove_cvb(iastr(1+nalf_fr(nestlevel)*(iter-1)+iastr_off(nestlevel)-1),iastr_acc(1,nestlevel),nalf_fr(nestlevel))
    iphase(nestlevel) = 1
  else
    iphase(nestlevel) = iphase(nestlevel-1)* &
                        ioemrg2_cvb(iastr_acc(1,nestlevel-1),nalf_acc(nestlevel-1), &
                                    iastr(1+nalf_fr(nestlevel)*(iter-1)+iastr_off(nestlevel)-1),nalf_fr(nestlevel), &
                                    iastr_acc(1,nestlevel))
    if (iphase(nestlevel) == 0) cycle outer
  end if
  ncombindex(nestlevel) = ncombindex(nestlevel-1)+nc_fac(nestlevel)*(iter-1)
  if (nestlevel == nfrag) then
    iatotindx = minind_cvb(iastr_acc(1,nestlevel),nalf,norb,xalf)
    ia12ind(ncombindex(nestlevel)) = iatotindx*iphase(nestlevel)
  end if

end do outer

! This is the end ...
return

end subroutine detsort2_cvb
