************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine detsort2_cvb(
     >  xalf,norb,nalf,nfrag,nda_fr,nalf_fr,
     >  nalf_acc,ia12ind,iphase,nc_fac,ncombindex,
     >  iastr,iastr_off,iastr_acc,istack,mxstack)
      implicit real*8 (a-h,o-w,y-z),integer(x)
      dimension xalf(0:norb,0:nalf)
      dimension istack(mxstack)
      dimension nalf_fr(nfrag)
      dimension nalf_acc(nfrag)
      dimension nda_fr(nfrag)
      dimension iastr(*),iastr_off(nfrag)
      dimension iastr_acc(norb,nfrag)
      dimension ncombindex(0:nfrag)
      dimension iphase(nfrag),nc_fac(nfrag)
      dimension ia12ind(*)

      call weightfl_cvb(xalf,nalf,norb)

      ncombindex(0)=1
      do i=1,nfrag
      if(i.eq.1)then
        nalf_acc(i)=nalf_fr(i)
        nc_fac(i)=1
      else
        nalf_acc(i)=nalf_acc(i-1)+nalf_fr(i)
        nc_fac(i)=nc_fac(i-1)*nda_fr(i-1)
      endif
      enddo

      nloop=nfrag
c  MXITERS -> NDA_FR

c  Following is code for a set of nested loops. To deal with the
c  complication that the number of nested loops is not known at
c  compile time, a simple integer stack is used.
c  NESTLEVEL=1 signifies we are doing outermost loop and so on.

      nestlevel=0
      call istkinit_cvb(istack,mxstack)

c  Here we go to the beginning of the next loop in the sequence :
1000  continue
      if(nestlevel.lt.nloop)then
        nestlevel=nestlevel+1
        iter=0
        mxiter=nda_fr(nestlevel)
        call istkpush_cvb(istack,iter)
        call istkpush_cvb(istack,mxiter)
      endif

c  Here we do the next loop iteration of the current loop :
2000  if(nestlevel.eq.0)goto 3000
      call istkpop_cvb(istack,mxiter)
      call istkpop_cvb(istack,iter)
      iter=iter+1
      if(iter.gt.mxiter)then
        nestlevel=nestlevel-1
        goto 2000
      else
        call istkpush_cvb(istack,iter)
        call istkpush_cvb(istack,mxiter)
      endif

c  Here goes the code specific to this loop level.
      if(nestlevel.eq.1)then
        call imove_cvb(
     >    iastr(1+nalf_fr(nestlevel)*(iter-1)+iastr_off(nestlevel)-1),
     >    iastr_acc(1,nestlevel),nalf_fr(nestlevel))
        iphase(nestlevel)=1
      else
        iphase(nestlevel)=iphase(nestlevel-1)*
     >    ioemrg2_cvb(iastr_acc(1,nestlevel-1),nalf_acc(nestlevel-1),
     >    iastr(1+nalf_fr(nestlevel)*(iter-1)+iastr_off(nestlevel)-1),
     >    nalf_fr(nestlevel),iastr_acc(1,nestlevel))
        if(iphase(nestlevel).eq.0)goto 1000
      endif
      ncombindex(nestlevel)=ncombindex(nestlevel-1)
     >  +nc_fac(nestlevel)*(iter-1)
      if(nestlevel.eq.nfrag)then
        iatotindx=minind_cvb(iastr_acc(1,nestlevel),nalf,norb,xalf)
        ia12ind(ncombindex(nestlevel))=iatotindx*iphase(nestlevel)
      endif

      goto 1000

c  This is the end ...
3000  continue

      return
      end
