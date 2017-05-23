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
* Copyright (C) 2004, Alexander Wolf                                   *
*               2004,2007, Markus Reiher                               *
************************************************************************
      subroutine initialize2 (nbas,isize,vv,nn,dd,yy,ff,gg,xx,ii,
     *                        jj,kk,ll,mm,revt,tran,sinv,aa,rr,tt,pp,e,
     *                        ew,snumber,tnumber,unumber,
     *                        scrno1,scrno2,nbasp,scr1,scr2,scr3,
     *                        scr5,scr8,dkhorder,xorder)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 14.01.2007 (M. Reiher, ETH Zurich)
c            * a couple of matrices are no longer set equal to zero
c             as this is not really necessary
c            * in addition, the last dimension of some arrays is now
c             '*' so that the call to this routine remains unchanged
c             if the no_prop option is switched on
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer nbas,isize,snumber,tnumber,unumber,scrno1,nbasp,
     *        scrno2,dkhorder,xorder
      integer i,j,k
c
      REAL*8 vv(nbas,*),nn(nbas,*),dd(nbas,*),
     *                 yy(nbas,*),ff(nbas,*),gg(nbas,*),
     *                 xx(nbasp,*),
     *                 ii(nbasp,*),jj(nbasp,*),kk(nbasp,*),
     *                 ll(nbasp,*),mm(nbasp,*)
c
      REAL*8 revt(nbas,nbas),tran(nbas,nbas),sinv(nbas,nbas),
     *                 aa(nbas),rr(nbas),tt(nbas),pp(nbas),e(nbas),
     *                 ew(nbas)
c
      REAL*8 scr1(nbas,nbas,snumber),scr2(nbas,nbas,tnumber),
     *                 scr3(nbas,nbas,unumber),scr5(nbas,nbas,scrno1),
     *                 scr8(isize,scrno2)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
CMR          vv(i,j)=0.0d0
CMR          nn(i,j)=0.0d0
CMR          dd(i,j)=0.0d0
CMR          yy(i,j)=0.0d0
CMR          ff(i,j)=0.0d0
CMR          gg(i,j)=0.0d0
          revt(i,j)=0.0d0
          tran(i,j)=0.0d0
          sinv(i,j)=0.0d0
CMR The following matrices are of dimension nbasp X nbasp
CMR          xx(i,j)=0.0d0
CMR          ii(i,j)=0.0d0
CMR          jj(i,j)=0.0d0
CMR          kk(i,j)=0.0d0
CMR          ll(i,j)=0.0d0
CMR          mm(i,j)=0.0d0
  20    continue
  10  continue
      do 40 i=1,nbas
        aa(i)=0.0d0
        rr(i)=0.0d0
        tt(i)=0.0d0
        pp(i)=0.0d0
        e(i)=0.0d0
        ew(i)=0.0d0
  40  continue
c
      do 110 k=1,snumber
        do 112 j=1,nbas
          do 114 i=1,nbas
            scr1(i,j,k)=0.0d0
 114      continue
 112    continue
 110  continue
c
      do 120 k=1,tnumber
        do 122 j=1,nbas
          do 124 i=1,nbas
            scr2(i,j,k)=0.0d0
 124      continue
 122    continue
 120  continue
c
      do 130 k=1,unumber
        do 132 j=1,nbas
          do 134 i=1,nbas
            scr3(i,j,k)=0.0d0
 134      continue
 132    continue
 130  continue
c
      do 150 k=1,scrno1
        do 152 j=1,nbas
          do 154 i=1,nbas
            scr5(i,j,k)=0.0d0
 154      continue
 152    continue
 150  continue
c
CMR      do 180 k=1,scrno2
CMR        do 182 j=1,isize
CMR          scr8(j,k)=0.0d0
CMR 182    continue
CMR 180  continue
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(vv)
        call Unused_real_array(nn)
        call Unused_real_array(dd)
        call Unused_real_array(yy)
        call Unused_real_array(ff)
        call Unused_real_array(gg)
        call Unused_real_array(xx)
        call Unused_real_array(ii)
        call Unused_real_array(jj)
        call Unused_real_array(kk)
        call Unused_real_array(ll)
        call Unused_real_array(mm)
        call Unused_real_array(scr8)
        call Unused_integer(dkhorder)
        call Unused_integer(xorder)
      end if
      end
