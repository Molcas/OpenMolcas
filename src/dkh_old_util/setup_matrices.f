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
      subroutine setup_matrices (nbas,isize,vv,nn,dd,yy,ff,gg,pp,xx,
     *                           ii,jj,kk,ll,mm,aa,rr,e,scrno2,
     *                           scr8,no_prop,nbasp,nbaso,dkhadr,
     *                           adrmem,dkh_48,adrnext)
c
************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     It completes what has been started in calc_E1().
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.2.0
c
c   last modified: 13.03.2007 (M. Reiher, ETH Zurich)
c          * only tiny changes but one may optimize this routine
c            by replacing all hand-made matrix routines by library routines
c          * no_prop and nbasp introduced
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer nbas,isize,scrno2,nbasp,nbaso,iscr
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),aa(nbas),
     *                 rr(nbas),pp(nbas),e(nbas),scr8(isize,scrno2)
      logical no_prop
c
      integer dkh_48,adrmem
      integer one,two,dkhadr(adrmem),adr,adrnext
      parameter(one=1,two=2)
c
c-----------------------------------------------------------------------
c
c    V    vv(nbas,nbas)  :  external electron--nucleus potential
c    N    nn(nbas,nbas)  :  damped potential V~ = vv(i,j)/(e(i)+e(j))
c    D    dd(nbas,nbas)  :  PVP
c    Y    yy(nbas,nbas)  :  [PVP] = PVP(i,j)/(e(i)+e(j))
c    F    ff(nbas,nbas)  :  E01 = E_1 = even term of first order
c    G    gg(nbas,nbas)  :  PE01P = P*E_1*P
c
c    X    xx(nbas,nbas)  :  nr property matrix X
c    I    ii(nbas,nbas)  :  damped property X~ = xx(i,j)/(e(i)+e(j))
c    J    jj(nbas,nbas)  :  PXP
c    K    kk(nbas,nbas)  :  [PXP] = PXP(i,j)/(e(i)+e(j))
c    L    ll(nbas,nbas)  :  CEO=Xeven(0)
c    M    mm(nbas,nbas)  :  P(CE0)P
c    ... if existing, otherwise nbasp=1
c
c  Note: store vv=V, dd=pVp, xx=X, jj=pXp  for set up of tilde terms
c
      if (Out_Of_Core) then
CMR
CMR     allocation can be skipped if aa() or rr(), which are no longer
CMR     needed in this revised routine, are destroyed
CMR
        call GetMem('SetMat1 ','ALLO','REAL',iscr,isize+4)
        adr=dkhadr(1)
        call ddafile(dkh_48,two,work(iscr),isize,adr)
        call mat_1_over_h_tri (work(iscr),nbas,e)
        dkhadr(5)=adrnext
        adr=dkhadr(5)
        call ddafile(dkh_48,one,work(iscr),isize,adr)
        dkhadr(6)=adr
c
        adr=dkhadr(2)
        call ddafile(dkh_48,two,work(iscr),isize,adr)
        call mat_1_over_h_tri (work(iscr),nbas,e)
        adr=dkhadr(6)
        call ddafile(dkh_48,one,work(iscr),isize,adr)
c
        dkhadr(7)=adr
        call ddafile(dkh_48,one,scr8(1,1),isize,adr)
        dkhadr(8)=adr
        call ddafile(dkh_48,one,scr8(1,2),isize,adr)
        adrnext=adr
      else
        call mat_copy (nn,nbas,nbas,vv)
        call mat_copy (yy,nbas,nbas,dd)
        call mat_1_over_h (nn,nbas,e)
        call mat_1_over_h (yy,nbas,e)
        call mat_sq_from_t (ff,nbas,scr8(1,1))
        call mat_sq_from_t (gg,nbas,scr8(1,2))
      end if
c
      if (.not.no_prop) then
        if (Out_Of_Core) then
          adr=dkhadr(3)
          call ddafile(dkh_48,two,work(iscr),isize,adr)
          call mat_1_over_h_tri (work(iscr),nbas,e)
          dkhadr(9)=adrnext
          adr=dkhadr(9)
          call ddafile(dkh_48,one,work(iscr),isize,adr)
          dkhadr(10)=adr
c
          adr=dkhadr(4)
          call ddafile(dkh_48,two,work(iscr),isize,adr)
          call mat_1_over_h_tri (work(iscr),nbas,e)
          adr=dkhadr(10)
          call ddafile(dkh_48,one,work(iscr),isize,adr)
          dkhadr(11)=adr
c
          call ddafile(dkh_48,one,scr8(1,3),isize,adr)
          dkhadr(12)=adr
          call ddafile(dkh_48,one,scr8(1,4),isize,adr)
          adrnext=adr
c
        else
          call mat_copy (ii,nbasp,nbasp,xx)
          call mat_copy (kk,nbasp,nbasp,jj)
          call mat_1_over_h (ii,nbasp,e)
          call mat_1_over_h (kk,nbasp,e)
          call mat_sq_from_t (ll,nbasp,scr8(1,3))
          call mat_sq_from_t (mm,nbasp,scr8(1,4))
        end if
      end if
c
      if (Out_Of_Core)
     *  call GetMem('SetMat1 ','FREE','REAL',iscr,isize+4)
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_real_array(pp)
        call Unused_real_array(aa)
        call Unused_real_array(rr)
      end if
      end
