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
* Copyright (C) 2004,2006,2007, Markus Reiher                          *
*               2004,2006, Alexander Wolf                              *
************************************************************************
      subroutine dkhparser_numeric (s,h,v,pvp,x,pxp,dkhorder,xorder,
     *               paramtype,dkhscfflg,nbas,isize,clight,vv,nn,dd,yy,
     *               ff,gg,xx,ii,jj,kk,ll,mm,revt,tran,sinv,aa,rr,
     *               tt,pp,e,ew,snumber,tnumber,unumber,scrno1,scrno2,
     *               scr1,scr2,scr3,scr5,scr8,no_hamil,no_prop,nbasp,
     *               nbaso,LDKroll,indx2,nAtom,maxsiz,nblock)
c
c*****************************************************************************************
c
c                  Implementation of the scalar-relativistic
c
c                       I N F I N I T E - O R D E R
c
c      G E N E R A L I Z E D   D O U G L A S - K R O L L - H E S S  (DKH)
c
c                         T R A N S F O R M A T I O N
c
c                                    F O R
c
c             H A M I L T O N I A N   A N D   P R O P E R T Y
c
c                 2006 Markus Reiher and Alexander Wolf, ETH Zurich
c                      Infinite-order DKH
c                      {markus.reiher,alexander.wolf}@phys.chem.ethz.ch
c
c
c   Reference to the infinite-order DKH method:
c
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 2037   part I   [Theory]
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 10945  part II  [Algorithm]
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064102 part III (Properties-Theory)
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064103 part IV  (Properties-Implem.
c
c   Reference to the generalized higher-order DKH method:
c
c     A. Wolf, M. Reiher, B.A. Hess, J. Chem. Phys. 117 (2002) 9215
c
c-----------------------------------------------------------------------------------------
c
c
c   This SR belongs to dkhparser_numeric (dkhparser2)
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena / ETH Zurich)
c
c   note: the package uses MOLCAS specific routines for
c           a) dynamical memory allocation via GetMem() in evalstring
c           b) IO to disk (ddafile etc.)
c         both may be replaced by a more general interface.
c
c   version:  2.2.0
c
c   last modified: 13.03.2007  M. Reiher (ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      logical dkhscfflg,no_hamil,no_prop,LDKroll
      integer dkhorder,xorder,nbas,isize,snumber,tnumber,nbasp,nbaso,
     *        unumber,scrno1,scrno2,nAtom,indx2(nAtom,4),maxsiz,nblock
      REAL*8 clight,s(isize),h(isize),v(isize),pvp(isize),
     *                 x(isize),pxp(isize)
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
c
      REAL*8 revt(nbas,nbas),tran(nbas,nbas),sinv(nbas,nbas),
     *                 aa(nbas),rr(nbas),tt(nbas),pp(nbas),e(nbas),
     *                 ew(nbas)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),scr5(nbas,nbas,scrno1),
     *                 scr8(isize,scrno2)
c
      character*(3) paramtype
C     REAL*8 DBLE
c
      integer poss(maxsnumber),post(maxsnumber),posu(maxunumber)
c
      integer i,iex,icontr
      REAL*8 det
c
*** addresses for matrices to disk in case of nbaso=1
***   ( isfreeunit() is an integer function of MOLCAS)
      integer adrmem,isfreeunit,dkh_48,adrnext
      parameter (adrmem=4000)
      integer dkhadr(adrmem)
#ifdef MOLPRO
      character*32 filnam
#endif
      dkh_48=0
c
c*******************************************************************************
      if (out_of_core) then
        dkh_48=48
#ifdef _MOLCAS_
        dkh_48=isfreeunit(dkh_48)
        call daname_mf(dkh_48,'DKHMX')
#else
        filnam='DKHMX'
        call tmp_filename(filnam)
        call daname_nolib(dkh_48,filnam)
#endif
      endif
c
      if (dkhorder.gt.maxorder) then
        write(stdout,1010) dkhorder,maxorder,maxorder
1010    format (/,'Warning: dkhorder = ',I2,', but maxorder = ',I2,'.',
     *          /'--> dkhorder will be reduced to maxorder = ',I2,'.'/)
        dkhorder=maxorder
      endif
c
CMR150107      call initialize2 (nbas,isize,vv,nn,dd,yy,ff,gg,xx,ii,jj,kk,
CMR150107     *                  ll,mm,revt,tran,sinv,aa,rr,tt,pp,e,ew,snumber,
CMR150107     *                  tnumber,unumber,scrno1,scrno2,nbasp,
CMR150107     *                  scr1,scr2,scr3,scr5,scr8,dkhorder,xorder)
CMR260207 "nbaso" has not yet been included in initialize2
      do 10 i=1,maxsnumber
        poss(i)=0
        post(i)=0
  10  continue
      do 11 i=1,maxunumber
        posu(i)=0
  11  continue
c
c     Verify that the overlap matrix is not singular
c
c  scr5(,,1) : overlap matrix (rectangular)
      call mat_sq_from_t (scr5(1,1,1),nbas,s)
      icontr=-1
      call mat_copy (scr5(1,1,2),nbas,nbas,scr5(1,1,1))
      call dcopiv (scr5(1,1,2),scr5(1,1,2),nbas,1,nbas,dkhzero,
     *             det,iex,icontr,scr8(1,3))
      if (icontr.ne.0) then
        write (stdout,1020) dkhzero
1020    format('  D K H |****** '/,
     *  '        |****** WARNING - OVERLAP MATRIX SINGULAR '/,
     *  '        |****** PIVOTAL ELEMENT LESS THAN ',D20.4,' FOUND'/,
     *  '        |******'//)
        call errex_rel(' D K H | singular overlap matrix')
      endif
*
*                                                                      *
************************************************************************
*     Local Douglas-Kroll with 2nd order perturbation                  *
*
      If (LDKroll) Then
         Call diag_ldkh(nAtom,nblock,nbas,isize,maxsiz,indx2,
     *                  scr5(1,1,1),scr8(1,1),s,h,scr5(1,1,2),e,
     *                  scr5(1,1,3),ew,tran,v,scr8(1,2),revt,scr8(1,3))
         Call ldkhpert(nbas,isize,h,s,sinv,scr5,scr8,nblock,
     *                 indx2,tran,ew,nAtom,revt)
         call mat_sq_from_t (scr5(1,1,1),nbas,s)
      Else
         Do i=1,isize
           scr8(i,1)=s(i)
         End Do
*** note: x is used as scratch in the following
      call sog(nbas,scr8(1,1),sinv,scr8(1,scrno2),scr8(1,2),scr5(1,1,3))
*** vorsicht: der vorhergehende Aufruf macht was falsch, wil da stand:
*** weshalb scrno2=10 verwendet wurde
***      call sog(nbas,scr8(1,1),sinv,scr8(1,scrno2),scr8(1,3),scr5(1,1,3scr5(1,1,3)))
*
         call diag_dkh(h,nbas,tran,ew,sinv,scr5(1,1,2),1)
      End If
      call calc_prefactors (nbas,isize,clight,aa,rr,tt,pp,e,ew)
c
c   e) Calculate reverse transformation 'revt'
c
      call calc_revt (nbas,revt,tran,sinv,scr5(1,1,1),scr5(1,1,2))
c
      call calc_E0 (nbas,isize,h,revt,ew)
      call calc_E1 (nbas,isize,xorder,dkhscfflg,h,v,pvp,x,pxp,vv,dd,xx,
     *              jj,revt,tran,sinv,aa,rr,tt,scrno1,scrno2,scr5,
     *              scr8,no_prop,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *              adrnext)
      call setup_matrices (nbas,isize,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,
     *                     kk,ll,mm,aa,rr,e,scrno2,scr8,no_prop,
     *                     nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
c      call trunc_dkh (nbas,isize,znuc,maxexp,minexp,dkhorder,
c     *                clight,scrno1,scrno2,scrno3,vv,nn,e,scr5,
c     *                scr8,scr9,tran,sinv)
c
      call calc_Uxxx (nbas,dkhorder,xorder,dkhscfflg,posu,post,poss,vv,
     *                nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     *                tnumber,unumber,scrno1,scr1,scr2,scr3,scr5,nbasp,
     *                nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
      call calc_Txxx (nbas,dkhorder,xorder,dkhscfflg,posu,post,poss,vv,
     *                nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     *                tnumber,unumber,scrno1,scr1,scr2,scr3,scr5,nbasp,
     *                nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
      call calc_Sxxx (nbas,dkhorder,xorder,dkhscfflg,posu,post,poss,vv,
     *                nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     *                tnumber,unumber,scrno1,scr1,scr2,scr3,scr5,nbasp,
     *                nbaso,dkhadr,adrmem,dkh_48,adrnext)
c
      if (.not.no_hamil)
     *  call calc_operators(nbas,isize,dkhorder,xorder,dkhscfflg,posu,
     *                      post,poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,
     *                      kk,ll,mm,e,snumber,tnumber,unumber,
     *                      scrno1,scrno2,scr1,scr2,scr3,scr5,scr8,
     *                      revt,h,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *                      adrnext)
c
      if (.not.dkhscfflg.and..not.no_prop)
     *  call calc_xoperators(nbas,isize,dkhorder,xorder,
     *                       dkhscfflg,posu,post,poss,vv,nn,dd,yy,
     *                       ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     *                       tnumber,unumber,scrno1,scrno2,scr1,
     *                       scr2,scr3,scr5,scr8,revt,x,nbasp,nbaso,
     *                       dkhadr,adrmem,dkh_48,adrnext)
c
      if (Out_Of_Core) then
#ifdef _MOLCAS_
        call DaEras(dkh_48)
#else
        call daclos(dkh_48)
#endif
      endif
c
CMR The following routine does not work any longer since all individual
CMR   operator contributions formerly stored in scr6 for H and scr7 for X
CMR   were removed.
CMR      call write_operators (nbas,isize,dkhorder,xorder,dkhscfflg,scrno2,
CMR     *                      scr6,scr7,scr8,paramtype,clight)
c
CMR      write (stdout,1150)
CMR1150  format ('Numeric parser routines "dkhparser2" successfully ',
CMR     *        'completed.',//90('-'),/2X)
c
      return
c Avoid unuser argument warnings
      if (.false.) Call Unused_character(paramtype)
      end
