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
c
      subroutine dkhparser_driver (nbas,isize,dkhscfflg,dkhorder,xorder,
     *                             s,t,v,pvp,x,pxp,clight,paramtype,
     *                             dkhmemmax,mem,no_hamil,no_prop,nbasp,
     *                             nbaso,LDKroll,indx2,nAtom,maxsiz,
     *                             nblock,LDKpert)
c
c*************************************************************************************************
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
c   Reference to the arbitrary/infinite-order DKH method:
c
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 2037-2047   part I   [Theory]
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 10945-10956 part II  [Algorithm]
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064102      part III (Properties-Theory)
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064103      part IV  (Properties-Implementation)
c
c   Reference to the generalized DKH method:
c
c     A. Wolf, M. Reiher, B.A. Hess, J. Chem. Phys. 117 (2002) 9215-9226
c
c--------------------------------------------------------------------------------------------------
c
c   This SR is the memory manager for the numeric part of DKHPACK,
c     i.e., the infinite-order DKH package.
c     (you may replace this routine by C or F90 memory allocation)
c
c   This is the first routine to be called by dkhparser2.
c
c   Input:
c
c       nbas        No. of basisfcts.
c       isize       nbas*(nbas+1)/2
c       clight      velocity of light (137.035989500D0)
c                     (if this is zero, it is set here)
c       t(isize)    non-relativistic kinetic energy: p**2/2
c       v(isize)    external electron-nucleus potential
c       pvp(isize)  matrix representation of <pxVpx>+<pyVpy>+<pzVpz>
c                     in the given original basis set
c       s(isize)    overlap matrix
c       x(isize)    property matrix X
c       pxp(isize)  matrix representation of <pxXpx>+<pyXpy>+<pzXpz>
c                     in the given original basis set
c
c       no_hamil    flag to switch of the one-electron Hamiltonian transformation
c       no_prop     flag to switch of the one-electron property transformation
c       nbasp       one dimension of matrices for property trafo (=nbas, if props
c                     are to be calculated, =1 if no_prop)
c       nbaso       one dimension of matrices for Hamiltonian trafo (=nbas, if suff.
c                     RAM available, =1 if matrices are stored to disk)
c       LDKroll      flag for the Local DKH procedure
c
c   Output:
c
c       t(isize)    DKH Hamiltonian (pos.space) up to dkhorder
c       x(isize)    DKH transformed property X (pos.space) up to xorder
c
c
c   This SR belongs to dkhparser_numeric (dkhparser2)
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.2.0
c
c   last modified: 26.02.2007  M. Reiher (ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
      Integer DKHMemMax
c
      logical dkhscfflg,no_hamil,no_prop,LDKroll
      integer nbas,isize,dkhorder,xorder,nbasp,nbaso
c
      REAL*8 s(isize),t(isize),v(isize),pvp(isize),x(isize),
     *                 pxp(isize),clight
c
      integer nexts,snumber,tnumber,unumber,scrno1,scrno2,
     *        nbassq,i,nbassqp,nbassqo
      parameter (scrno1=3)
c
c  Note: scrno1 must be 4.
c        scrno2 must be 5 for hamil and prop or 3 for hamil-only.
c          hence, scrno2 is set below
c
      integer iv,in,id,iy,if,ig,ix,ii,ij,ik,il,im,irevt,itran,
     *        isinv,iaa,irr,itt,ipp,ie,iew,iscr1,iscr2,iscr3,
     *        iscr5,iscr8,iscrt,iscrsq,nAtom,maxsiz,nblock,
     *        ihsm,ieigsm,issm,indx2(nAtom),iewsm,exnbas
      Logical LDKpert
c
c  Allocate super vector 'mem'
      REAL*8 mem(dkhmemmax)
c
c
      character*(3) paramtype
      character*(maxlength) scrtxt
c
c
c------------------------------------------------------------------------------
c
c          Meaning of parameters, array sizes, and arrays:
c
c   dkhmemax   :  length of global memory array (set in dkhparameters.h)
c   maxorder :  max. DKH order 'parser2' can deal with (set in dkhparameters.h)
c
c   eigunit  :  unit for output of 1s-eigenvalues for Z=1,..,137
c                 ('eigenvalues.dat')
c   outunit1 :  unit for output of final even DKH operators stored as lower
c                 triangular matrices  ('dkhops.21')
c   outunit2 :  unit for output of final even property operators X_i
c                 (pos-space) stored as lower triangular matrices ('dkhops.22')
c  dkhunit1  :  unit for final DKH operators ('dkhops.11')
c  dkhunit2  :  unit for final property operators X ('dkhops.12')
c  dkhunit3  :  unit for auxiliary matrices Sxxx ('dkhops.13')
c  dkhunit4  :  unit for auxiliary matrices Txxx ('dkhops.14')
c  dkhunit5  :  unit for auxiliary matrices Uxxx ('dkhops.15')
c
c   snumber  :  no. of Sxxx matrices
c   tnumber  :  no. of Txxx matrices
c   unumber  :  no. of Uxxx matrices
c   scrno1   :  no. of temp. rectangular scratch matrices
c   scrno2   :  no. of temp. triangular scratch matrices
c
c   scr1(1,1,snumber) :  Sxxx matrices               (rectangular)
c   scr2(1,1,tnumber) :  Txxx matrices               (rectangular)
c   scr3(1,1,unumber) :  Uxxx matrices               (rectangular)
c   scr5(1,1,scrno1)  :  temp. scratch matrices      (rectangular)
c   scr8(1,scrno2)    :  temp. scratch matrices of
c                          length 'isize'            (triangular)
c
c
c
      if (no_prop) then
        scrno2=3
      else
        scrno2=5
      end if
c
c-----------------------------------------------------------
c  1. Read in number of auxiliary matrices Sxxx, Txxx, Uxxx,
c       and check clight, dkhorder, xorder
c-----------------------------------------------------------
c
c  Initialize scrtxt, clight
      do 10 i=1,maxlength
        scrtxt(i:i)=' '
  10  continue
      if (clight.eq.0.0d0) clight=137.035989500d0
c
c  Read in number 'snumber' of required scratch Sxxx-matrices, (i.e., the number
c    of employed Sxxx-matrices of symbolic PARSER procedure 'dkhparser1') and
c    the chosen unitary parametrization 'paramtype'
#ifndef MOLPRO
#ifdef _MOLCAS_
      Call Molcas_Open(dkhunit3,'dkhops.13')
#else
      open (dkhunit3, file='dkhops.13', status='OLD',
     *      form='FORMATTED')
#endif
#endif
      rewind(dkhunit3)
      read (dkhunit3,'(A50)') scrtxt(1:50)
      read (dkhunit3,'(A50)') scrtxt(1:50)
c  Check dkhorder and determine paramtype
      call ordertest (dkhorder,paramtype,scrtxt,1)
      read (dkhunit3,'(A50)') scrtxt(1:50)
c  Check xorder
      call ordertest (xorder,paramtype,scrtxt,0)
      read (dkhunit3,'(A50)') scrtxt(1:50)
      if (scrtxt(15:15).eq.'F' .and. dkhscfflg) then
        write (stdout,1040) dkhscfflg
1040    format (/2X,'ERROR in SR "dkhparser_driver": dkhscfflg = ',L1,
     *          ' does not correspond to the settings of "dkhops.13"',
     *          //2X,'STOP.',/2X)
        Call Abend()
      endif
      if (scrtxt(15:15).eq.'T' .and. .not.dkhscfflg) then
        write (stdout,1041) dkhscfflg
1041    format (/2X,'ERROR in SR "dkhparser_driver": dkhscfflg = ',L1,
     *          ' does not correspond to the settings of "dkhops.13"',
     *          //2X,'STOP.',/2X)
        call Abend
      endif
1050  read (dkhunit3,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1050
      read (dkhunit3,'(I3)') snumber
#ifndef MOLPRO
      close (dkhunit3)
#endif
c
c  Read in number 'tnumber' of required scratch Txxx-matrices, (i.e., the number
c    of employed Txxx-matrices of symbolic PARSER procedure 'dkhparser1')
#ifndef MOLPRO
#ifdef _MOLCAS_
      Call Molcas_Open(dkhunit4,'dkhops.14')
#else
      open (dkhunit4, file='dkhops.14', status='OLD',
     *      form='FORMATTED')
#endif
#endif
      rewind(dkhunit4)
1051  read (dkhunit4,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1051
      read (dkhunit4,'(I3)') tnumber
#ifndef MOLPRO
      close (dkhunit4)
#endif
c
c  Read in number 'unumber' of required scratch Uxxx-matrices, (i.e., the number
c    of employed Uxxx-matrices of symbolic PARSER procedure 'dkhparser1')
#ifndef MOLPRO
#ifdef _MOLCAS_
      Call Molcas_Open(dkhunit5,'dkhops.15')
#else
      open (dkhunit5, file='dkhops.15', status='OLD',
     *      form='FORMATTED')
#endif
#endif
      rewind(dkhunit5)
1052  read (dkhunit5,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1052
      read (dkhunit5,'(I3)') unumber
#ifndef MOLPRO
      close (dkhunit5)
#endif
      If (LDKroll.and.(.not.LDKpert)) Then
        exnbas=nbas
        nbas=maxsiz
        if (nbaso.eq.exnbas) nbaso=maxsiz
        if (nbasp.eq.exnbas) nbasp=maxsiz
      End If
c
c
c-----------------------------------
c  2. Allocate memory vector mem
c-----------------------------------
c
      nbassq=nbas*nbas
      nbassqp=nbasp*nbasp
      nbassqo=nbaso*nbaso
      nexts=1
c
c  Note: Full square matrices are (consequently) used throughout, though it
c        should be possible to use symmetry, i.e., lower triangular matrices
c        only.
c        --> This is left to future work!
c
c  Allocate basic matrices (V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M) first
      iv=nexts
      in=iv+nbassqo
      id=in+nbassqo
      iy=id+nbassqo
      if=iy+nbassqo
      ig=if+nbassqo
      ix=ig+nbassqo
      ii=ix+nbassqp
      ij=ii+nbassqp
      ik=ij+nbassqp
      il=ik+nbassqp
      im=il+nbassqp
      nexts=im+nbassqp
c
c  Allocate transformation matrices
      irevt=nexts
      itran=irevt+nbassq
      If (LDKroll.and.(.not.LDKpert)) Then
         isinv=itran+exnbas*maxsiz
         ieigsm=isinv+maxsiz*maxsiz
         ihsm=ieigsm+maxsiz*maxsiz
         issm=ihsm+maxsiz*(maxsiz+1)/2
         nexts=issm+maxsiz*maxsiz
      Else
         ieigsm=ip_Dummy
         ihsm=ip_Dummy
         issm=ip_Dummy
         isinv=itran+nbassq
         nexts=isinv+nbassq
      End If
c
c  Allocate kinematical factors
      iaa=nexts
      irr=iaa+nbas
      itt=irr+nbas
      If (LDKroll.and.(.not.LDKpert)) Then
        iewsm=itt+maxsiz
        ipp=iewsm+maxsiz
      Else
        iewsm=ip_Dummy
        ipp=itt+nbas
      End If
      ie=ipp+nbas
      iew=ie+nbas
      If (LDKroll.and.(.not.LDKpert)) Then
         nexts=iew+exnbas
      Else
         nexts=iew+nbas
      End If
c
c  Allocate (permanent) auxiliary matrices: Sxxx,  Txxx,  Uxxx
c                                           (scr1) (scr2) (scr3)
      iscr1=nexts
      iscr2=iscr1+nbassqo*snumber
      iscr3=iscr2+nbassqo*tnumber
      nexts=iscr3+nbassqo*unumber
c
c  Allocate rectangular (scr5), triangular (scr8)
c    scratch structures
      iscr5=nexts
      If (LDKroll.and.(.not.LDKpert)) Then
         iscr8=iscr5+maxsiz*exnbas*2
         iscrt=iscr8+isize*3
         iscrsq=iscrt+scrno2*maxsiz*(maxsiz+1)/2
         nexts=iscrsq+maxsiz*maxsiz*scrno1
      Else
         iscrt=ip_Dummy
         iscrsq=ip_Dummy
         iscr8=iscr5+nbassq*scrno1
         nexts=iscr8+isize*scrno2
      End If
c
c  Check memory
      if (nexts-1.ne.dkhmemmax) then
        write (stdout,1070) dkhmemmax,nexts
1070    format (2X,'SR dkhparser_driver: There is not enough memory ',
     *          'for DKH!',//2X,'Now:      dkhmemmax = ',I9,/2X,
     *          'Required: dkhmemmax = ',I9,//2X,'Increase parameter ',
     *          'dkhmemmax in "dkhparameters.h".',//2X,'STOP.',/)
        call Abend
      endif
c
c
c-------------------------------------------------------
c  3. Call dkh_main, where DKH Hamiltonian is evaluated
c-------------------------------------------------------
c
      If (LDKroll.and.(.not.LDKpert)) Then
         Call ldkh (s,t,v,pvp,x,pxp,dkhorder,xorder,paramtype,dkhscfflg,
     *           exnbas,isize,clight,mem(iv),mem(in),mem(id),mem(iy),
     *           mem(if),mem(ig),mem(ix),mem(ii),mem(ij),mem(ik),
     *           mem(im),mem(irevt),mem(itran),mem(isinv),mem(iaa),
     *           mem(irr),mem(itt),mem(ipp),mem(ie),mem(iew),snumber,
     *           tnumber,unumber,scrno1,scrno2,mem(iscr1),mem(iscr2),
     *           mem(iscr3),mem(iscr5),mem(iscr8),no_hamil,no_prop,
     *           nbasp,nbaso,indx2,maxsiz,nblock,nAtom,
     *           mem(iscrsq),mem(iscrt),mem(issm),mem(ihsm),
     *           mem(iewsm),mem(ieigsm))
      Else
      call dkhparser_numeric (s,t,v,pvp,x,pxp,dkhorder,xorder,paramtype,
     *           dkhscfflg,nbas,isize,clight,mem(iv),mem(in),mem(id),
     *           mem(iy),mem(if),mem(ig),mem(ix),
     *           mem(ii),mem(ij),mem(ik),mem(il),mem(im),mem(irevt),
     *           mem(itran),mem(isinv),mem(iaa),mem(irr),mem(itt),
     *           mem(ipp),mem(ie),mem(iew),snumber,tnumber,unumber,
     *           scrno1,scrno2,mem(iscr1),mem(iscr2),mem(iscr3),
     *           mem(iscr5),mem(iscr8),no_hamil,no_prop,nbasp,nbaso,
     *           LDKroll,indx2,nAtom,maxsiz,nblock)
      End If
c
      return
      end
c
