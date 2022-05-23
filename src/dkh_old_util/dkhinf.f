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
* Copyright (C) 2004-2007, Markus Reiher                               *
*               2004-2006, Alexander Wolf                              *
************************************************************************
c
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
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 2037    part I   [Theory]
c     M. Reiher, A. Wolf, J.Chem.Phys. 121 (2004) 10945   part II  [Algorithm]
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064102  part III (Properties: theory)
c     A. Wolf, M. Reiher, J.Chem.Phys. 124 (2006) 064103  part IV  (Properties: implem.)
c
c   Reference to the generalized higher-order DKH method:
c
c     A. Wolf, M. Reiher, B.A. Hess, J. Chem. Phys. 117 (2002) 9215-9226
c
c
c*****************************************************************************************
c
c
c
c
      subroutine adjust_param (dkhorder,xorder)
c
c***************************************************************************
c
c   Adjust parameters "maxoperators" and "maxuops" approriately,
c   i.e., depending on the chosen value of dkhorder.
c
c   written by:  M. Reiher and A. Wolf  (Univ. Jena)
c
c   last modified:  14.11.2006 (MR, ETH Zurich)
c
c***************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder
c
      maxoperators=0
      maxuops=0
c
      if (dkhorder.le.6) then
         maxoperators=4000+xorder*2500
         maxuops=500+xorder*100
      else if (dkhorder.eq.7) then
         maxoperators=7000+xorder*2500
         maxuops=500+xorder*100
      else if (dkhorder.eq.8) then
        maxoperators=9500+xorder*8000
        maxuops=500+xorder*250
      else if (dkhorder.eq.9) then
        maxoperators=12500+xorder*10000
        maxuops=500+xorder*250
      else if (dkhorder.eq.10) then
        maxoperators=19000+xorder*17000
        maxuops=1000+xorder*300
      else if (dkhorder.eq.11) then
        maxoperators=30000+xorder*30000
        maxuops=1500+xorder*200
      else if (dkhorder.eq.12) then
        maxoperators=100000+xorder*100000
        maxuops=2000+xorder*200
      else
        write (stdout,1001)
 1001   format (/2X,'ERROR: In Subroutine adjust_param (dkhinf.f):',
     *    //11X,'Increase values for parameters ',
     *    'maxoperators and maxuops and recompile the code.',//2X,
     *    'STOP.',/2X)
         CALL Abend
      endif
c
      return
      end
c

c
      subroutine build_U (i,uniorder,dkhorder,xorder,dkhscfflg,ucounter,
     *               uopsleng,uoporder,eouops,duops,puop,uops,opcounter,
     *               operleng,oporder,evenodd,doperators,operators,
     *               wordercounter,wopscounter,wopsleng,woporder,eowops,
     *               dwops,wops,oddcounter,oddleng,oddorder,eoodd,dodd,
     *               odd,ducoeffs)
c
c*************************************************************************
c
c   Construct actual high-level U_i-operator, i.e., W_i and W_i^n.
c
c   Simplify it by neglect of all terms, which are irrelevant for
c     the chosen order, and eliminate all terms with 'coefficient=0'.
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 18.07.2005
c
c   first version: 22.07.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer i,uniorder,dkhorder,xorder
      integer ucounter,uopsleng(maxuops),uoporder(maxuops),
     *        eouops(maxuops),puop(maxuops)
      real*8 duops(maxuops)
      character*(maxlength) uops(maxuops)
c
      integer opcounter,operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        wordercounter(maxorder),wopscounter,wopsleng(maxuops),
     *        woporder(maxuops,3),eowops(maxuops),
     *        oddcounter,oddleng(maxuops),oddorder(maxuops,3),
     *        eoodd(maxuops)
      real*8 doperators(maxoperators),dwops(maxuops),
     *                 dodd(maxuops)
      character*(maxlength) wops(maxuops),odd(maxuops)
#if defined(_MOLCAS_) || defined(MOLPRO)
      character operators(*)
      character*(maxlength) opstring
#else
      character*(maxlength) operators(maxoperators)
#endif
c
      real*8 ducoeffs(maxorder)
c
      character*(3) scrchar, dkh_int2char
      integer maxpuop,j,k,idum
      intrinsic DBLE,INT
      do 10 j=1,maxuops
        uopsleng(j)=0
        uoporder(j)=0
        eouops(j)=0
        duops(j)=0.0d0
        puop(j)=0
        do 20 k=1,maxlength
          uops(j)(k:k)=' '
  20    continue
  10  continue
c   wordercounter(i):  number of terms (high-level) contributing to W_i.
c   NB: wordercounter(i)=oddcounter for actual 'i'
      wordercounter(i)=0
      oddcounter=0
      idum=0
      do 30 j=1,i-1
        idum=idum+wordercounter(j)
  30  continue
c
      do 50 j=1,opcounter
        if (evenodd(j).eq.-1 .and. oporder(j,3).eq.i) then
          wordercounter(i)=wordercounter(i)+1
          woporder(idum+wordercounter(i),1)=oporder(j,1)
          woporder(idum+wordercounter(i),2)=oporder(j,2)
          woporder(idum+wordercounter(i),3)=oporder(j,3)
          eowops(idum+wordercounter(i))=evenodd(j)
          dwops(idum+wordercounter(i))=doperators(j)
          wopsleng(idum+wordercounter(i))=operleng(j)+3
          if (wopsleng(idum+wordercounter(i)).gt.maxlength) then
            write (stdout,1001)
1001        format (/2X,'ERROR1 in SR "build_U":  wopsleng > maxlength',
     *              //2X,'STOP.',/2X)
            CALL Abend
          endif
#if defined(_MOLCAS_) || defined(MOLPRO)
          call get_dkoperators(j,opstring,operators)
          wops(idum+wordercounter(i))='B['//opstring(1:operleng(j))//']'
#else
          wops(idum+wordercounter(i))(1:wopsleng(idum+wordercounter(i)))
     *             ='B['//operators(j)(1:operleng(j))//']'
#endif
          oddcounter=oddcounter+1
          oddorder(oddcounter,1)=oporder(j,1)
          oddorder(oddcounter,2)=oporder(j,2)
          oddorder(oddcounter,3)=oporder(j,3)
          eoodd(oddcounter)=evenodd(j)
          dodd(oddcounter)=doperators(j)
          oddleng(oddcounter)=operleng(j)
#if defined(_MOLCAS_) || defined(MOLPRO)
          odd(oddcounter)=opstring(1:operleng(j))
#else
          odd(oddcounter)(1:oddleng(oddcounter))=
     *             operators(j)(1:operleng(j))
#endif
        endif
  50  continue
c
      wopscounter=wopscounter+wordercounter(i)
      maxpuop=INT(DBLE(max(dkhorder,xorder))/DBLE(i))
      if(i.ge.1.and.i.le.99) then
          scrchar=dkh_int2char(i)
          scrchar(1:1)='W'
      endif
cvv      do 13 k=1,99
cvv         if (i.eq.k) then
cvv          scrchar=dkh_int2char(k)
cvv          scrchar(1:1)='W'
cvv        endif
cvv 13   continue
      if (uniorder.gt.99) then
        write (stdout,1002)
1002    format (/2X,'Note: So far only 99 unitary transformations U_i',
     *          ' are possible in SR "build_U".')
        CALL Abend
      endif
c
      ucounter=1
      uoporder(1)=i
      eouops(1)=-1
      duops(1)=ducoeffs(1)
      puop(1)=1
      uopsleng(1)=3
      uops(1)(1:uopsleng(1))=scrchar
c
      do 100 j=2,maxpuop
        ucounter=ucounter+1
        uoporder(ucounter)=i*j
        eouops(ucounter)=(-1)**j
        duops(ucounter)=ducoeffs(j)
        puop(j)=j
        uopsleng(ucounter)=uopsleng(ucounter-1)+3
        if (uopsleng(ucounter).gt.maxlength) then
CMR040307        if (wopsleng(idum+wordercounter(i)).gt.maxlength) then
          write (stdout,1003)
CMR0403071003      format (/2X,'ERROR2 in SR "build_U:  wopsleng > maxlength'
1003      format (/2X,'ERROR2 in SR "build_U:  uleng > maxlength',
     *            //2X,'STOP.',/2X)
          CALL Abend
        endif
        uops(ucounter)(1:uopsleng(ucounter))=scrchar//uops(ucounter-1)
 100  continue
2000  continue
      do 110 j=1,ucounter
       if (abs(duops(j)).lt.dkhzero) then
         do 120 k=j+1,ucounter
            uoporder(k-1)=uoporder(k)
            eouops(k-1)=eouops(k)
            duops(k-1)=duops(k)
            puop(k-1)=puop(k)
            uopsleng(k-1)=uopsleng(k)
            uops(k-1)=uops(k)
 120      continue
          ucounter=ucounter-1
          goto 2000
        endif
 110  continue
CMR      if (i.le.9) then
CMR        write (stdout,1005) i,wordercounter(i)
CMR1005    format (2X,'Number of terms contributing to W',I1,':',2X,I6)
CMR      endif
CMR      if (i.ge.10) then
CMR        write (stdout,1006) i,wordercounter(i)
CMR1006    format (2X,'Number of terms contributing to W',I2,':',1X,I6)
CMR      endif
c
      if (ucounter.gt.maxuops) then
        write (stdout,1007)
1007    format (/2X,'WARNING: Too many U-operators (ucounter >',
     *          ' maxuops)!',/15X,'Increase maxuops ',
     *          '(e.g., by one order of magnitude)!')
        CALL Abend
      endif
c
CMR      if (i.le.9) then
CMR        write (stdout,1009) i,ucounter
CMR1009    format (2X,'Number of terms contributing to U',I1,':',2X,I6)
CMR      endif
CMR      if (i.ge.10) then
CMR        write (stdout,1010) i,ucounter
CMR1010    format (2X,'Number of terms contributing to U',I2,':',1X,I6)
CMR      endif
      if (dbgflg.ge.1 .and. i.le.9) write (dbgunit,1021) i,i,
     *                                                  wordercounter(i)
1021  format (//2X,'Transformation U',I1,' constructed. Operator W',I1,
     *        ' (',I4,' terms) stored in wops.')
      if (dbgflg.ge.1 .and.i.ge.10) write (dbgunit,1022) i,i,
     *                                                  wordercounter(i)
1022  format (//2X,'Transformation U',I2,' constructed. Operator W',I2,
     *        ' (',I4,' terms) stored in wops.')
c
      if (dbgflg.ne.0.and.(dbgflg.ge.2 .or. i.eq.uniorder)) then
        if (i.le.9) write (dbgunit,1023) i
1023    format (/2X,'Form of stored W-operators (wops) after ',
     *        'transformation U',I1,':',/2X,58('-'))
        if (i.ge.10) write (dbgunit,1024) i
1024    format (/2X,'Form of stored W-operators (wops) after ',
     *        'transformation U',I2,':',/2X,59('-'))
        call output2 (dbgunit,wopscounter,wopsleng,woporder,eowops,
     *                dwops,wops)
      endif
c
      if (dbgflg.ge.2) then
        write (dbgunit,1025)
1025    format(2X/)
        if (i.le.9.and.dbgflg.gt.0) then
          write (dbgunit,1027) i
1027      format (2X,'Terms contributing to transformation U',I1,':',
     *            /2X,40('-'))
        endif
        if (i.ge.10.and.dbgflg.gt.0) then
          write (dbgunit,1028) i
1028      format (2X,'Terms contributing to transformation U',I2,':',
     *            /2X,41('-'))
        endif
        call output2b (dbgunit,dkhscfflg,ucounter,uopsleng,uoporder,
     *                 eouops,duops,uops)
      endif
c
      return
      end
c
c
c
      subroutine calc_E0 (nbas,isize,h,revt,ew)
c
************************************************************************
c
c   Calculate  E0=h  (position space)
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 15.01.2007 (M. Reiher, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer nbas,isize,i,j,ij,k
      real*8 h(isize),revt(nbas,nbas),ew(nbas)
      ij=0
      do 100 i=1,nbas
        do 110 j=1,i
          ij=ij+1
          h(ij)=0.0D0
          do 120 k=1,nbas
            h(ij)=h(ij)+revt(i,k)*revt(j,k)*ew(k)
  120     continue
  110   continue
  100 continue
c
      return
      end
c
c
      subroutine calc_E1 (nbas,isize,xorder,dkhscfflg,h,v,pvp,x,pxp,vv,
     *                    dd,xx,jj,revt,tran,sinv,aa,rr,tt,scrno1,
     *                    scrno2,scr5,scr8,no_prop,nbasp,nbaso,dkhadr,
     *                    adrmem,dkh_48,adrnext)
c
************************************************************************
c
c   a) Calculate  E1=scr8(,1)  and  PE1P=scr8(,2)
c                 X0=scr8(,3)  and  PX0P=scr8(,4)
c
c   b) Set up  vv=V,  dd=pVp,  xx=X,  jj=pXp
c
c   c) CAUTION: x() is used a scratch here, so x() must be defined even
c               if no_prop=.true. !
c               (it contains the final result at the end but until then
c                it is considered scratch)
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena / ETH Zurich)
c
c   version:  2.2.0
c
c   last modified: 26.02.2007  (M. Reiher, ETH Zurich)
c            * completely rewritten to save memory
c            * required memory now: scrno2=5
c              (or only scrno2=3 if no_prop=.true.)
c            * note scrno1=3 here
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg,no_prop
      integer nbas,isize,xorder,scrno1,scrno2,nbasp,nbaso
      real*8 h(isize),v(isize),pvp(isize),x(isize),pxp(isize),
     *                 vv(nbaso,nbaso),dd(nbaso,nbaso),
     *                 xx(nbasp,nbasp),jj(nbasp,nbasp)
      real*8 revt(nbas,nbas),tran(nbas,nbas),sinv(nbas,nbas),
     *                 aa(nbas),rr(nbas),tt(nbas)
      real*8 scr5(nbas,nbas,scrno1),scr8(isize,scrno2)
c
      integer i,j,ij
c
      integer dkh_48,adrmem
      integer one,dkhadr(adrmem),adr,adrnext
      parameter(one=1)
c
c     consider Hamiltonian first; transform V to p**2-basis
c
      call trsm_dkh(v,sinv,scr5(1,1,3),nbas,scr5(1,1,1),scr5(1,1,2))
      call trsm_dkh(scr5(1,1,3),tran,scr8(1,3),nbas,scr5(1,1,1),
     *              scr5(1,1,2))
***   now: scr8(,3)=V (momentum-space)
c
c     Store V in momentum-space as square matrix vv if Out_Of_Core=.false.
c
      call mat_axa_tri (scr8(1,3),nbas,aa)
      if (Out_Of_Core) then
        dkhadr(1)=0
        adr=dkhadr(1)
        call ddafile(dkh_48,one,scr8(1,3),isize,adr)
      else
        call mat_sq_from_t (vv,nbas,scr8(1,3))
***     now:  vv=AVA (momentum space) (rectangular)
      end if
c
c     Multiply V with aa and rr factors
c
      ij=0
      do 10 i=1,nbas
        do 20 j=1,i
          ij=ij+1
          scr8(ij,1)=scr8(ij,3)
          scr8(ij,2)=4.0D0*rr(i)*rr(i)*tt(i)*
     *               scr8(ij,3)*tt(j)*rr(j)*rr(j)
  20   continue
  10  continue
***   now:  scr8(,1) = AVA (momentum-space)
c
c     Transform pVp integrals to p**2-basis
c
      call trsm_dkh(pvp,sinv,scr5(1,1,3),nbas,scr5(1,1,1),scr5(1,1,2))
      call trsm_dkh(scr5(1,1,3),tran,scr8(1,3),nbas,scr5(1,1,1),
     *              scr5(1,1,2))
***   now:  scr8(,3) = pVp (momentum-space)
c
c     Store pVp in momentum-space as square matrix dd
c
      call mat_arxra_tri (scr8(1,3),nbas,aa,rr)
      if (Out_Of_Core) then
        dkhadr(2)=adr
        call ddafile(dkh_48,one,scr8(1,3),isize,adr)
        adrnext=adr
      else
        call mat_sq_from_t (dd,nbas,scr8(1,3))
***     now:  dd=ApVpA (momentum space) (rectangular)
      end if
c
c     Multiply pVp  with aa and rr factors
c
      ij=0
      do 30 i=1,nbas
        do 40 j=1,i
          ij=ij+1
          scr8(ij,1)=scr8(ij,1)+scr8(ij,3)
          scr8(ij,2)=scr8(ij,2)+scr8(ij,3)
***   now: scr8(,2) = ARpVpRa = APVPA (momentum space)
***   now: scr8(,1) = E1 = AVA + APVPA (momentum-space)
c
  40    continue
  30  continue
c
      if (.not.no_prop) then
c
c     next: consider property operator; transform X to p**2-basis
c
        call trsm_dkh(x,sinv,scr5(1,1,3),nbas,scr5(1,1,1),scr5(1,1,2))
        call trsm_dkh(scr5(1,1,3),tran,scr8(1,5),nbas,scr5(1,1,1),
     *                scr5(1,1,2))
***   now: scr8(,5)=X (momentum-space)
c
c     Store X in momentum-space as square matrix xx
c
        call mat_axa_tri (scr8(1,5),nbas,aa)
        if (Out_Of_Core) then
          dkhadr(3)=adr
          call ddafile(dkh_48,one,scr8(1,5),isize,adr)
        else
          call mat_sq_from_t (xx,nbas,scr8(1,5))
***       now: xx=AXA (momentum space) (rectangular)
        end if
c
c     Multiply X with aa and rr factors
c
        ij=0
        do 11 i=1,nbas
          do 21 j=1,i
            ij=ij+1
            scr8(ij,3)=scr8(ij,5)
            scr8(ij,4)=4.0D0*rr(i)*rr(i)*tt(i)*
     *                 scr8(ij,5)*tt(j)*rr(j)*rr(j)
  21     continue
  11    continue
***   scr8(,3) = AXA (momentum-space)
c
c     Transform pXp integrals to p**2-basis
c
        call trsm_dkh(pxp,sinv,scr5(1,1,3),nbas,scr5(1,1,1),scr5(1,1,2))
        call trsm_dkh(scr5(1,1,3),tran,scr8(1,5),nbas,scr5(1,1,1),
     *                scr5(1,1,2))
***   now:  scr8(,5) = pXp (momentum-space)
c
c     Store pXp in momentum-space as square matrix jj
c
        call mat_arxra_tri (scr8(1,5),nbas,aa,rr)
        if (Out_Of_Core) then
          dkhadr(4)=adr
          call ddafile(dkh_48,one,scr8(1,5),isize,adr)
          adrnext=adr
        else
          call mat_sq_from_t (jj,nbas,scr8(1,5))
***       now: jj=ApXpA (momentum space) (rectangular)
        end if
c
c     Multiply pXp  with aa and rr factors
c
        ij=0
        do 31 i=1,nbas
          do 41 j=1,i
            ij=ij+1
            scr8(ij,3)=scr8(ij,3)+scr8(ij,5)
            scr8(ij,4)=scr8(ij,4)+scr8(ij,5)
***   now: scr8(,4) = ARpXpRa = APXPA (momentum space)
***   now: scr8(,3) = X0 = AXA + APXPA (momentum-space)
c
  41      continue
  31    continue
c
      end if
c
c  --> triangular E1, PE1P and X0, PX0P (momentum space) set up finished
c
c
c   the next steps are different for variational and perturbative treatment of X
c
      if (dkhscfflg.and..not.no_prop) then
***   scr8(,5) = E1 = E1 + X0 = E1(lambda) (momentum space)
        call mat_copy2 (scr8(1,5),isize,scr8(1,1))
        if (xorder.gt.0) call mat_tadd (1.0d0,scr8(1,5),isize,1.0d0,
     *                                  scr8(1,3))
c     transform E1 back to position-space
        call trsmt(scr8(1,5),revt,x,nbas,scr5(1,1,1),
     *               scr5(1,1,2))
***   now:  x = E1(X) (position space)
c
c     add E1 to final Hamiltonian
        call mat_tadd (1.0d0,h,isize,1.0d0,x)
c
      else
c
c     transform E1, X0 back to position-space
        call trsmt(scr8(1,1),revt,scr5(1,1,3),nbas,scr5(1,1,1),
     *               scr5(1,1,2))
c     add E1 to final Hamiltonian
        call mat_tadd (1.0d0,h,isize,1.0d0,scr5(1,1,3))
***   now:  h = E0+E1  (position-space)
        if (.not.no_prop) call trsmt(scr8(1,3),revt,x,nbas,scr5(1,1,1),
     *                               scr5(1,1,2))
***   now:  x        = X0 (position space)  --- if no_prop=.false.
c
      endif
c
      return
      end
c
c
c
c
      subroutine calc_operators (nbas,isize,dkhorder,xorder,dkhscfflg,
     *                       posu,post,poss,vv,nn,dd,yy,ff,gg,pp,xx,
     *                       ii,jj,kk,ll,mm,e,snumber,tnumber,unumber,
     *                       scrno1,scrno2,scr1,scr2,scr3,
     *                       scr5,scr8,revt,h,nbasp,nbaso,
     *                       dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Calculate Hamiltonian up to desired order
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 13.03.2007 (MR, ETH Zurich)
c           * matrices can be read from disk
c           * scr4 and scr6 were eliminated from this routine and have been
c             replaced by scr5 and scr8, resp.;
c             hence, scrno1 has to be at least 3 (scrno1=3)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer nbas,isize,dkhorder,xorder,posu(maxunumber),nbasp,nbaso,
     *        post(maxsnumber),poss(maxsnumber),snumber,tnumber,unumber,
     *        scrno1,scrno2
c
      real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
c
      real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 scr5(nbas,nbas,scrno1),
     *                 scr8(isize,scrno2),revt(nbas,nbas),h(isize)
c
      integer length,order(3),termcounter
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
      character*(maxlength) scrtxt,term
CMR      character*(2) dkh_int2char2
      real*8 coeff
c
      integer i,j,k,idum1
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit1=11
      dkhunit1=IsFreeUnit(dkhunit1)
      Call Molcas_Open(dkhunit1,'dkhops.11')
#else
      open (dkhunit1, file='dkhops.11', status='OLD',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit1)
1009  read (dkhunit1,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1009
c
CMR      if (dkhorder.gt.1 .and. dkhorder.le.9) then
CMR        write (stdout,1010) dkhorder
CMR1010    format (/,'SR calc_operators:  Start evaluation of the DKH',I1,
CMR     *          ' Hamiltonian.',/2X)
CMR      endif
CMR      if (dkhorder.ge.10) then
CMR        write (stdout,1011) dkhorder
CMR1011    format (/,'SR calc_operators:  Start evaluation of the DKH',I2,
CMR     *          ' Hamiltonian.',/2X)
CMR      endif
CMR      if (dkhorder.le.1) then
CMR        write (stdout,1012)
CMR1012    format (/,'SR calc_operators:  Nothing to do for DKH0 or DKH1.',
CMR     *          /2X)
CMR      endif
      do 10 i=2,dkhorder
CMR        write(*,*) "dkhorder=",i
c
1022    read (dkhunit1,'(A3)') scrtxt(1:3)
        if (scrtxt(1:3).ne.'***') goto 1022
        call mat_zero (scr5(1,1,1),nbas)
        do 20 j=1,maxlength
          scrtxt(j:j)=' '
  20    continue
        idum1=0
        do 23 j=1,3
          order(j)=0
  23    continue
        read (dkhunit1,'(5X,I2,19X,I6)') order(3),termcounter
        do 30 j=1,termcounter
          do 40 k=1,maxlength
            term(k:k)=' '
  40      continue
          idum1=0
c
          read (dkhunit1,
     *         '(I7,1X,I3,4X,I2,1X,I2,1X,I2,1X,A90,4X,F17.14)')
     *          idum1,length,order(1),order(2),order(3),
     *                       term(1:90),coeff
          do 50 k=1,length
            term(k:k)=term(90-length+k:90-length+k)
            term(90-length+k:90-length+k)=' '
  50      continue
CMR       write(*,*) "termcounter,term(1:length)=",j,term(1:length)
          call evalstring (length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                     dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                     snumber,tnumber,unumber,scrno1,scr1,
     *                     scr2,scr3,scr5,nbasp,nbaso,dkhadr,adrmem,
     *                     dkh_48,adrnext)
          if (.not.(dkhscfflg .and. order(2).gt.xorder))
     *      call mat_plainadd (scr5(1,1,1),nbas,scr5(1,1,scrno1))
c
  30    continue
CMR      write(*,*) "calc_operators2: scr5(1,1,1)=",scr5(1,1,1)
        call mat_triang (scr8(1,1),nbas,scr5(1,1,1))
        call TrSmt(scr8(1,1),revt,scr8(1,2),nbas,scr5(1,1,scrno1-1),
     *             scr5(1,1,scrno1))
        call mat_tadd (1.0d0,h,isize,1.0d0,scr8(1,2))
CMR        write (stdout,1060) dkh_int2char2(i)
CMR1060    format (15X,'All terms belonging to E',A2,' have been ',
CMR     *          'evaluated.')
c
  10  continue
c
#ifndef MOLPRO
      close (dkhunit1)
#endif
c
      return
#ifdef _WARNING_WORKAROUND_
      if (.false.) call Unused_integer(idum1)
#endif
      end
c
c
      subroutine calc_orders (dkhscfflg,length,orderv,orderx,ordertot,
     *                        term,sorder,uorder)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 18.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer length,orderv,orderx, ordertot,sorder(maxsnumber,3),
     *        uorder(maxunumber,3)
      character*(maxlength) term
c
      integer poss,post,posw,i,idum,dkh_char2int
c
      posw=index(term(1:length),'W')
      if (posw.gt.0) then
        write (stdout,1002)
 1002   format (/2X,'ERROR in SR "calc_orders": W-operator occurs in ',
     *          'expression. This is not possible at this stage.',//2X,
     *          'STOP.',/2X)
        Call Abend
      endif
      poss=index(term(1:length),'S')
      post=index(term(1:length),'T')
c
      if (dkhscfflg) then
        if (poss.eq.0 .and. post.eq.0) then
          orderv=0
          orderx=0
          do 10 i=1,length
            if (term(i:i).eq.'V'.or.term(i:i).eq.'N'.or.term(i:i).eq.'D'
     *          .or. term(i:i).eq.'Y'.or.term(i:i).eq.'F'
     *          .or. term(i:i).eq.'G')  orderv=orderv+1
            if (term(i:i).eq.'X'.or.term(i:i).eq.'I'.or.term(i:i).eq.'J'
     *          .or. term(i:i).eq.'K'.or.term(i:i).eq.'L'
     *          .or. term(i:i).eq.'M')  orderx=orderx+1
            if (term(i:i).eq.'U') then
              idum=dkh_char2int(3,term(i+1:i+3))
              orderv=orderv+uorder(idum,1)
              orderx=orderx+uorder(idum,2)
            endif
  10      continue
        else
          orderv=0
          orderx=0
        endif
c
      else
c
        orderv=0
        orderx=0
        do 20 i=1,length
          if (term(i:i).eq.'V'.or.term(i:i).eq.'N'.or.term(i:i).eq.'D'
     *        .or. term(i:i).eq.'Y'.or.term(i:i).eq.'F'
     *        .or. term(i:i).eq.'G')  orderv=orderv+1
          if (term(i:i).eq.'X'.or.term(i:i).eq.'I'.or.term(i:i).eq.'J'
     *        .or. term(i:i).eq.'K'.or.term(i:i).eq.'L'
     *        .or. term(i:i).eq.'M')  orderx=orderx+1
          if (term(i:i).eq.'S' .or. term(i:i).eq.'T') then
            idum=dkh_char2int(3,term(i+1:i+3))
            orderv=orderv+sorder(idum,1)
            orderx=orderx+sorder(idum,2)
          endif
          if (term(i:i).eq.'U') then
            idum=dkh_char2int(3,term(i+1:i+3))
            orderv=orderv+uorder(idum,1)
            orderx=orderx+uorder(idum,2)
          endif
  20    continue
c
      endif
        if (.not.dkhscfflg .and. (orderv+orderx .ne. ordertot)) then
          write (stdout,1023) orderv,orderx,ordertot
1023      format (2X,'ERROR in SR "calc_orders":  orderv = ',I2,',',3X,
     *            'orderx = ',I2,',',5X,'ordertot = ',I2,//2X,'STOP.',
     *            //2X)
          CALL Abend
        endif
c
      return
      end
c
      subroutine calc_prefactors (nbas,isize,clight,aa,rr,tt,pp,e,ew)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 25.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer nbas,isize
      real*8 clight,aa(nbas),rr(nbas),tt(nbas),pp(nbas),
     *                 e(nbas),ew(nbas)
      real*8 cl2,invcl2,ratio,tv1,tv2,tv3,tv4
      integer i
c
      cl2=clight*clight
      invcl2=1d0/(clight*clight)
c
      do 10 i=1,nbas
        if (ew(i).lt.0.0D0) then
          write(stdout,1000) i,ew(i)
1000      format (2X,'ERROR in SR "calc_prefactors":  ew(',I4,') = ',
     *            F16.8,' is less than zero.',//2X,'STOP.',/)
          CALL Abend
        endif
c
        ratio=ew(i)/clight
        tt(i)=ew(i)
c
        if (ratio.gt.0.02d0) then
          ew(i)=cl2*(sqrt(1.0D0+(invcl2+invcl2)*ew(i))-1.0d0)
        else
          tv1=ew(i)
          tv2=-tv1*ew(i)*invcl2/2.D0
          tv3=-tv2*ew(i)*invcl2
          tv4=-tv3*ew(i)*invcl2*1.25D0
          ew(i)=tv1+tv2+tv3+tv4
        endif
        e(i)=ew(i)+cl2
c
        aa(i)=sqrt((cl2+e(i))/(2.0D0*e(i)))
        rr(i)=clight/(cl2+e(i))
        pp(i)=2.0D0*tt(i)*rr(i)*rr(i)
c
  10  continue
c
c
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(isize)
      end
c
c
      subroutine calc_revt (nbas,revt,tran,sinv,scrmat1,scrmat2)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 25.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer nbas
      real*8 revt(nbas,nbas),tran(nbas,nbas),sinv(nbas,nbas),
     *                 scrmat1(nbas,nbas),scrmat2(nbas,nbas)
      integer i,j,k
      call mat_zero (scrmat2,nbas)
      call mat_zero (revt,nbas)
c
      do 10 i=1,nbas
        do 20 j=1,nbas
          do 30 k=i,nbas
            scrmat2(i,j)=scrmat2(i,j)+sinv(i,k)*tran(k,j)
  30       continue
  20     continue
  10   continue
c
      do 40 i=1,nbas
        do 50 j=1,nbas
          do 60 k=1,nbas
            revt(i,j)=revt(i,j)+scrmat1(i,k)*scrmat2(k,j)
  60      continue
  50    continue
  40  continue
c
      return
      end
c
c
      subroutine calc_smult (jnumber,smult,termcounter,termleng,term)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 18.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer jnumber,smult,termcounter,termleng(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
      integer term(*)
      character*(maxlength) termstr
#else
      character*(maxlength) term(maxoperators)
#endif
c
      integer k,l
c
      smult=0
      do 100 k=1,termcounter
#if defined(_MOLCAS_) || defined(MOLPRO)
        call get_dkoperators_i(k,termstr,term)
#endif
        do 200 l=1,termleng(k)
#if defined(_MOLCAS_) || defined(MOLPRO)
          if ( termstr(l:l).eq.'V' .or. termstr(l:l).eq.'N' .or.
     *         termstr(l:l).eq.'D' .or. termstr(l:l).eq.'Y' .or.
     *         termstr(l:l).eq.'F' .or. termstr(l:l).eq.'G' .or.
     *         termstr(l:l).eq.'Z' .or. termstr(l:l).eq.'Q' .or.
     *         termstr(l:l).eq.'X' .or. termstr(l:l).eq.'I' .or.
     *         termstr(l:l).eq.'J' .or. termstr(l:l).eq.'K' .or.
     *         termstr(l:l).eq.'L' .or. termstr(l:l).eq.'M' .or.
     *         termstr(l:l).eq.'S' .or. termstr(l:l).eq.'T' .or.
     *         termstr(l:l).eq.'U' )   smult=smult+1
#else
          if ( term(k)(l:l).eq.'V' .or. term(k)(l:l).eq.'N' .or.
     *         term(k)(l:l).eq.'D' .or. term(k)(l:l).eq.'Y' .or.
     *         term(k)(l:l).eq.'F' .or. term(k)(l:l).eq.'G' .or.
     *         term(k)(l:l).eq.'Z' .or. term(k)(l:l).eq.'Q' .or.
     *         term(k)(l:l).eq.'X' .or. term(k)(l:l).eq.'I' .or.
     *         term(k)(l:l).eq.'J' .or. term(k)(l:l).eq.'K' .or.
     *         term(k)(l:l).eq.'L' .or. term(k)(l:l).eq.'M' .or.
     *         term(k)(l:l).eq.'S' .or. term(k)(l:l).eq.'T' .or.
     *         term(k)(l:l).eq.'U' )   smult=smult+1
#endif
 200    continue
 100  continue
      smult=smult-termcounter
c
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(jnumber)
      end
c
c
c
      subroutine calc_stimes1 (knumber,wstart,sused,stimes,wstimes,
     *                         ttimes,reslengl,reslengr,rescharl,
     *                         rescharr)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 18.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer knumber,wstart,sused,reslengl,reslengr,stimes(maxsnumber),
     *        m,wstimes(maxuops,maxsnumber),ttimes(maxsnumber)
      character*(maxlength) rescharl,rescharr
      integer dummyleng,idum,istart,iposs,dkh_char2int
c
      if (knumber.gt.1) then
        istart=1
1083    continue
        iposs=0
        iposs=index(rescharl(istart:reslengl),'S')
        if (iposs.gt.0) then
          dummyleng=3
          idum=dkh_char2int(dummyleng,
     *                        rescharl(istart+iposs:istart+iposs+2))
          stimes(idum)=stimes(idum)+1
          istart=istart+iposs+3
          goto 1083
        endif
        istart=1
1084    continue
        iposs=0
        if (istart.le.reslengr) then
          iposs=index(rescharr(istart:reslengr),'S')
          if (iposs.gt.0) then
            dummyleng=3
            idum=dkh_char2int(dummyleng,
     *                          rescharr(istart+iposs:istart+iposs+2))
            stimes(idum)=stimes(idum)+1
            istart=istart+iposs+3
            goto 1084
          endif
        endif
c
      endif
      if (knumber.gt.1) then
        istart=1
1093    continue
        iposs=0
        iposs=index(rescharl(istart:reslengl),'T')
        if (iposs.gt.0) then
          dummyleng=3
          idum=dkh_char2int(dummyleng,
     *                        rescharl(istart+iposs:istart+iposs+2))
          ttimes(idum)=ttimes(idum)+1
          istart=istart+iposs+3
          goto 1093
        endif
        istart=1
1094    continue
        iposs=0
        if (istart.le.reslengr) then
          iposs=index(rescharr(istart:reslengr),'T')
          if (iposs.gt.0) then
            dummyleng=3
            idum=dkh_char2int(dummyleng,
     *                          rescharr(istart+iposs:istart+iposs+2))
            ttimes(idum)=ttimes(idum)+1
            istart=istart+iposs+3
            goto 1094
          endif
        endif
c
      endif
      do 352 m=1,sused
        stimes(m)=stimes(m)+wstimes(knumber+wstart-1,m)
 352  continue
c
      return
      end
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine calc_stimes2 (stimes,ttimes,reslengl,reslengr,
     *                         rescharl,rescharr)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 11.10.2006 (MR, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer reslengl,reslengr,stimes(maxsnumber),
     *        ttimes(maxsnumber)
      character*(maxlength) rescharl,rescharr
      integer dummyleng,idum,istart,iposs,dkh_char2int
c
      istart=1
1083  continue
      iposs=0
      iposs=index(rescharl(istart:reslengl),'S')
      if (iposs.gt.0) then
        dummyleng=3
        idum=dkh_char2int(dummyleng,
     *                      rescharl(istart+iposs:istart+iposs+2))
        stimes(idum)=stimes(idum)+1
        istart=istart+iposs+3
        goto 1083
      endif
      istart=1
1084  continue
      iposs=0
      if (istart.le.reslengr) then
        iposs=index(rescharr(istart:reslengr),'S')
        if (iposs.gt.0) then
          dummyleng=3
          idum=dkh_char2int(dummyleng,
     *                        rescharr(istart+iposs:istart+iposs+2))
          stimes(idum)=stimes(idum)+1
          istart=istart+iposs+3
          goto 1084
        endif
      endif
      istart=1
1093  continue
      iposs=0
      iposs=index(rescharl(istart:reslengl),'T')
      if (iposs.gt.0) then
        dummyleng=3
        idum=dkh_char2int(dummyleng,
     *                      rescharl(istart+iposs:istart+iposs+2))
        ttimes(idum)=ttimes(idum)+1
        istart=istart+iposs+3
        goto 1093
      endif
      istart=1
1094  continue
      iposs=0
      if (istart.le.reslengr) then
        iposs=index(rescharr(istart:reslengr),'T')
        if (iposs.gt.0) then
          dummyleng=3
          idum=dkh_char2int(dummyleng,
     *                        rescharr(istart+iposs:istart+iposs+2))
          ttimes(idum)=ttimes(idum)+1
          istart=istart+iposs+3
          goto 1094
        endif
      endif
c
      return
      end
c
c
      subroutine calc_Sxxx (nbas,dkhorder,xorder,dkhscfflg,posu,post,
     *                      poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,
     *                      mm,e,snumber,tnumber,unumber,scrno1,
     *                      scr1,scr2,scr3,scr5,nbasp,nbaso,
     *                      dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 13.03.2007 (M. Reiher, ETH Zurich)
c            * result matrix can now be written to disk
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      logical dkhscfflg
      integer nbas,dkhorder,xorder,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,scrno1,nbasp,
     *        nbaso
c
      real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
c
      real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 scr5(nbas,nbas,scrno1)
c
      integer length,order(3),termcounter
      character*(maxlength) scrtxt,term
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
#ifdef _DEBUGPRINT_
      character*(3) dkh_int2char
#endif
      real*8 coeff
c
      integer i,j,k,idum1,idum2
c
      integer dkh_48,adrmem,idum3,iscr,adrnext,adr
      integer one,dkhadr(adrmem)
      parameter(one=1)
      idum3=0
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit3=13
      dkhunit3=IsFreeUnit(dkhunit3)
      Call Molcas_Open(dkhunit3,'dkhops.13')
#else
      open (dkhunit3, file='dkhops.13', status='OLD',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit3)
1009  read (dkhunit3,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1009
      read (dkhunit3, '(I3)') idum1
c
CMR      if (snumber.gt.0) then
CMR        write (stdout,1010) snumber
CMR1010    format (/,'SR calc_Sxxx:  Start evaluation of the ',I3,
CMR     *          ' Sxxx matrices.',/2X)
CMR      endif
CMR      if (snumber.eq.0) then
CMR        write (stdout,1011)
CMR1011    format (/,'SR calc_Sxxx:  Nothing to do here.',/2X)
CMR      endif
c
      if (Out_Of_Core)
     *  call GetMem('calc_Un ','ALLO','REAL',iscr,nbas*nbas+4)
c
      do 10 i=1,snumber
c
1022    read (dkhunit3,'(A3)') scrtxt(1:3)
        if (scrtxt(1:3).ne.'***') goto 1022
        do 20 j=1,maxlength
          scrtxt(j:j)=' '
  20    continue
        idum1=0
        do 23 j=1,3
          order(j)=0
  23    continue
        read (dkhunit3,'(A1,I3,1X,A9,1X,I2,1X,I2,1X,I2,2X,I7)')
     *          scrtxt(1:1),idum1,scrtxt(1:9),order(1),order(2),
     *          order(3),termcounter
        poss(idum1)=i
        if (Out_Of_Core) then
          idum3=1000+idum1
          call mat_zero (work(iscr),nbas)
        else
          call mat_zero (scr1(1,1,poss(idum1)),nbas)
        end if
        if (order(3).gt.Max(dkhorder,xorder+1)) then
#ifdef _DEBUGPRINT_
          write (stdout,1040) i,dkh_int2char(idum1)
1040      format (15X,I3,2X,'S',A3,' is not required here; thus skip ',
     *            'it.')
#endif
          goto 10
        endif
        do 30 j=1,termcounter
          do 40 k=1,maxlength
            term(k:k)=' '
  40      continue
          idum2=0
          do 43 k=1,3
            order(k)=0
  43      continue
c
          read (dkhunit3,
     *         '(I7,2X,I3,3X,I2,1X,I2,1X,I2,1X,A90,4X,F17.14)')
     *          idum2,length,order(1),order(2),order(3),
     *                      term(1:90),coeff
          do 50 k=1,length
            term(k:k)=term(90-length+k:90-length+k)
            term(90-length+k:90-length+k)=' '
  50      continue
          call evalstring (length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                     dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                     snumber,tnumber,unumber,scrno1,scr1,
     *                     scr2,scr3,scr5,nbasp,nbaso,dkhadr,adrmem,
     *                     dkh_48,adrnext)
          if (.not.(dkhscfflg .and. order(2).gt.xorder)) then
            if (Out_Of_Core) then
              call mat_plainadd (work(iscr),nbas,scr5(1,1,scrno1))
            else
              call mat_plainadd (scr1(1,1,poss(idum1)),nbas,
     *                          scr5(1,1,scrno1))
            end if
          end if
c
  30    continue
c
        if (Out_Of_Core) then
          if (i.eq.1) adr=adrnext
          dkhadr(idum3)=adr
          call ddafile(dkh_48,one,work(iscr),nbas*nbas,adr)
          if (i.eq.snumber) adrnext=adr
CMR        write(*,*) "S : scr1(1,1,poss(idum1))=",work(iscr)
        else
CMR        write(*,*) "S : scr1(1,1,poss(idum1))=",scr1(1,1,poss(idum1))
        end if
c
#ifdef _DEBUGPRINT_
        write (stdout,1060) i,dkh_int2char(idum1)
1060    format (15X,I3,2X,'S',A3,' has now been completed.')
#endif
c
  10  continue
c
      if (Out_Of_Core)
     *  call GetMem('calc_Un ','FREE','REAL',iscr,nbas*nbas+4)
c
#ifndef MOLPRO
      close (dkhunit3)
#endif
c
CMR      write (stdout,1090)
CMR1090  format (/15X,'All Sxxx matrices have been stored in scr1().')
c
      return
#if _WARNING_WORKAROUND_
      if (.false.) call Unused_integer(idum2)
#endif
      end
c
c
      subroutine calc_Txxx (nbas,dkhorder,xorder,dkhscfflg,posu,post,
     *                      poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,
     *                      mm,e,snumber,tnumber,unumber,scrno1,
     *                      scr1,scr2,scr3,scr5,nbasp,nbaso,
     *                      dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 13.03.2007 (M. Reiher, ETH Zurich)
c            * result matrix can now be written to disk
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      logical dkhscfflg
      integer nbas,dkhorder,xorder,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,scrno1,nbasp,
     *        nbaso
c
      real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
c
      real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 scr5(nbas,nbas,scrno1)
c
      integer length,order(3),termcounter
      character*(maxlength) scrtxt,term
CMR      character*(3) dkh_int2char
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
      real*8 coeff
c
      integer i,j,k,idum1,idum2
c
      integer dkh_48,adrmem,idum3,iscr,adrnext,adr
      integer one,dkhadr(adrmem)
      parameter(one=1)
      idum3=0
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit4=14
      dkhunit4=IsFreeUnit(dkhunit4)
      Call Molcas_Open(dkhunit4,'dkhops.14')
#else
      open (dkhunit4, file='dkhops.14', status='OLD',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit4)
1009  read (dkhunit4,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1009
      read (dkhunit4, '(I3)') idum1
c
CMR      if (tnumber.gt.0) then
CMR        write (stdout,1010) tnumber
CMR1010    format (/,'SR calc_Txxx:  Start evaluation of the ',I3,
CMR     *          ' Txxx matrices.',/2X)
CMR      endif
CMR      if (tnumber.eq.0) then
CMR        write (stdout,1011)
CMR1011    format (/,'SR calc_Txxx:  Nothing to do here.')
CMR      endif
c
      if (Out_Of_Core)
     *  call GetMem('calc_Un ','ALLO','REAL',iscr,nbas*nbas+4)
c
      do 10 i=1,tnumber
c
1022    read (dkhunit4,'(A3)') scrtxt(1:3)
        if (scrtxt(1:3).ne.'***') goto 1022
        do 20 j=1,maxlength
          scrtxt(j:j)=' '
  20    continue
        idum1=0
        do 23 j=1,3
          order(j)=0
  23    continue
        read (dkhunit4,'(A1,I3,1X,A11,1X,I2,1X,I2,1X,I2,2X,I7)')
     *          scrtxt(1:1),idum1,scrtxt(1:11),order(1),order(2),
     *          order(3),termcounter
        post(idum1)=i
        if (Out_Of_Core) then
          idum3=2000+idum1
          call mat_zero (work(iscr),nbas)
        else
          call mat_zero (scr2(1,1,post(idum1)),nbas)
        end if
        if (order(3).gt.Max(dkhorder,xorder+1)) then
CMR          write (stdout,1040) i,dkh_int2char(idum1)
CMR1040      format (15X,I3,2X,'T',A3,' is not required here; thus skip ',
CMR     *            'it.')
          goto 10
        endif
        do 30 j=1,termcounter
          do 40 k=1,maxlength
            term(k:k)=' '
  40      continue
          idum2=0
          do 43 k=1,3
            order(k)=0
  43      continue
c
          read (dkhunit4,
     *         '(I7,2X,I3,5X,I2,1X,I2,1X,I2,1X,A90,4X,F17.14)')
     *          idum2,length,order(1),order(2),order(3),
     *                       term(1:90),coeff
          do 50 k=1,length
            term(k:k)=term(90-length+k:90-length+k)
            term(90-length+k:90-length+k)=' '
  50      continue
          call evalstring (length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                     dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                     snumber,tnumber,unumber,scrno1,scr1,
     *                     scr2,scr3,scr5,nbasp,nbaso,dkhadr,adrmem,
     *                     dkh_48,adrnext)
          if (.not.(dkhscfflg .and. order(2).gt.xorder)) then
            if (Out_Of_Core) then
              call mat_plainadd (work(iscr),nbas,scr5(1,1,scrno1))
            else
              call mat_plainadd (scr2(1,1,post(idum1)),nbas,
     *                           scr5(1,1,scrno1))
            end if
          end if
c
  30    continue
c
        if (Out_Of_Core) then
          if (i.eq.1) adr=adrnext
          dkhadr(idum3)=adr
          call ddafile(dkh_48,one,work(iscr),nbas*nbas,adr)
          if (i.eq.tnumber) adrnext=adr
        end if
c
CMR        write (stdout,1060) i,dkh_int2char(idum1)
CMR1060    format (15X,I3,2X,'T',A3,' has now been completed.')
c
  10  continue
c
      if (Out_Of_Core)
     *  call GetMem('calc_Un ','FREE','REAL',iscr,nbas*nbas+4)
c
#ifndef MOLPRO
      close (dkhunit4)
#endif
c
CMR      if (tnumber.gt.0) write (stdout,1090)
CMR1090  format (/15X,'All Txxx matrices have been stored in scr2().')
c
      return
#if _WARNING_WORKAROUND_
      if (.false.) call Unused_integer(idum2)
#endif
      end
c
c
      subroutine calc_Uxxx (nbas,dkhorder,xorder,dkhscfflg,posu,post,
     *                      poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,
     *                      mm,e,snumber,tnumber,unumber,scrno1,
     *                      scr1,scr2,scr3,scr5,nbasp,nbaso,
     *                      dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 13.03.2007 (M. Reiher, ETH Zurich)
c            * result matrix may be written to disk
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer nbas,dkhorder,xorder,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,scrno1,nbasp,
     *        nbaso
c
      real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
c
      real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 scr5(nbas,nbas,scrno1)
c
      character*(maxlength) scrtxt
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
      integer length,i,j,idum1,order(3)
      real*8 coeff
c
      integer dkh_48,adrmem,idum2,adrnext,adr
      integer one,dkhadr(adrmem)
      parameter(one=1)
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit5=15
      dkhunit5=IsFreeUnit(dkhunit5)
      Call Molcas_Open(dkhunit5,'dkhops.15')
#else
      open (dkhunit5, file='dkhops.15', status='OLD',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit5)
1009  read (dkhunit5,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1009
      read (dkhunit5, '(I3)') idum1
c
      do 10 i=1,unumber
        do 20 j=1,maxlength
          scrtxt(j:j)=' '
  20    continue
        length=0
        idum1=0
        do 23 j=1,3
          order(j)=0
  23    continue
        read (dkhunit5,'(I1,3X,A1,I3,3X,A3,8X,I2,8X,I2,9X,I2)')
     *          length,scrtxt(1:1),idum1,scrtxt(1:3),order(1),order(2),
     *          order(3)
        posu(idum1)=i
        coeff=1.0d0
        call evalstring (length,scrtxt,coeff,nbas,posu,post,poss,vv,nn,
     *                   dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,snumber,
     *                   tnumber,unumber,scrno1,scr1,scr2,scr3,scr5,
     *                   nbasp,nbaso,dkhadr,adrmem,dkh_48,adrnext)
        if (dkhscfflg .and. order(2).gt.xorder) then
          call mat_zero (scr5(1,1,scrno1),nbas)
        else
          if (Out_Of_Core) then
            if (i.eq.1) adr=adrnext
            idum2=idum1+3000
            dkhadr(idum2)=adr
            call ddafile(dkh_48,one,scr5(1,1,scrno1),nbas*nbas,
     *                   adr)
            if (i.eq.unumber) adrnext=adr
          else
            call mat_copy (scr3(1,1,posu(idum1)),nbas,nbas,
     *                     scr5(1,1,scrno1))
          endif
CMR       write(*,*) "U : scr5(1,1,scrno1)=",scr5(1,1,scrno1)
        endif
c
  10  continue
CMR      if (unumber.gt.0) then
CMR        write (stdout,1072) unumber
CMR1072    format (/,'SR calc_Uxxx:  All ',I3,' Uxxx matrices ',
CMR     *          'have successfully been evaluated.',/15X,'They have ',
CMR     *          'been stored in scr3().')
CMR      endif
CMR      if (unumber.eq.0) then
CMR        write (stdout,1073)
CMR1073    format (/,'SR calc_Uxxx:  Nothing to do here.')
CMR      endif
c
#ifndef MOLPRO
      close (dkhunit5)
#endif
c
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dkhorder)
      end
c
c
      subroutine calc_xoperators (nbas,isize,dkhorder,xorder,dkhscfflg,
     *                      posu,post,poss,vv,nn,dd,yy,ff,gg,pp,xx,
     *                      ii,jj,kk,ll,mm,e,snumber,tnumber,unumber,
     *                      scrno1,scrno2,scr1,scr2,scr3,
     *                      scr5,scr8,revt,x,nbasp,nbaso,
     *                      dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 13.03.2007 (MR, ETH Zurich)
c             * scr4 and scr7 were eliminated and replaced by scr5 and
c               scr7, resp.; hence,
c               scrno1 now must be at least 3 (scrno1=3)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer nbas,isize,dkhorder,xorder,posu(maxunumber),nbasp,nbaso,
     *        post(maxsnumber),poss(maxsnumber),snumber,tnumber,unumber,
     *        scrno1,scrno2
c
      real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
c
      real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),
     *                 scr5(nbas,nbas,scrno1),
     *                 scr8(isize,scrno2),revt(nbas,nbas),x(isize)
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
      integer length,order,termcounter
      character*(maxlength) scrtxt,term
#ifdef _MOLCAS_
      Integer IsFreeUnit
#endif
CMR      character*(2) dkh_int2char2
      real*8 coeff
c
      integer i,j,k,idum1
c
#ifndef MOLPRO
#ifdef _MOLCAS_
      dkhunit2=12
      dkhunit2=IsFreeUnit(dkhunit2)
      Call Molcas_Open(dkhunit2,'dkhops.12')
#else
      open (dkhunit2, file='dkhops.12', status='OLD',
     *        form='FORMATTED')
#endif
#endif
      rewind(dkhunit2)
1009  read (dkhunit2,'(A3)') scrtxt(1:3)
      if (scrtxt(1:3).ne.'+++') goto 1009
c
CMR      if (xorder.ge.1 .and. xorder.le.9) then
CMR        write (stdout,1010) xorder
CMR1010    format (/,'SR calc_xoperators:  Start evaluation of the DKH',I1,
CMR     *          '-transformed property X.',/2X)
CMR      endif
CMR      if (xorder.ge.10) then
CMR        write (stdout,1011) xorder
CMR1011    format (/,'SR calc_xoperators:  Start evaluation of the DKH',I2,
CMR     *          '-transformed property X.',/2X)
CMR      endif
CMR      if (xorder.le.0) then
CMR        write (stdout,1012)
CMR1012    format (/,'SR calc_xoperators:  Nothing to do for X(,0).',
CMR     *          /2X)
CMR      endif
c
      do 10 i=1,xorder
c
1022    read (dkhunit2,'(A3)') scrtxt(1:3)
        if (scrtxt(1:3).ne.'***') goto 1022
        call mat_zero (scr5(1,1,1),nbas)
        do 20 j=1,maxlength
          scrtxt(j:j)=' '
  20    continue
        idum1=0
        order=0
        read (dkhunit2,'(5X,I2,17X,I6)') order,termcounter
        do 30 j=1,termcounter
          do 40 k=1,maxlength
            term(k:k)=' '
  40      continue
          idum1=0
          order=0
c
          read (dkhunit2,'(I7,1X,I3,5X,I2,5X,A90,4X,F17.14)')
     *            idum1,length,order,term(1:90),coeff
          do 50 k=1,length
            term(k:k)=term(90-length+k:90-length+k)
            term(90-length+k:90-length+k)=' '
  50      continue
          call evalstring (length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                     dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                     snumber,tnumber,unumber,scrno1,scr1,scr2,
     *                     scr3,scr5,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *                     adrnext)
          call mat_plainadd (scr5(1,1,1),nbas,scr5(1,1,scrno1))
c
  30    continue
        call mat_triang (scr8(1,1),nbas,scr5(1,1,1))
        call TrSmt(scr8(1,1),revt,scr8(1,2),nbas,scr5(1,1,scrno1-1),
     *                                               scr5(1,1,scrno1))
        call mat_tadd (1.0d0,x(1),isize,1.0d0,scr8(1,2))
CMR        write (stdout,1060) dkh_int2char2(i)
CMR1060    format (15X,'All terms belonging to X',A2,' have been ',
CMR     *          'evaluated.')
c
  10  continue
c
#ifndef MOLPRO
      close (dkhunit2)
#endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(dkhorder)
        call Unused_logical(dkhscfflg)
#ifdef _WARNING_WORKAROUND_
        call Unused_integer(idum1)
        call Unused_integer(order)
#endif
      end if
      end
c
c
      subroutine concatenate (operleng,operator,reslengl,rescharl,
     *                        termleng,term,reslengr,rescharr)
c
c***************************************************************************
c
c   Concatenate:  operator = rescharl // term // rescharr
c
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 12.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer operleng,reslengl,termleng,reslengr
      character*(maxlength) operator,rescharl,rescharr
      character*(*) term
c
c  Consistency check
      if (operleng.ne.reslengl+termleng+reslengr) then
        write (stdout,3000) operleng,reslengl,termleng,reslengr
3000    format (/2X,'ERROR in subroutine concatenate: operleng =',
     *          I4,' is not equal to reslengl+termleng+reslengr (=',
     *          I3,'+',I3,'+',I3,') as it should.',//2X,'STOP.',/)
        CALL Abend
      endif
c
#if defined(_MOLCAS_) || defined(MOLPRO)
      if (reslengl.ne.0.and.reslengr.ne.0)
     *    operator=rescharl(1:reslengl)//term(1:termleng)
     *       //rescharr(1:reslengr)
      if (reslengl.ne.0.and.reslengr.eq.0)
     *    operator=rescharl(1:reslengl)//term(1:termleng)
      if (reslengl.eq.0.and.reslengr.ne.0)
     *    operator=term(1:termleng)//rescharr(1:reslengr)
      if (reslengl.eq.0.and.reslengr.eq.0)
     *    operator=term(1:termleng)
#else
      if (reslengl.ne.0.and.reslengr.ne.0)
     *    operator(1:operleng)=rescharl(1:reslengl)//term(1:termleng)
     *       //rescharr(1:reslengr)
      if (reslengl.ne.0.and.reslengr.eq.0)
     *    operator(1:operleng)=rescharl(1:reslengl)//term(1:termleng)
      if (reslengl.eq.0.and.reslengr.ne.0)
     *    operator(1:operleng)=term(1:termleng)//rescharr(1:reslengr)
      if (reslengl.eq.0.and.reslengr.eq.0)
     *    operator(1:operleng)=term(1:termleng)
#endif
c
      return
      end
c
c
      subroutine determine_factor (length,term,iact,nbas,posu,post,poss,
     *                       vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                       snumber,tnumber,unumber,scr1,scr2,scr3,
     *                       factor,nbasp,nbaso,dkhadr,adrmem,dkh_48,
     *                       adrnext)
c
c****************************************************************************
c
c   Determine which factor stands at position 'iact' of 'term'.
c   Store its matrix representation in 'factor'.
c
c   Adjust value of 'iact', such that it points at beginning of next factor.
c
c   Note: There are no brackets [] occurring in term.
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.2.0
c
c   last modified: 13.03.2007 (M. Reiher, ETH Zurich)
c                  * reading from disk was included
c                  * a couple of loops were re-arranged for speed up
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c****************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,iact,nbas,posu(maxunumber),post(maxsnumber),
     *        poss(maxsnumber),snumber,tnumber,unumber,nbasp,nbaso
      character*(maxlength) term
c
      real*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp)
      real*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),factor(nbas,nbas)
c
      integer idum1,dkh_char2int,i,j
c
      integer dkh_48,adrmem,adr,iscr,isize,adrnext
      integer two,dkhadr(adrmem)
      logical flag
      parameter(two=2)
      flag=.true.
      isize=(nbas*(nbas+1)/2)
c
c------------------------------------------------------------------------
c
      if (Out_Of_Core) then
       if (term(iact:iact).eq.'S') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=1000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'T') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=2000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'U') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         idum1=3000+idum1
         iact=iact+4
         adr=dkhadr(idum1)
         isize=nbas*nbas
       else if (term(iact:iact).eq.'V') then
         iact=iact+1
         adr=dkhadr(1)
       else if (term(iact:iact).eq.'N') then
         iact=iact+1
         adr=dkhadr(5)
       else if (term(iact:iact).eq.'D') then
         iact=iact+1
         adr=dkhadr(2)
       else if (term(iact:iact).eq.'Y') then
         iact=iact+1
         adr=dkhadr(6)
       else if (term(iact:iact).eq.'F') then
         iact=iact+1
         adr=dkhadr(7)
       else if (term(iact:iact).eq.'G') then
         iact=iact+1
         adr=dkhadr(8)
       else if (term(iact:iact).eq.'Z') then
         call mat_sq_from_d (factor,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'Q') then
         call mat_sq_dev_d (factor,nbas,pp)
         iact=iact+1
         flag=.false.
       else if (term(iact:iact).eq.'X') then
         iact=iact+1
         adr=dkhadr(3)
       else if (term(iact:iact).eq.'I') then
         iact=iact+1
         adr=dkhadr(9)
       else if (term(iact:iact).eq.'J') then
         iact=iact+1
         adr=dkhadr(4)
       else if (term(iact:iact).eq.'K') then
         iact=iact+1
         adr=dkhadr(10)
       else if (term(iact:iact).eq.'L') then
         iact=iact+1
         adr=dkhadr(11)
       else if (term(iact:iact).eq.'M') then
         iact=iact+1
         adr=dkhadr(12)
       end if
       if (flag) then
        if (isize.eq.nbas*nbas) then
         call ddafile(dkh_48,two,factor,nbas*nbas,adr)
        else
CMR it would be advantageous to have a triangular scratch
CMR   matrix permanently available
         call GetMem('DetFac  ','ALLO','REAL',iscr,isize+4)
         call ddafile(dkh_48,two,work(iscr),isize,adr)
         call mat_sq_from_t(factor,nbas,work(iscr))
         call GetMem('DetFac  ','FREE','REAL',iscr,isize+4)
        end if
       end if
      else
       if (term(iact:iact).eq.'S') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         do 10 i=1,nbas
           do 20 j=1,nbas
             factor(j,i)=scr1(j,i,poss(idum1))
  20       continue
  10     continue
         iact=iact+4
       else if (term(iact:iact).eq.'T') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         do 30 i=1,nbas
           do 40 j=1,nbas
             factor(j,i)=scr2(j,i,post(idum1))
  40       continue
  30     continue
         iact=iact+4
       else if (term(iact:iact).eq.'U') then
         idum1=dkh_char2int(3,term(iact+1:iact+3))
         do 50 i=1,nbas
           do 60 j=1,nbas
             factor(j,i)=scr3(j,i,posu(idum1))
  60       continue
  50     continue
         iact=iact+4
       else if (term(iact:iact).eq.'V') then
         call mat_copy (factor,nbas,nbas,vv)
         iact=iact+1
       else if (term(iact:iact).eq.'N') then
         call mat_copy (factor,nbas,nbas,nn)
         iact=iact+1
       else if (term(iact:iact).eq.'D') then
         call mat_copy (factor,nbas,nbas,dd)
         iact=iact+1
       else if (term(iact:iact).eq.'Y') then
         call mat_copy (factor,nbas,nbas,yy)
         iact=iact+1
       else if (term(iact:iact).eq.'F') then
         call mat_copy (factor,nbas,nbas,ff)
         iact=iact+1
       else if (term(iact:iact).eq.'G') then
         call mat_copy (factor,nbas,nbas,gg)
         iact=iact+1
       else if (term(iact:iact).eq.'Z') then
         call mat_sq_from_d (factor,nbas,pp)
         iact=iact+1
       else if (term(iact:iact).eq.'Q') then
         call mat_sq_dev_d (factor,nbas,pp)
         iact=iact+1
       else if (term(iact:iact).eq.'X') then
CMR         if (nbasp.eq.1) stop "PROGRAM ERROR"
         call mat_copy (factor,nbasp,nbasp,xx)
         iact=iact+1
       else if (term(iact:iact).eq.'I') then
CMR         if (nbasp.eq.1) stop "PROGRAM ERROR"
         call mat_copy (factor,nbasp,nbasp,ii)
         iact=iact+1
       else if (term(iact:iact).eq.'J') then
CMR         if (nbasp.eq.1) stop "PROGRAM ERROR"
         call mat_copy (factor,nbasp,nbasp,jj)
         iact=iact+1
       else if (term(iact:iact).eq.'K') then
CMR         if (nbasp.eq.1) stop "PROGRAM ERROR"
         call mat_copy (factor,nbasp,nbasp,kk)
         iact=iact+1
       else if (term(iact:iact).eq.'L') then
CMR         if (nbasp.eq.1) stop "PROGRAM ERROR"
         call mat_copy (factor,nbasp,nbasp,ll)
         iact=iact+1
       else if (term(iact:iact).eq.'M') then
CMR         if (nbasp.eq.1) stop "PROGRAM ERROR"
         call mat_copy (factor,nbasp,nbasp,mm)
         iact=iact+1
       else
         write (stdout,1083)
1083        format (2X,'ERROR in determine_factor(): could not ',
     *        'determine factor!'//2X,'STOP.',/2X)
         write(stdout,*) "term(iact:iact)=",term(iact:iact)
            CALL Abend
        endif
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(length)
        call Unused_integer(adrnext)
      end if
      end
c
c
c***********************************************************************
c
c   'dkhmatut.f' :  SRs for matrix operations for DKHPACK
c   -----------------------------------------------------
c
c   Contains all subroutines (SR) for matrix operations needed by DKHPACK.
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_mulP (p,nbas,d,q)
c
c***********************************************************************
c
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),d(nbas),q(nbas,nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          q(i,j)=p(i,j)*d(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_mulQ (p,nbas,d,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),d(nbas),q(nbas,nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          q(i,j)=p(i,j)/d(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_axa (p,nbas,aa)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),aa(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=p(i,j)*aa(i)*aa(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_axa_tri (p,nbas,aa)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j,ii
      real*8 p(nbas*(nbas+1)/2),aa(nbas)
c
      ii=0
      do 10 i=1,nbas
        do 20 j=1,i-1
          ii=ii+1
          p(ii)=p(ii)*aa(i)*aa(j)
  20    continue
        ii=ii+1
        p(ii)=p(ii)*aa(i)*aa(i)
  10  continue
c
      return
      end

c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_arxra (p,nbas,aa,rr)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),aa(nbas),rr(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=p(i,j)*aa(i)*aa(j)*rr(i)*rr(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_arxra_tri (p,nbas,aa,rr)

c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j,ii
      real*8 p(nbas*(nbas+1)/2),aa(nbas),rr(nbas)
c
      ii=0
      do 10 i=1,nbas
        do 20 j=1,i-1
          ii=ii+1
          p(ii)=p(ii)*aa(i)*aa(j)*rr(i)*rr(j)
  20    continue
        ii=ii+1
        p(ii)=p(ii)*aa(i)*aa(i)*rr(i)*rr(i)
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_copy (p,nbas1,nbas2,q)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas1,nbas2,i,j
      real*8 p(nbas1,nbas2),q(nbas1,nbas2)
c
      do 10 i=1,nbas2
        do 20 j=1,nbas1
          p(j,i)=q(j,i)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_copy_c (p,nbas1,nbas2,q,c)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  2.1.0
c
c   first version: 14.01.2007  (Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas1,nbas2,i,j
      real*8 p(nbas1,nbas2),q(nbas1,nbas2),c
c
      if (c.ne.1.0d0) then
        do 10 i=1,nbas2
          do 20 j=1,nbas1
            p(j,i)=q(j,i)*c
  20      continue
  10    continue
      else
        do 30 i=1,nbas2
          do 40 j=1,nbas1
            p(j,i)=q(j,i)
  40      continue
  30    continue
      end if
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_copy2 (p,isize,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer isize,i
      real*8 p(isize),q(isize)
c
      do 10 i=1,isize
        p(i)=q(i)
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_copy3 (p,nbas,nbasmax,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,nbasmax,i,j,ij
      real*8 p(nbasmax,nbasmax),q(nbasmax*nbasmax)
c
      call mat_zero (p,nbasmax)
      ij=0
      do 10 i=1,nbas
        do 20 j=1,nbas
          ij=ij+1
          p(j,i)=q(ij)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_zero (p,nbas)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas)
c
      do 10 i=1,nbas
        do 20 j=1,nbas
          p(j,i)=0.0D0
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_zero2 (p,isize)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer isize,i
      real*8 p(isize)
c
      do 10 i=1,isize
        p(i)=0.0D0
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_triang (p,nbas,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j,ii
      real*8 p(nbas*(nbas+1)/2),q(nbas,nbas)
c
      ii=0
      do 10 i=1,nbas
        do 20 j=1,i
          ii=ii+1
          p(ii)=q(j,i)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_sq_from_t (p,nbas,q)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer ii,i,j,nbas
      real*8 p(nbas,nbas),q(nbas*(nbas+1)/2)
c
      ii=0
      do 10 i=1,nbas
        do 20 j=1,i-1
          ii=ii+1
          p(i,j)=q(ii)
          p(j,i)=q(ii)
  20    continue
        ii=ii+1
        p(i,i)=q(ii)
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_sq_from_d (p,nbas,q)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer i,nbas
      real*8 p(nbas,nbas),q(nbas)
c
      call mat_zero (p,nbas)
      do 10 i=1,nbas
        p(i,i)=q(i)
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_sq_dev_d (p,nbas,q)

c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer i,nbas
      real*8 p(nbas,nbas),q(nbas)
c
      call mat_zero (p,nbas)
      do 10 i=1,nbas
        p(i,i)=1.0d0/q(i)
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_plainadd (p,nbas,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),q(nbas,nbas)
c
      do 10 i=1,nbas
        do 20 j=1,nbas
          p(j,i)=p(j,i)+q(j,i)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_tadd (alpha,p,isize,beta,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer isize,i
      real*8 alpha,beta,p(isize),q(isize)
c
      do 10 i=1,isize
          p(i) = alpha*p(i) + beta*q(i)
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_1_over_h (p,nbas,e)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),e(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=p(i,j)/(e(i)+e(j))
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_1_over_h_tri (p,nbas,e)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j,ii
      real*8 p(nbas*(nbas+1)/2),e(nbas)
c
      ii=0
      do 10 i=1,nbas
        do 20 j=1,i-1
          ii=ii+1
          p(ii)=p(ii)/(e(i)+e(j))
  20    continue
        ii=ii+1
        p(ii)=p(ii)/(e(i)+e(i))
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_1_over_h2 (p,nbas,e,q)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),e(nbas),q(nbas,nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=q(i,j)/(e(i)+e(j))
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_times_p2 (p,nbas,pp)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     multiplication with diag(pp(j)) from the right
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),pp(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=p(i,j)*pp(j)
  20    continue
  10  continue
c
      return
      end
c
c
      subroutine mat_times_p2b (p,q,nbas,pp)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     multiplication with diag(pp(j)) from the right
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),q(nbas,nbas),pp(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=q(i,j)*pp(j)
  20    continue
  10  continue
c
      return
      end
c
c
      subroutine mat_times_p2c (p,q,nbas,pp)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     multiplication with diag(pp(j)) from the left
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),q(nbas,nbas),pp(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(j,i)=q(j,i)*pp(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_div_p2 (p,nbas,pp)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     multiplication with diag(1/pp(j)) from the right
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),pp(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=p(i,j)/pp(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_div_p2b (p,q,nbas,pp)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     multiplication with diag(1/pp(j)) from the right
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),q(nbas,nbas),pp(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(i,j)=q(i,j)/pp(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_div_p2c (p,q,nbas,pp)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c     multiplication with diag(1/pp(j)) from the left
c
c   version:  1.0.0
c
c   first version: 13.03.2007  (M. Reiher, Theoretical Chemistry, ETH Zurich)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j
      real*8 p(nbas,nbas),q(nbas,nbas),pp(nbas)
c
      do 10 j=1,nbas
        do 20 i=1,nbas
          p(j,i)=q(j,i)/pp(j)
  20    continue
  10  continue
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_Strans_op (p,nbas,isize,sinv,q,scrmat,p_mat,q_mat)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 22.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,isize,i,j,k
      real*8 p(isize),q(isize),sinv(nbas,nbas),
     *                 scrmat(nbas,nbas),p_mat(nbas,nbas),
     *                 q_mat(nbas,nbas)
c
c   Initialize all matrices
      call mat_zero2 (p,isize)
      call mat_zero (scrmat,nbas)
      call mat_zero (q_mat,nbas)
      call mat_zero (p_mat,nbas)
c
c   Build q_mat
      call mat_sq_from_t (q_mat,nbas,q)
c
c   Calculate scrmat = q * sinv
      do 50 i=1,nbas
        do 60 j=1,nbas
          scrmat(i,j)=0.0d0
          do 70 k=1,j
c  note: sinv is upper triangular matrix!
            scrmat(i,j)=scrmat(i,j)+q_mat(i,k)*sinv(k,j)
  70      continue
  60    continue
  50  continue
c
c   Calculate p_mat = sinv^dagger * scrmat = sinv^dagger * q * sinv
      do 150 i=1,nbas
        do 160 j=1,nbas
          p_mat(i,j)=0.0d0
          do 170 k=1,i
            p_mat(i,j)=p_mat(i,j)+sinv(k,i)*scrmat(k,j)
 170      continue
 160    continue
 150  continue
c
c   Copy p_mat back to array of length isize
      call mat_triang (p,nbas,p_mat)
c
      return
      end
c
c
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine mat_Strans_vec (nbas,sinv,eigvect1,eigvect2)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
c
      integer nbas,i,j,k
      real*8 eigvect1(nbas,nbas),eigvect2(nbas,nbas),
     *                 sinv(nbas,nbas)
c
      do 10 i=1,nbas
        call mat_zero2 (eigvect1(1,i),nbas)
        do 20 j=1,nbas
          do 30 k=j,nbas
            eigvect1(j,i)=eigvect1(j,i)+sinv(j,k)*eigvect2(k,i)
  30     continue
  20    continue
  10  continue
c
      return
      end
