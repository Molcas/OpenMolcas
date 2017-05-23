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
* Copyright (C) 2007, Markus Reiher                                    *
************************************************************************
      subroutine evalstring (length,term,coeff,nbas,posu,post,poss,vv,
     *                       nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                       snumber,tnumber,unumber,scrno1,
     *                       scr1,scr2,scr3,scr5,nbasp,nbaso,
     *                       dkhadr,adrmem,dkh_48,adrnext)
c
c************************************************************************
c
c   Evaluate the low-level character string 'term' of length 'length'
c     with coefficient 'coeff'. The result is stored in scr5(,,scrno1).
c
c   Note:  'term' might in general contain Sxxx, Txxx, and Uxxx expressions
c           as well as brackets [,], and V,N,D,Y,F,G,Z,Q,X,I,J,K,L,M  !!
c
c   The logic in this sub works only for auxiliary matrices of the
c     type AXX, i.e. with 2 numbers XX --- at most 99 matrices (A00 is
c     used elsewhere)
c
c   This routine analyzes term with respect to the occurrence
c     of brackets (which may show a complicated structure of nested
c     and contiguous brackets especially at high orders).
c
c   NOTE: This sub includes a program package dependent dynamic
c         memory allocation of square matrices (number is depending
c         on the bracket structure of term).
c
c   Possibility for Improvement:
c   Unfortunately, the memory manager GetMem() does not allow to have
c     "unbraced" memory allocations. Hence, the auxiliary matrices are
c     not allocated per level but per term. If one would switch to
c     a different memory manager (like available in F90) one might
c     activate the ievalold and countold lines to sae more memory.
c
c   This SR belongs to dkhparser_numeric.
c
c   written by:  Markus Reiher  (ETH Zurich)
c
c   version:  1.1.0
c
c   modified: 13.03.2007 MR@ETH
c   first version: 19.01.2007  (Theoretical Chemistry, ETH Zurich)
c
c************************************************************************
c
      implicit none
#include "dkhparameters.fh"
#include "WrkSpc.fh"
c
      integer length,nbas,posu(maxunumber),post(maxsnumber),nbasp,
     *        poss(maxsnumber),snumber,tnumber,unumber,scrno1,nbaso
      character*(maxlength) term
      character*(2) dkh_int2char2
      REAL*8 coeff,alpha
c
      REAL*8 vv(nbaso,nbaso),nn(nbaso,nbaso),dd(nbaso,nbaso),
     *                 yy(nbaso,nbaso),ff(nbaso,nbaso),gg(nbaso,nbaso),
     *                 pp(nbas),xx(nbasp,nbasp),
     *                 ii(nbasp,nbasp),jj(nbasp,nbasp),kk(nbasp,nbasp),
     *                 ll(nbasp,nbasp),mm(nbasp,nbasp),e(nbas)
      REAL*8 scr1(nbaso,nbaso,snumber),
     *                 scr2(nbaso,nbaso,tnumber),
     *                 scr3(nbaso,nbaso,unumber),scr5(nbas,nbas,scrno1)
c
      integer ilevel(maxlength),maxlevel,levelcounter,imarker,count,
     *        iscr,ieval(99),minlevel
c
      integer i,idum1,ileft,iright,scrleng,j
      character*(maxlength) scrtxt
c
      integer dkh_48,adrmem,adrnext
      integer dkhadr(adrmem)
c
*** get one additional scratch matrix
      iscr=0
      call GetMem('EvalSt1 ','ALLO','REAL',iscr,nbas*nbas+4)
c
c     analyze term with respect to number and level of brackets
c     i.e., a string like ABC[DE]FGH[[IJ]K]
c             is coded as 11122221112333322
c
      maxlevel=1
      levelcounter=1
      do i=1,length
        if (term(i:i).eq.'[') levelcounter=levelcounter+1
        if (maxlevel.lt.levelcounter) maxlevel=levelcounter
        ilevel(i)=levelcounter
        if (term(i:i).eq.']') levelcounter=levelcounter-1
      end do
      minlevel=maxlevel
      do i=1,length
        if (ilevel(i).lt.minlevel) minlevel=ilevel(i)
      end do
c
      if (maxlevel.gt.2) then
c
*** find innermost brackets if present
        count=0
        do levelcounter=maxlevel,minlevel,-1
          if (levelcounter.gt.minlevel) then
            alpha=1.0d0
          else
*** set multiplicative scalar only the last step
            alpha=coeff
          end if
          ileft=1
          imarker=1
CMR          count=0
2000      continue
          iright=length
          do j=imarker,length
            if (ilevel(j).eq.levelcounter) then
              ileft=j
              do i=ileft+1,length
                if (ilevel(i).eq.levelcounter) then
                  iright=i
                else if (ilevel(i).eq.levelcounter-1) then
*** found; leave loop
                  goto 3000
                end if
              end do
3000          continue
*** pass over string in scrtxt copied from term(ileft:iright) to evalstring2()
              scrleng=iright-ileft+1
              scrtxt(1:scrleng)=term(ileft:iright)
*** allocate a single scratch matrix
              if (levelcounter.gt.1) then
                count=count+1
                if (count.gt.99.and.levelcounter.gt.1) then
                  write(6,*) "ERROR in evalstring(): Not enough memory"
     *                       //" for auxiliary matrices !"
                  call Abend
                end if
                call GetMem('EvalSt  ','ALLO','REAL',ieval(count),
     *                    nbas*nbas+4)
              end if
              if (levelcounter.lt.maxlevel) then
                call evalstring2b(scrleng,scrtxt,alpha,nbas,posu,post,
     *                   poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                   e,snumber,tnumber,unumber,scrno1,scr1,
     *                   scr2,scr3,scr5,nbasp,nbaso,work(iscr),
     *                   ieval,count,dkhadr,adrmem,dkh_48,adrnext)
CMR     *                   ievalold,countold)
CMR                call mat_copy(work(ievalold(countold)),nbas,nbas,
                call mat_copy(work(ieval(count)),nbas,nbas,
     *                        scr5(1,1,scrno1))
              else
*** there are no auxiliary matrices yet as we are at the lowest level
                call evalstring2(scrleng,scrtxt,alpha,nbas,posu,post,
     *                   poss,vv,nn,dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,
     *                   e,snumber,tnumber,unumber,scrno1,scr1,
     *                   scr2,scr3,scr5,nbasp,nbaso,work(iscr),
     *                   dkhadr,adrmem,dkh_48,adrnext)
                call mat_copy(work(ieval(count)),nbas,nbas,
     *                        scr5(1,1,scrno1))
              end if
              if (levelcounter.gt.1) then
*** introduce new AXX matrix symbol and adjust term as well as ilevel
                term(ileft:ileft+2)='A'//dkh_int2char2(count)
                idum1=(iright-ileft+1)-3
*** adjust indices as well
                ilevel(ileft)=ilevel(ileft)-1
                ilevel(ileft+1)=ilevel(ileft+1)-1
                ilevel(ileft+2)=ilevel(ileft+2)-1
                do 10 i=ileft+3,length-idum1
                  term(i:i)=term(i+idum1:i+idum1)
                  ilevel(i)=ilevel(i+idum1)
10              continue
                do 15 i=length-idum1+1,length
                  term(i:i)=' '
15              continue
                length=length-idum1
*** adust iright as well for the reduced character string size
                iright=iright-idum1
              end if
*** check whether there are additional contributions from this level
              if (iright.ge.length) goto 2001
CMR              if (iright.ne.length.and.levelcounter.ne.1) then
              imarker=iright+1
              goto 2000
CMR              end if
            end if
          end do
*** now all brackets have been resolved for the next level but they may be
***   discontinuously distributed
*
*** free all memory required for the level i-1 as these are now no longer needed
2001       continue
        end do
            do i=count,1,-1
              CALL GetMem('EvalSt  ','FREE','REAL',ieval(i),
     *                         nbas*nbas+4)
            end do
c
      else
*** this string contains at most single brackets -> one call to evalstring2()
***   note: scr5 is used as scr5(1,1,scrno1) and scr5(1,1,scrno1-1)
        call evalstring2(length,term,coeff,nbas,posu,post,poss,vv,nn,
     *                   dd,yy,ff,gg,pp,xx,ii,jj,kk,ll,mm,e,
     *                   snumber,tnumber,unumber,scrno1,scr1,
     *                   scr2,scr3,scr5,nbasp,nbaso,work(iscr),
     *                   dkhadr,adrmem,dkh_48,adrnext)
      end if
      call GetMem('EvalSt1 ','FREE','REAL',iscr,nbas*nbas+4)
c
      return
      end
c
c
