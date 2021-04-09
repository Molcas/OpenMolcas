************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine expmpq (wrk,wrksize,
     & syma,typv3,typm,typp,typq,directyes,
     &                    inverseyes)
c
c     this routine realize reorganization to
c     #3 <m,a p q> with given typv3 and typm,p,q  <-  #2 <a,m|pp,qq>
c     #2 is in shape <a,m|pp,qq> for symm,sympp,symqq with types
c     _a  m,pp,qq - 1,5,5
c     #3 <m a p q> may be antisymetrized or not, two parameters (directyes,
c     inverseyes) can be deduced trom typv3 and typm,p,q and syma
c     but for simplicity these are as input parameters
c     this routine allow to use typv3=0 and 2
c
c     syma      - irrep of a
c     typv3     - typ of final #2 (I)
c     typm,p,q  - types of ind. m,p,q (I)
c     directyes - 1 if direct <pqij> integrals are included (I)
c     inverseyes- 1 if inverse <qpij> integrals are included (I)
c
c     foreingh routines used: grc0
c     ccsort_mv0zero
c
c     it also defines new mapd2,mapi2 corresponding to #2
c
#include "wrk.fh"
#include "reorg.fh"
#include "ccsort.fh"
       integer typv3,typp,typq,typm,syma,directyes,inverseyes
c
c     help variables
c
       integer symp,symq,symm,possv3,length
       integer ii,iiv2d,iiv2i,possv2d,possv2i
       integer posst
c
c*    get mapd mapi of <m,a|p,q> as _a(m,p q) into mapd3,mapi3
c
       call ccsort_grc0 (3,typv3,typm,typp,typq,0,syma,
     & poss30,posst,mapd3,mapi3)
c
c*    realize reorganization psb
c
       do 100 ii=1,mapd3(0,5)
c
c*    def parameters of #3
       possv3=mapd3(ii,1)
       length=mapd3(ii,2)
       symm=mapd3(ii,3)
       symp=mapd3(ii,4)
       symq=mapd3(ii,5)
c
c*    vanish #3
       call ccsort_mv0zero (length,length,wrk(possv3))
c
       if (directyes.eq.1) then
c
c**   def possition #2 direct (i.e.
       iiv2d=mapi2(symm,symq,1)
       possv2d=mapd2(iiv2d,1)
c
c**   do #3 <m a p q> <- #2 <a m q p> (i.e. direct)
       call  mreorg (wrk,wrksize,
     & symm,symp,symq,typm,typp,typq,
     & 1,3,2,1,5,5,
     & typv3,possv2d,possv3,1.0d0)
c
       end if
c
       if (inverseyes.eq.1) then
c
c**   def possition #2 inverse (i.e. #2 <symq,symp| symi,symj>)
       iiv2i=mapi2(symm,symp,1)
       possv2i=mapd2(iiv2i,1)
c
c**   do #3 <m a q p> <- - #2 <a m p q> (i.e. inverse)
       call  mreorg (wrk,wrksize,
     & symm,symp,symq,typm,typp,typq,
     & 1,2,3,1,5,5,
     & typv3,possv2i,possv3,-1.0d0)
c
       end if
c
 100    continue
c
       return
       end
