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
       subroutine exppqij (wrk,wrksize,
     & typv2,typp,typq,typr,typs,directyes,
     &                     inverseyes)
c
c     this routine realize reorganization to
c     #2 <p q i j> with given typv2 and typp-typs  <-  #1 <pp,qq|i,j>
c     #1 is in shape <pp,qq|i,j> for sympp,symqq,symi>=symj with types
c     pp,qq,i,j -5,5,1,1
c     #2 <p q i j> may be antisymetrized or not, two parameters (directyes,
c     inverseyes) can be deduced trom typv2 and typp-s, but for simplicity
c     thes are as input parameters
c     this routine do not allow to use typv2=2
c
c     typv2     - typ of final #2 (I)
c     typp-s    - types of ind. p-s (I)
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
       integer typv2,typp,typq,typr,typs,directyes,inverseyes
c
c     help variables
c
       integer symp,symq,symi,symj,possv2,length
       integer ii,iiv1d,iiv1i,possv1d,possv1i
       integer posst
c
c*    get mapd mapi of <p,q r,s> into mapd2,mapi2
c
       call ccsort_grc0 (4,typv2,typp,typq,typr,typs,1,
     & poss20,posst,mapd2,mapi2)
c
c*    realize reorganization psb
c
       do 100 ii=1,mapd2(0,5)
c
c*    def parameters of #2
       possv2=mapd2(ii,1)
       length=mapd2(ii,2)
       symp=mapd2(ii,3)
       symq=mapd2(ii,4)
       symi=mapd2(ii,5)
       symj=mapd2(ii,6)
c
c*    skip this step if length=0
       if (length.eq.0) then
       goto 100
       end if
c
c*    vanish #2
       call ccsort_mv0zero (length,length,wrk(possv2))
c
       if (symi.ge.symj) then
c*    case symi>=symj - integrals in #1 are in that shape
c
       if (directyes.eq.1) then
c**   def possition #1 direct (i.e. #1 <symp,symq| symi,symj>)
c     direct integrals are always used
       iiv1d=mapi1(symp,symq,symi)
       possv1d=mapd1(iiv1d,1)
c
c**   do #2 <p q i j> <- #1 <p,q|i,j> (i.e. direct)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 is the same
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 1,2,3,4,5,5,1,1,
     & typv2,possv1d,possv2,1.0d0)
c
       end if
c
       if (inverseyes.eq.1) then
c
c**   def possition #1 inverse (i.e. #1 <symq,symp| symi,symj>)
c     inverse integrals are used only if antysymetry is required
       iiv1i=mapi1(symq,symp,symi)
       possv1i=mapd1(iiv1i,1)
c
c**   do #2 <p q i j> <- - #1 <q,p|i,j> (i.e. inverse)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 are inversed 1<->2
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 2,1,3,4,5,5,1,1,
     & typv2,possv1i,possv2,-1.0d0)
c
       end if
c
       else
c*    case symi<symj - integrals in #1 are in inverse (symi>=symj) shape
c
       if (directyes.eq.1) then
c
c**   def possition #1 direct (i.e. #1 <symq,symp| symj,symi>)
c     direct integrals are always used
       iiv1d=mapi1(symq,symp,symj)
       possv1d=mapd1(iiv1d,1)
c
c**   do #2 <p q i j> <- #1 <q,p|j,i> (i.e. direct)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 is inversed 1<->2, 3<->4
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 2,1,4,3,5,5,1,1,
     & typv2,possv1d,possv2,1.0d0)
c
       end if
c
       if (inverseyes.eq.1) then
c
c**   def possition #1 inverse (i.e. #1 <symp,symq| symj,symi>)
c     inverse integrals are used only if antysymetry is required
       iiv1i=mapi1(symp,symq,symj)
       possv1i=mapd1(iiv1i,1)
c
c**   do #2 <p q i j> <- - #1 <p,q|j,i> (i.e. inverse)
c     N.B. Since #1 is  always >= symj
c     so in this case orede of indices in #1 and #2 are inversed 3<->4
       call  ireorg (wrk,wrksize,
     & symp,symq,symi,symj,typp,typq,typr,typs,
     & 1,2,4,3,5,5,1,1,
     & typv2,possv1i,possv2,-1.0d0)
c
       end if
c
       end if
c
 100    continue
c
       return
       end
