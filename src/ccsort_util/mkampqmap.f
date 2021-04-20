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
       subroutine mkampqmap (ammap,syma,rc)
c
c     this routine prepair ammap
c
#include "reorg.fh"
#include "ccsort.fh"
c
       integer syma,rc
       integer ammap(1:mbas,1:8,1:8)
c
c     help variables
c
       integer a,symp,symq,symm,symam
       integer lengthmpq,nrecc,nrest,irec
c
cT    test, if there are any a in this symmtry
c
       if (nvb(syma).eq.0) then
       rc=1
c     RC=1 : there are no a in this symmetry
       return
       else
       rc=0
       end if
c
c*    def initial address
c
       irec=1
c
c*    loop over all combinations
c
       do 100 symm=1,nsym
       symam=mul(syma,symm)
       do 101 symp=1,nsym
       symq=mul(symam,symp)
c
c*    define number of records, required to store this block
c     and determine shift in initial possitions
c
       lengthmpq=noa(symm)*norb(symp)*norb(symq)
       nrecc=int(lengthmpq/recl)
       nrest=lengthmpq-nrecc*recl
       if (nrest.gt.0) then
       nrecc=nrecc+1
       end if
c
       do 102 a=1,nvb(syma)
c
       ammap(a,symm,symp)=irec
       irec=irec+nrecc
c
 102    continue
 101    continue
 100    continue
c
       return
       end
