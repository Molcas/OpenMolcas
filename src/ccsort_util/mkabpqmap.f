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
       subroutine mkabpqmap (abmap,syma,symb,rc)
c
c     this routine prepair abmap
c
#include "reorg.fh"
#include "ccsort.fh"
c
       integer abmap(1:mbas,1:mbas,1:8)
       integer syma,symb,rc
c
c     help variables
c
       integer a,b,bup,symp,symq,symab
       integer lengthpq,nrecc,nrest,irec
c
cT    test, if there are any ab pair
c
       if (nvb(syma)*nvb(symb).eq.0) then
       rc=1
c     RC=1 : there are no ab pair in this symmetry
       return
       else
       rc=0
       end if
c
c*    def initial address
c
       irec=1
       symab=mul(syma,symb)
c
c*    loop over all combinations
c
       do 100 symp=1,nsym
       symq=mul(symab,symp)
c
c*    define number of records, required to store this block
c     and determine shift in initial possitions
c
       lengthpq=norb(symp)*norb(symq)
       nrecc=int(lengthpq/recl)
       nrest=lengthpq-nrecc*recl
       if (nrest.gt.0) then
       nrecc=nrecc+1
       end if
c
       do 101 a=1,nvb(syma)
c
       if (syma.eq.symb) then
       bup=a
       else
       bup=nvb(symb)
       end if
c
       do 102 b=1,bup
c
       abmap(a,b,symp)=irec
       irec=irec+nrecc
c
 102    continue
 101    continue
 100    continue
c
       return
       end
