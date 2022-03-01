!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
       subroutine mkampqmap (ammap,syma,rc)
!
!     this routine prepair ammap
!
#include "reorg.fh"
#include "ccsort.fh"
!
       integer syma,rc
       integer ammap(1:mbas,1:8,1:8)
!
!     help variables
!
       integer a,symp,symq,symm,symam
       integer lengthmpq,nrecc,nrest,irec
!
!T    test, if there are any a in this symmtry
!
       if (nvb(syma).eq.0) then
       rc=1
!     RC=1 : there are no a in this symmetry
       return
       else
       rc=0
       end if
!
!*    def initial address
!
       irec=1
!
!*    loop over all combinations
!
       do 100 symm=1,nsym
       symam=mul(syma,symm)
       do 101 symp=1,nsym
       symq=mul(symam,symp)
!
!*    define number of records, required to store this block
!     and determine shift in initial possitions
!
       lengthmpq=noa(symm)*norb(symp)*norb(symq)
       nrecc=int(lengthmpq/recl)
       nrest=lengthmpq-nrecc*recl
       if (nrest.gt.0) then
       nrecc=nrecc+1
       end if
!
       do 102 a=1,nvb(syma)
!
       ammap(a,symm,symp)=irec
       irec=irec+nrecc
!
 102    continue
 101    continue
 100    continue
!
       return
       end
