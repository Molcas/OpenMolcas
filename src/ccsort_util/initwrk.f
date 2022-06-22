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
       subroutine initwrk (length)
c
c      this routine calculate required size of work space and
c      definie initial possitions of work vectors
c
#include "ccsort.fh"
#include "reorg.fh"
c
c     help variables
c
       integer n
       integer sizevint,sizev1,sizev2,sizempq,length,norbmax,sizeri
       integer symp,symq,symi,symj,sympq,sympqi,symm,symmp,syma
c
c1*   def maxzie of vint
c
       norbmax=norb(1)
       do 10 n=1,nsym
       if (norb(n).gt.norbmax) then
       norbmax=norb(n)
       end if
 10     continue
c
       sizevint=norbmax*norbmax*norbmax
c
c2*   def size of <pq|i>=j>, <pq|i,j>
c
       sizev1=0
       sizev2=0
       do 20 symp=1,nsym
       do 21 symq=1,nsym
       sympq=mul(symp,symq)
       do 22 symi=1,nsym
       sympqi=mul(sympq,symi)
       symj=sympqi
c      calc. length
       if (symj.gt.symi) then
         sizev2=sizev2+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
       else
         sizev1=sizev1+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
         sizev2=sizev2+noa(symi)*noa(symj)*NORB(symp)*NORB(symq)
       end if
22     continue
21     continue
20     continue
c
c3*   def maxsize of <_am|pq>
c

       sizempq=0
       do 50 syma=1,nsym
c
       length=0
       do 30 symm=1,nsym
       do 31 symp=1,nsym
       symmp=mul(symm,symp)
       symq=mul(syma,symmp)
c     calc. length
       length=length+noa(symm)*NORB(symp)*NORB(symq)
 31     continue
 30     continue
c
       if (sizempq.lt.length) then
       sizempq=length
       end if
c
 50     continue
c
c4*   def maxsize of R_i if needed
c
       sizeri=0
c
       if (t3key.eq.1) then
       do 60 symi=1,nsym
       call ccsort_t3grc0 (3,8,4,4,4,0,symi,1,length,mapdri,mapiri)
       length=length-1
       if (length.gt.sizeri) then
       sizeri=length
       end if
 60     continue
       end if

c     ******* distribution of memory ******
c
       poss10=1+sizevint
       poss20=poss10+sizev1
       poss30=poss20+sizev2
       possri0=poss30+sizempq
       length=possri0+sizeri-1
c
       if (fullprint.gt.1) then
       write(6,*)
       write(6,'(6X,A)') 'size of help (work) vectors:'
       write(6,'(6X,A)') '----------------------------'
       write(6,*)
       write(6,'(6X,A,I8)') 'Vints     V0 required : ',sizevint
       write(6,'(6X,A,I8)') 'PQIJ ints V1 required : ',sizev1
       write(6,'(6X,A,I8)') '          V2 required : ',sizev2
       write(6,'(6X,A,I8)') 'AMIJ ints V3 required : ',sizempq
       write(6,'(6X,A,I8)') 'R_i mtx   Ri required : ',sizeri
       end if
c
       if (fullprint.ge.0)
     &     write(6,'(6X,A,I20)') 'Required WRK size-sum : ',length
c
       return
       end
