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
       subroutine multc0 (wrk,wrksize,
     & mvec,ix,mapdc,key)
c
c     This routine realize multiplying according mvec
c     for C=A*B
c     N.B. if key=0, C file is not vanished (ie can be used for
c     adding to some existing file)
c
c     If C=A*B process is faster or comparambe with C=AT*B then mchntyp should be set to 1.
c     If C=AT*B is significantly faster than C=A*T (more than 20%), than mchntyp should be set
c     to 2. (default is 1)
c     if mchntyp is 2, than
c     1) proceses with scale(A)/scale(B) > scalelim will be calculated as C=A*B
c     2) processes with scale(A)/scale(B) < scalelim will be calculated as C=AT*B
c     Note, that for mchntyp =2 more memory is required, due to requirement of
c     aditional o2v2 help file possd0        (parameter possd0 is transported through ccsd1.fh, not
c     through ccsd2.fh)
c
c
#include "ccsd1.fh"
#include "wrk.fh"
       integer mvec(1:4096,1:7)
       integer ix,key
       integer mapdc(0:512,1:6)

c     help variables
c
       integer nhelp1,nhelp2,nhelp3,nhelp4,nhelp5,nhelp6
       integer iix,ic
       real*8 scale
c
c1    set C=0
c
       if (key.eq.1) then
c
c     C matrix must be vanished
c
       do 100 ic=1,mapdc(0,5)
       nhelp1=mapdc(ic,1)
       nhelp2=mapdc(ic,2)
       call mv0zero (nhelp2,nhelp2,wrk(nhelp1))
 100    continue
c
       end if
c
c2    C=C+A*B
c
       if (ix.eq.0) then
       return
       end if
c
       do 200 iix=1,ix
c
c     skip this sumation if yes/no=0
       if (mvec(iix,1).eq.0) goto 200
c
c     realize individial sumation
c
c     def possitions of A,B,C
       nhelp1=mvec(iix,2)
       nhelp2=mvec(iix,3)
       nhelp3=mvec(iix,4)
c
c     def rowA(rowC), colA(rowB,sum), colB(colC)
       nhelp4=mvec(iix,5)
       nhelp5=mvec(iix,6)
       nhelp6=mvec(iix,7)
c
       if (mchntyp.eq.1) then
c
c*    Typ 1
       call mc0c1a3b (nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,
     & nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))
c
       else
c
c*    Typ2
       scale=(1.0d0*nhelp4)/(1.0d0*nhelp6)
       if (scale.gt.slim) then
       call mc0c1a3b (nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,
     & nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

       else
c     map D=AT
       call map21 (wrk(nhelp1),wrk(possd0),nhelp4,nhelp5,2,1,1)
c     calc C=DT*B
       call mc0c1at3b (nhelp5,nhelp4,nhelp5,nhelp6,nhelp4,nhelp6,
     & nhelp4,nhelp5,nhelp6,wrk(possd0),wrk(nhelp2),wrk(nhelp3))
       end if
c
       end if
c
 200    continue
c
       return
       end
