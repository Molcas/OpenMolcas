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
       subroutine unpckhelp7 (a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof)
c
c     this routine do:
c     b(e,f) =  -a(pf,qe)
c
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf)
c
c     help variables
       integer qe,pf,f
c
       do 100 pf=fadd+1,fadd+nof
       f=pf-fadd
       do 101 qe=eadd+1,eadd+noe
       b(qe-eadd,f)=-a(pf,qe)
 101    continue
 100    continue
c
       return
       end
