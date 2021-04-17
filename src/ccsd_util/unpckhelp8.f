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
       subroutine unpckhelp8 (a,b,dimp,dimef,eadd,noe,bb,dimb)
c
c     this routine do:
c     b(ef,_Bb) = a(pe,qf)-a(qf,pe) for symp=symq
c
       integer dimp,dimef,eadd,noe,bb,dimb
       real*8 a(1:dimp,1:dimp)
       real*8 b(1:dimef,1:dimb)
c
c     help variables
       integer pe,qf,ef
c
       ef=0
       do 100 pe=eadd+2,eadd+noe
       do 101 qf=eadd+1,pe-1
       ef=ef+1
       b(ef,bb)=a(pe,qf)-a(qf,pe)
 101    continue
 100    continue
c
       return
       end
