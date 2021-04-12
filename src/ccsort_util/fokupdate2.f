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
       subroutine fokupdate2 (foka,symp,i,vint,ndimv1,ndimv2,ndimv3)
c
c     this routine realize update
c     foka(p,q) = foka(p,q) - <ip|qi>
c
c     N.B. integrals are of type <symi, symp| symp, symi>
c
c     foka    - packed Fokaa matrix (I,O)
c     symp    - irrep or p (and also q) index (I)
c     i       - value of i, (I)
c     vint    - array of integrals <ip|iq> for given i (I)
c     ndimv1  - first dimension (norb(symp)) (I)
c     ndimv2  - second dimension (norb(symi)) (I)
c     ndimv3  - third dimension (norb(symp)) (I)
c
#include "ccsort.fh"
       real*8 foka(*)
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
       integer symp,i,ndimv1,ndimv2,ndimv3
c
c     help variables
c
       integer nhelp1,nhelp2,p,q,pq
c
c*    calculate shift
c
       nhelp1=0
       if (symp.gt.1) then
       do 100 nhelp2=1,symp-1
       nhelp1=nhelp1+(norb(nhelp2)**2+norb(nhelp2))/2
 100    continue
       end if
c
c*    add integral
c
       pq=nhelp1
       do 200 p=1,norb(symp)
       do 201 q=1,p
       pq=pq+1
       foka(pq)=foka(pq)-vint(p,q,i)
 201    continue
 200    continue
c
       return
       end
