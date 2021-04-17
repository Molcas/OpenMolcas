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
       subroutine fokunpck5 (symp,foka,fokb,dpa,dpb,dimfok,rc)
c
c     this routine produce dpa,dpb from foka,fokb
c     for some cases
c     shifto,shiftv will be also added
c
c     symp   - symmtry of this block
c     foka   - Fok aa matrix (I)
c     fokb   - Fok bb matrix (I)
c     dpa    - Diagonal part alfa vector (O)
c     dpa    - Diagonal part beta vector (O)
c     dimfok - dimension for Fok matrix - norb (I)
c     rc     - return (error) code
c
#include "ccsd1.fh"
       integer symp,dimfok,rc
c
       real*8 foka(1:dimfok,1:dimfok)
       real*8 fokb(1:dimfok,1:dimfok)
       real*8 dpa(1:dimfok)
       real*8 dpb(1:dimfok)
c
c     help variables
c
       integer p,nhelp1,nhelp2
c
       rc=0
c
       if (typden.eq.0) then
c1    diagonal elements are required
c
       do 100 p=1,dimfok
       dpa(p)=foka(p,p)
       dpb(p)=fokb(p,p)
 100    continue
c
       else if (typden.eq.1) then
c2    (faa+fbb)/2 are required
c
       do 200 p=1,dimfok
       dpa(p)=(foka(p,p)+fokb(p,p))/2
       dpb(p)=dpa(p)
 200    continue
c
       else if (typden.eq.2) then
c3    orbital energies are required
c
c3.1  def shift
       if (symp.eq.1) then
       nhelp1=0
       else
       nhelp1=0
       do 300 nhelp2=1,symp-1
       nhelp1=nhelp1+norb(nhelp2)
 300    continue
       end if
c
c3.2  map oe to dp
       do 400 p=1,dimfok
       dpa(p)=eps(nhelp1+p)
       dpb(p)=eps(nhelp1+p)
 400    continue
c
       else
c     RC=1 : invalid key (NCI/Stup)
       rc=1
       end if
c
       if ((keysa.eq.3).or.(keysa.eq.4)) then
c     for full adaptation scheme only D and V orbitals are shifted
c
       do 501 p=1,nob(symp)
       dpa(p)=dpa(p)-shifto
       dpb(p)=dpb(p)-shifto
 501    continue
c
       do 502 p=1+noa(symp),norb(symp)
       dpa(p)=dpa(p)+shiftv
       dpb(p)=dpb(p)+shiftv
 502    continue
c
       else
c     for other schemes all orbitals are shifted
c
       do 511 p=1,noa(symp)
       dpa(p)=dpa(p)-shifto
 511    continue
c
       do 512 p=1,nob(symp)
       dpb(p)=dpb(p)-shifto
 512    continue
c
       do 513 p=1+noa(symp),norb(symp)
       dpa(p)=dpa(p)+shiftv
 513    continue
c
       do 514 p=1+nob(symp),norb(symp)
       dpb(p)=dpb(p)+shiftv
 514    continue
c
       end if
c
        if (fullprint.ge.2) then
        write (6,*) ' Diagonal part Dp aa, bb for irrep: ',symp
        do p=1,norb(symp)
        write (6,99) p,dpa(p),dpb(p)
99      format (2x,i4,2(f20.14,2x))
        end do
        end if
c
       return
       end
