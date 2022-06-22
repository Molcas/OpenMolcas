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
       subroutine unpackk_pck (i,vint,ndimv1,ndimv2,ndimv3,key)
c
c     this routine expand integrals packed in i-th TEMP file
c     to vint(j,k,l) = <i,j|k,l>
c
c     i      - value of pivot index (I)
c     vint   - array of integrals (O)
c     ndimv1 - first dimension of vint (norb(symj)) (I)
c     ndimv2 - second dimension of vint (norb(symk)) (I)
c     ndimv3 - third dimension of vint (norb(syml)) (I)
c     key    - reduced storing key (I)
c     = 0 if symj is not syml
c     = 1 if symj = syml
c
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimv1,ndimv2,ndimv3,key
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
c
c     help variables
c
       integer nhelp,length,daddr,nrec
c
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
c
       character*(RtoB+ItoB) pp(1:nsize),pphelp
       real*8 rhelp
       integer ihelp,ires
c
c*    set vint=0
c
       nhelp=ndimv1*ndimv2*ndimv3
       call ccsort_mv0zero (nhelp,nhelp,vint)
c
c*     open corresponding TEMP file
c
       if (iokey.eq.1) then
c      Fortran IO
       call molcas_binaryopen_vanilla(lunpublic,tmpnam(i))
c       open (unit=lunpublic,file=tmpnam(i),form='unformatted')
       else
c      MOLCAS IO
       call daname (lunpublic,tmpnam(i))
       daddr=0
       end if
c
      do nrec=1,nrectemp(i)
c
        if (nrec.ne.nrectemp(i)) then
        length=nsize
        else
        length=lrectemp(i)
        end if
c
      if (iokey.eq.1) then
c     Fortran IO
      call getpp_pck (lunpublic,pp,length)
      else
c     MOLCAS IO
      call cdafile (lunpublic,2,pp,(RtoB+ItoB)*length,daddr)
      end if

c
c*     get indexes jh,kh,lh and value valh from packed form
c
      do nhelp=1,length
      pphelp=pp(nhelp)
      rhelp=transfer(pphelp(1:RtoB),rhelp)
      ihelp=transfer(pphelp(RtoB+1:),ihelp)
      valh(nhelp)=rhelp
      jh(nhelp)=int(ihelp/constj)
      ires=ihelp-constj*jh(nhelp)
      kh(nhelp)=int(ires/constk)
      lh(nhelp)=ires-constk*kh(nhelp)
      end do
c
      if (key.eq.0) then
      do 100 nhelp=1,length
      vint(jh(nhelp),kh(nhelp),lh(nhelp))=valh(nhelp)
 100  continue
      else
      do 200 nhelp=1,length
      vint(jh(nhelp),kh(nhelp),lh(nhelp))=valh(nhelp)
      vint(lh(nhelp),kh(nhelp),jh(nhelp))=valh(nhelp)
 200  continue
      end if
c
       end do
c
      if (iokey.eq.1) then
c     Fortran IO
       close (lunpublic)
      else
c     Molcas IO
       call daclos (lunpublic)
      end if
c
       return
       end
