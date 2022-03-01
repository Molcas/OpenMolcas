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
       subroutine unpackk_pck (i,vint,ndimv1,ndimv2,ndimv3,key)
!
!     this routine expand integrals packed in i-th TEMP file
!     to vint(j,k,l) = <i,j|k,l>
!
!     i      - value of pivot index (I)
!     vint   - array of integrals (O)
!     ndimv1 - first dimension of vint (norb(symj)) (I)
!     ndimv2 - second dimension of vint (norb(symk)) (I)
!     ndimv3 - third dimension of vint (norb(syml)) (I)
!     key    - reduced storing key (I)
!     = 0 if symj is not syml
!     = 1 if symj = syml
!
#include "reorg.fh"

#include "SysDef.fh"
       integer i,ndimv1,ndimv2,ndimv3,key
       real*8 vint(1:ndimv1,1:ndimv2,1:ndimv3)
!
!     help variables
!
       integer nhelp,length,daddr,nrec
!
       integer constj
       parameter (constj=1048576)
       integer constk
       parameter (constk=1024)
!
       character*(RtoB+ItoB) pp(1:nsize),pphelp
       real*8 rhelp
       integer ihelp,ires
!
!*    set vint=0
!
       nhelp=ndimv1*ndimv2*ndimv3
       call ccsort_mv0zero (nhelp,nhelp,vint)
!
!*     open corresponding TEMP file
!
       if (iokey.eq.1) then
!      Fortran IO
       call molcas_binaryopen_vanilla(lunpublic,tmpnam(i))
!       open (unit=lunpublic,file=tmpnam(i),form='unformatted')
       else
!      MOLCAS IO
       call daname (lunpublic,tmpnam(i))
       daddr=0
       end if
!
      do nrec=1,nrectemp(i)
!
        if (nrec.ne.nrectemp(i)) then
        length=nsize
        else
        length=lrectemp(i)
        end if
!
      if (iokey.eq.1) then
!     Fortran IO
      call getpp_pck (lunpublic,pp,length)
      else
!     MOLCAS IO
      call cdafile (lunpublic,2,pp,(RtoB+ItoB)*length,daddr)
      end if

!
!*     get indexes jh,kh,lh and value valh from packed form
!
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
!
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
!
       end do
!
      if (iokey.eq.1) then
!     Fortran IO
       close (lunpublic)
      else
!     Molcas IO
       call daclos (lunpublic)
      end if
!
       return
       end
