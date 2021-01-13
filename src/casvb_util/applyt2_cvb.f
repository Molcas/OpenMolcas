************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine applyt2_cvb(vec,gjorb,igjorb,
     > i1alf,i1bet,iato,ibto,phato,phbto)
c  Apply T(O) to the vector VEC. O is defined in terms of GJORB.
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension vec(nda,ndb)
      dimension gjorb(norb*norb),igjorb(2,norb*norb)
      dimension i1alf(n1a,norb),i1bet(n1b,norb)
      dimension iato(norb,0:nam1),ibto(norb,0:nbm1)
      dimension phato(norb,nam1),phbto(norb,nbm1)
      save thresh
      data thresh/1.d-10/

      do 1000 ij=1,norb*norb
      iorb=igjorb(2,ij)
      jorb=igjorb(1,ij)
      scale=gjorb(ij)
      if(iorb.ne.jorb.and.abs(scale).gt.thresh)then
c  a) Alpha excitation
        if(absym(2))then
          do 1100 ia=1,n1a
          iaxtmp=i1alf(ia,iorb)
          jax=iato(jorb,iaxtmp)
          if(jax.ne.0)then
            iax=iato(iorb,iaxtmp)
            tcof=scale*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
            if(jax.gt.iax)then
              call daxpy_(ndb-jax+1,tcof,vec(iax,jax),nda,
     >          vec(jax,jax),nda)
            else
              call daxpy_(iax-jax,tcof,vec(jax,iax),1,vec(jax,jax),nda)
              vec(jax,iax)=vec(jax,iax)+tcof*vec(iax,iax)
              if(ndb-iax.gt.0) call daxpy_(ndb-iax,tcof,
     >          vec(iax,iax+1),nda,vec(jax,iax+1),nda)
            endif
          endif
1100      continue
        else
          do 2100 ia=1,n1a
          iaxtmp=i1alf(ia,iorb)
          jax=iato(jorb,iaxtmp)
          if(jax.ne.0)then
            iax=iato(iorb,iaxtmp)
            tcof=scale*phato(iorb,iaxtmp)*phato(jorb,iaxtmp)
            call daxpy_(ndb,tcof,vec(iax,1),nda,vec(jax,1),nda)
          endif
2100      continue
        endif
c  b) Beta excitation
        if(absym(2))then
          do 3100 ib=1,n1b
          ibxtmp=i1bet(ib,iorb)
          jbx=ibto(jorb,ibxtmp)
          if(jbx.ne.0)then
            ibx=ibto(iorb,ibxtmp)
            tcof=scale*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            if(jbx.gt.ibx)then
              call daxpy_(ibx-1,tcof,vec(1,ibx),1,vec(1,jbx),1)
              vec(ibx,jbx)=vec(ibx,jbx)+tcof*vec(ibx,ibx)
              if(jbx-ibx.gt.0) call daxpy_(jbx-ibx,tcof,
     >          vec(ibx,ibx+1),nda,vec(ibx+1,jbx),1)
            else
              call daxpy_(jbx,tcof,vec(1,ibx),1,vec(1,jbx),1)
            endif
          endif
3100      continue
        else
          do 4100 ib=1,n1b
          ibxtmp=i1bet(ib,iorb)
          jbx=ibto(jorb,ibxtmp)
          if(jbx.ne.0)then
            ibx=ibto(iorb,ibxtmp)
            tcof=scale*phbto(iorb,ibxtmp)*phbto(jorb,ibxtmp)
            call daxpy_(nda,tcof,vec(1,ibx),1,vec(1,jbx),1)
          endif
4100      continue
        endif
      elseif(iorb.eq.jorb.and.abs(scale-one).gt.thresh)then
c Alpha singly occupied
        if(absym(2))then
          do 5100 ia=1,n1a
          iak=iato(iorb,i1alf(ia,iorb))
          call dscal_(ndb-iak+1,scale,vec(iak,iak),nda)
5100      continue
        else
          do 5200 ia=1,n1a
          iak=iato(iorb,i1alf(ia,iorb))
          call dscal_(ndb,scale,vec(iak,1),nda)
5200      continue
        endif
c Beta singly occupied
        if(absym(2))then
          do 5300 ib=1,n1b
          ibk=ibto(iorb,i1bet(ib,iorb))
          call dscal_(ibk,scale,vec(1,ibk),1)
5300      continue
        else
          do 5400 ib=1,n1b
          ibk=ibto(iorb,i1bet(ib,iorb))
          call dscal_(nda,scale,vec(1,ibk),1)
5400      continue
        endif
      endif
1000  continue
      if(absym(2))then
        do 6000 ia=1,nda
        do 6001 ib=ia+1,ndb
        vec(ib,ia)=vec(ia,ib)
6001    continue
6000    continue
      endif
      return
      end
