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
      subroutine mkciinfo2_cvb(i1alf,i1bet,iafrm,ibfrm,iato,ibto,
     >  phato,phbto,
     >  xalf,xbet,xalf2,xbet2,mingrph,maxgrph,
     >  nk,locc,lunocc,inewocc,iaccm)
      implicit real*8 (a-h,o-w,y-z),integer(x)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

      dimension i1alf(n1a,norb),i1bet(n1b,norb)
      dimension iafrm(norb,nda),ibfrm(norb,ndb)
      dimension iato(norb,0:nam1),ibto(norb,0:nbm1)
      dimension phato(norb,nam1),phbto(norb,nbm1)
      dimension xalf(0:norb,0:nalf),xbet(0:norb,0:nbet)
      dimension xalf2(0:norb,0:nalf-1),xbet2(0:norb,0:nbet-1)
      dimension mingrph(0:norb),maxgrph(0:norb)
      dimension nk(0:norb),locc(norb+1),lunocc(norb+1)
      dimension inewocc(norb),iaccm(norb)

      call izero(iafrm,nda*norb)
      call izero(ibfrm,ndb*norb)
      call izero(iato,(nam1+1)*norb)
      call izero(ibto,(nbm1+1)*norb)
      call fzero(phato,nam1*norb)
      call fzero(phbto,nbm1*norb)
c Alpha loop:
      call izero(iaccm,norb)
      do 100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf,0)
      maxgrph(iorb)=min(iorb,nalf)
100   continue
      call weight_cvb(xalf,mingrph,maxgrph,nalf,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
200   continue
      call izero(inewocc,norb)
      do 300 i=1,nalf
      inewocc(locc(i))=1
300   continue
      do 325 iel=1,norb
      if(inewocc(iel).eq.1)then
        iaccm(iel)=iaccm(iel)+1
        i1alf(iaccm(iel),iel)=index
      endif
325   continue
c Indexing arrays :
      do 350 i=1,nalf
      inewocc(locc(i))=0
      iafrm(locc(i),index)=indget_cvb(inewocc,nalf,norb,xalf)
      inewocc(locc(i))=1
350   continue
      call loind_cvb(norb,nalf,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xalf,*200)
c Beta loop:
      call izero(iaccm,norb)
      do 400 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet,0)
      maxgrph(iorb)=min(iorb,nbet)
400   continue
      call weight_cvb(xbet,mingrph,maxgrph,nbet,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
500   continue
      call izero(inewocc,norb)
      do 600 i=1,nbet
      inewocc(locc(i))=1
600   continue
      do 625 iel=1,norb
      if(inewocc(iel).eq.1)then
        iaccm(iel)=iaccm(iel)+1
        i1bet(iaccm(iel),iel)=index
      endif
625   continue
c Indexing arrays :
      do 650 i=1,nbet
      inewocc(locc(i))=0
      ibfrm(locc(i),index)=indget_cvb(inewocc,nbet,norb,xbet)
      inewocc(locc(i))=1
650   continue
      call loind_cvb(norb,nbet,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xbet,*500)

c  Altered definitions of I1ALF & I1BET :
      do 700 iorb=1,norb
      do 701 ia=1,n1a
      iax=i1alf(ia,iorb)
      iaxtmp=iafrm(iorb,iax)
      i1alf(ia,iorb)=iaxtmp
701   continue
700   continue
      if(absym(4))then
c  I1ALF & I1BET may share memory:
        if(iiloc(i1alf).ne.iiloc(i1bet))
     >    call imove_cvb(i1alf,i1bet,norb*n1a)
      else
        do 800 iorb=1,norb
        do 801 ib=1,n1b
        ibx=i1bet(ib,iorb)
        ibxtmp=ibfrm(iorb,ibx)
        i1bet(ib,iorb)=ibxtmp
801     continue
800     continue
      endif
c More indexing arrays :
c Second alpha loop (NALF-1):
      do 1100 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nalf-1,0)
      maxgrph(iorb)=min(iorb,nalf-1)
1100  continue
      call weight_cvb(xalf2,mingrph,maxgrph,nalf-1,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
1200  continue
      call izero(inewocc,norb)
      do 1300 i=1,nalf-1
      inewocc(locc(i))=1
1300  continue
      do 1350 i=1,norb-nalf+1
      inewocc(lunocc(i))=1
      iato(lunocc(i),index)=indget_cvb(inewocc,nalf,norb,xalf)
      inewocc(lunocc(i))=0
      if(mod(nalf+i+1+lunocc(i),2).eq.0)then
        phato(lunocc(i),index)= one
      else
        phato(lunocc(i),index)=-one
      endif
1350  continue
      call loind_cvb(norb,nalf-1,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xalf2,*1200)

      if(nbet.le.0)goto 1700
c Second beta loop (NBET-1):
      do 1400 iorb=0,norb
      mingrph(iorb)=max(iorb-norb+nbet-1,0)
      maxgrph(iorb)=min(iorb,nbet-1)
1400  continue
      call weight_cvb(xbet2,mingrph,maxgrph,nbet-1,norb)
      call imove_cvb(maxgrph,nk,norb+1)
      call occupy_cvb(nk,norb,locc,lunocc)
      index=1
1500  continue
      call izero(inewocc,norb)
      do 1600 i=1,nbet-1
      inewocc(locc(i))=1
1600  continue
      do 1650 i=1,norb-nbet+1
      inewocc(lunocc(i))=1
      ibto(lunocc(i),index)=indget_cvb(inewocc,nbet,norb,xbet)
      inewocc(lunocc(i))=0
      if(mod(nbet+i+1+lunocc(i),2).eq.0)then
        phbto(lunocc(i),index)= one
      else
        phbto(lunocc(i),index)=-one
      endif
1650  continue
      call loind_cvb(norb,nbet-1,nk,mingrph,maxgrph,
     >                       locc,lunocc,index,xbet2,*1500)
1700  continue
      return
      end
c  ****************************************
c  ** Routines involving CI, ORBS and VB **
c  ****************************************
c  *********************************************************************
c  *                                                                   *
c  *  EXC1   := apply combination of single-excitations                *
c  *  DENS1  := calculate one-electron (transition) density            *
c  *                                                                   *
c  *  :> Icfrom:    right-hand vector                                  *
c  *  :> Icto:      left-hand vector (dens1: input exc1: output)       *
c  *  :> vij:       Orbital matrix (exc1: input dens1: output)         *
c  *  :> diag:      Include diagonal of vij?                           *
c  *  :> iPvb:      Include Pvb (0: no, 1: on left, 2: on right)       *
c  *                                                                   *
c  *  :< Icto:      left-hand vector (dens1: input exc1: output)       *
c  *  :< vij:       Orbital matrix (exc1: input dens1: output)         *
c  *                                                                   *
c  *********************************************************************
