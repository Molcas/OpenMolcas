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
C
C----------------------------------------------------------------------|
C
      subroutine dkh_prop(n,s,t,v,w,X,pXp,clight,dkord,xord,dkparam)
C
C Apply the arbitrary order DKH transformation to property integral
C
      implicit none
#include "WrkSpc.fh"
C Input
C X     matrix of property operator
C pXp   matrix representation of <pxXpx>+<pyXpy>+<pzXpz>
C w     aka pVp
      integer n,dkord,xord,dkparam
      Real*8 s(n,n),t(n,n),v(n,n),w(n,n),clight
      Real*8 X(n,n),pXp(n,n)
C Output : X store the transformed property integral
C Temp
      integer nn,vord,nz,m
      integer iTr,iBk,iEL,iES,iOL,iOS,iEp,iE0,iKC,iCo,iSco
      integer iM,iZ,iW
C
C Transform to free-particle FW picture
C
      nn=n*n+4
      call getmem('Tr  ','ALLOC','REAL',iTr ,nn)
      call getmem('Back','ALLOC','REAL',iBk ,nn)
      call getmem('mEL ','ALLOC','REAL',iEL ,nn)
      call getmem('mES ','ALLOC','REAL',iES ,nn)
      call getmem('mOL ','ALLOC','REAL',iOL ,nn)
      call getmem('mOS ','ALLOC','REAL',iOS ,nn)
      call getmem('Ep  ','ALLOC','REAL',iEp ,n+4)
      call getmem('E0  ','ALLOC','REAL',iE0 ,n+4)
      call getmem('KC  ','ALLOC','REAL',iKC ,n*3+4)
      call XDR_fpFW(n,s,t,v,w,Work(iTr),Work(iBk),Work(iEL),Work(iES),
     &              Work(iOL),Work(iOS),Work(iEp),work(iE0),Work(iKC),
     &              Work(iKC+n),Work(iKC+2*n),clight)
C
C Calculate the DKH unitary transformation ( in terms of a series of W matrices 1..xord )
C
      vord = xord * 2
      m = n*n
      nz = m*vord
      call getmem('Wsav','ALLOC','REAL', iW, m*xord*2+4)
      call getmem('Cof ','ALLOC','REAL',iCo, vord+8)
      call dkh_cofu(vord,dkparam,Work(iCo))
C
      call getmem('Cof2','ALLOC','REAL',iSCo,vord+8)
      call getmem('Mat ','ALLOC','REAL', iM, m*6+4)
      call getmem('Mat2','ALLOC','REAL', iZ, nz*10+4)
      call dkh_ham(n,vord,xord,vord,Work(iEL),Work(iES),Work(iOL),
     &             Work(iOS),Work(iEp),Work(iE0),Work(iCo),Work(iSco),
     &             Work(iM),Work(iM+m),Work(iM+m*2),Work(iM+m*3),
     &             Work(iM+m*4),Work(iM+m*5),
     &             Work(iZ),Work(iZ+nz),Work(iZ+nz*2),Work(iZ+nz*3),
     &             Work(iZ+nz*4),Work(iZ+nz*5),Work(iZ+nz*6),
     &             Work(iZ+nz*7),Work(iZ+nz*8),Work(iZ+nz*9),
     &             Work(iW) )
C
C Apply W[1..xord] determined transformation to property operator X
C
C     ! convert X to fpFW picture
      call XDR_fpFWprop(n,Work(iTr),X,pXp,Work(iKC),Work(iKC+n),
     &                  Work(iKC+2*n),Work(iEL),Work(iES),Work(iOL),
     &                  Work(iOS),Work(iM) )
      call dkh_xpx(n,vord,xord,vord,Work(iEL),Work(iES),Work(iOL),
     &             Work(iOS),Work(iEp),Work(iE0),Work(iCo),Work(iSco),
     &             Work(iM),Work(iM+m),Work(iM+m*2),Work(iM+m*3),
     &             Work(iM+m*4),Work(iM+m*5),
     &             Work(iZ),Work(iZ+nz),Work(iZ+nz*2),Work(iZ+nz*3),
     &             Work(iZ+nz*4),Work(iZ+nz*5),Work(iZ+nz*6),
     &             Work(iZ+nz*7),Work(iZ+nz*8),Work(iZ+nz*9),
     &             Work(iW) )
      call getmem('Cof2','FREE','REAL',iSCo,vord+8)
      call getmem('Mat ','FREE','REAL',iM, m*6+4)
      call getmem('Mat2','FREE','REAL',iZ, nz*10+4)
C
C Back transform to original non-orthogonal basis picture
C
      call dmxma(n,'C','N',Work(iBk),Work(iEL),Work(iES),1.d0)
      call dmxma(n,'N','N',Work(iES),Work(iBk),X,1.d0)
C
C Free temp memories
C
      call getmem('Cof ','FREE','REAL',iCo, vord+8)
      call getmem('Wsav','FREE','REAL',iW, m*xord*2+4)
      call getmem('Tr  ','FREE','REAL',iTr ,nn)
      call getmem('Back','FREE','REAL',iBk ,nn)
      call getmem('mEL ','FREE','REAL',iEL ,nn)
      call getmem('mES ','FREE','REAL',iES ,nn)
      call getmem('mOL ','FREE','REAL',iOL ,nn)
      call getmem('mOS ','FREE','REAL',iOS ,nn)
      call getmem('Ep  ','FREE','REAL',iEP ,n+4)
      call getmem('E0  ','FREE','REAL',iE0 ,n+4)
      call getmem('KC  ','FREE','REAL',iKC ,n*3+4)
C
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(dkord)
      end
