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
      subroutine XDR_Prop(nbas,isize,jsize,imethod,paratyp,dkhorder,
     &                    xorder,inS,inK,inV,inpVp,inX,inpXp,
     &                    inUL,inUS,clight)
C
C Driver for relativistic transformation of property integrals
C
C   also called the "picture change correction" since the physical operator
C   X is defined in four-component picture, one need to transform them
C   to two/one- component picture as well as the Hamiltonian in the two/one-
C   component relativistic calculations
C
      implicit none
#include "WrkSpc.fh"
C Input variables
      integer nbas,isize,jsize,imethod,paratyp,dkhorder,xorder
      Real*8 clight
      Real*8 inS(isize),inK(isize),inpVp(isize),inpXp(isize)
      Real*8 inV(isize),inUL(jsize),inUS(jsize)
C Input/Output variables
      Real*8 inX(isize)
C Local variables
      integer nn,i,j,k
      integer jS,jK,jV,jpVp,jX,jpXp,itmp
C
C Convert to square matrices
C
      nn = nbas*nbas+4
      call getmem('skin ','ALLOC','REAL',jK   ,nn)
      call getmem('sSS  ','ALLOC','REAL',jS   ,nn)
      call getmem('sV   ','ALLOC','REAL',jV   ,nn)
      call getmem('spVp ','ALLOC','REAL',jpVp ,nn)
      call getmem('sX   ','ALLOC','REAL',jX   ,nn)
      call getmem('spXp ','ALLOC','REAL',jpXp ,nn)
#ifdef MOLPRO
      call square(inK  ,Work(jK),nbas,nbas)
      call square(inS  ,Work(jS),nbas,nbas)
      call square(inV  ,Work(jV),nbas,nbas)
      call square(inpVp,Work(jpVp),nbas,nbas)
      call square(inX  ,Work(jX),nbas,nbas)
      call square(inpXp,Work(jpXp),nbas,nbas)
#else
      call square(inK  ,Work(jK),nbas,1,nbas)
      call square(inS  ,Work(jS),nbas,1,nbas)
      call square(inV  ,Work(jV),nbas,1,nbas)
      call square(inpVp,Work(jpVp),nbas,1,nbas)
      call square(inX  ,Work(jX),nbas,1,nbas)
      call square(inpXp,Work(jpXp),nbas,1,nbas)
#endif
C
C Calculate the relativistic transformed property integrals
C
      if(imethod.eq.2.or.imethod.eq.3.or.
     &   (imethod.eq.1.and.xorder.ge.15) )then
C
C X2C/BSS transformation
C
C   because the transformation matrix in non-orthogonal basis picture has
C   obtained via the Hamiltonian drivers, here we just need to simply apply
C   it to property integrals ( X, pXp in four-component picture )
C
C   high order DKH can also employ this formulation, only negligible
C   contribution from higher orders is included
C
        call getmem('TMP ','ALLOC','REAL',itmp ,nn)
C        ! eval U_L^{\dag} X U_L
        call dmxma(nbas,'C','N',inUL,Work(jX),Work(itmp),1.d0)
        call dmxma(nbas,'N','N',Work(itmp),inUL,Work(jX),1.d0)
C        ! eval U_S^{\dag}pXp U_S
        call dmxma(nbas,'C','N',inUS,Work(jpXp),Work(itmp),1.d0)
        call dmxma(nbas,'N','N',Work(itmp),inUS,Work(jpXp),1.d0)
C        ! sum
        do k=0,nbas*nbas-1
          Work(jX+k) = Work(jX+k) + Work(jpXp+k)
        end do
        call getmem('TMP ','FREE','REAL',itmp ,nn)
C
      else if(imethod.eq.1)then
C
C Arbitrary order DKH transformation
C
        call dkh_prop(nbas,Work(jS),Work(jK),Work(jV),Work(jpVp),
     &              Work(jX),Work(jpXp),clight,dkhorder,xorder,paratyp)
      endif
C
C Copy transformed property integral back to inX
C
      k=0
      do i=1,nbas
        do j=1,i
          k=k+1
          inX(k)=Work(jX+j-1+(i-1)*nbas)
        enddo
      enddo
C
C Free temp memories
C
      call getmem('skin ','FREE','REAL',jK   ,nn)
      call getmem('sSS  ','FREE','REAL',jS   ,nn)
      call getmem('sV   ','FREE','REAL',jV   ,nn)
      call getmem('spVp ','FREE','REAL',jpVp ,nn)
      call getmem('sX   ','FREE','REAL',jX   ,nn)
      call getmem('spXp ','FREE','REAL',jpXp ,nn)
C
      return
      end
