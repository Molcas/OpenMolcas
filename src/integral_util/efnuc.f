!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine EFNuc(CoOP,Chrg,Coor,nAtm,ESIT,nOrdOp)
!***********************************************************************
!                                                                      *
! Object: to compute the electricstatic interaction tensor contribution*
!         from the nuclei. In the case that the test charge coincide   *
!         with a nucleau we will remove that center.                   *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, April '95.                                      *
!***********************************************************************
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nAtm, nOrdOp
      Real*8 Chrg(nAtm), Coor(3,nAtm), ESIT((nOrdOp+1)*(nOrdOp+2)/2)

      Integer, Allocatable:: C_ESIT(:)
      Real*8 CoOp(3)
      Integer nTot, iPowr, iAtom, ix, iy, iz
      Real*8 Fact, x, y, z, r2, Thr, r, eix, eiy, eiz, temp
#ifdef _DEBUGPRINT_
      Integer n, nElem
!
!---- Statement function
!
      nElem(n)=(n+1)*(n+2)/2
#endif
!
!     Compute the nuclear contribution to the electrostatic interation
!     tensor, ESIT.
!
      ESIT(:)=Zero
!
      nTot=(nOrdOp+1)**6
      Call mma_allocate(C_ESIT,nTot,Label='ESIT')
      Call InitIA(C_ESIT,nOrdOp)
!
      iPowR=2*nOrdOp+1
      Fact=One
      If (nOrdOp.ge.1) Fact=-One
      Do iAtom = 1, nAtm
         x = CoOp(1) - Coor(1,iAtom)
         y = CoOp(2) - Coor(2,iAtom)
         z = CoOp(3) - Coor(3,iAtom)
         r2 = x**2 + y**2 + z**2
         Thr=1.0D-12
         If (r2.gt.Thr) Then
            r  = Chrg(iAtom)/Sqrt(r2)**iPowR
            Do ix = nOrdOp, 0, -1
               Do iy = nOrdOp-ix, 0, -1
                  iz = nOrdOp - ix - iy
                  If (ix.eq.0) Then
                     EIx=One
                  Else
                     EIx=x**ix
                  End If
                  If (iy.eq.0) Then
                     EIy=One
                  Else
                     EIy=y**iy
                  End If
                  If (iz.eq.0) Then
                     EIz=One
                  Else
                     EIz=z**iz
                  End If
                  temp=Fact*EIx*EIy*EIz*r
!
                  Call ContEI(C_ESIT,nOrdOp,ESIT,ix,iy,iz,temp)
!
               End Do
            End Do       ! End loop over cartesian combinations
         End If
      End Do             ! End loop over atoms
!
      Call mma_deallocate(C_ESIT)
!
#ifdef _DEBUGPRINT_
      Call RecPrt(' The Electrostatic Interaction'
     &                 //' Tensor',' ',ESIT,nElem(nOrdOp),1)
#endif
      Return
      End SubRoutine EFNuc

      Subroutine InitIA(I,mDeg)
      implicit None
      Integer mDeg
      Integer I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg)

      Integer n, a, b, c, p, q, r, new
!
! Purpose: Express the interaction tensor, defined by the
! quantities T(a,b,c) as functions of the vector R=(x,y,z),
! where a,b, and c are nonnegative integers and
! T(a,b,c)=((d/dx)**a)((d/dy)**b)((d/dz)**c) 1/R, in terms
! of a polynomial:
! T(a,b,c)=
!  (sum over p,q,r of I(a,b,c,p,q,r) x**p y**q z**r)/(R**(2*n+1)),
! where n=a+b+c.
! The polynomial coefficients are integers, and are 0 unless
! p+q+r=n.
! Author: PAM
!
!----- Statement function
!
!      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
! initialize:
      I(:,:,:,:,:,:)=0
      I(0,0,0,0,0,0)=1
      If (mDeg.gt.0) Then
         I(1,0,0,1,0,0)=-1
         I(0,1,0,0,1,0)=-1
         I(0,0,1,0,0,1)=-1
      End If
      do 100 n=2,mDeg
      do 101 a=0,n
      do 102 b=0,n-a
      c=n-a-b
      do 103 p=0,n
      do 104 q=0,n-p
      r=n-p-q
      new=0
      if(a.gt.0) then
        if(p.gt.0) new=(p-(2*n))*I(a-1,b,c,p-1,q,r)
        if(q.gt.1) new=new+(p+1)*I(a-1,b,c,p+1,q-2,r)
        if(r.gt.1) new=new+(p+1)*I(a-1,b,c,p+1,q,r-2)
      else if(b.gt.0) then
        if(q.gt.0) new=(q-(2*n))*I(a,b-1,c,p,q-1,r)
        if(r.gt.1) new=new+(q+1)*I(a,b-1,c,p,q+1,r-2)
        if(p.gt.1) new=new+(q+1)*I(a,b-1,c,p-2,q+1,r)
      else
        if(r.gt.0) new=(r-(2*n))*I(a,b,c-1,p,q,r-1)
        if(p.gt.1) new=new+(r+1)*I(a,b,c-1,p-2,q,r+1)
        if(q.gt.1) new=new+(r+1)*I(a,b,c-1,p,q-2,r+1)
      end if
      I(a,b,c,p,q,r)=new
 104  continue
 103  continue
 102  continue
 101  continue
 100  continue
!
! write out only elements with a>=b>=c. The others are obtained
! by index permutation.
! This restriction has been removed! (Roland Lindh)
!     n=mDeg
!     do 200 a=n,0,-1
!     do 200 b=n-a,0,-1
!     c=n-a-b
!     write(*,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,
!    &     Ind(n,a,c)
!     do 150 p=n,0,-1
!     do 150 q=n-p,0,-1
!     r=n-p-q
!     coef=I(a,b,c,p,q,r)
!     if(coef.eq.0) goto 150
!     write(*,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)')
!    &  coef,p,q,r,Ind(n,p,r)
!150  continue
!200  continue
!
      Return
      End Subroutine InitIA

      Subroutine ContEI(I,mDeg,ESIT,ix,iy,iz,temp)
      implicit None
      Integer mDeg, ix, iy, iz
      Integer I(0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg,0:mDeg)
      Real*8 ESIT((mDeg+1)*(mDeg+2)/2), Temp

      Integer n, ip, a, b, c
!
! Purpose: Express the interaction tensor, defined by the
! quantities T(a,b,c) as functions of the vector R=(x,y,z),
! where a,b, and c are nonnegative integers and
! T(a,b,c)=((d/dx)**a)((d/dy)**b)((d/dz)**c) 1/R, in terms
! of a polynomial:
! T(a,b,c)=
!  (sum over p,q,r of I(a,b,c,p,q,r) x**p y**q z**r)/(R**(2*n+1)),
! where n=a+b+c.
! The polynomial coefficients are integers, and are 0 unless
! p+q+r=n.
!
!
!----- Statement function
!
!      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
!
!
! write out only elements with a>=b>=c. The others are obtained
! by index permutation.
! This restriction has been removed! (Roland Lindh)
!     n=mDeg
!     do 200 a=n,0,-1
!     do 200 b=n-a,0,-1
!     c=n-a-b
!     write(*,'(5x,''T('',i1,'','',i1,'','',i1,'')='',i5)')a,b,c,
!    &     Ind(n,a,c)
!     do 150 p=n,0,-1
!     do 150 q=n-p,0,-1
!     r=n-p-q
!     coef=I(a,b,c,p,q,r)
!     if(coef.eq.0) goto 150
!     write(*,'(10x,i8,''*x**'',i1,'' *y**'',i1,'' *z**'',i1,i5)')
!    &  coef,p,q,r,Ind(n,p,r)
!150  continue
!200  continue
!
!     Write (*,*) ' Temp,ix,iy,iz=',temp,ix,iy,iz
      n=mDeg
      ip = 0
      Do a=n,0,-1
         Do b=n-a,0,-1
            c=n-a-b
            ip = ip + 1
!           Write (*,*) ip, I(a,b,c,ix,iy,iz)
            If (I(a,b,c,ix,iy,iz).ne.0)
     &      ESIT(ip)=ESIT(ip)+DBLE(I(a,b,c,ix,iy,iz))*temp
         End Do
      End Do
!
      Return
      End Subroutine ContEI
