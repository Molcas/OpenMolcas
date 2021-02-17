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
* Copyright (C) 1996,1999, Niclas Forsberg                             *
*               1996,1999, Anders Bernhardsson                         *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module PotKin
C!
C!  Contains:
C!    PotEnergy      (A,nMat,energy,grad,Hess,D3,D4,max_term)
C!    KinEnergy      (A,nMat,G,Gprime,Gdbleprime,max_term)
C!
C!  Written by:
C!    Niclas Forsberg & Anders Bernhardsson,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!    Dept. of Theoretical Chemistry, Lund University, 1999.
C!
C!-----------------------------------------------------------------------!
C!
Cvv       Private
C!
c       Contains


C!-----------------------------------------------------------------------!
C!
      Subroutine KinEnergy_drv(A,nMat,iCre,iAnn,G,Gprime,
     &       Gdbleprime,max_term,
     &       C,W,alpha1,alpha2,beta,r_diff,
     &   max_Ord, nOsc,nOscOld)
C!
C! Out : A
C! In  :
C!      Real G,Gprime,Gdbleprime,W,C,alpha1,alpha2,beta,r_diff
C!      integer nMat,iCre,iAnn,max_term
C!
C!  Purpose:
C!    Calculate matrix elements of kinetic energy terms.
C!
C!  Written by:
C!    Niclas Forsberg,Anders Bernhardsson
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!    Dept. of Theoretical Chemistry, Lund University, 1999.
C!
c       Use TabMod
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Parameter             ( Thrs = 1.0d-15)
      Real*8 A  (0:mdim1,0:ndim1)
      Integer nMat  (0:ndim1,ndim2)
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
c       Real*8 rdx(4)
      Real*8 alpha1  (nosc,nosc)
      Real*8 alpha2  (nosc,nosc)
      Real*8 beta  (nosc,nosc)
      Real*8 C   (nosc,nosc)
      Real*8 W   (nosc,nosc)
      Real*8 G     (nosc,nosc)
      Real*8 r_diff  (noscold)
      Real*8 Gprime  (nosc,nosc,nosc)
      Real*8 Gdbleprime  (nosc,nosc,nosc,nosc)
#include "WrkSpc.fh"

      nOsc2=nOsc*nOsc
      nOsc3=nOsc2*nOsc

      Call GetMem('Tempa','Allo','Real',ipTempa,nOsc2)
      Call GetMem('Tempb','Allo','Real',ipTempb,nOsc2)
      Call GetMem('r_temp','Allo','Real',ipr_temp,nOsc)
      Call GetMem('G_2','Allo','Real',ipG_2,nOsc2)
      Call GetMem('Temp','Allo','Real',ipTemp,nOsc2)
      Call GetMem('T1','Allo','Real',ipT1,nOsc)
      Call GetMem('T2','Allo','Real',ipT2,nOsc)
      Call GetMem('Temp1','Allo','Real',ipTemp1,nOsc3)
      Call GetMem('Temp2','Allo','Real',ipTemp2,nOsc3)
      Call GetMem('Temp3','Allo','Real',ipTemp3,nOsc3)
      Call GetMem('Temp4','Allo','Real',ipTemp4,nOsc3)


      Call KinEnergy(A,nMat,iCre,iAnn,G,Gprime,
     &       Gdbleprime,max_term,
     &       C,W,alpha1,alpha2,beta,r_diff,
     &   max_Ord, nOsc,nOscOld,
     &   Work(ipTempa), Work(ipTempb),Work(ipr_temp),
     &   Work(ipG_2),Work(ipTemp), Work(ipT1),Work(ipT2),
     &   Work(ipTemp1),Work(ipTemp2),Work(ipTemp3),Work(ipTemp4))

      Call GetMem('Tempa','Free','Real',ipTempa,nOsc2)
      Call GetMem('Tempb','Free','Real',ipTempb,nOsc2)
      Call GetMem('r_temp','Free','Real',ipr_temp,nOsc)
      Call GetMem('G_2','Free','Real',ipG_2,nOsc2)
      Call GetMem('Temp','Free','Real',ipTemp,nOsc2)
      Call GetMem('T1','Free','Real',ipT1,nOsc)
      Call GetMem('T2','Free','Real',ipT2,nOsc)
      Call GetMem('Temp1','Free','Real',ipTemp1,nOsc3)
      Call GetMem('Temp2','Free','Real',ipTemp2,nOsc3)
      Call GetMem('Temp3','Free','Real',ipTemp3,nOsc3)
      Call GetMem('Temp4','Free','Real',ipTemp4,nOsc3)

      End


C!-----------------------------------------------------------------------!
C!
      Subroutine KinEnergy(A,nMat,iCre,iAnn,G,Gprime,
     &       Gdbleprime,max_term,
     &       C,W,alpha1,alpha2,beta,r_diff,
     &   max_Ord, nOsc,nOscOld,Tempa,Tempb,r_temp,G_2,Temp,
     &   T1,T2,Temp1,Temp2,Temp3,Temp4)
C!
C! Out : A
C! In  :
C!      Real G,Gprime,Gdbleprime,W,C,alpha1,alpha2,beta,r_diff
C!      integer nMat,iCre,iAnn,max_term
C!
C!  Purpose:
C!    Calculate matrix elements of kinetic energy terms.
C!
C!  Written by:
C!    Niclas Forsberg,Anders Bernhardsson
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!    Dept. of Theoretical Chemistry, Lund University, 1999.
C!
c       Use TabMod
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Parameter             ( Thrs = 1.0d-15)
      Real*8 A  (0:mdim1,0:ndim1)
      Integer nMat  (0:ndim1,ndim2)
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
      Real*8 rdx(4)
      Real*8 alpha1  (nosc,nosc)
      Real*8 alpha2  (nosc,nosc)
      Real*8 beta  (nosc,nosc)
      Real*8 C   (nosc,nosc)
      Real*8 W   (nosc,nosc)
      Real*8 G     (nosc,nosc)
      Real*8 r_diff  (noscold)
      Real*8 Gprime  (nosc,nosc,nosc)
      Real*8 Gdbleprime  (nosc,nosc,nosc,nosc)
      Real*8  Temp1 (nosc,nosc,nosc)
      Real*8  Temp (nosc,nosc)
      Real*8  r_Temp (nosc)
      Real*8  Tempb (nosc,nosc)
      Real*8  Tempa  (nosc,nosc)
      Real*8  Temp2 (nosc,nosc,nosc)
      Real*8  G_2 (nosc,nosc)
      Real*8  Temp3 (nosc,nosc,nosc,nosc)
      Real*8  Temp4 (nosc,nosc,nosc,nosc)
      Real*8  T1(nosc),T2(nosc)
C!
C!---- Initialize.
      mPlus = max_Ord+1
      nOscSqr = nOsc**2
      r_norm = Dnrm2_(nOsc,r_diff,1)
      ran=1.0d0
      Call Dgesub(alpha1,nOsc,'N',alpha2,nOsc,'N',Tempa,
     &       nOsc,nOsc,nOsc)
      alpha_norm = Dnrm2_(nOscSqr,Tempa,1)
      If (( r_norm.gt.Thrs ).or.( alpha_norm.gt.Thrs )) Then
      Call Dgesub(alpha1,nOsc,'N',alpha2,nOsc,'N',Tempa,
     &    nOsc,nOsc,nOsc)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,Tempa,nOsc,
     &            W,nOsc,
     &            0.0d0,Tempb,nOsc)
      Call DGEMM_('N','N',
     &            nOsc,1,nOsc,
     &            1.0d0,beta,nOsc,
     &            r_diff,nOsc,
     &            0.0d0,r_temp,nOsc)
      call dscal_(nOsc,-2.0d0,r_temp,1)
      rdx(1)=+1.0d0
      rdx(2)=-1.0d0
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,G,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp,nosc)
      Call DGEMM_('T','T',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            C,nosc,
     &            0.0d0,G_2,nosc)
      call dscal_(nosc**2,-ran,G_2,1)
      Call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
      rdx(1)=-1.0d0
      rdx(2)=+1.0d0
      Call DGEMM_('T','T',
     &            nosc,nosc,nosc,
     &            1.0d0,G,nosc,
     &            C,nosc,
     &            0.0d0,Temp,nosc)
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            Tempb,nosc,
     &            0.0d0,G_2,nosc)
      call dscal_(nosc**2,-ran,G_2,1)
      Call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
      rdx(1)=+1.0d0
      rdx(2)=+1.0d0
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,G,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp,nosc)
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            Tempb,nosc,
     &            0.0d0,G_2,nosc)
      call dscal_(nosc**2,-1.0d0,G_2,1)
      Call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
      Call dGeMV_('N',nosc,nosc,1.0d0,G,nosc,r_temp,1,0.0d0,T1,1)
      Call dGeMV_('T',nosc,nosc,1.0d0,G,nosc,r_temp,1,1.0d0,T1,1)
      Call dGeMV_('N',nosc,nosc,1.0d0,C,nosc,T1,1,0.0d0,T2,1)
      call dscal_(nosc,-0.5d0*ran,T2,1)
      rdx(1)=-1.0d0
      Call Mul1(nmat,A,iCre,iAnn,T2,max_ord,nosc,rdx)
      Call dGeMV_('N',nosc,nosc,1.0d0,G,nosc,r_temp,1,0.0d0,T1,1)
      Call dGeMV_('T',nosc,nosc,1.0d0,G,nosc,r_temp,1,1.0d0,T1,1)
      Call dGeMV_('T',nosc,nosc,-0.5d0,tempb,nosc,T1,1,0.0d0,T2,1)
      rdx(1)=+1.0d0
      Call Mul1(nmat,A,iCre,iAnn,T2,max_ord,nosc,rdx)
      Call dGeMV_('N',nosc,nosc,1.0d0,G,nosc,r_temp,1,0.0d0,T1,1)
      r = Ddot_(nosc,T1,1,r_temp,1)
      call daxpy_(mplus,r,-0.5d0,0,A,mplus+1)
      If ( max_term.gt.2 ) Then
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Gprime,nosc,
     &            tempb,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','T',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            C,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=-1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp2(i,k,j)=-3.0d0*ran*temp1(i,j,k)
      End Do
      End Do
      End Do
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
      Call DGEMM_('T','T',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Gprime,nosc,
     &            C,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      rdx(1)=-1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp2(i,k,j)=-3.0d0*ran*temp1(i,j,k)
      End Do
      End Do
      End Do
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Gprime,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp2(i,k,j)=-3.0d0*temp1(i,j,k)
      End Do
      End Do
      End Do
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

c            temp=0.0d0
      call dcopy_(nOsc*nOsc,[0.0d0],0,temp,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp(k,j)=-ran*Gprime(i,j,k)*R_temp(i)+temp(k,j)
      End Do
      End Do
      End Do
      rdx(1)=1.0d0
      rdx(2)=-1.0d0
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            W,nosc,
     &            0.0d0,G_2,nosc)
      Call DGEMM_('T','T',
     &            nosc,nosc,nosc,
     &            1.0d0,G_2,nosc,
     &            C,nosc,
     &            0.0d0,Temp,nosc)
      Call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)
C!
c            temp=0.0d0
      call dcopy_(nOsc*nOsc,[0.0d0],0,temp,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp(k,j)=-Gprime(i,j,k)*R_temp(i)+temp(k,j)
      End Do
      End Do
      End Do
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            W,nosc,
     &            0.0d0,G_2,nosc)
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,G_2,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp,nosc)
      Call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)

c            temp=0.0d0
      call dcopy_(nOsc*nOsc,[0.0d0],0,temp,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp(i,k)=-ran*Gprime(i,j,k)*R_temp(j)+temp(i,k)
      End Do
      End Do
      End Do
      rdx(1)=-1.0d0
      rdx(2)=1.0d0
      Call DGEMM_('T','T',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            C,nosc,
     &            0.0d0,G_2,nosc)
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,G_2,nosc,
     &            W,nosc,
     &            0.0d0,Temp,nosc)
      Call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)

c            temp=0.0d0
      call dcopy_(nOsc*nOsc,[0.0d0],0,temp,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp(i,k)=-Gprime(i,j,k)*R_temp(j)+temp(i,k) !!
      End Do
      End Do
      End Do
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      Call DGEMM_('T','n',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            Tempb,nosc,
     &            0.0d0,G_2,nosc)
      Call DGEMM_('T','n',
     &            nosc,nosc,nosc,
     &            1.0d0,G_2,nosc,
     &            W,nosc,
     &            0.0d0,Temp,nosc)
      Call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)
C!
      call dcopy_(nOsc,[0.0d0],0,t1,1)
c            t1=0.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      T1(k)=-0.5d0*Gprime(i,j,k)*r_temp(i)*R_temp(j)+t1(k)
      End Do
      End Do
      End Do
      Call DGEMM_('T','n',
     &            1,nosc,nosc,
     &            1.0d0,T1,nosc,
     &            W,nosc,
     &            0.0d0,T2,1)
      rdx(1)=1.0d0
      Call Mul1(nmat,A,iCre,iAnn,T2,max_ord,nosc,rdx)
      End If
      If ( max_term.gt.3 ) Then
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Gdbleprime,nosc,
     &            tempb,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','T',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            C,nosc,
     &            0.0d0,Temp4,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp4,nosc,
     &            W,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            W,nosc,
     &            0.0d0,Temp4,nosc**3)
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=+1.0d0
      rdx(4)=-1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp3(i,k,l,j)=-ran*6.0d0*temp4(i,j,k,l)
      End Do
      End Do
      End Do
      End Do
      Call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)
      Call DGEMM_('T','T',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Gdbleprime,nosc,
     &            C,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp4,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp4,nosc,
     &            W,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            W,nosc,
     &            0.0d0,Temp4,nosc**3)
      rdx(1)=-1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      rdx(4)=1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp3(i,k,l,j)=-ran*6.0d0*temp4(i,j,k,l)
      End Do
      End Do
      End Do
      End Do
      Call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)
      Call DGEMM_('T','n',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Gdbleprime,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp4,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp4,nosc,
     &            W,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            W,nosc,
     &            0.0d0,Temp4,nosc**3)
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      rdx(4)=1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp3(i,k,l,j)=-6.0d0*temp4(i,j,k,l)
      End Do
      End Do
      End Do
      End Do
      Call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)

c            temp1=0.0d0
      call dcopy_(nOsc*nOsc*nOsc,[0.0d0],0,temp1,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp1(k,l,j)=-ran*1.5d0*Gdbleprime(i,j,k,l)*
     &          R_temp(i)+temp1(k,l,j)
      End Do
      End Do
      End Do
      End Do
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=-1.0d0
      Call DGEMM_('T','n',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            W,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','n',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','T',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            C,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
c            temp1=0.0d0
      call dcopy_(nOsc*nOsc*nOsc,[0.0d0],0,temp1,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp1(k,l,j)=-1.5d0*Gdbleprime(i,j,k,l)*
     &          R_temp(i)+temp1(k,l,j)
      End Do
      End Do
      End Do
      End Do
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      Call DGEMM_('T','n',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            W,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','n',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','n',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

c            temp1=0.0d0
      call dcopy_(nOsc*nOsc*nOsc,[0.0d0],0,temp1,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp1(i,k,l)=-ran*1.5d0*Gdbleprime(i,j,k,l)*
     &          R_temp(j)+temp1(i,k,l)
      End Do
      End Do
      End Do
      End Do
      rdx(1)=-1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      Call DGEMM_('T','T',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            C,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            W,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

c            temp1=0.0d0
      call dcopy_(nOsc*nOsc*nOsc,[0.0d0],0,temp1,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp1(i,k,l)=-1.5d0*Gdbleprime(i,j,k,l)*
     &          R_temp(j)+temp1(i,k,l)
      End Do
      End Do
      End Do
      End Do
      rdx(1)=1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      Call DGEMM_('T','n',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            Tempb,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            W,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)

      rdx(1)=1.0d0
      rdx(2)=1.0d0
c            temp=0.0d0
      call dcopy_(nOsc*nOsc,[0.0d0],0,temp,1)
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp(k,l)=-0.5d0*Gdbleprime(i,j,k,l)*r_temp(i)*
     &          R_temp(j)+temp(k,l)
      End Do
      End Do
      End Do
      End Do
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            W,nosc,
     &            0.0d0,G_2,nosc)
      Call DGEMM_('T','N',
     &            nosc,nosc,nosc,
     &            1.0d0,G_2,nosc,
     &            W,nosc,
     &            0.0d0,Temp,nosc)
      Call Mul2(nmat,A,iCre,iAnn,Temp,max_ord,nosc,rdx)
      End If
      End If
C!
C!
C!---- If higher terms than quadratic are used in the polynomial
C!     fit of the potential surface, then we have to use a
C!     Taylor expansion of the inverse mass tensor.
      rdx(1)=-1.0d0
      rdx(2)=-1.0d0
      Call DGEMM_('T','T',
     &            nosc,nosc,nosc,
     &            1.0d0,G,nosc,
     &            C,nosc,
     &            0.0d0,Temp,nosc)
      Call DGEMM_('T','T',
     &            nosc,nosc,nosc,
     &            1.0d0,Temp,nosc,
     &            C,nosc,
     &            0.0d0,G_2,nosc)
      Call DSCAL_(nosc**2,-1.0d0,G_2,1)
      Call Mul2(nmat,A,iCre,iAnn,G_2,max_ord,nosc,rdx)
      If ( max_term.gt.2 ) Then
      Call DGEMM_('T','T',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Gprime,nosc,
     &            C,nosc,
     &            0.0d0,Temp1,nosc**2)
      Call DGEMM_('T','T',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp1,nosc,
     &            C,nosc,
     &            0.0d0,Temp2,nosc**2)
      Call DGEMM_('T','N',
     &            nosc**2,nosc,nosc,
     &            1.0d0,Temp2,nosc,
     &            W,nosc,
     &            0.0d0,Temp1,nosc**2)
      rdx(1)=-1.0d0
      rdx(2)=1.0d0
      rdx(3)=-1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Temp2(i,k,j)=-3.0d0*temp1(i,j,k)
      End Do
      End Do
      End Do
      Call Mul3(nmat,A,iCre,iAnn,Temp2,max_ord,nosc,rdx)
      End If
      If ( max_term.gt.3 ) Then
      Call DGEMM_('T','T',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Gdbleprime,nosc,
     &            C,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','T',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            C,nosc,
     &            0.0d0,Temp4,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp4,nosc,
     &            W,nosc,
     &            0.0d0,Temp3,nosc**3)
      Call DGEMM_('T','N',
     &            nosc**3,nosc,nosc,
     &            1.0d0,Temp3,nosc,
     &            W,nosc,
     &            0.0d0,Temp4,nosc**3)

      rdx(1)=-1.0d0
      rdx(2)=1.0d0
      rdx(3)=1.0d0
      rdx(4)=-1.0d0
      Do i=1,nosc
      Do j=1,nosc
      Do k=1,nosc
      Do l=1,nosc
      Temp3(i,k,l,j)=-6.0d0*Temp4(i,j,k,l)
      End Do
      End Do
      End Do
      End Do
      Call Mul4(nmat,A,iCre,iAnn,Temp3,max_ord,nosc,rdx)
      End If
C!
      End
