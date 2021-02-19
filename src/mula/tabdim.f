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
* Copyright (C) 1995, Niclas Forsberg                                  *
*               1999, Anders Bernhardsson                              *
************************************************************************
      Subroutine TabDim_drv(nDim,nOsc,nTabDim)
C!
      Integer   nDim,nOsc
      Integer   nTabDim

#include "WrkSpc.fh"
      call GetMem('binomCoef','Allo','INTE',
     &  ipbinomCoef,(nDim+1)*nOsc)
      call TabDim(nDim,nOsc,nTabDim,iWork(ipbinomCoef))
      call GetMem('binomCoef','Free','INTE',
     &  ipbinomCoef,(nDim+1)*nOsc)
      Return
      End
C!-----------------------------------------------------------------------!
C!
      Subroutine TabDim(nDim,nOsc,nTabDim,binomCoef)
C!
C!  Purpose:
C!    Fill a nDim*nOsc matrix with binomial coefficients.
C!    This matrix is then used to calculate the dimension
C!    of a table containing excitations for nOsc oscillators.
C!    nDim is the maximum sum of the oscillator quanta.
C!
C!  Input:
C!    nDim    : Integer variable - maximun excitation.
C!    nOsc    : Integer variable - number of dimensions.
C!
C!  Result:
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Implicit None
      Integer   iRow,jCol,nDim,nOsc
      Integer   nTabDim
      Integer binomCoef(0:nDim,nOsc)
C!
C!---- Initialize.
C!
C!---- Calculate binomial coefficients.
      If ( nDim .gt. 0 ) Then
      Do jCol = 1,nOsc
      binomCoef(0,jCol) = 1
      End Do
      Do iRow = 0,nDim
      binomCoef(iRow,1) = 1
      End Do
      Do jCol = 2,nOsc
      Do iRow = 1,nDim
      binomCoef(iRow,jCol) = binomCoef(iRow-1,jCol)+
     &          binomCoef(iRow,jCol-1)
      End Do
      End Do
C!
C!---- Sum all elements of the nOsc'th column.
      nTabDim = 0
      Do iRow = 0,nDim
      nTabDim = nTabDim+binomCoef(iRow,nOsc)
      End Do
      Else
      nTabDim = 1
      End If
C!
C!
      End
C!
C!-----------------------------------------------------------------------!
C!
      Subroutine TabDim2_drv(nDim,nOsc,nTabDim)
C!
      Integer   nDim,nOsc
      Integer   nTabDim

#include "WrkSpc.fh"
      call GetMem('binomCoef','Allo','INTE',
     &  ipbinomCoef,(nDim+1)*nOsc)
      call TabDim(nDim,nOsc,nTabDim,iWork(ipbinomCoef))
      call GetMem('binomCoef','Free','INTE',
     &  ipbinomCoef,(nDim+1)*nOsc)
      Return
      End

C!-----------------------------------------------------------------------!
C!
      Integer Function iDetNr(iocc,graph2,nosc,m)
      implicit integer(a-z)
      integer iocc(nosc),graph2(0:m,0:m,nosc)
C!
C! iocc:occupation vector
C! graph2:vertex graph table
C! nosc number of nodes
C! m  number of quantas
C! number of determinants in with lower number of quantas.
C!
C! Calculate the index of occupation string iocc
C!
      iqnew=0
      iqold=0
      n=0
      Do i=1,nosc
      iqnew=iqnew+iocc(i)
      n=n+graph2(iqnew,iqold,i)
      iqold=iqnew
      End Do
C!
      idetnr=n
C!
      End

C! Muln  is a set of subroutines that calculates
C! <i|H|j> where H is a operator described
C! M_1,2..,n(a_1+a^t_1)*(a_2+a^t_2)...(a_n+a^t_n)
C! you can find n=1,2,3,4 in this file.
C!
C! nmat : Occupation of slater det.
C! F    : Output <i|A|j>
C! iCre,iAnn : Gives the resulting slater determinant if a^t (a) is acticting on SD
C! mat  : Matrix describing the operator expanded in Normal modes
C! m_ord: Number of slater determinants.
C!
C! Anders Bernhardsson Friday the 13th august 1999
C!
C!-----------------------------------------------------------------------!
C!
      Subroutine Mul1(nMat,F,iCre,iAnn,mat,m_ord,nosc,rdx)
C!
#include "dims.fh"
      Integer nMat (0:m_ord,nosc)
      Real*8 Mat(nOsc)
      Real*8 F( 0:mdim1,0:ndim1)
      Real*8 rdx ( 1)
      Real*8 sqr( 0:50 )
      Real*8 fact
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
C!
      Do i = 0,50
      sqr(i) = sqrt(dble(i)/2.0d0)
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Fact=sqr(nmat(iord,iosc))*rdx(1)
      F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*Mat(iOsc)
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iCre(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Fact=sqr(nmat(jord,iosc))
      F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*Mat(iOsc)
      End If
      End Do
      End Do


      End
C!
C!-----------------------------------------------------------------------!

C!-----------------------------------------------------------------------!
C!
      Subroutine Mul2(nMat,F,iCre,iAnn,mat,m_ord,nosc,rdx)
C!
#include "dims.fh"
      Integer nMat (0:m_ord,nosc)
      Real*8 Mat(nOsc,nOsc)
      Real*8 F(0:mdim1,0:ndim1)
      Real*8 sqr( 0:50 )
      Real*8 fact,r,rsym
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
      Real*8 rdx(2)
C!
      rsym=1/2.0d0
      Do i = 0,50
      sqr(i) = sqrt(dble(i)/2.0d0)
      End Do
      r=rdx(1)*rdx(2)
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iAnn(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*r*rsym
      F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*Mat(iOsc,josc)
      End If
      End Do
      End If
      End Do
      End Do
      If (r .ne. 0.0d0) Then
      Do iOrd = 0,m_Ord
      Do iOsc = 1,nOsc
      jOrd = iann(iord,iosc)
      If (jOrd .ge. 0) Then
      Do jOsc = 1,nOsc
      kOrd = iCre(jord,josc)
      If ( kOrd .ge. 0 ) Then
      r=rdx(1)*Mat(iosc,jOsc)+rdx(2)*Mat(josc,iOsc)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*
     &              r*rsym
      F(iOrd,kOrd) = F(iOrd,kOrd)+Fact
      End If
      End Do
      End If
      End Do
      End Do
      End If
      Do iOrd = 0,m_Ord
      Do iOsc = 1,nOsc
      jOrd = iCre(iord,iosc)
      If (jOrd .ge. 0) Then
      Do jOsc = 1,nOsc
      kOrd = iCre(jord,josc)
      If ( kOrd .ge. 0 ) Then
      Fact=sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*rsym
      F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*Mat(iosc,jOsc)
      End If
      End Do
      End If
      End Do
      End Do
      r=rdx(2)
      Do iOrd = 0,m_Ord
      Do iOsc = 1,nOsc
      F(iOrd,iOrd) = F(iOrd,iOrd)+Mat(iosc,iOsc)*r*rsym/2.0d0
      End Do
      End Do

C!
      End
C!
C!-----------------------------------------------------------------------!
C!-----------------------------------------------------------------------!
C!
      Subroutine Mul3(nMat,F,iCre,iAnn,mat,m_ord,nosc,rdx)
C!
#include "dims.fh"
      Integer nMat (0:m_ord,nosc)
      Real*8 Mat(nOsc,nOsc,nOsc)
      Real*8 F(0:mdim1,0:ndim1)
      Real*8 sqr( 0:50 )
      Real*8 rdx( 3 )
      Real*8 fact,r,rsym,relem
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
C!
      rsym=1.0d0/6.0d0
      Do i = 0,50
      sqr(i) = sqrt(dble(i)/2.0d0)
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iAnn(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iAnn(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      r=rdx(1)*rdx(2)*rdx(3)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*
     &          sqr(nmat(kord,kosc))*r*rsym
      F(iOrd,lOrd) = F(iOrd,lOrd)+Fact*Mat(iOsc,josc,kosc)
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do

      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iAnn(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iCre(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      r=rdx(1)*rdx(2)*Mat(iOsc,josc,kosc)+rdx(2)*rdx(3)*
     &          Mat(kOsc,iosc,josc)+rdx(3)*rdx(1)*Mat(iOsc,kosc,josc)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*
     &          sqr(nmat(lord,kosc))*r*rsym
      F(iOrd,lOrd) = F(iOrd,lOrd)+Fact
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do

      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iCre(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iCre(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      r=rdx(1)*mat(iosc,josc,kosc)+rdx(2)*
     &          mat(josc,iosc,kosc)+rdx(3)*mat(josc,kosc,iosc)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*
     &          sqr(nmat(lord,kosc))*rsym*r
      F(iOrd,lOrd) = F(iOrd,lOrd)+Fact
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do

      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iCre(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iCre(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iCre(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      Fact=sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*
     &          sqr(nmat(lord,kosc))*rsym
      F(iOrd,lOrd) = F(iOrd,lOrd)+Fact*Mat(iOsc,josc,kosc)
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do

      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iCre(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nosc
      Fact=sqr(nmat(jord,iosc))
      relem=Mat(josc,josc,iosc)*rdx(2)+rdx(3)*
     &         (mat(iosc,josc,josc)+mat(josc,iosc,josc))
      F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*relem*rsym*0.5d0
      End Do
      End If
      End Do
      End Do

      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iann(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nosc
      Fact=sqr(nmat(iord,iosc))
      relem=Mat(iosc,josc,josc)*rdx(1)*rdx(3)+rdx(2)*rdx(3)*
     &         (mat(josc,josc,iosc)+mat(josc,iosc,josc))
      F(iOrd,jOrd) = F(iOrd,jOrd)+Fact*relem*rsym*0.5d0
      End Do
      End If
      End Do
      End Do



      End
C!
C!-----------------------------------------------------------------------!
C!-----------------------------------------------------------------------!
C!
      Subroutine Mul4(nMat,F,iCre,iAnn,mat,m_ord,nosc,rdx)
C!
#include "dims.fh"
      Integer nMat (0:m_ord,nosc)
      Real*8 Mat(nOsc,nOsc,nOsc,nOsc)
      Real*8 F(0:mdim1,0:ndim1)
      Real*8 sqr( 0:50 )
      Real*8 rdx( 4 )
      Real*8 fact,r,rsym,relem
      Integer iAnn  (0:ndim1,ndim2)
      Integer iCre  (0:ndim1,ndim2)
C!
      Do i = 0,50
      sqr(i) = sqrt(dble(i)/2.0d0)
      End Do
      rsym=1.0d0/24.0d0
      r=Rdx(1)*rdx(2)*rdx(3)*rdx(4)
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iAnn(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iAnn(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      Do lOsc=1,nOsc
      mOrd = iAnn(lOrd,lOsc)
      If (mOrd .ge. 0) Then
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*
     &            sqr(nmat(kord,kosc))*sqr(nmat(lord,losc))*r*rsym
      F(iOrd,mOrd) = F(iOrd,mOrd)+Fact*
     &            Mat(iOsc,josc,kosc,losc)
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iAnn(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iAnn(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      Do lOsc=1,nOsc
      mOrd = iCre(lOrd,lOsc)
      If (mOrd .ge. 0) Then
      r=rdx(1)*rdx(2)*rdx(3)*mat(iosc,josc,kosc,losc)+
     &            rdx(2)*rdx(3)*rdx(4)*mat(losc,iosc,josc,kosc)+
     &                      rdx(3)*rdx(4)*rdx(1)*
     &     mat(iosc,losc,josc,kosc)+rdx(4)*rdx(1)*
     &     rdx(2)*mat(iosc,josc,losc,kosc)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*
     &            sqr(nmat(kord,kosc))*sqr(nmat(mord,losc))*r*rsym
      F(iOrd,mOrd) = F(iOrd,mOrd)+Fact
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iAnn(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iCre(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      Do lOsc=1,nOsc
      mOrd = iCre(lOrd,lOsc)
      If (mOrd .ge. 0) Then
      r=rdx(1)*rdx(2)*mat(iosc,josc,kosc,losc)+rdx(2)*rdx(3)*
     &   mat(kosc,iosc,josc,losc)+
     &             rdx(3)*rdx(4)*mat(kosc,losc,iosc,josc)+rdx(4)*
     &     rdx(1)*mat(iosc,kosc,losc,josc)+
     &                   rdx(1)*rdx(3)*mat(iosc,kosc,josc,losc)+
     &     rdx(4)*rdx(2)*mat(kosc,iosc,losc,josc)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*
     &            sqr(nmat(lord,kosc))*sqr(nmat(mord,losc))*r*rsym
      F(iOrd,mOrd) = F(iOrd,mOrd)+Fact
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iAnn(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iCre(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = icre(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      Do lOsc=1,nOsc
      mOrd = iCre(lOrd,lOsc)
      If (mOrd .ge. 0) Then
      r=rdx(1)*mat(iosc,josc,kosc,losc)+rdx(2)*
     &            mat(josc,iosc,kosc,losc)+rdx(3)*
     &     mat(josc,kosc,iosc,losc)
     &                       +rdx(4)*mat(josc,kosc,losc,iosc)
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*
     &            sqr(nmat(lord,kosc))*sqr(nmat(mord,losc))*r*rsym
      F(iOrd,mOrd) = F(iOrd,mOrd)+Fact
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iCre(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iCre(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      lOrd = iCre(kOrd,kOsc)
      If (lOrd .ge. 0) Then
      Do lOsc=1,nOsc
      mOrd = iCre(lOrd,lOsc)
      If (mOrd .ge. 0) Then
      Fact=sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*
     &            sqr(nmat(lord,kosc))*sqr(nmat(mord,losc))*rsym
      F(iOrd,mOrd) = F(iOrd,mOrd)+Fact*
     &            Mat(iOsc,josc,kosc,losc)
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iann(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iann(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(jord,josc))*rsym
      relem=rdx(1)*rdx(2)*rdx(4)*mat(iosc,josc,kosc,kosc)
      relem=relem+rdx(1)*rdx(3)*rdx(4)*
     &          (mat(iosc,kosc,kosc,josc)+mat(iosc,kosc,josc,kosc))
      relem=relem+rdx(2)*rdx(3)*rdx(4)*
     &          (mat(kosc,kosc,josc,iosc)+mat(kosc,iosc,kosc,josc)+
     &     mat(kosc,josc,iosc,kosc))
      F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*relem*0.5d0
      End Do
      End If
      End Do
      End If
      End Do
      End Do

      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = icre(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = icre(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      Fact=sqr(nmat(jord,iosc))*sqr(nmat(kord,josc))*rsym
      relem=rdx(2)*mat(kosc,kosc,iosc,josc)
      relem=relem+rdx(3)*(mat(iosc,kosc,kosc,josc)+
     &          mat(kosc,iosc,kosc,josc))
      relem=relem+rdx(4)*(mat(iosc,josc,kosc,kosc)+
     &          mat(iosc,kosc,josc,kosc)+mat(kosc,iosc,josc,kosc))
      F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*relem*0.5d0
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      jOrd = iann(iOrd,iOsc)
      If (jOrd .ge. 0) Then
      Do jOsc=1,nOsc
      kOrd = iCre(jOrd,jOsc)
      If ( kOrd .ge. 0 ) Then
      Do kOsc=1,nOsc
      Fact=sqr(nmat(iord,iosc))*sqr(nmat(kord,josc))*rsym
      relem=(Mat(jOsc,kosc,kosc,iosc)+
     &          Mat(kOsc,josc,kosc,iosc)+Mat(jOsc,kosc,iosc,kosc)+
     &     Mat(kOsc,josc,iosc,kosc))*rdx(3)*rdx(4)
      relem=relem+(Mat(kOsc,kosc,josc,iosc)+
     &          Mat(jOsc,iosc,kosc,kosc)+Mat(kOsc,iosc,josc,kosc))*
     &     rdx(2)*rdx(4)
      relem=relem+(Mat(kOsc,kosc,iosc,josc)+
     &          Mat(kOsc,iosc,kosc,josc))*rdx(2)*rdx(3)
      relem=relem+(Mat(iOsc,kosc,josc,kosc)+
     &          Mat(iOsc,josc,kosc,kosc))*rdx(1)*rdx(4)
      relem=relem+Mat(iOsc,kosc,kosc,josc)*rdx(1)*rdx(3)
      F(iOrd,kOrd) = F(iOrd,kOrd)+Fact*relem*0.5d0
      End Do
      End If
      End Do
      End If
      End Do
      End Do
      Do iOrd = 0,m_Ord
      Do iOsc=1,nOsc
      Do jOsc=1,nOsc
      relem=rdx(4)*rdx(3)*(mat(iosc,josc,josc,iosc)+
     &          mat(iosc,josc,iosc,josc))
      relem=relem+rdx(2)*rdx(4)*mat(iosc,iosc,josc,josc)
      F(iOrd,iOrd) = F(iOrd,iOrd)+relem*rsym*0.25d0
      End Do
      End Do
      End Do


      End
