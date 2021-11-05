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
! Copyright (C) 1996,1999, Niclas Forsberg                             *
!               1996,1999, Anders Bernhardsson                         *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
!       Module MatElMod
!!
!!  Contains:
!!    MatrixElements (L,U,FC00,Hmat,C,W,r_diff,mMat,nMat,
!!                    max_nOrd,energy,grad,Hess,D3,D4,G,
!!                    Gprime,Gdbleprime,alpha1,alpha2,beta,max_term)
!!    LSPotFit       (r01,energy1,grad1,Hess1,D3_1,D4_1,
!!                    r02,energy2,grad2,Hess2,D3_2,D4_2,
!!                    r00,energy0,r_min,FitCoef,mMat,stand_dev,max_err,
!!                    use_weight,max_term,pot)
!!    SetUpHmat      (energy0,r_min,ipow,var,yin,coef,r00,trfName,max_term,
!!                    C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,
!!                    max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,
!!                    nInc,mDec,nDec,L,U,H,S,G1,G2,G0,Gprime1,Gprime2,
!!                    Gprime0,Gdbleprime1,Gdbleprime2,Gdbleprime0,
!!                    C0,W0,det0,Mass,rOrigin)
!!
!!  Written by:
!!    Niclas Forsberg & Anders Bernhardsson,
!!    Dept. of Theoretical Chemistry, Lund University, 1996.
!!    Dept. of Theoretical Chemistry, Lund University, 1999.
!!
!!-----------------------------------------------------------------------!
!!

!vv       Private
!!
!       Contains

!!-----------------------------------------------------------------------!
!!
      Subroutine SetUpHmat(energy0,r_min,ipow,var,yin,r00,              &
     &       trfName,max_term,                                          &
     &        C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,          &
     &        max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,     &
     &        nInc,mDec,nDec,H,S,G1,G2,G0,Gprime1,Gprime2,Gprime0,      &
     &        Gdbleprime1,Gdbleprime2,Gdbleprime0,C0,W0,det0,Mass,      &
     &        rOrigin,Base,r0,r1,r2,nnsiz,                              &
     &  nterm,nvar,ndata,nosc,ndimtot,numofat)
!!
!!  Purpose:
!!    Set up Hamilton matrix.
!!
!!  Input:
!!
!!   Energy
!!   r_min
!!   ipow
!!   var
!!   yin
!!  coeff
!!   r00
!!   trfname
!!   Max_term
!!   C1,W1,det1,r01
!!   C2,W2,det2,r02
!!   H S
!!   G1 G2 G0 1' G2' g0' g0'' g1'' g2''
!!   C0 W0 det0
!!   mass rorigin
!!
!!
!!  The expansion point for lspotfit is r00!!!!!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Integer ipow  (nterm,nvar)
      Real*8 var (ndata,nvar)
      Real*8 yin (ndata)
      Real*8 r01(nOsc),r02(nOsc),r00(nOsc),r_min(nOsc)
      Real*8 r1(nOsc),r2(nOsc),r0(nOsc)
      Real*8 rOrigin(nOsc)
      Character*80 trfName (nvar)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),                               &
     &  W1(nOsc,nOsc),W2(nOsc,nOsc),C0(nOsc,nOsc),W0(nOsc,nOsc)
      Real*8 G1(nOsc,nOsc),G2(nOsc,nOsc),G0(nOsc,nOsc)
      Real*8 Gprime1(ngdim,ngdim,ngdim)
      Real*8 Gprime2(ngdim,ngdim,ngdim)
      Real*8 Gprime0(ngdim,ngdim,ngdim)
      Real*8 Gdbleprime1(ngdim,ngdim,ngdim,ngdim)
      Real*8 Gdbleprime2(ngdim,ngdim,ngdim,ngdim)
      Real*8 Gdbleprime0(ngdim,ngdim,ngdim,ngdim)
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),                  &
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),                  &
     &  nDec(0:ndim1,ndim2)
      Real*8 H(ndimtot,ndimtot),S (ndimtot,ndimtot)
      Real*8 Mass (numOfAt)
      Real*8 Base (nosc,nosc)
!       Logical   find_minimum,use_weight
!       Real*8    stand_dev,max_err
!       Logical   pot

#include "WrkSpc.fh"

      nOscOld = nOsc

      Call GetMem('Hij','Allo','Real',                                  &
     &  ipHij,(max_mOrd+1)*(max_nOrd+1))
      Call GetMem('Hijtrans','Allo','Real',                             &
     &  ipHijTrans,(max_nOrd+1)*(max_mOrd+1))
      Call GetMem('Sij','Allo','Real',                                  &
     &  ipSij,(max_mOrd+1)*(max_nOrd+1))
      Call GetMem('Sijtrans','Allo','Real',                             &
     &  ipSijTrans,(max_nOrd+1)*(max_mOrd+1))
      Call GetMem('r0vec','Allo','Real',ipr0vec,nOscOld)
      Call GetMem('r_diff','Allo','Real',ipr_diff,nOscOld)
      Call GetMem('C','Allo','Real',ipC,nOsc*nOsc)
      Call GetMem('W','Allo','Real',ipW,nOsc*nOsc)
      Call GetMem('grad','Allo','Real',ipgrad,nOscOld)
      Call GetMem('grad1','Allo','Real',ipgrad1,nOscOld)
      Call GetMem('grad2','Allo','Real',ipgrad2,nOscOld)
      Call GetMem('Hess','Allo','Real',ipHess,nOscOld*nOscOld)

      Call GetMem('Hess1','Allo','Real',ipHess1,nOscOld*nOscOld)
      Call GetMem('Hess2','Allo','Real',ipHess2,nOscOld*nOscOld)

      Call GetMem('D3','Allo','Real',ipD3,nOscOld*nOscOld*nOscOld)

      Call GetMem('D3_1','Allo','Real',                                 &
     &  ipD3_1,nOscOld*nOscOld*nOscOld)
      Call GetMem('D3_2','Allo','Real',                                 &
     &  ipD3_2,nOscOld*nOscOld*nOscOld)
      Call GetMem('D4','Allo','Real',                                   &
     &  ipD4,nOscOld*nOscOld*nOscOld*nOscOld)
      Call GetMem('D4_1','Allo','Real',                                 &
     &  ipD4_1,nOscOld*nOscOld*nOscOld*nOscOld)
      Call GetMem('D4_2','Allo','Real',                                 &
     &  ipD4_2,nOscOld*nOscOld*nOscOld*nOscOld)
      Call GetMem('GTemp','Allo','Real',ipGtemp,nOsc*nOsc)
      Call GetMem('GprimeTemp','Allo','Real',                           &
     &  ipGprimetemp,nOsc*nOsc*nOsc)
      Call GetMem('GdbleprimeTemp','Allo','Real',                       &
     &  ipGdbleprimetemp,nOsc*nOsc*nOsc*nOsc)
      Call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
      Call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)
      Call GetMem('L','Allo','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('U','Allo','Real',ipU,(max_nOrd+1)*(max_nOrd+1))
!!

      Call SetUpHmat_a(energy0,r_min,ipow,var,yin,r00,                  &
     &       trfName,max_term,                                          &
     &        C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,          &
     &        max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,     &
     &        nInc,mDec,nDec,H,S,G1,G2,G0,Gprime1,Gprime2,Gprime0,      &
     &        Gdbleprime1,Gdbleprime2,Gdbleprime0,C0,W0,det0,Mass,      &
     &        rOrigin,Base,r0,r1,r2,nnsiz,                              &
     &  nterm,nvar,ndata,nosc,ndimtot,numofat,nOscOld,                  &
     & Work(ipL),Work(ipU),Work(ipr_diff),Work(ipr0vec),                &
     & Work(ipC),Work(ipW),                                             &
     & Work(ipgrad),Work(ipHess),Work(ipD3),Work(ipD4),                 &
     & Work(ipgrad1),Work(ipHess1),Work(ipD3_1),Work(ipD4_1),           &
     & Work(ipgrad2),Work(ipHess2),Work(ipD3_2),Work(ipD4_2),           &
     & Work(ipHij),Work(ipHijTrans),Work(ipSij),Work(ipSijTrans),       &
     & Work(ipGtemp),Work(ipGprimeTemp),Work(ipGdbleprimeTemp),         &
     & Work(ipalpha1),Work(ipalpha2),Work(ipbeta))


      Call GetMem('Hij','Free','Real',                                  &
     &  ipHij,(max_mOrd+1)*(max_nOrd+1))
      Call GetMem('Hijtrans','Free','Real',                             &
     &  ipHijTrans,(max_nOrd+1)*(max_mOrd+1))
      Call GetMem('Sij','Free','Real',                                  &
     &  ipSij,(max_mOrd+1)*(max_nOrd+1))
      Call GetMem('Sijtrans','Free','Real',                             &
     &  ipSijTrans,(max_nOrd+1)*(max_mOrd+1))
      Call GetMem('r0vec','Free','Real',ipr0vec,nOscOld)
      Call GetMem('r_diff','Free','Real',ipr_diff,nOscOld)
      Call GetMem('C','Free','Real',ipC,nOsc*nOsc)
      Call GetMem('W','Free','Real',ipW,nOsc*nOsc)
      Call GetMem('grad','Free','Real',ipgrad,nOscOld)
      Call GetMem('grad1','Free','Real',ipgrad1,nOscOld)
      Call GetMem('grad2','Free','Real',ipgrad2,nOscOld)
      Call GetMem('Hess','Free','Real',ipHess,nOscOld*nOscOld)

      Call GetMem('Hess1','Free','Real',ipHess1,nOscOld*nOscOld)
      Call GetMem('Hess2','Free','Real',ipHess2,nOscOld*nOscOld)

      Call GetMem('D3','Free','Real',ipD3,nOscOld*nOscOld*nOscOld)

      Call GetMem('D3_1','Free','Real',                                 &
     &  ipD3_1,nOscOld*nOscOld*nOscOld)
      Call GetMem('D3_2','Free','Real',                                 &
     &  ipD3_2,nOscOld*nOscOld*nOscOld)
      Call GetMem('D4','Free','Real',                                   &
     &  ipD4,nOscOld*nOscOld*nOscOld*nOscOld)
      Call GetMem('D4_1','Free','Real',                                 &
     &  ipD4_1,nOscOld*nOscOld*nOscOld*nOscOld)
      Call GetMem('D4_2','Free','Real',                                 &
     &  ipD4_2,nOscOld*nOscOld*nOscOld*nOscOld)
      Call GetMem('GTemp','Free','Real',ipGtemp,nOsc*nOsc)
      Call GetMem('GprimeTemp','Free','Real',                           &
     &  ipGprimetemp,nOsc*nOsc*nOsc)
      Call GetMem('GdbleprimeTemp','Free','Real',                       &
     &  ipGdbleprimetemp,nOsc*nOsc*nOsc*nOsc)
      Call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
      Call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)
      Call GetMem('L','Free','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('U','Free','Real',ipU,(max_nOrd+1)*(max_nOrd+1))
!!
      End


!!-----------------------------------------------------------------------!
!!
      Subroutine SetUpHmat_a(energy0,r_min,ipow,var,yin,r00,            &
     &       trfName,max_term,                                          &
     &        C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,          &
     &        max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,     &
     &        nInc,mDec,nDec,H,S,G1,G2,G0,Gprime1,Gprime2,Gprime0,      &
     &        Gdbleprime1,Gdbleprime2,Gdbleprime0,C0,W0,det0,Mass,      &
     &        rOrigin,Base,r0,r1,r2,nnsiz,                              &
     &  nterm,nvar,ndata,nosc,ndimtot,numofat,nOscOld,                  &
     & L,U,r_diff,r0vec,C,W,                                            &
     & grad,Hess,D3,D4,grad1,Hess1,D3_1,D4_1,                           &
     & grad2,Hess2,D3_2,D4_2,Hij,HijTrans,Sij,SijTrans,                 &
     & Gtemp,GprimeTemp,GdbleprimeTemp,alpha1,alpha2,beta)
!!
!!  Purpose:
!!    Set up Hamilton matrix.
!!
!!  Input:
!!
!!   Energy
!!   r_min
!!   ipow
!!   var
!!   yin
!!  coeff
!!   r00
!!   trfname
!!   Max_term
!!   C1,W1,det1,r01
!!   C2,W2,det2,r02
!!   H S
!!   G1 G2 G0 1' G2' g0' g0'' g1'' g2''
!!   C0 W0 det0
!!   mass rorigin
!!
!!
!!  The expansion point for lspotfit is r00!!!!!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
      Integer ipow  (nterm,nvar)
      Real*8 var (ndata,nvar)
      Real*8 yin (ndata)
      Real*8 r01(nOsc),r02(nOsc),r00(nOsc),r_min(nOsc)
      Real*8 r1(nOsc),r2(nOsc),r0(nOsc)
      Real*8 rOrigin(nOsc)
      Character*80 trfName (nvar)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),                               &
     &  W1(nOsc,nOsc),W2(nOsc,nOsc),C0(nOsc,nOsc),W0(nOsc,nOsc)
      Real*8 G1(nOsc,nOsc),G2(nOsc,nOsc),G0(nOsc,nOsc)
      Real*8 Gprime1(ngdim,ngdim,ngdim)
      Real*8 Gprime2(ngdim,ngdim,ngdim)
      Real*8 Gprime0(ngdim,ngdim,ngdim)
      Real*8 Gdbleprime1(ngdim,ngdim,ngdim,ngdim)
      Real*8 Gdbleprime2(ngdim,ngdim,ngdim,ngdim)
      Real*8 Gdbleprime0(ngdim,ngdim,ngdim,ngdim)
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),                  &
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),                  &
     &  nDec(0:ndim1,ndim2)
      Real*8 H(ndimtot,ndimtot),S (ndimtot,ndimtot)
      Real*8 Mass (numOfAt)
      Real*8 Base (nosc,nosc)
      Logical   find_minimum,use_weight
      Real*8    stand_dev,max_err
      Logical   pot
      Real*8  Hij(0:max_mOrd,0:max_nOrd)
      Real*8  HijTrans(0:max_nOrd,0:max_mOrd)
      Real*8  Sij(0:max_mOrd,0:max_nOrd)
      Real*8  SijTrans(0:max_nOrd,0:max_mOrd)
      Real*8  r0vec(nOscOld)
      Real*8  r_diff(nOscOld)
      Real*8  C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8  grad(nOscOld)
      Real*8  grad1(nOscOld),grad2(nOscOld)
      Real*8  Hess(nOscOld,nOscOld)
      Real*8  Hess1(nOscOld,nOscOld),Hess2(nOscOld,nOscOld)
      Dimension  D3(nOscOld,nOscOld,nOscOld)
      Dimension  D3_1(nOscOld,nOscOld,nOscOld)
      Dimension  D3_2(nOscOld,nOscOld,nOscOld)
      Real*8  D4(nOscOld,nOscOld,nOscOld,nOscOld)
      Dimension  D4_1(nOscOld,nOscOld,nOscOld,nOscOld)
      Dimension  D4_2(nOscOld,nOscOld,nOscOld,nOscOld)
      Real*8  Gtemp(nOsc,nOsc)
      Real*8  GprimeTemp(nOsc,nOsc,nOsc)
      Real*8  GdbleprimeTemp(nOsc,nOsc,nOsc,nOsc)
      Real*8  alpha1(nOsc,nOsc),alpha2(nOsc,nOsc),beta(nOsc,nOsc)
      Real*8  L(0:max_mOrd,0:max_mOrd)
      Real*8  U(0:max_nOrd,0:max_nOrd)

#include "WrkSpc.fh"
!!
!!---- Initialize.
      l_r2=nOsc
      find_minimum = .false.
      use_weight = .false.

      call GetMem('Coef','Allo','Real',ipCoef,nPolyTerm)
!!
!!---- Fit polynomial and calculate energy, gradient, Hessian and third
!!     and fourth order force constants around r01 and r02.
      Call PotFit(nterm,nvar,ndata,ipow,var,yin,                        &
     &  Work(ipcoef),r1,nOscOld,energy1,                                &
     &       grad1,Hess1,D3_1,D4_1,trfName,                             &
     &       stand_dev,max_err,find_minimum,max_term,use_weight,        &
     &     nOscOld,nOscOld,nOscOld)
      Call PotFit(nterm,nvar,ndata,ipow,var,yin,                        &
     &  Work(ipcoef),r2,l_r2,energy2,                                   &
     &       grad2,Hess2,D3_2,D4_2,trfName,                             &
     &       stand_dev,max_err,find_minimum,max_term,use_weight,        &
     &     nOscOld,nOscOld,nOscOld)
      call GetMem('Coef','Allo','Real',ipCoef,nPolyTerm)
!!
!!---- Perform a least squares fit.
      Call TabDim_drv(max_term,nOsc,numCoef)
      Call GetMem('FitCoef','Allo','Real',ipFitCoef,numCoef)

      Call GetMem('jPow','Allo','Inte',ipjPow,numCoef*nOsc)
      pot = .true.
      Call LSPotFit(r1,energy1,grad1,Hess1,D3_1,D4_1,                   &
     &    r2,energy2,grad2,Hess2,D3_2,D4_2,                             &
     &    r0,energy0,r_min,Work(ipFitCoef),iWork(ipjPow),               &
     &    stand_dev,max_err,                                            &
     &    use_weight,max_term,pot,nosc,numcoef)
      Call GetMem('kPow','Allo','Inte',ipkPow,numCoef*nOsc)
!       kPow = jPow
      do iv=0,numCoef*nOsc
      iWork(ipkPow+iv)=iWork(ipjPow+iv)
      enddo
!       call dcopy_(numCoef*nOsc,iWork(ipjPow),1,iWork(ipkPow),1)
!!
!!
!!---- For each of the centers:
!!     - Call Franck-Condon routine.
!!     - Calculate matrix elements.
!!
!!---- Block 11.
      l_C1=nOsc
      Call Calc_r00(C1,C1,W1,W1,C,W,alpha1,alpha2,r0vec,r01,            &
     &       r01,det0,det1,det1,FC00,l_C1)
      Call FCval(C1,W1,det1,r01,C1,W1,det1,r01,Sij,                     &
     &       max_mOrd,max_nOrd,max_nOrd2,                               &
     &       max_mInc,max_nInc,max_nInc2,mMat,                          &
     &       nMat,mInc,nInc,mDec,nDec,                                  &
     &       C,W,det1,r0vec,L,U,FC00,alpha1,alpha2,beta,                &
     &       l_C1,nnsiz)
      do iv=1,nOsc
      r0vec(iv)  = r1(iv)-r0(iv)
      enddo
      call funcval(r0vec,Work(ipFitCoef),iWork(ipkPow),                 &
     &  energy,nterm,nvar)
      call gradient(r0vec,Work(ipFitCoef),iWork(ipkPow),                &
     &  grad,nterm,nvar)
      call Hessian(r0vec,Work(ipFitCoef),iWork(ipkPow),                 &
     &  Hess,nterm,nvar)
      call thirdDer(r0vec,Work(ipFitCoef),iWork(ipkPow),                &
     &  D3,nterm,nvar)
      call fourthDer(r0vec,Work(ipFitCoef),iWork(ipkPow),               &
     &  D4,nterm,nvar)
      energy = energy-energy0
      call dcopy_(nOsc,[0.0d0],0,r0vec,1)
!       r0vec = 0.0d0
      call dcopy_(nOsc,[0.0d0],0,r_diff,1)
!              r_diff = 0.0d0
!       Gtemp = G1
      call dcopy_(nOsc*nOsc,G1,1,Gtemp,1)
!              GprimeTemp = Gprime1
      call dcopy_(nOsc*nOsc*nOsc,Gprime1,1,GprimeTemp,1)

!              GdbleprimeTemp = Gdbleprime1
      call dcopy_(nOsc**4,Gdbleprime1,1,GdbleprimeTemp,1)
      Call MatrixElements(L,U,FC00,Hij,C,W,r_diff,mMat,nMat,            &
     &       ninc,ndec,max_nOrd, max_mOrd,nOsc,                         &
     &       energy,grad,Hess,D3,D4,Gtemp,GprimeTemp,                   &
     &       GdbleprimeTemp,alpha1,alpha2,beta,max_term,Base)

!       H(1:max_mOrd+1,1:max_nOrd+1) = Hij
!       S(1:max_mOrd+1,1:max_nOrd+1) = Sij
      call dcopy_((max_mOrd+1)*(max_nOrd+1),Hij,1,H,1)
      call dcopy_((max_mOrd+1)*(max_nOrd+1),Sij,1,S,1)
!!
!!---- Block 22.
      l_C2=nOsc
      Call Calc_r00(C2,C2,W2,W2,C,W,alpha1,alpha2,r0vec,r02,            &
     &       r02,det0,det2,det2,FC00,l_C2)
      Call FCval(C2,W2,det2,r02,C2,W2,det2,r02,Sij,max_mOrd,            &
     &       max_nOrd,max_nOrd2,                                        &
     &       max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,           &
     &       mDec,nDec,                                                 &
     &       C,W,det2,r0vec,L,U,FC00,alpha1,alpha2,beta,                &
     &       l_C2,nnsiz)
      do iv=1,nOsc
      r0vec(iv)  = r2(iv)-r0(iv)
      enddo
      call funcval(r0vec,Work(ipFitCoef),iWork(ipkPow),                 &
     &  energy,nterm,nvar)
      call gradient(r0vec,Work(ipFitCoef),iWork(ipkPow),                &
     &  grad,nterm,nvar)
      call Hessian(r0vec,Work(ipFitCoef),iWork(ipkPow),                 &
     &  Hess,nterm,nvar)
      call thirdDer(r0vec,Work(ipFitCoef),iWork(ipkPow),                &
     &  D3,nterm,nvar)
      call fourthDer(r0vec,Work(ipFitCoef),iWork(ipkPow),               &
     &  D4,nterm,nvar)
      energy = energy-energy0
!       r0vec = 0.0d0
!              r_diff = 0.0d0
      call dcopy_(nOsc,[0.0d0],0,r0vec,1)
      call dcopy_(nOsc,[0.0d0],0,r_diff,1)
!       Gtemp = G2
      call dcopy_(nOsc*nOsc,G2,1,Gtemp,1)
!              GprimeTemp = Gprime2
      call dcopy_(nOsc*nOsc*nOsc,Gprime2,1,GprimeTemp,1)
!              GdbleprimeTemp = Gdbleprime2
      call dcopy_(nOsc*nOsc*nOsc*nOsc,Gprime2,1,GprimeTemp,1)
      Call MatrixElements(L,U,FC00,Hij,C,W,r_diff,mMat,nMat,ninc,       &
     &       ndec,max_nOrd, max_mOrd,nOsc,                              &
     &        energy,grad,Hess,D3,D4,Gtemp,GprimeTemp,                  &
     &        GdbleprimeTemp,alpha1,alpha2,beta,max_term,Base)
!       H(max_mOrd+2:2*max_mOrd+2,max_nOrd+2:2*max_nOrd+2) = Hij
!       S(max_mOrd+2:2*max_mOrd+2,max_nOrd+2:2*max_nOrd+2) = Sij
      do iv=max_mOrd+2,2*max_mOrd+2
      do jv=max_nOrd+2,2*max_nOrd+2
      H(iv,jv)=Hij(iv-max_mOrd-2,jv-max_nOrd-2)
      S(iv,jv)=Sij(iv-max_mOrd-2,jv-max_nOrd-2)
      enddo
      enddo
!!
!!---- Block 12 and 21.
      l_C1=nOsc
      Call Calc_r00(C1,C2,W1,W2,C,W,alpha1,alpha2,r0vec,r01,r02,        &
     &       det0,det1,det2,FC00,l_C1)
      Call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,max_mOrd,            &
     &       max_nOrd,max_nOrd2,                                        &
     &       max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,                &
     &       nInc,mDec,nDec,                                            &
     &         C,W,det0,r0vec,L,U,FC00,alpha1,alpha2,beta,              &
     &       l_C1,nnsiz)
!       r0vec  = r0-r0
      call dcopy_(nOsc,[0.0d0],0, r0vec,1)
      call funcval(r0vec,Work(ipFitCoef),iWork(ipkPow),                 &
     &  energy,nterm,nvar)
      call gradient(r0vec,Work(ipFitCoef),iWork(ipkPow),                &
     &  grad,nterm,nvar)
      Call Hessian(r0vec,Work(ipFitCoef),iWork(ipkPow),                 &
     &  Hess,nterm,nvar)
      call thirdDer(r0vec,Work(ipFitCoef),iWork(ipkPow),                &
     &  D3,nterm,nvar)
      call fourthDer(r0vec,Work(ipFitCoef),iWork(ipkPow),               &
     &  D4,nterm,nvar)
      energy = energy-energy0
      call dcopy_(nOsc,[0.0d0],0, r0vec,1)
!       r0vec = 0.0d0
      do iv=1,nOsc
      r_diff(iv) = r01(iv)-r02(iv)
      enddo
!       Gtemp = G0
      call dcopy_(nOsc*nOsc,G0,1, Gtemp,1)
!              GprimeTemp = Gprime0
      call dcopy_(nOsc*nOsc*nOsc,Gprime0,1, GprimeTemp,1)
!              GdbleprimeTemp = Gdbleprime0
      call dcopy_(nOsc*nOsc*nOsc*nOsc,Gdbleprime0,1,GdbleprimeTemp,1)
      Call MatrixElements(L,U,FC00,Hij,C,W,r_diff,mMat,nMat,nInc,       &
     &       nDec,max_nOrd, max_mOrd,nOsc,                              &
     &       energy,grad,Hess,D3,D4,Gtemp,GprimeTemp,                   &
     &       GdbleprimeTemp,alpha1,alpha2,beta,max_term,Base)
!       H(1:max_mOrd+1,max_nOrd+2:2*max_nOrd+2) = Hij
!       S(1:max_mOrd+1,max_nOrd+2:2*max_nOrd+2) = Sij

      do iv=1,max_mOrd+1
      do jv=max_nOrd+2,2*max_nOrd+2
      H(iv,jv)=Hij(iv-1,jv-max_nOrd-2)
      S(iv,jv)=Sij(iv-1,jv-max_nOrd-2)
      enddo
      enddo


      Do i = 0,max_mOrd
      Do j = 0,max_nOrd
      HijTrans(j,i) = Hij(i,j)
      SijTrans(j,i) = Sij(i,j)
      End Do
      End Do
!       H(max_mOrd+2:2*max_mOrd+2,1:max_nOrd+1) = HijTrans
!       S(max_mOrd+2:2*max_mOrd+2,1:max_nOrd+1) = SijTrans

      do iv=max_mOrd+2,2*max_mOrd+2
      do jv=1,max_nOrd+1
      H(iv,jv)=HijTrans(iv-max_mOrd-2,jv-1)
      S(iv,jv)=SijTrans(iv-max_mOrd-2,jv-1)
      enddo
      enddo


!!
      Call GetMem('jPow','Free','Inte',ipjPow,numCoef*nOsc)
      Call GetMem('kPow','Free','Inte',ipkPow,numCoef*nOsc)
      Call GetMem('FitCoef','Free','Real',ipFitCoef,numCoef)
!!
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(r00)
         Call Unused_real_array(Gdbleprime2)
         Call Unused_real_array(C0)
         Call Unused_real_array(W0)
         Call Unused_real_array(Mass)
         Call Unused_real_array(rOrigin)
      End If
      End
