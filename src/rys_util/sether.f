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
* Copyright (C) 1992, Per Ake Malmqvist                                *
*               1992, Roland Lindh                                     *
************************************************************************
       Subroutine SetHer(nDiff)
************************************************************************
*                                                                      *
* Object: to setup the roots and weights of the Hermite polynomials    *
*         for the evaluation of one electron integrals.                *
*                                                                      *
*    Authors: Per-Ake Malmqvist and Roland Lindh,                      *
*             March 1992.                                              *
************************************************************************
      use Her_RW
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "status.fh"
      Real*8, Dimension(:), Allocatable :: Beta, BInv, Herm
*
      If (nPrp.gt.nPrpMx) Then
         Write (6,*) 'nPrp, nPrpMx=',nPrp, nPrpMx
         Call WarningMessage(2,'SetHer: nPrp too large!')
         Call Abend()
      End If
*
*
*     1) Hermite-Gauss
*     2) Rys-Gauss (asymtotic formula)
*
      n_1111 = (2*iAngMx+nPrp+2+nDiff)/2
      n_2222 = 4*iAngMx+2+nDiff
*
      If (Allocated(HerR) .and. Max(n_1111,n_2222).le.MaxHer) Then
         Return
      Else If (Allocated(HerR)) Then
         Call Free_HerRW()
      End If
      MaxHer = Max(n_1111,n_2222)
      Call mma_allocate(iHerR,MaxHer,label='iHerR')
      Call mma_allocate(iHerW,MaxHer,label='iHerW')
*
*     Set up square of roots and weights for Hermite polynomials
*
      nMem = (MaxHer*MaxHer+MaxHer)/2
      Call mma_Allocate(HerR,nMem,label='HerR')
      iHerR(1)=1
      Call dCopy_(nMem,0.0d0,0,HerR,1)
      Call mma_allocate(HerW,nMem,label='HerW')
      iHerW(1)=1
      Call dCopy_(nMem,0.0d0,0,HerW,1)
      Call mma_allocate(Beta,MaxHer,label='Beta')
      Call dCopy_(MaxHer,0.0d0,0,Beta,1)
      Call mma_allocate(BInv,MaxHer,label='BInv')
      Call dCopy_(MaxHer,0.0d0,0,BInv,1)
      Call mma_allocate(Herm,MaxHer+1,label='Herm')
      Call dCopy_(MaxHer+1,0.0d0,0,Herm,1)
      DO 10 K=1,MaxHer
        b_1111 = HALF*DBLE(K)
        B=SQRT(b_1111)
        Beta(K)=B
        BInv(K)=1.0d0/B
  10  CONTINUE
      HerR(iHerR(1))=0.0d0
      HerR(iHerR(1)+2)=SQRT(HALF)
      HerR(iHerR(1)+1)=-HerR(iHerR(1)+2)
      HerW(iHerW(1))=SQRT(PI)
      HerW(iHerW(1)+1)=HerW(iHerW(1))*HALF
      HerW(iHerW(1)+2)=HerW(iHerW(1)+1)
      Herm(1)=1.0d0/SQRT( HerW(iHerW(1)) )
      Do 11 iHer = 2, MaxHer
        i_1111 = (iHer*iHer-iHer)/2
        iHerR(iHer) = iHerR(1) + i_1111
        iHerW(iHer) = iHerW(1) + i_1111
 11   Continue
*
      Alpha = BInv(1)
      DO 2000 IDEG=3,MaxHer
        i_0000 = (IDEG*IDEG-IDEG)/2
        IR=iHerR(1)-1+i_0000
        IW=iHerW(1)-1+i_0000
        IDH=IDEG/2
        i_1111 = IR+IDH+1
        i_3333 = i_1111-IDEG
        w_3333 = HerR(i_3333)
        i_2222 = i_3333+1
        w_2222 = HerR(i_2222)
        X=(w_2222-w_3333)/2.0d0
        HerR(i_1111)=0.0d0
        DO 20 IROOT=2,IDEG,2
          j_0000 = IROOT/2
          j_1111 = IR+j_0000
          j_2222 = IR-IDEG+1+j_0000
          j_3333 = IR+IDEG+1-j_0000
          R =  HerR(j_2222)-X
          HerR(j_1111) = R
          HerR(j_3333) = -R
  20    CONTINUE
        DO 1000 IROOT=1,IDH
          j_0000 = IR+IROOT
          Z = HerR(j_0000)
          CORR = 0.0d0
          Do j = 1,ideg
            If ( j.ne.iroot ) then
              c_0000 = Z-HerR(IR+J)
              CORR=CORR+(1.0d0/c_0000)
            End If
          End Do
  99      CONTINUE
          Herm(2)=Z*Herm(1)*Alpha
          DO 110 K=1,IDEG-1
            w_1111 = Herm(K+1)
            w_3333 = Herm(K)
            w_4444 = Beta(K)
            w_5555 = BInv(K+1)
            w_2222 = (Z*w_1111-w_4444*w_3333)*w_5555
            Herm(K+2) = w_2222
 110      CONTINUE
          HDER=2.0d0*Beta(IDEG)*Herm(IDEG)
          DELTA=-Herm(IDEG+1)/(HDER-CORR*Herm(IDEG+1))
          Z=Z+DELTA
          IF(ABS(DELTA).GT.1.0d-8) then
             if(abs(DELTA).gt.1.0d8) then
                Call WarningMessage(1,'Warning: large value in sether')
c               write(6,*) delta
              endif
              goto 99
             endif
          HerR(IR+IROOT)=Z
          HerR(IR+IDEG+1-IROOT)=-Z
 1000   CONTINUE
        DO 3010 IROOT=1,IDH+1
          j_0000 = IR+IROOT
          Z = HerR(j_0000)
          Herm(2)=Z*Herm(1)*Alpha
          SUM=Herm(1)**2
          SUM=SUM+Herm(2)**2
          DO 3020 K=1,IDEG-2
            w_1111 = Herm(K+1)
            w_3333 = Herm(K)
            w_4444 = Beta(K)
            w_5555 = BInv(K+1)
            w_2222 = (Z*w_1111-w_4444*w_3333)*w_5555
            Herm(K+2) = w_2222
            SUM=SUM+w_2222*w_2222
 3020     CONTINUE
          W = 1.0d0/SUM
          HerW(IW+IROOT) = W
          HerW(IW+IDEG+1-IROOT) = W
 3010   CONTINUE
 2000 CONTINUE
      Call mma_deallocate(Beta)
      Call mma_deallocate(BInv)
      Call mma_deallocate(Herm)
*
*define _DEBUG_
#ifdef _DEBUG_
      Call TriPrt(' Hermite roots',' ',HerR(iHerR(1)),MaxHer)
      Call TriPrt(' Hermite weights',' ',HerW(iHerW(1)),MaxHer)
      Write (6,*) ' MaxHer=',MaxHer,nPrp,iAngMx
#endif
      Return
      End
