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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2003, Valera Veryazov                                  *
*               2016, Roland Lindh                                     *
************************************************************************
      SubRoutine DIIS_i(CInter,nCI,TrDh,TrDP,TrDD,nTr,nD,iOpt_DIIS,Ind)
************************************************************************
*                                                                      *
*     purpose: density matrix optimization                             *
*                                                                      *
*              EDIIS optimization:                                     *
*                   K. N. Kudin, G. E. Scuseria, and E. Cances         *
*                   JCP, 116, , 8255 (2002)                            *
*                   doi: 10.1063/1.1470195                             *
*                                                                      *
*     input:                                                           *
*       TrDh    : Traces of D(i)*h of size (nTr)                       *
*       TrDP    : Traces of D(i)*P(j) of size (nTr,nTr)                *
*       TrDD    : Traces of D(i)*D(j) of size (nTr,nTr)                *
*                                                                      *
*     output:                                                          *
*       CInter  : Interpolation coefficients of length nCI             *
*                                                                      *
*     called from: WfCtl                                               *
*                                                                      *
*     calls to: Optim                                                  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: UHF - V.Veryazov, 2003                                  *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 CInter(nCI,nD),TrDh(nTr,nTr,nD),TrDP(nTr,nTr,nD),
     &                      TrDD(nTr,nTr,nD)
      Logical Do_SpinAV
      COMMON  / SPAVE_L  / Do_SpinAV
*
*---- Define local variables
      Real*8 Eline(MxOptm,2),Equad(MxOptm**2,2),DD(MxOptm**2,2)
      Integer Ind(MxOptm)
*define _NEW_CODE_
#ifdef _NEW_CODE_
      Logical Ignore
#endif
*     Save Eline,Equad
      Real*8 EPred(MxIter+1), h, r_SO
      Save h, EPred
      Data h/0.35D0/
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*define _DEBUG_
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) 'E_pred=',(EPred(i),i=1,iter)
      If (nD.eq.1) Then
         Write (6,*) 'E_actu=',(Elst(i,1),i=1,iter)
      Else
         Write (6,*) 'E_actu=',(Elst(i,1)+Elst(i,2),i=1,iter)
      End If
#endif
      If (iOpt_DIIS.eq.1) Then
         AccCon = 'EDIIS    '
      Else
         AccCon = 'ADIIS    '
      End If
*
#ifdef _NEW_CODE_
      Do i = kOptim, 1, -1
         Ind(i)=0
*
         tmp0= Zero
         Do j = 1, iter
*
            Ignore=.False.
            Do k = kOptim, i+1, -1
               Ignore = Ignore .or. Ind(k).eq.j
            End Do
            If (Ignore) Cycle
*
            tmp = Zero
            Do iD = 1, nD
               tmp = tmp + Elst(j,iD)
            End Do
*
            If (tmp.lt.tmp0) Then
               tmp=tmp0
               Ind(i)=j
            End If
*
         End Do
      End Do
#else
      Do i = 1, kOptim
         Ind(i) = iter-kOptim+i
      End Do
#endif
*
      Do iD = 1, nD
*
*        EDIIS optimization, doi:10.1063/1.1470195, Eq. (8)
*
*        Noticed the change in sign - in optim the quadratic terms,
*        however, are added in the evaluation of Eq. (8). Additionally,
*        the paper is written in a spin-orbital notation. A factor of
*        1/2 for the two-electron term is assumed and does not need to
*        be included.
*
         If (iOpt_DIIS.eq.1) Then
*
            Do i=1,kOptim
               ii = Ind(i)
               Eline(i,iD)=Elst(ii,iD)
               Do j=1,kOptim
                  jj= Ind(j)
                  DiFi=TrDh(ii,ii,iD)+TrDP(ii,ii,iD)
                  DiFj=TrDh(ii,ii,iD)+TrDP(ii,jj,iD)
                  DjFi=TrDh(jj,jj,iD)+TrDP(jj,ii,iD)
                  DjFj=TrDh(jj,jj,iD)+TrDP(jj,jj,iD)
                  Equad(kOptim*(i-1)+j,iD)=
     &             -Half*(DiFi-DiFj-DjFi+DjFj)
                  DD(kOptim*(i-1)+j,iD)=TrDD(ii,jj,iD)
               End Do
            End Do
*
         Else
*
*        ADIIS optimization (Eq. 7). We arbitrarily set E(Dn) to zero.
*        The option is not tested yet.
*
            kk = iter ! This needs to be properly set!
            Do i=1,kOptim
               ii = Ind(i)
               DiFn=TrDh(ii,ii,iD)+TrDP(ii,kk,iD)
               DnFn=TrDh(kk,kk,iD)+TrDP(kk,kk,iD)
               Eline(i,iD)= DiFn-DnFn
               Do j=1,kOptim
                  jj = Ind(j)
                  DiFj=TrDh(ii,ii,iD)+TrDp(ii,jj,iD)
                  DnFj=TrDh(kk,kk,iD)+TrDp(kk,jj,iD)
                  Equad(kOptim*(i-1)+j,iD)=  DiFj-DiFn-DnFj+DnFn
                  DD(kOptim*(i-1)+j,iD)=TrDD(ii,jj,iD)
               End Do
            End Do
*
         End If
*
      End Do
*
*     Tweak for UHF
*
      If (nD.eq.2) Then
         Do i=1,kOptim
            tmp_a =Eline(i,1)
            tmp_b =Eline(i,2)
            Eline(i,1) = tmp_a + tmp_b
            Eline(i,2) = tmp_a
            Do j=1,kOptim
               tmp_a = Equad(kOptim*(i-1)+j,1)
               tmp_b = Equad(kOptim*(i-1)+j,2)
               Equad(kOptim*(i-1)+j,1) = tmp_a + tmp_b
               Equad(kOptim*(i-1)+j,2) = tmp_a
            End Do
         End Do
      End If
*
************************************************************************
*
      If(kOptim.ge.3) Then
         BigOne=Zero
         BigTwo=Zero
         Do iD = 1, nD
            Do i=2,kOptim
               Do j=1,i-1
                  BigOne=Max(BigOne,Abs(Eline(i,iD)-Eline(j,iD)))
               End Do
            End Do
            Do i=1,kOptim
               Do j=1,kOptim
                  BigTwo=Max(BigTwo,Abs(Equad(i+(j-1)*kOptim,iD)))
               End Do
            End Do
         End Do
         Big=Max(BigOne,BigTwo)
         If(Big.lt.1.0d-8) Then
            EmConv=.true.
            WarnPocc=.true.
         Else
            EmConv=.false.
            WarnPocc=.false.
         End If
         If (Do_SpinAV) Then ! it is not a diagnostic in this case
            WarnPocc=.false.
         EndIf
      End If
*                                                                      *
************************************************************************
*                                                                      *
*---- Find interpolation coefficients with the constaints
*
*     Sum_i c_i = 1
*
*     0 =< c_i =< 1
*
      Call Optim(E_Pred,Eline,Equad,CInter(1,1),kOptim,kOptim)
      EPred(iter+1)=E_Pred
*
#ifdef _DEBUG_
         Write(6,*)' Interpolation coefficients:'
         Write(6,'(5f16.8)')(CInter(i,1), i = 1, kOptim)
#endif
*
*     Temporary fix for UHF
*
      If (nD.eq.2) Call DCopy_(nCI,CInter(1,1),1,CInter(1,2),1)
*
#ifdef _DEBUG_
      Write(6,*)' Interpolation coefficients:'
      Write(6,'(5f16.8)')(CInter(i,1), i = 1, kOptim)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Update the trust radius, h, by checking to which degree previous
*     interpolations estimated the next energy correctly.
*
*     The first two energies have not been predicted!
*
      If (iter.gt.3) Then
         E_Pred=EPred(iter)
         E_n1 = Zero
         E_n  = Zero
         Do iD = 1, nD
            E_n1 = E_n1 + Elst(iter,iD)
            E_n  = E_n  + Elst(iter-1,iD)
         End Do
*
#ifdef _DEBUG_
         Write (6,*) 'iter=',iter
         Write (6,*) 'Energy of iter  =',E_n1
         Write (6,*) 'E_pred of iter  =',E_Pred
         Write (6,*) 'Energy of iter-1=',E_n
         Write (6,*)
#endif
*
*        In some cases the DIIS will come out with a set to
*        coefficient which corresponds to the last density. In this
*        case the E_Pred(i+1)=E(i). We add a small number to avoid
*        dividing with zero.
*
         r_SO = (E_n1 - E_n)/(E_Pred - E_n + 1.0D-12 )
      Else
         r_SO = One
      End If
*
*     Update the trust radius according to this ad hoc scheme.
*
      If (r_SO.ge.0.75D0) Then
         h = 1.2D0 * h
*     Else If (r_SO.ge.0.25D0) Then
*        h = 1.0D0 * h
      Else If (r_SO.lt.0.25D0) Then
         h = 0.7D0 * h
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Find the reference density for the density with the
*     lowest energy.
*
      E_Min=Zero
      n_min = 0
      n1    = 0
      Do i=1,kOptim
         E_tot=Zero
         Do iD = 1, nD
            E_tot = E_tot + Elst(Ind(i),iD)
         End Do
         If (E_tot.lt.E_Min) Then
            n1    = n_min
            n_min = i
            E_Min = E_tot
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the square of the trace of the difference between the
*     reference density and the new interpolated density.
*
      If (.false.) Then
      nn = Ind(n_min)
      r2 = Zero
      Do i = 1, kOptim
         ii = Ind(i)
         Do j = 1, kOptim
            jj = Ind(j)
            Do iD = 1, nD
               r2 = r2 + CInter(i,iD)*CInter(j,iD) *
     &              (TrDD(ii,jj,iD)-TrDD(nn,jj,iD)
     &              -TrDD(ii,nn,iD)+TrDD(nn,nn,iD))
            End Do
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Perform a RS-DIIS interpolation if required
*
      If (Sqrt(r2).gt.h) Then
         Write (6,*) 'Apply optimization with step restriction'
         Write (6,*) 'r,h =', Sqrt(r2),h
         Call Abend()
         Call Optim2(E_Pred,Eline,Equad,DD,CInter(1,1),kOptim,kOptim,
     &               n_min,n1,r2)
         EPred(iter+1)=E_Pred
*
*        Temporary fix for UHF
*
         If (nD.eq.2) Call DCopy_(nCI,CInter(1,1),1,CInter(1,2),1)
*
      End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Check that the coefficients sum up correctly.
*
      Do iD = 1, nD
         CSum = Zero
         Do i = 1, kOptim
            CSum = CSum + CInter(i,iD)
         End Do
         If (Abs(CSum - One).gt.1.0D-5) Then
            Write (6,*) 'diis_i: Abs(CSum - One).gt.1.0D-5'
            Write (6,*) 'CSum=',CSum
            Call QTrace
            Call Abend()
         End If
*
*        If the coefficient for the last density is zero we are in
*        problem.
*
         If (CInter(kOptim,iD).eq.Zero.and.kOptim.gt.2.and.
     &       CInter(kOptim-1,iD).eq.One) Then
            Write (6,*) 'DIIS_I optimization failed!'
*           Write (6,*) 'iD=',iD
*           Write (6,*) (CInter(i,iD),i=1,kOptim)
            CInter(kOptim,iD)=1.0D-6
            CInter(kOptim,iD-1)=One-1.0D-6
*           EmConv=.True.
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld( 6) = TimFld( 6) + (Cpu2 - Cpu1)
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
