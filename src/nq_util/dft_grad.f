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
* Copyright (C) 2002, Roland Lindh                                     *
************************************************************************
      Subroutine DFT_Grad(Grad,nGrad,nD,Grid,mGrid,dRho_dR,ndRho_dR,
     &                    nGrad_Eff,Weights,iNQ)
************************************************************************
*                                                                      *
*     Object: to trace the correct parts to get the contributions to   *
*             the gradient due to the DFT energy.                      *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics, University of   *
*             Lund, Sweden.  May 2002 in Bologna, Italy.               *
************************************************************************
      use nq_Grid, only: F_xc, GradRho, vRho, vSigma, vTau, vLapl
      use nq_Grid, only: Pax
      use nq_Grid, only: IndGrd, iTab, Temp, dW_dR
      use nq_Structure, only: NQ_data
      use nq_Info
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "debug.fh"
#include "Molcas.fh"
#include "itmax.fh"
#include "ksdft.fh"
#include "stdalloc.fh"
      Parameter (Mxdc=MxAtom)
#include "disp.fh"
      Real*8 Grad(nGrad), Grid(3,mGrid),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff), OV(3,3), V(3,3),
     &       R_Grid(3), Weights(mGrid), OVT(3)
      Real*8, Allocatable:: Aux(:,:)
*                                                                      *
************************************************************************
*                                                                      *
      R_Grid(:)=NQ_Data(iNQ)%Coor(:)
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Debug=.True.
      If (Debug) Then
         Call RecPrt('R_Grid',' ',R_Grid,1,3)
         Call RecPrt('Grid',' ',Grid,3,mGrid)
         Call RecPrt('Weights',' ',Weights,1,mGrid)
         Call RecPrt('dW_dR',' ',dW_dR,nGrad_Eff,mGrid)
         Call RecPrt('dRho_dR(1)',' ',dRho_dR,ndRho_dR,mGrid)
         Call RecPrt('dF_dRho',' ',dF_dRho,ndF_dRho,mGrid)
         Do iEff = 1, nGrad_Eff
            Write (6,*) 'iTab=',iTab(1,iEff),iTab(2,iEff),iTab(3,iEff),
     &                          iTab(3,iEff)
            Write (6,*) 'IndGrd=',IndGrd(iEff)
         End Do
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     We have that the DFT energy is expressed as
*
*     E_DFT = Sum_Gg  w(r_g(G))  f(G,r_g(G))
*
*     r_g = R_G + O s_g
*
*     The first derivative is computed as
*
*     E_DFT^x = Sum w^x f  + w f^x
*
*     where
*
*     f^x = f^(x) + <nabla_r  f * r^x >
*
*     where
*
*     nabla_r f is the functional differentiated with respect to a
*     displacement of a grid point.
*
*     and
*
*     r^x = delta_AG e_i + O^x s
*                                                                      *
************************************************************************
*                                                                      *
*     Add the contributions
*
*     w * f^x to centers other than the origin of the current
*             set of grid points.
*
*     Note that for x being one of the cartesian components of the
*     center of the atomic grid we do not have f^x but rather f^(x).
*     The correct contribution will be added below.
*
*     Here we also accumulate contributions for the rotational
*     invariance.
*
*
*                                                                      *
************************************************************************
*                                                                      *
      Select Case (Functional_type)
*                                                                      *
************************************************************************
*                                                                      *
      Case (LDA_type)
*                                                                      *
************************************************************************
*                                                                      *

         Call mma_Allocate(Aux,1*nD,mGrid,Label='Aux')
         If (nD.eq.1) Then
            Do j = 1, mGrid
               Aux(1,j)=vRho(1,j)
            End Do
         Else
            Do j = 1, mGrid
               Aux(1,j)=vRho(1,j)
               Aux(2,j)=vRho(2,j)
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Case (GGA_type)
*                                                                      *
************************************************************************
*
         Call mma_Allocate(Aux,4*nD,mGrid,Label='Aux')
         If (nD.eq.1) Then
            Do j = 1, mGrid
               Aux(1,j)=vRho(1,j)
               Aux(2,j)=2.0d0*vSigma(1,j)*Gradrho(1,j)
               Aux(3,j)=2.0d0*vSigma(1,j)*Gradrho(2,j)
               Aux(4,j)=2.0d0*vSigma(1,j)*Gradrho(3,j)
            End Do
         Else
            Do j = 1, mGrid
               gxa=Gradrho(1,j)
               gya=Gradrho(2,j)
               gza=Gradrho(3,j)
               gxb=Gradrho(4,j)
               gyb=Gradrho(5,j)
               gzb=Gradrho(6,j)

               Aux(1,j)=vRho(1,j)
               Aux(2,j)=vRho(2,j)
               Aux(3,j)=( 2.0d0*vSigma(1,j)*gxa
     &                         +vSigma(2,j)*gxb )
               Aux(4,j)=( 2.0d0*vSigma(1,j)*gya
     &                         +vSigma(2,j)*gyb )
               Aux(5,j)=( 2.0d0*vSigma(1,j)*gza
     &                         +vSigma(2,j)*gzb )
               Aux(6,j)=( 2.0d0*vSigma(3,j)*gxb
     &                         +vSigma(2,j)*gxa )
               Aux(7,j)=( 2.0d0*vSigma(3,j)*gyb
     &                         +vSigma(2,j)*gya )
               Aux(8,j)=( 2.0d0*vSigma(3,j)*gzb
     &                         +vSigma(2,j)*gza )
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Case (meta_GGA_type1)
*                                                                      *
************************************************************************
*                                                                      *
         Call mma_Allocate(Aux,5*nD,mGrid,Label='Aux')
         If (nD.eq.1) Then
            Do j = 1, mGrid
               Aux(1,j)=vRho(1,j)
               Aux(2,j)=2.0d0*vSigma(1,j)*Gradrho(1,j)
               Aux(3,j)=2.0d0*vSigma(1,j)*Gradrho(2,j)
               Aux(4,j)=2.0d0*vSigma(1,j)*Gradrho(3,j)
               Aux(5,j)=0.25D0*vTau(1,j)
            End Do
         Else
            Do j = 1, mGrid
               gxa=Gradrho(1,j)
               gya=Gradrho(2,j)
               gza=Gradrho(3,j)
               gxb=Gradrho(4,j)
               gyb=Gradrho(5,j)
               gzb=Gradrho(6,j)

               Aux(1,j)=vRho(1,j)
               Aux(2,j)=vRho(2,j)
               Aux(3,j)=( 2.0d0*vSigma(1,j)*gxa
     &                         +vSigma(2,j)*gxb )
               Aux(4,j)=( 2.0d0*vSigma(1,j)*gya
     &                         +vSigma(2,j)*gyb )
               Aux(5,j)=( 2.0d0*vSigma(1,j)*gza
     &                         +vSigma(2,j)*gzb )
               Aux(6,j)=( 2.0d0*vSigma(3,j)*gxb
     &                         +vSigma(2,j)*gxa )
               Aux(7,j)=( 2.0d0*vSigma(3,j)*gyb
     &                         +vSigma(2,j)*gya )
               Aux(8,j)=( 2.0d0*vSigma(3,j)*gzb
     &                         +vSigma(2,j)*gza )
               Aux(9,j)=0.5D0*vTau(1,j)
               Aux(10,j)=0.5D0*vTau(2,j)
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Case (meta_GGA_type2)
*                                                                      *
************************************************************************
*                                                                      *
         Call mma_Allocate(Aux,6*nD,mGrid,Label='Aux')
         If (nD.eq.1) Then
            Do j = 1, mGrid
               Aux(1,j)=vRho(1,j)
               Aux(2,j)=2.0d0*vSigma(1,j)*Gradrho(1,j)
               Aux(3,j)=2.0d0*vSigma(1,j)*Gradrho(2,j)
               Aux(4,j)=2.0d0*vSigma(1,j)*Gradrho(3,j)
               Aux(5,j)=0.25D0*vTau(1,j)
               Aux(6,j)=vLapl(1,j)
            End Do
         Else
            Do j = 1, mGrid
               gxa=Gradrho(1,j)
               gya=Gradrho(2,j)
               gza=Gradrho(3,j)
               gxb=Gradrho(4,j)
               gyb=Gradrho(5,j)
               gzb=Gradrho(6,j)

               Aux(1,j)=vRho(1,j)
               Aux(2,j)=vRho(2,j)
               Aux(3,j)=( 2.0d0*vSigma(1,j)*gxa
     &                         +vSigma(2,j)*gxb )
               Aux(4,j)=( 2.0d0*vSigma(1,j)*gya
     &                         +vSigma(2,j)*gyb )
               Aux(5,j)=( 2.0d0*vSigma(1,j)*gza
     &                         +vSigma(2,j)*gzb )
               Aux(6,j)=( 2.0d0*vSigma(3,j)*gxb
     &                         +vSigma(2,j)*gxa )
               Aux(7,j)=( 2.0d0*vSigma(3,j)*gyb
     &                         +vSigma(2,j)*gya )
               Aux(8,j)=( 2.0d0*vSigma(3,j)*gzb
     &                         +vSigma(2,j)*gza )
               Aux(9,j)=0.5D0*vTau(1,j)
               Aux(10,j)=0.5D0*vTau(2,j)
               Aux(11,j)=vLapl(1,j)
               Aux(12,j)=vLapl(2,j)
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Case Default
*                                                                      *
************************************************************************
*                                                                      *
         Call WarningMessage(2,'Do_Grad: wrong functional type!')
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
      End Select
*                                                                      *
************************************************************************
*                                                                      *
      OV(:,:)=Zero
      Do i_Eff=1, nGrad_Eff
         tmp=Zero
         OVT(:)=Zero
         Do j = 1, mGrid
            dF_dr = Weights(j)*DOT_Product(Aux(:,j),dRho_dR(:,j,i_Eff))
            tmp = tmp + dF_dr
*
*           Accumulate stuff for rotational invariance
*
            OVT(:) = OVT(:) + dF_dr * Grid(:,j)
         End Do
         ixyz=iTab(1,i_Eff)
         OV(ixyz,:) = OV(ixyz,:) + OVT(:) - tmp * R_Grid(:)
         Temp(i_Eff)=-tmp
      End Do

      Call mma_deAllocate(Aux)

      Do i_Eff=1, nGrad_Eff
         If (iTab(2,i_Eff)==Off) Temp(i_Eff)=Zero
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (Debug) Then
         Call RecPrt('w * f^x before translational contributions',
     &               ' ',Temp,1,nGrad_Eff)
         Call RecPrt('OV',' ',OV,3,3)
      End If
#endif
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Here we compute the term
*
*     w * f^x
*
*     for x being a cartesian component of the center of the atomic
*     grid. This is done using the translational invariance condition.
*
      If (Grid_Type.eq.Moving_Grid) Then
         Do ixyz = 1, 3
            iGrad=0
            Do jGrad = 1, nGrad_Eff
               If (iTab(1,jGrad).eq.ixyz .and.
     &             iTab(2,jGrad).eq.Off  .and.
     &             IndGrd(jGrad).gt.0           ) iGrad=jGrad
            End Do
#ifdef _DEBUGPRINT_
            If (Debug) Write (6,*) 'iGrad=',iGrad
#endif
            If (iGrad.ne.0) Then
*
*              Evaluate indirectly via the translational invariance
*              the sum of the direct and indirect term
*
               Do jGrad = 1, nGrad_Eff
                  If (jGrad.ne.iGrad .and. iTab(1,jGrad).eq.ixyz ) Then
*
                      Temp(iGrad)=Temp(iGrad)-Temp(jGrad)
*
#ifdef _DEBUGPRINT_
            If (Debug) Write (6,*) 'jGrad,Temp(jGrad)=',
     &                              jGrad,Temp(jGrad)
#endif
                  End If
               End Do
*
            End If
         End Do
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('w * f^x',' ',Temp,1,nGrad_Eff)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     For a "moving" grid add contributions due to the derivative with
*     respect to the partitioning.
*
*     w^x * f
*
         Call DGEMM_('N','N',nGrad_Eff,1,mGrid,
     &              One,dW_dR,nGrad_Eff,
     &                  F_xc,mGrid,
     &              One,Temp,nGrad_Eff)
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('w * f^x + w^x * f',' ',Temp,1,nGrad_Eff)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*        Add the rotational invariance term
*
*        First transform back to the cartesian coordinates system.
*
         Fact=DBLE(2-(nD/2))
         Call DGEMM_('N','N',
     &               3,3,3,
     &               Fact,OV,3,
     &               Pax,3,
     &               0.0d0,V,3)
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('V',' ',V,3,3)
#endif
*
         Do i_Eff = 1, nGrad_Eff
            iCar = iTab(1,i_Eff)
            jNQ  = iTab(3,i_Eff)
*
*           Compute < nabla_r f * r^x > as Tr (O^x V)
*
            Tmp = DDot_(9,NQ_Data(jNQ)%dOdx(:,:,iCar),1,V,1) * Half
#ifdef _DEBUGPRINT_
            If (Debug) Then
               Write (6,*)
               Write (6,*) 'iCar,jNQ=',iCar,jNQ
               Call RecPrt('dOdx',' ',NQ_Data(jNQ)%dOdx(:,:,iCar),3,3)
               Write (6,*) 'Tmp=',Tmp
            End If
#endif
            Temp(i_Eff) = Temp(i_Eff) - Tmp
         End Do
*
      End If !moving grid
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('Gradient contribution from this block',
     &                       ' ',Temp,1,nGrad_Eff)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Accumulate and symmetry adapt.
*
      Do i_Eff = 1, nGrad_Eff
         i = IndGrd(i_Eff)
         If (i.ge.1) Then
            Fact=DBLE(iTab(4,i_Eff))
            Grad(i) = Grad(i) + Fact*Temp(i_Eff)
         End If
      End Do
#ifdef _DEBUGPRINT_
      If (Debug) Call RecPrt('Gradient accumulated so far',
     &                     ' ',Grad,1,nGrad)
      Debug=.False.
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
