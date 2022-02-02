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
      Subroutine DFT_Grad(Grad,nGrad,iSpin,Grid,mGrid,dRho_dR,ndRho_dR,
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
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "nq_info.fh"
#include "debug.fh"
#include "Molcas.fh"
#include "itmax.fh"
#include "ksdft.fh"
      Parameter (Mxdc=MxAtom)
#include "disp.fh"
      Real*8 Grad(nGrad), Grid(3,mGrid),
     &       dRho_dR(ndRho_dR,mGrid,nGrad_Eff), OV(3,3), V(3,3),
     &       R_Grid(3), Weights(mGrid)
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
      Call FZero(Temp,nGrad_Eff)
      Call FZero(OV,9)
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
      If (Functional_type.eq.LDA_type) Then
*
         If (iSpin.eq.1) Then
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  dF_dr = vRho(1,j)    *dRho_dR(1,j,i_Eff)
*
                  tmp = tmp + Weights(j) * dF_dr
*
*                 For rotational invariance accumulate
*
*                 (nabla_r f_g)^T O s_g
*
                  OV(ixyz,1) = OV(ixyz,1) + Two* Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) + Two* Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) + Two* Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         Else
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  dF_dr = vRho(1,j)    *dRho_dR(1,j,i_Eff)
     &                   +vRho(2,j)    *dRho_dR(2,j,i_Eff)
                  tmp = tmp + Weights(j) * dF_dr
*
*                 Accumulate stuff for rotational invariance
*
                  OV(ixyz,1) = OV(ixyz,1) +      Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) +      Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) +      Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do

*****************************************************************************************************************
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.GGA_type) Then
*
         If (iSpin.eq.1) Then
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  gx=Gradrho(1,j)
                  gy=Gradrho(2,j)
                  gz=Gradrho(3,j)
                  Temp0=vRho(1,j)
                  Temp1=2.0d0*vSigma(1,j)*gx
                  Temp2=2.0d0*vSigma(1,j)*gy
                  Temp3=2.0d0*vSigma(1,j)*gz
*
                  dF_dr = Temp0*dRho_dR(1,j,i_Eff)
     &                  + Temp1*dRho_dR(2,j,i_Eff)
     &                  + Temp2*dRho_dR(3,j,i_Eff)
     &                  + Temp3*dRho_dR(4,j,i_Eff)
                  tmp = tmp  + Weights(j) * dF_dr
*
*                 Accumulate stuff for rotational invariance

*****************************************************************************************************************
*
                  OV(ixyz,1) = OV(ixyz,1) + Two* Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) + Two* Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) + Two* Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         Else
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  gxa=Gradrho(1,j)
                  gya=Gradrho(2,j)
                  gza=Gradrho(3,j)
                  gxb=Gradrho(4,j)
                  gyb=Gradrho(5,j)
                  gzb=Gradrho(6,j)

                  Temp0a=vRho(1,j)
                  Temp0b=vRho(2,j)
                  Temp1a=( 2.0d0*vSigma(1,j)*gxa
     &                          +vSigma(2,j)*gxb )
                  Temp1b=( 2.0d0*vSigma(3,j)*gxb
     &                          +vSigma(2,j)*gxa )
                  Temp2a=( 2.0d0*vSigma(1,j)*gya
     &                          +vSigma(2,j)*gyb )
                  Temp2b=( 2.0d0*vSigma(3,j)*gyb
     &                          +vSigma(2,j)*gya )
                  Temp3a=( 2.0d0*vSigma(1,j)*gza
     &                          +vSigma(2,j)*gzb )
                  Temp3b=( 2.0d0*vSigma(3,j)*gzb
     &                          +vSigma(2,j)*gza )
*
                  dF_dr = Temp0a*dRho_dR(1,j,i_Eff)
     &                  + Temp0b*dRho_dR(2,j,i_Eff)
     &                  + Temp1a*dRho_dR(3,j,i_Eff)
     &                  + Temp2a*dRho_dR(4,j,i_Eff)
     &                  + Temp3a*dRho_dR(5,j,i_Eff)
     &                  + Temp1b*dRho_dR(6,j,i_Eff)
     &                  + Temp2b*dRho_dR(7,j,i_Eff)
     &                  + Temp3b*dRho_dR(8,j,i_Eff)
                  tmp = tmp + Weights(j) * dF_dR
*
*                 Accumulate stuff for rotational invariance
*
                  OV(ixyz,1) = OV(ixyz,1) +      Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) +      Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) +      Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type1) Then
         If (iSpin.eq.1) Then
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  gx=Gradrho(1,j)
                  gy=Gradrho(2,j)
                  gz=Gradrho(3,j)
                  Temp0=vRho(1,j)
                  Temp1=2.0d0*vSigma(1,j)*gx
                  Temp2=2.0d0*vSigma(1,j)*gy
                  Temp3=2.0d0*vSigma(1,j)*gz
                  Temp4=0.25D0*vTau(1,j)
*
                  dF_dr = Temp0*dRho_dR(1,j,i_Eff)
     &                  + Temp1*dRho_dR(2,j,i_Eff)
     &                  + Temp2*dRho_dR(3,j,i_Eff)
     &                  + Temp3*dRho_dR(4,j,i_Eff)
     &                  + Temp4*dRho_dR(5,j,i_Eff)
                  tmp = tmp  + Weights(j) * dF_dr
*
*                 Accumulate stuff for rotational invariance
*
                  OV(ixyz,1) = OV(ixyz,1) + Two* Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) + Two* Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) + Two* Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         Else
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  gxa=Gradrho(1,j)
                  gya=Gradrho(2,j)
                  gza=Gradrho(3,j)
                  gxb=Gradrho(4,j)
                  gyb=Gradrho(5,j)
                  gzb=Gradrho(6,j)

                  Temp0a=vRho(1,j)
                  Temp0b=vRho(2,j)
                  Temp1a=( 2.0d0*vSigma(1,j)*gxa
     &                          +vSigma(2,j)*gxb )
                  Temp1b=( 2.0d0*vSigma(3,j)*gxb
     &                          +vSigma(2,j)*gxa )
                  Temp2a=( 2.0d0*vSigma(1,j)*gya
     &                          +vSigma(2,j)*gyb )
                  Temp2b=( 2.0d0*vSigma(3,j)*gyb
     &                          +vSigma(2,j)*gya )
                  Temp3a=( 2.0d0*vSigma(1,j)*gza
     &                          +vSigma(2,j)*gzb )
                  Temp3b=( 2.0d0*vSigma(3,j)*gzb
     &                          +vSigma(2,j)*gza )
                  Temp4a=0.5D0*vTau(1,j)
                  Temp4b=0.5D0*vTau(2,j)
*
                  dF_dr = Temp0a*dRho_dR(1,j,i_Eff)
     &                  + Temp0b*dRho_dR(2,j,i_Eff)
     &                  + Temp1a*dRho_dR(3,j,i_Eff)
     &                  + Temp2a*dRho_dR(4,j,i_Eff)
     &                  + Temp3a*dRho_dR(5,j,i_Eff)
     &                  + Temp1b*dRho_dR(6,j,i_Eff)
     &                  + Temp2b*dRho_dR(7,j,i_Eff)
     &                  + Temp3b*dRho_dR(8,j,i_Eff)
     &                  + Temp4a*dRho_dR(9,j,i_Eff)
     &                  + Temp4b*dRho_dR(10,j,i_Eff)
                  tmp = tmp + Weights(j) * dF_dR
*
*                 Accumulate stuff for rotational invariance
*
                  OV(ixyz,1) = OV(ixyz,1) +      Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) +      Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) +      Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Functional_type.eq.meta_GGA_type2) Then
         If (iSpin.eq.1) Then
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  gx=Gradrho(1,j)
                  gy=Gradrho(2,j)
                  gz=Gradrho(3,j)
                  Temp0=vRho(1,j)
                  Temp1=2.0d0*vSigma(1,j)*gx
                  Temp2=2.0d0*vSigma(1,j)*gy
                  Temp3=2.0d0*vSigma(1,j)*gz
                  Temp4=0.25D0*vTau(1,j)
                  Temp5=vLapl(1,j)
*
                  dF_dr = Temp0*dRho_dR(1,j,i_Eff)
     &                  + Temp1*dRho_dR(2,j,i_Eff)
     &                  + Temp2*dRho_dR(3,j,i_Eff)
     &                  + Temp3*dRho_dR(4,j,i_Eff)
     &                  + Temp4*dRho_dR(5,j,i_Eff)
     &                  + Temp5*dRho_dR(6,j,i_Eff)
                  tmp = tmp  + Weights(j) * dF_dr
*
*                 Accumulate stuff for rotational invariance
*
                  OV(ixyz,1) = OV(ixyz,1) + Two* Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) + Two* Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) + Two* Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         Else
            Do i_Eff=1, nGrad_Eff
               tmp=Zero
               ixyz=iTab(1,i_Eff)
               Do j = 1, mGrid
                  gxa=Gradrho(1,j)
                  gya=Gradrho(2,j)
                  gza=Gradrho(3,j)
                  gxb=Gradrho(4,j)
                  gyb=Gradrho(5,j)
                  gzb=Gradrho(6,j)

                  Temp0a=vRho(1,j)
                  Temp0b=vRho(2,j)
                  Temp1a=( 2.0d0*vSigma(1,j)*gxa
     &                          +vSigma(2,j)*gxb )
                  Temp1b=( 2.0d0*vSigma(3,j)*gxb
     &                          +vSigma(2,j)*gxa )
                  Temp2a=( 2.0d0*vSigma(1,j)*gya
     &                          +vSigma(2,j)*gyb )
                  Temp2b=( 2.0d0*vSigma(3,j)*gyb
     &                          +vSigma(2,j)*gya )
                  Temp3a=( 2.0d0*vSigma(1,j)*gza
     &                          +vSigma(2,j)*gzb )
                  Temp3b=( 2.0d0*vSigma(3,j)*gzb
     &                          +vSigma(2,j)*gza )
                  Temp4a=0.5D0*vTau(1,j)
                  Temp4b=0.5D0*vTau(2,j)
                  Temp5a=vLapl(1,j)
                  Temp5b=vLapl(2,j)
*
                  dF_dr = Temp0a*dRho_dR(1,j,i_Eff)
     &                  + Temp0b*dRho_dR(2,j,i_Eff)
     &                  + Temp1a*dRho_dR(3,j,i_Eff)
     &                  + Temp2a*dRho_dR(4,j,i_Eff)
     &                  + Temp3a*dRho_dR(5,j,i_Eff)
     &                  + Temp1b*dRho_dR(6,j,i_Eff)
     &                  + Temp2b*dRho_dR(7,j,i_Eff)
     &                  + Temp3b*dRho_dR(8,j,i_Eff)
     &                  + Temp4a*dRho_dR(9,j,i_Eff)
     &                  + Temp4b*dRho_dR(10,j,i_Eff)
     &                  + Temp5a*dRho_dR(11,j,i_Eff)
     &                  + Temp5b*dRho_dR(12,j,i_Eff)
                  tmp = tmp + Weights(j) * dF_dR
*
*                 Accumulate stuff for rotational invariance
*
                  OV(ixyz,1) = OV(ixyz,1) +      Weights(j) *
     &                        dF_dr * (Grid(1,j)-R_Grid(1))
                  OV(ixyz,2) = OV(ixyz,2) +      Weights(j) *
     &                        dF_dr * (Grid(2,j)-R_Grid(2))
                  OV(ixyz,3) = OV(ixyz,3) +      Weights(j) *
     &                        dF_dr * (Grid(3,j)-R_Grid(3))
               End Do
               If (iTab(2,i_Eff).ne.Off)
     &            Temp(i_Eff)=Temp(i_Eff)-tmp
            End Do
         End If
*                                                                      *
************************************************************************
*                                                                      *
      Else
*
         Call WarningMessage(2,'Do_Grad: wrong functional type!')
         Call Abend()
      End If
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
         Call DGEMM_('N','N',
     &               3,3,3,
     &               1.0d0,OV,3,
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
c Avoid unused argument warnings
      End
