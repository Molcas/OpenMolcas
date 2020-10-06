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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
Subroutine Start_Kriging(nPoints_In,nD,nInter,x_,dy_,y_)
  use kriging_mod
  Implicit None
#include "stdalloc.fh"
!
!    nPoints: the number of sample points, n
!    nInter: the dimensionality of the function, d
!    y_: the values of the function at the sample points
!    dy_: the gradient of the function at the sample points
!    x_: the coordinates of the sample points
!
  Integer nInter,nPoints_In,nD
  Real*8 x_(nInter,nPoints_In)
  Real*8 y_(nPoints_In)
  Real*8 dy_(nInter,nPoints_In-nd)
!
!
!#define _DEBUG_
#ifdef _DEBUG_
  Call RecPrt('Start_Kriging: x',' ',x_,nInter,nPoints_In)
  Call RecPrt('Start_Kriging: y',' ',y_,     1,nPoints_In)
  Call RecPrt('Start_Kriging: dy',' ',dy_,nInter,nPoints_In-nD)
#endif
!
! Call Setup_Kriging to store the data in some internally protected arrays and scalars.
!
  Call Setup_Kriging(nPoints_In,nD,nInter,x_,dy_,y_)
!
! Allocate x0, which is a n-dimensional vector of the coordinats of the last iteration computed
!
  Call mma_Allocate(x0,nInter,Label="nx")
!
! m_t is the dimentionality of the square correlation matrix Gradient-Psi
! (equation (2) on:
!-------- ref. = doi:10.1007/s00366-015-0397-y)-------
!
!
! In the case the nPoint_v last energies and nPoint_g last gradients were use
! m_t would be computed as
!
  m_t=nPoints_v + nInter*nPoints_g
!
! npx is the number of new points (Energy and Gradient) to be predict
! according to the iteration that was computed in update_sl subroutine
!
  npx = 1
!
! full_R correspond to the gradient of Psi (eq. (2) ref.)
!
  Call mma_Allocate(full_R,m_t,m_t,Label="full_R")
  Call mma_Allocate(full_RInv,m_t,m_t,Label="full_RInv")
!
  If (mblAI) sbmev = y(maxloc(y,dim=1))
!
! rl and dl are temporary matrices for the contruction of Psi which is inside of
! Grad-Psi (eq.(2) ref.) dl=rl^2=Sum[i] [(x_i-x0_i)/l)^2]
! more inoformation is given in subsequen files.
! Mat is the final matrix after the distance (between source data rl and dl) has
! passed through the Matern correlation function (ISBN 0-486-61272-4 & eq. (11-12)
! ref.).
!Iden is just an identity matrix necesary to avoid that the Grad-Psi becomes
! Singular after been multiplied by EPS factor
!
  Call mma_Allocate(rl,nPoints_v,npx,nInter,Label="rl")
  Call mma_Allocate(dl,nPoints_v,npx,Label="dl")
  Call mma_Allocate(Rones,m_t,Label="Rones")
!
!kv is the vector that contains the dot product of the inverse of Grad-Psi and
!Grad-y minus the dot product of the inverse of Grad-Psi and f-ones multiplied
! by the constant Grad-Trend function (eq. (3), (5), (6) and (7)
!ref.).
!Pred, gpred and hpred are the predicted energy, gradient and hessian according
! to the coordinates given by update_sl
!var is part of the GEK variance and the concentrated Likelihood function (eq.
! (8) and (9) ref.).
!Sigma is a dispersion function +-sigma=1.96*sqrt(abs(var*variance)) (need reference).
!cv is the covariant vector that contains the correlation function values and
! gradients between the sample data and the prediction
! (eq. 4 ref.).
!l is a n-dimensional vector of the width of the Matern function.
!ll is the likelihood function.
!
! Allocate additional variables needed for the kriging.
!
  Call mma_allocate(kv,m_t,Label="kv")
  Call mma_allocate(gpred,nInter,Label="gpred")
  Call mma_allocate(hpred,nInter,nInter,Label="hpred")
  Call mma_allocate(var,npx,Label="var")
  Call mma_allocate(sigma,npx,Label="sigma")
  Call mma_allocate(l,nInter,Label="l")
  Call mma_allocate(ll,int(lb(3)),Label="ll")
!
  Call mma_allocate(cv,m_t,npx,nInter,nInter,Label="cv")
  Call mma_allocate(cvMatFder,nPoints_v,npx,Label="cvMatFder")
  Call mma_allocate(cvMatSder,nPoints_v,npx,Label="cvMatSder")
  Call mma_allocate(cvMatTder,nPoints_v,npx,Label="cvMatTder")
!
  return
End Subroutine Start_Kriging
