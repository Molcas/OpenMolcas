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

      Subroutine Start_Kriging(nPoints,nInter,x_,dy_,y_)
      use globvar
#include "stdalloc.fh"
!
        Integer nInter,nPoints
        Real*8 x_(nInter,nPoints),dy_(nInter,nPoints),y_(nPoints)
        nInter_save=nInter
        nPoints_save=nPoints
!
        Call mma_Allocate(x,nInter,nPoints,Label="x")
        Call mma_Allocate(dy,nInter*nPoints,Label="dy")
        Call mma_Allocate(y,nInter,Label="y")
        Call mma_Allocate(nx,nInter,1,Label="nx")
        Call mma_Allocate(nx_SaveG,nInter,1,Label="nx_SaveG")
        Call mma_Allocate(nx_SaveH,nInter,1,Label="nx_SaveH")
!m_t is the dimentionality of the square correlation matrix Gradient-Psi
! (equation (2) on:
!-------- ref. = DOI 10.1007/s00366-015-0397-y)-------
          m_t=nPoints*(1+nInter)
!npx is the number of new points (Energy and Gradient) to be predict
! according to the iteration that was computed in update_sl subroutine
          npx = 1
!full_R correspond to the gradient of Psi (eq. (2) ref.)
        Call mma_Allocate(full_R,m_t,m_t,Label="full_R")
!nx is the n-dimensional vector of the last iteration cumputed in update_sl
! subroutine
!x is the n-dimensional internal coordinates
          x = x_
          ! write(6,*) 'x',x
!y is the energy
          y = y_
          ! write(6,*) 'y',y
!dy it's a vector of Grad-y (eq. (5)  ref. gradients of
! the energy with respect to the internal coordinates
          do i=1,nInter
            do j=1,nPoints
              dy((i-1)*nPoints+j) = dy_(i,j)
            enddo
          enddo
          ! write(6,*) 'dy',dy
!
          If (mblAI) sbmev = y(maxloc(y,dim=1))
!rl and dl are temporary matrices for the contruction of Psi which is inside of
! Grad-Psi (eq.(2) ref.) dl=rl^2=Sum[i] [(x_i-x0_i)/l)^2]
! more inoformation is given in subsequen files.
!Mat is the final matrix after the distance (between source data rl and dl) has
! passed through the Mat'ern correlation function (ISBN 0-486-61272-4 & eq. (11-12)
! ref.).
!Iden is just an identity matrix necesary to avoid that the Grad-Psi becomes
! Singular after been multiplied by EPS factor
          Call mma_Allocate(rl,nPoints,npx,nInter,Label="rl")
          Call mma_Allocate(dl,nPoints,npx,Label="dl")
          Call mma_Allocate(Rones,m_t,Label="Rones")
!kv is the vector that contains the dot product of the inverse of Grad-Psi and
!Grad-y minus the dot product of the inverse of Grad-Psi and f-ones multiplied
! by the constant Grad-Trend function (eq. (3), (5), (6) and (7)
!ref.).
!Pred, gpred and hpred are the predicted energy, gradient and hessian according
! to the coordinates given by update_sl
!var is part of the GEK variance and the concentrated Likelihood function (eq.
! (8) and (9) ref.).
!Sigma is a dispersion function +-sigma=1.96*sqrt(abs(var*variance)) (need reference).
!cv is the corvariant vector that contains the correlation function values and
! gradients between the sample data and the prediction
! (eq. 4 ref.).
!l is a n-dimensional vector of the width of the Mat'ern function.
!ll is the likelihood function.
!
          Call mma_allocate(kv,m_t,Label="kv")
          Call mma_allocate(pred,npx,Label="pred")
          Call mma_allocate(gpred,npx,nInter,Label="gpred")
          Call mma_allocate(hpred,npx,nInter,nInter,Label="hpred")
          Call mma_allocate(var,npx,Label="var")
          Call mma_allocate(sigma,npx,Label="sigma")
          Call mma_allocate(l,nInter,Label="l")
          Call mma_allocate(ll,int(lb(3)),Label="ll")
          Call mma_allocate(cv,m_t,npx,nInter,nInter,Label="cv")
          Call mma_allocate(cvMatFder,nPoints,npx,Label="cvMatFder")
          Call mma_allocate(cvMatSder,nPoints,npx,Label="cvMatSder")
          Call mma_allocate(cvMatTder,nPoints,npx,Label="cvMatTder")
!
        return
      end
