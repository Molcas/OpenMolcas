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

subroutine Start_Kriging(nPoints_In,nInter_In,x_,dy_,y_)

use kriging_mod

implicit none
#include "stdalloc.fh"

!    nPoints: the number of sample points, n
!    nInter: the dimensionality of the function, d
!    y_: the values of the function at the sample points
!    dy_: the gradient of the function at the sample points
!    x_: the coordinates of the sample points

integer nInter_In, nPoints_In
real*8 x_(nInter_In,nPoints_In)
real*8 y_(nPoints_In)
real*8 dy_(nInter_In,nPoints_In)

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('Start_Kriging: x',' ',x_,nInter_In,nPoints_In)
call RecPrt('Start_Kriging: y',' ',y_,1,nPoints_In)
call RecPrt('Start_Kriging: dy',' ',dy_,nInter_In,nPoints_In)
#endif

! Call Setup_Kriging to store the data in some internally protected arrays and scalars.

call Prep_Kriging(nPoints_In,nInter_In,x_,dy_,y_)

! Development code for partial gradient enhanced Kriging (PGEK) based on Mutual Information between
! the coordinates and the energy.

if (PGEK_On .and. nPoints >= 2) call PGEK()

! m_t is the dimentionality of the square correlation matrix Gradient-Psi
! (equation (2) on:
!-------- ref. = doi:10.1007/s00366-015-0397-y)-------
!
! In the case the nPoint_v last energies and nPoint_g last gradients were use
! m_t would be computed as

#ifdef _DEBUGPRINT_
write(6,*) 'nD=',nD
write(6,*) 'nPoints_v,nPoints_g=',nPoints,nPoints-nD
#endif
m_t = nPoints+nInter_Eff*(nPoints-nD)

! full_R correspond to the gradient of Psi (eq. (2) ref.)

call mma_Allocate(full_R,m_t,m_t,Label='full_R')
call mma_Allocate(full_RInv,m_t,m_t,Label='full_RInv')

if (mblAI) sbmev = y(maxloc(y,dim=1))

! Allocate x0, which is a n-dimensional vector of the coordinats of the last iteration computed

call mma_Allocate(x0,nInter,Label='nx')

! rl and dl are temporary matrices for the contruction of Psi which is inside of
! Grad-Psi (eq.(2) ref.) dl=rl^2=Sum[i] [(x_i-x0_i)/l)^2]
! more inoformation is given in subsequen files.
! Mat is the final matrix after the distance (between source data rl and dl) has
! passed through the Matern correlation function (ISBN 0-486-61272-4 & eq. (11-12)
! ref.).
!Iden is just an identity matrix necesary to avoid that the Grad-Psi becomes
! Singular after been multiplied by EPS factor

call mma_Allocate(rl,nPoints,nInter,Label='rl')
call mma_Allocate(dl,nPoints,Label='dl')
call mma_Allocate(Rones,m_t,Label='Rones')

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

! Allocate additional variables needed for the kriging.

call mma_allocate(kv,m_t,Label='kv')
call mma_allocate(gpred,nInter,Label='gpred')
call mma_allocate(hpred,nInter,nInter,Label='hpred')
call mma_allocate(l,nInter,Label='l')
call mma_allocate(ll,int(lb(3)),Label='ll')

call mma_allocate(cv,m_t,nInter,nInter,Label='cv')
call mma_allocate(cvMatFder,nPoints,Label='cvMatFder')
call mma_allocate(cvMatSder,nPoints,Label='cvMatSder')
call mma_allocate(cvMatTder,nPoints,Label='cvMatTder')

return

end subroutine Start_Kriging
