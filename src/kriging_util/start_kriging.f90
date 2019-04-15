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

      Subroutine Start_Kriging(iter,nInter,x,temp_dy,y)
        use globvar
        Integer nInter,iter
        Real*8 x(nInter,iter),temp_dy(nInter,iter),y(iter),dy(nInter*iter)
!
!       if (.NOT. allocated(x)) then
          allocate (x(nInter,iter),y(iter),dy(nInter*iter), &
                    nx(nInter,1))
!m_t is the dimentionality of the square correlation matrix Gradient-Psi
! (equation (2) DOI 10.1007/s00366-015-0397-y)
          m_t=iter*(1+nInter)
!npx is the number of new points (Energy and Gradient) to be predict
! according to the iteration that was computed in update_sl subroutine
          npx = 1
!full_R correspond to the gradient of Psi (eq. (2) DOI 10.1007/s00366-015-0397-y)
          allocate (full_R(m_t,m_t))
!nx is the n-dimensional vector of the last iteration cumputed in update_sl
! subroutine
          nx = qInt(:,iter+1:iter+1)
!x is the n-dimensional internal coordinates
          x = qInt(1:nInter,1:iter)
!y is the energy
          y = Energy
!dy it's a vector of Grad-y (eq. (5)  DOI 10.1007/s00366-015-0397-y) gradients of
! the energy with respect to the internal coordinates
          do i=1,nInter
            do j=1,iter
              dy(j+(i-1)*iter) = temp_dy(i,j)
            enddo
          enddo
!rl and dl are temporary matrices for the contruction of Psi which is inside of
! Grad-Psi (eq.(2) DOI 10.1007/s00366-015-0397) dl=rl^2=Sum[i] [(x_i-x0_i)/l)^2]
! more inoformation is given in subsequen files.
!Mat is the final matrix after the distance (between source data rl and dl) has
! passed through the Mat'ern correlation function (ISBN 0-486-61272-4 & eq. (11-12)
! DOI 10.1007/s00366-015-0397).
!Iden is just an identity matrix necesary to avoid that the Grad-Psi becomes
! Singular
          allocate (rl(iter,iter),dl(iter,iter), &
                        mat(iter,iter),Iden(iter,iter))
!kv is the vector that contains the dot product of the inverse of Grad-Psi and
!Grad-y minus the dot product of the inverse of Grad-Psi and f-ones multiplied
! by the constant Grad-Trend function (eq. (3), (5), (6) and (7)
!DOI 10.1007/s00366-015-0397).
!Pred, gpred and hpred are the predicted energy, gradient and hessian according
! to the coordinates given by update_sl
!var is part of the GEK variance and the concentrated Likelihood function (eq.
! (8) and (9) DOI 10.1007/s00366-015-0397).
!Sigma is a dispersion function +-sigma=1.96*sqrt(abs(var*variance)) (need reference).
!cv is the corvariant vector that contains the correlation of unction values and
! gradients between the sample data and the prediction
! (eq. 4 DOI 10.1007/s00366-015-0397).
!l is a n-dimensional vector of the width of the Mat'ern function.
!ll is the likelihood function.
          allocate (kv(m_t),pred(npx),gpred(npx),hpred(npx),var(npx), &
          sigma(npx),cv(m_t,npx,nInter), l(nInter),ll(3,int(lb(3))))
          make_parameters=.True.
        ! else
        !   nx = qInt(:,iter+1:iter+1)
        !   make_parameters=.False.
        ! endif
!
        call kernels(iter,nInter)
!
!         Energy(iter+1)=pred(npx)
!         Grad(:,iter+1)=gpred
!         write(6,*) 'New values of Energy and grad', pred(npx), gpred
! !
        return
      end
