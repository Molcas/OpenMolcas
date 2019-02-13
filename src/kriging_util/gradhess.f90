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
subroutine gradhess
   use globvar
!   use ogpf
!   type(gpf):: gp
   !Calculating the derivatives of the Kriging function stored as pred
   allocate (gpred(size(pred)),hpred(size(pred)))
    Write(6,*) 'Kriging',pred
    call deriv(pred,gpred,1,size(pred,1),0)
    Write(6,*) 'Grad Kriging',gpred
    call deriv(gpred,hpred,1,size(pred,1),0)
    Write(6,*) 'Hessian Kriging',hpred
!    call gp%multiplot(1,2)
!    call gp%plot(nx, pred, 't "Pred" w lines lc "blue" lt 1 lw 2','', &
!            nx, ny, 't "Original" w lines lc "red" lt 0 lw 2','', &
!            x, y, 't "Data Points" w points pt 6','', &
!            nx, gpred, 't "Gradient" w lines lc "green" lt 0 lw 2')
!    call gp%plot(nx, pred, 't "Pred" w lines lc "blue" lt 1 lw 2','', &
!            nx, ny, 't "Original" w lines lc "red" lt 0 lw 2','', &
!            x, y, 't "Data Points" w points pt 6','', &
!            nx, hpred, 't "Hessian" w lines lc "green" lt 0 lw 2')
end
