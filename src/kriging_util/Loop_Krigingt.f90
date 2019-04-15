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

      Subroutine Loop_Kriging(LastqInt)
        use globvar
        Integer nInter,iter
        Real*8 LastqInt(nInter,1),Grad(nInter,iter),Energy(iter)
!
!nx is the n-dimensional vector of the last iteration cumputed in update_sl
! subroutine
        nx = LastqInt(:,1)
!
        ! make_parameters=.False.
        call covarvector(0,iter,nInter) ! for: 0-GEK, 1-Gradient of GEK, 2-Hessian of GEK
        call predict(0,iter,nInter)
!
        ! Energy(iter+1)=pred(npx)
        ! Grad(:,iter+1)=gpred
        write(6,*) 'New values of Energy and grad', pred(npx), gpred
!
        return
      end
