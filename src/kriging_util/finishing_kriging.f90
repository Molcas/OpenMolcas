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

      Subroutine Finish_Kriging()
        use globvar
!
!       write(6,*) 'Deallocating all kriging variables'
        deallocate (x,y,dy,nx,l)
        deallocate (Rones)
        deallocate (full_R,rl)
        deallocate (kv,pred,gpred,hpred,var,sigma,cv,ll)
        deallocate (cvMatFder,cvMatSder,cvMatTder)
!
        return
      end
