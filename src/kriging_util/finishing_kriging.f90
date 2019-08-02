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
#include "stdalloc.fh"
!
!       write(6,*) 'Deallocating all kriging variables'
        Call mma_deallocate(x)
        Call mma_deallocate(y)
        Call mma_deallocate(nx)
        Call mma_deallocate(dy)
        Call mma_deallocate(full_R)
        Call mma_deallocate(Ys)
        Call mma_deallocate(dl)
        Call mma_deallocate(rl)
        Call mma_deallocate(Rones)
        Call mma_deallocate(kv)
        Call mma_deallocate(pred)
        Call mma_deallocate(gpred)
        Call mma_deallocate(hpred)
        Call mma_deallocate(var)
        Call mma_deallocate(sigma)
        Call mma_deallocate(l)
        Call mma_deallocate(ll)
        Call mma_deallocate(cv)
        Call mma_deallocate(cvMatFder)
        Call mma_deallocate(cvMatSder)
        Call mma_deallocate(cvMatTder)
!
        return
      end
