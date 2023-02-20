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
! Copyright (C) 2021, Roland Lindh                                     *
!***********************************************************************
Subroutine Kriging_Update(nQQ,iter,qInt,E_Disp)
Use Slapaf_Info, only: Energy, dqInt
Use Kriging_Mod, only: nSet
Implicit None
Integer nQQ, iter
Real*8  qInt(nQQ), E_Disp

#include "real.fh"
#include "stdalloc.fh"
Real*8, Allocatable :: Aux(:,:), Demp(:), Temp(:)

Call mma_allocate(Temp,nSet,Label='Temp')
Call mma_allocate(Demp,nSet,Label='Demp')
Call mma_allocate(Aux,nQQ,nSet,Label='Aux')

#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: qInt',' ',qInt,nQQ,1)
#endif

Call Energy_Kriging_layer(qInt,Temp,nQQ)

Call Dispersion_Kriging_Layer(qInt,Demp,nQQ)

Call Gradient_Kriging_layer(qInt,Aux,nQQ)

#ifdef _DEBUGPRINT_
Call RecPrt('Kriging_Update: Temp',' ',Temp,1,nSet)
Call RecPrt('Kriging_Update: Demp',' ',Demp,1,nSet)
Call RecPrt('Kriging_Update: Aux',' ',Aux,nQQ,nSet)
#endif

Energy(iter) = Temp(1)

E_Disp = Demp(1)

dqInt(:,iter) = -Aux(:,1)

Call mma_deallocate(Temp)
Call mma_deallocate(Demp)
Call mma_deallocate(Aux)
!
End Subroutine Kriging_Update
