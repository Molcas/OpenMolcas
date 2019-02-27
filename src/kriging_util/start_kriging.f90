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

      Subroutine Start_Kriging(iter,nInter,qInt,Grad)
        use globvar
        Real*8 qInt(nInter,iter+1), Grad(nInter,iter), Energy(iter)
        Write (6,*) 'Kriging values in Start Kriging'
        Write (6,*) 'iter', iter
        Write (6,*) 'nInter', nInter
        Write (6,*) 'Grad: ',Grad
        Write (6,*) 'Grad size: ',size(Grad)
        Write (6,*) 'Grad shape: ',shape(Grad)
        Write (6,*) 'Energy: ',Energy
        Write (6,*) 'Coord: ',qInt
        Write (6,*) 'Coord size: ',size(qInt)
        Write (6,*) 'Coord shape: ',shape(qInt)
      end