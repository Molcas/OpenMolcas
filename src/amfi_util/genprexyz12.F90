!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      subroutine genprexyz12(preY)
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
      Dimension preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
!bs #####################################################################
!bs   additional (-) signs from the (-i) factors  in the
!bs   (-) linear combinations   (see tosigX(Y,Z).f)
!bs #####################################################################
!bs   - -  + +  >   -
      do M4=0,Lmax
      do M3=0,Lmax
      do M2=-Lmax,-1
!      do M1=-Lmax,-1
!              preY(m1,m2,m3,m4)=-preY(m1,m2,m3,m4)
            call dscal_(Lmax,-1d0,preY(-Lmax,m2,m3,m4),1)
!      enddo
      enddo
      enddo
      enddo
      return
      end
