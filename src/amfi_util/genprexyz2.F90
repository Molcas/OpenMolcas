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
      subroutine genprexyz2(preXZ)
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
      Dimension preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
      roottwo=sqrt(2d0)
!bs #####################################################################
!bs   prefactors preXZ und preY include the factors 1/root(2)
!bs   for the +/- linear combinations of spherical harmonics
!bs #####################################################################
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
              preXZ(m1,m2,m3,0)=preXZ(m1,m2,m3,0)*roottwo
      enddo
      enddo
      enddo
      return
      end
