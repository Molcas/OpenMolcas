************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine genprexyz6(preY,preXZ)
      implicit real*8(a-h,o-z)
#include "para.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
      Dimension preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax),
     *preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
cbs #####################################################################
cbs   prefactors preXZ und preY include the factors 1/root(2)
cbs   for the +/- linear combinations of spherical harmonics
cbs #####################################################################
      do M4=-Lmax,Lmax
      do M3=-Lmax,Lmax
      do M2=-Lmax,Lmax
      do M1=-Lmax,Lmax
              preY(m1,m2,m3,m4)=preXZ(m1,m2,m3,m4)
      enddo
      enddo
      enddo
      enddo
      return
      end
