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
      subroutine mult_sro(A,nA,C,nC,B,nB,Fact,Res,Tmp)
      Implicit Real*8(a-h,o-z)
      Real*8 A(*),B(*),C(*),Res(*),Tmp(*)

      Call DGEMM_('N','N', nA,nC,nC,
     &           1.0d0,A,nA, C,nC,
     &           0.0d0,Tmp,nA)
      Call DGEMM_('N','N', nA,nB,nC,
     &           Fact,Tmp,nA, B,nC,
     &           1.0d0,Res,nA)
      return
      end
