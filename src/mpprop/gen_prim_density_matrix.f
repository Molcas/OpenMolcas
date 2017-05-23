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
      Subroutine Gen_Prim_Density_Matrix(nBas,nPrim,ip_D_p,nOcOb,
     &           oNum,oCof)

      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Dimension oNum(nBas)
      Dimension oCof(nBas,nPrim)

      Call GetMem('D_p','ALLO','REAL',ip_D_p,nPrim*(nPrim+1)/2)
       Do K = 1 , nPrim
        Do L = 1 , K
         Work(ip_D_p+k*(k-1)/2+l-1) = 0.0d0
        EndDo
       EndDo
       Do K=1,nPrim
         Do L=1,K
           Do i=1,nOcOb
             Work(ip_D_p+k*(k-1)/2+l-1) = Work(ip_D_p+k*(k-1)/2+l-1) +
     &       oNum(I)*oCof(I,K)*oCof(I,L)
           EndDo
         EndDo
       EndDo


      Return
      End
*
