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
      Subroutine Get_Prim_Density_Matrix(ip_D,nBas,ip_D_p,nPrim,TM)

      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Dimension TM(nPrim,nBas)

      Do i=1,nBas
        Do k=1,i-1
          Work(ip_D+i*(i-1)/2+k-1)=Work(ip_D+i*(i-1)/2+k-1)/2.0d0
        EndDo
      EndDo
      Call GetMem('D_p','ALLO','REAL',ip_D_p,nPrim*(nPrim+1)/2)
      Do i=1,nPrim
        Do k=1,i
          TmpDensity=0.0d0
          Do j=1,nBas
            Do l=1,nBas
              TMij=TM(i,j)
              TMkl=TM(k,l)
              If(j.lt.l) Then
                Djl=Work(ip_D+l*(l-1)/2+j-1)
              Else
                Djl=Work(ip_D+j*(j-1)/2+l-1)
              EndIf
              TmpDensity=TmpDensity+TMij*TMkl*Djl
            EndDo
          EndDo
          Work(ip_D_p+i*(i-1)/2+k-1)=TmpDensity
        EndDo
      EndDo

      Return
      End
*
