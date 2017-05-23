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
      Subroutine Ttotal(T1,T2,T3,T4,Ttot,Ttot_Inv,nDim)
      Implicit ReaL*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Real*8 T1(nDim*nDim),T2(nDim*nDim),T3(nDim*nDim),
     &       T4(nDim,nDim),Ttot(nDim,nDim),
     &       Ttot_Inv(nDim,nDim)
*
      Call Allocate_Work(ipTemp,nDim**2)
      Call Allocate_Work(ipTemp2,nDim**2)
      Call Ttotal_(T1,T2,T3,T4,Ttot,Ttot_Inv,nDim,Work(ipTemp),
     &            Work(ipTemp2))
      Call Free_Work(ipTemp2)
      Call Free_Work(ipTemp)
*
      Return
      End
      Subroutine Ttotal_(T1,T2,T3,T4,Ttot,Ttot_Inv,nDim,Temp,Temp2)
      Implicit Real*8 (A-H,O-Z)
      Real*8 T1(nDim*nDim),T2(nDim*nDim),T3(nDim*nDim),
     &       T4(nDim,nDim),Ttot(nDim,nDim),
     &       Ttot_Inv(nDim,nDim),
     &       Temp(nDim,nDim), Temp2(nDim,nDim)
*                                                                      *
************************************************************************
*                                                                      *
*     Ttot=T1*T2*T3*T4
*
clg   write (*,*) 'Ttotal ', nDim
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,T1,nDim,
     &            T2,nDim,
     &            0.0d0,Temp,nDim)
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,Temp,nDim,
     &            T3,nDim,
     &            0.0d0,Temp2,nDim)
      Call DGEMM_('N','N',
     &            nDim,nDim,nDim,
     &            1.0d0,Temp2,nDim,
     &            T4,nDim,
     &            0.0d0,Ttot,nDim)
clg   Call RecPrt('T_TOT',' ',Ttot,nDim,nDim)
      Call MINV(Ttot,Ttot_Inv,ISING,DET,nDim)
*
      Return
      End
