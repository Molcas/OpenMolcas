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
      Subroutine Eq_Solver(Mode,M,N,NRHS,B,Curvilinear,Degen,dSS,DFC)
      Implicit Real*8 (a-h,o-z)
************************************************************************
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 B(M,*),Degen(M),dSS(*),DFC(*)
      Logical Curvilinear
      Character*1 Mode
      Dimension Temp(1)
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Solve the equation Ax=b
*
      LDA=M
      LDB=Max(1,M,N)
      If (Mode.eq.'T') Then
         Call Allocate_Work(ipA,M*M)
         Call FZero(Work(ipA),M*M)
         call dcopy_(M*N,B,1,Work(ipA),1)
         If (.NOT.Curvilinear) Then
            Do i = 0, M-1
               Call DScal_(M,Sqrt(Degen(i+1)),Work(ipA+i),M)
            End Do
         End If
#ifdef _DEBUG_
         Call RecPrt('A',' ',Work(ipA),M,M)
#endif
      Else
         Call Allocate_Work(ipA,M*N)
         call dcopy_(M*N,B,1,Work(ipA),1)
         If (.NOT.Curvilinear) Then
            Do i = 0, M-1
               Call DScal_(N,Sqrt(Degen(i+1)),Work(ipA+i),M)
            End Do
         End If
#ifdef _DEBUG_
         Call RecPrt('A',' ',Work(ipA),M,N)
#endif
      End If
*
      Call Allocate_Work(ipB,LDB*NRHS)
      Call FZero(Work(ipB),LDB*NRHS)
*
      jpB=ipB
      If (Mode.eq.'T') Then
         Do iRHS = 1, nRHS
            ij = (iRHS-1)*N + 1
            call dcopy_(N,dss(ij),1,Work(jpB),1)
            jpB=jpB+LDB
         End Do
      Else
#ifdef _DEBUG_
         Call RecPrt('B(raw)',' ',dss,M,nRHS)
#endif
         Do iRHS = 1, nRHS
            If (.Not.Curvilinear) Then
               Do i = 0, M-1
                  ij = (iRHS-1)*M + i+1
                  Work(jpB+i) = dss(ij)*Sqrt(Degen(i+1))
               End Do
            Else
               ij = (iRHS-1)*M + 1
               call dcopy_(M,dss(ij),1,Work(jpB),1)
            End If
            jpB=jpB+LDB
         End Do
      End If
#ifdef _DEBUG_
      Call RecPrt('B(in)',' ',Work(ipB),LDB,NRHS)
#endif
*
      LWork=-1
      call dgels_(Mode,M,N,NRHS,Work(ipA),LDA,Work(ipB),LDB,
     &           Temp,LWork,INFO)
*     Write (6,*) 'Temp,INFO=',Temp,INFO
      LWork=INT(Temp(1))
      Call Allocate_Work(ipWork,LWork)
      call dgels_(Mode,M,N,NRHS,Work(ipA),LDA,Work(ipB),LDB,
     &           Work(ipWork),LWork,INFO)
*
#ifdef _DEBUG_
      Call RecPrt('B(out)',' ',Work(ipB),LDB,NRHS)
#endif
      jpB=ipB
      If (Mode.eq.'T') Then
         Do iRHS = 1, nRHS
            If (.Not.Curvilinear) Then
               Do i = 0, M-1
                  Work(jpB+i) = Work(jpB+i)/Sqrt(Degen(i+1))
               End Do
            End If
            ij = (iRHS-1)*M + 1
            call dcopy_(M,Work(jpB),1,DFC(ij),1)
            jpB=jpB+LDB
         End Do
      Else
         Do iRHS = 1, nRHS
            ij = (iRHS-1)*N + 1
            call dcopy_(N,Work(jpB),1,DFC(ij),1)
            jpB=jpB+LDB
         End Do
      End If
*
      Call Free_Work(ipWork)
      Call Free_Work(ipB)
      Call Free_Work(ipA)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
