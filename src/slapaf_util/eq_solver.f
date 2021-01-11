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
#include "stdalloc.fh"
#include "warnings.fh"
      Real*8 B(M,*),Degen(M),dSS(*),DFC(*)
      Logical Curvilinear
      Character(LEN=1) Mode
      Real*8 Temp(1)
      Real*8, Allocatable:: A(:), Btmp(:), Work(:)
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Solve the equation Ax=b
*
      LDA=M
      LDB=Max(1,M,N)
      If (Mode.eq.'T') Then
         Call mma_allocate(A,M*M,Label='A')
         A(:)=Zero
         call dcopy_(M*N,B,1,A,1)
         If (.NOT.Curvilinear) Then
            Do i = 1, M
               Call DScal_(M,Sqrt(Degen(i)),A(i),M)
            End Do
         End If
#ifdef _DEBUGPRINT_
         Call RecPrt('A',' ',A,M,M)
#endif
      Else
         Call mma_allocate(A,M*N,Label='A')
         call dcopy_(M*N,B,1,A,1)
         If (.NOT.Curvilinear) Then
            Do i = 1, M
               Call DScal_(N,Sqrt(Degen(i)),A(i),M)
            End Do
         End If
#ifdef _DEBUGPRINT_
         Call RecPrt('A',' ',A,M,N)
#endif
      End If
*
      Call mma_allocate(Btmp,LDB*NRHS,Label='Btmp')
      Btmp(:)=Zero
*
      jpB=1
      If (Mode.eq.'T') Then
         Do iRHS = 1, nRHS
            ij = (iRHS-1)*N + 1
            call dcopy_(N,dss(ij),1,Btmp(jpB),1)
            jpB=jpB+LDB
         End Do
      Else
#ifdef _DEBUGPRINT_
         Call RecPrt('B(raw)',' ',dss,M,nRHS)
#endif
         Do iRHS = 1, nRHS
            If (.Not.Curvilinear) Then
               Do i = 0, M-1
                  ij = (iRHS-1)*M + i+1
                  Btmp(jpB+i) = dss(ij)*Sqrt(Degen(i+1))
               End Do
            Else
               ij = (iRHS-1)*M + 1
               call dcopy_(M,dss(ij),1,Btmp(jpB),1)
            End If
            jpB=jpB+LDB
         End Do
      End If
#ifdef _DEBUGPRINT_
      Call RecPrt('B(in)',' ',Btmp,LDB,NRHS)
#endif
*
      LWork=-1
      call dgels_(Mode,M,N,NRHS,A,LDA,Btmp,LDB,Temp,LWork,INFO)
*     Write (6,*) 'Temp,INFO=',Temp,INFO
      LWork=INT(Temp(1))
      Call mma_allocate(Work,LWork,Label='Work')
      call dgels_(Mode,M,N,NRHS,A,LDA,Btmp,LDB,Work,LWork,INFO)
      If (INFO.gt.0) Then
         Call WarningMessage(2,'Error in Eq_Solver')
         Write (6,*)
         Write (6,*) '***********************************************'
         Write (6,*) ' ERROR: Eq_Solver could not find a solution.   '
         Write (6,*) ' The matrix is rank deficient.                 '
         Write (6,*) '***********************************************'
         Call Quit(_RC_INTERNAL_ERROR_)
      End If
*
#ifdef _DEBUGPRINT_
      Call RecPrt('B(out)',' ',Btmp,LDB,NRHS)
#endif
      jpB=1
      If (Mode.eq.'T') Then
         Do iRHS = 1, nRHS
            If (.Not.Curvilinear) Then
               Do i = 0, M-1
                  Btmp(jpB+i) = Btmp(jpB+i)/Sqrt(Degen(i+1))
               End Do
            End If
            ij = (iRHS-1)*M + 1
            call dcopy_(M,Btmp(jpB),1,DFC(ij),1)
            jpB=jpB+LDB
         End Do
      Else
         Do iRHS = 1, nRHS
            ij = (iRHS-1)*N + 1
            call dcopy_(N,Btmp(jpB),1,DFC(ij),1)
            jpB=jpB+LDB
         End Do
      End If
*
      Call mma_deallocate(Work)
      Call mma_deallocate(Btmp)
      Call mma_deallocate(A)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
