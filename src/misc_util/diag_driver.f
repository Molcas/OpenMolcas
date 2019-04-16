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
      Subroutine Diag_Driver(JobZ, Range, UpLo, nDim, Triangular, Aux,
     &                       lDimAux, vLower, vUpper, iLower, iUpper,
     &                       EigVal, EigVec, lDimVec, iUnit_Matrix,
     &                       iSort, Method, nFound, Ierr)
      Implicit None
#include "WrkSpc.fh"
* Input/Output variables
      Real*8 EigVal, EigVec, Aux, vLower, vUpper, Triangular
      Integer nDim, iLower, iUpper, nFound, iErr, lDimAux, lDimVec
      Integer iUnit_Matrix, iSort
      Character JobZ, Range, UpLo, Method
      Dimension EigVal(*), EigVec(*), Aux(*), Triangular(*)
* Local Variables
      Real*8 Tollerance
      Integer iMethod, nDim2, iSize, lWork, liWork,liW(2),liErr(2)
* Scratch and pointers to scratch arrays
      Real*8 Work_L(1)
      Integer lScr, liScr, iSuppZ
* External Functions ..
      Logical            lSame
      Real*8             dLamCh_
      External           lSame, dLamCh_
*
* Sizes for arrays
*
      nDim2 = nDim * nDim
      iSize = nDim * (nDim + 1) / 2
*
* Determine which algorithm to use
*
      iMethod = 0
      If (lSame(Method, 'A')) Then
         iMethod = 0
      Else If (lSame(Method, 'Q')) Then
         iMethod = 1
      Else If (lSame(Method, 'J')) Then
         iMethod = 2
      Else
         Write(6,*) '!!! Diag_Driver called with an unknown method: ',
     &              Method
         Write(6,*) '!!! Supported methods: Q, J, and A'
         Write(6,*) '    Method = ''',Method,''''
         Call Abend()
      End If
*
      If (iMethod .le. 1) Then
*
* Use the QL algorithm (dSyevR)
*
         Call Square(Triangular,Aux,lDimAux,1,nDim)
         Call dCopy_(nDim2,[0.0D0],0,EigVec,1)
         Call dCopy_(nDim,[1.0D0],0,EigVec,nDim+1)
*
* Determine safe tollerance for dSyevR
*
         Tollerance = dLamCh_( 'Safe minimum' )
*
* Determine optimal sizes of scratch arrays
*
         CALL GetMem('ISUPPZ  ','ALLO','INTE',iSuppZ,2*nDim)
         lWork  = -1
         liWork = -1
CC AOM 03.08.2005 - Added LiW(2), otherwise on Opteron crashed
CC AOM 04.08        Also added liErr for the same reason
         call dsyevr_(JobZ, Range, UpLo, nDim, Aux, lDimAux, vLower,
     &               vUpper, iLower, iUpper, Tollerance, nFound, EigVal,
     &               EigVec, lDimVec, iWork(iSuppZ), Work_L, lWork,
     &               liW, liWork, liErr(1))
         lWork = Int(Work_L(1))
         liWork = liW(1)
*
* Allocate scratch arrays
*
         CALL GetMem('SCRATCH ','ALLO','REAL',lScr,lWork)
         CALL GetMem('ISCRATCH','ALLO','INTE',liScr,liWork)
*
* Run actual QL routine
*
         call dsyevr_(JobZ, Range, UpLo, nDim, Aux, lDimAux, vLower,
     &               vUpper, iLower, iUpper, Tollerance, nFound, EigVal,
     &               EigVec, lDimVec, iWork(iSuppZ), Work(lScr), lWork,
     &               iWork(liScr), liWork, liErr(1))
         iErr=liErr(1)
*
* Free scratch
*
         CALL GetMem('SCRATCH ','FREE','REAL',lScr,lWork)
         CALL GetMem('ISCRATCH','FREE','INTE',liScr,liWork)
         CALL GetMem('ISUPPZ  ','FREE','INTE',iSuppZ,2*nDim)
*
* Check for convergence
*
         If (iErr .ne. 0) Then
            Write(6,*) '!!! No Convergence in the QL algorithm.'
            If (lSame(Method, 'A')) Then
               Write(6,*) '!!! Trying Jacobi instead.'
               Write(6,*) '!!! Warning: This might be very slow.'
               iMethod = 2
            Else
               Call Abend()
            End If
         Else
            Call Chk4NAN(nDim**2,EigVec,Ierr)
            If (iErr .gt. 0) Then
               Write(6,*) 'At least one of the eigenvectors found with'
               Write(6,*) 'DSYEVR contained a NAN.'
               If (lSame(Method, 'A')) Then
                  Write(6,*) 'Trying Jacobi instead.'
                  Write(6,*) 'Warning: This might be very slow.'
                  iMethod = 2
               Else
                  Call Abend()
               End If
            End If
         End If
      End If
      If (iMethod .eq. 2) Then
*
* Use the Jacobi algorithm
*
         Call dCopy_(iSize,Triangular,1,Aux,1)
         If (iUnit_Matrix .eq. 1) Then
            Call dCopy_(nDim2,[0.0D0],0,EigVec,1)
            Call dCopy_(nDim,[1.0D0],0,EigVec,nDim+1)
         End If
         Call Jacob(Aux,EigVec,nDim,lDimVec)
         Call vEig(nDim,Aux,EigVal)
      End If
*
* Sort the eigenvalues and eigenvectors?
*
      If (iSort .eq. 1) Then
         Call JacOrd2(EigVal, EigVec, nDim, lDimVec)
      Else If (iSort .eq. -1) Then
         Call SortEig(EigVal, EigVec, nDim, lDimVec)
      End If
      Return
      End
