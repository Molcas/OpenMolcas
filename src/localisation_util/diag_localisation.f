************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2005, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Diag_Localisation(A,EVR,EVI,n,iGetVecs)
C
C     Thomas Bondo Pedersen, Dec. 2005.
C
C     Purpose: diagonalize a real general matrix A and return
C              the real and imaginary part of eigenvalues in EVR and
C              EVI. If iGetVecs=0 no eigenvectors are computed; else,
C              eigenvectors are returned in A (such that column one
C              corresponds to eigenvalue 1 in EVR and EVI). Uses XEIGEN
C              from the linalg_util directory (which uses LAPACK).
C              See LAPACK for details about the storage of real and
C              complex eigenvectors.
C
      Implicit None
      Integer n, iGetVecs
      Real*8  A(n,n), EVR(n), EVI(n)
#include "WrkSpc.fh"

      Character*17 SecNam
      Parameter (SecNam = 'Diag_Localisation')

      Integer iErr
      Integer ip_Vecs, l_Vecs

      l_Vecs = n*n
      Call GetMem('Vecs','Allo','Real',ip_Vecs,l_Vecs)

      iErr = 0
      Call xEigen(iGetVecs,n,n,A,EVR,EVI,Work(ip_Vecs),iErr)
      If (iErr .ne. 0) Then
         Write(6,*) SecNam,': xEigen returned ',iErr
         Call SysAbendMsg(SecNam,'Error in xEigen',' ')
      End If

      If (iGetVecs .ne. 0) Call dCopy_(n*n,Work(ip_Vecs),1,A,1)

      Call GetMem('Vecs','Free','Real',ip_Vecs,l_Vecs)

      End
