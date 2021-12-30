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
      SubRoutine Ortho_Orb(Xmo,Smat,nBas,nOrb2Loc,nPass,                &
     &                     Test)
!
!
!     Purpose: orthonormalization of orbitals according to
!
!              V = X^T*S*X
!              X <- X*V^(-1/2)
!
!              where S is the AO overlap matrix.
!              The orthonormalization is carried out nPass times.
!              After this routine, X will satisfy X^T*S*X=1.
!
!
      Implicit Real*8 (a-h,o-z)
      Real*8  Xmo(*), Smat(*)
      Integer nBas, nOrb2Loc
      Logical Test
#include "WrkSpc.fh"

      Character*9 SecNam
      Parameter (SecNam = 'Ortho_Orb')

      Parameter (Tol = 1.0d-10)

      external ddot_

!     Check for quick return.
!     -----------------------

      If (nPass .lt. 1) Return

!     Allocations.
!     ------------

      l_V = nOrb2Loc**2
      l_VSqrt = l_V
      l_VISqrt = l_V
      l_Scr = 2*(nBas**2) + nBas*(nBas+1)/2  ! needed in SqrtMt
      Call GetMem('V','Allo','Real',ip_V,l_V)
      Call GetMem('VSqrt','Allo','Real',ip_VSqrt,l_VSqrt)
      Call GetMem('VISqrt','Allo','Real',ip_VISqrt,l_VISqrt)
      Call GetMem('Scr','Allo','Real',ip_Scr,l_Scr)

!     Orthonormalization passes.
!     --------------------------

      Do iPass = 1,nPass

!        Compute V = X^T*S*X.
!        --------------------

         Call GetUmat_Localisation(Work(ip_V),                          &
     &                             Xmo(1),Smat(1),Xmo(1),               &
     &                             Work(ip_Scr),l_Scr,                  &
     &                             nBas,nOrb2Loc)

!        Compute V^(-1/2).
!        -----------------

         iTask = 2 ! compute sqrt as well as inverse sqrt
         Call SqrtMt(Work(ip_V),nOrb2Loc,iTask,                         &
     &               Work(ip_VSqrt),Work(ip_VISqrt),Work(ip_Scr))

!        Compute orthonormal X <- X*V^(-1/2).
!        ------------------------------------

         nB = max(nBas,1)
         nO2L = max(nOrb2Loc,1)
         Call dCopy_(nBas*nOrb2Loc,Xmo(1),1,                            &
     &                            Work(ip_Scr),1)
         Call DGEMM_('N','N',nBas,nOrb2Loc,nOrb2Loc,                    &
     &                1.0d0,Work(ip_Scr),nB,Work(ip_VISqrt),nO2L,       &
     &                0.0d0,Xmo(1),nB)

      End Do

!     Test orthonormalization (i.e. V=1?).
!     ------------------------------------

      If (Test) Then
         nErr = 0
         Call GetUmat_Localisation(Work(ip_V),                          &
     &                             Xmo(1),Smat(1),Xmo(1),               &
     &                             Work(ip_Scr),l_Scr,                  &
     &                             nBas,nOrb2Loc)
         kOff = ip_V - 1
         Do i = 1,nOrb2Loc
            Work(kOff+nOrb2Loc*(i-1)+i) =                               &
     &      Work(kOff+nOrb2Loc*(i-1)+i) - 1.0d0
         End Do
         xNrm = sqrt(dDot_(nOrb2Loc**2,Work(ip_V),1,Work(ip_V),1))
         If (xNrm .gt. Tol) Then
            Write(6,'(A,A,D16.8,A,I2,A)')                               &
     &      SecNam,': ERROR: ||X^TSX - 1|| = ',xNrm
            nErr = nErr + 1
         End If
         If (nErr .ne. 0) Then
            Write(6,*) SecNam,': failure after ',nPass,' passes'
            Call SysAbendMsg(SecNam,'Orthonormalization failure!',' ')
         End If
      End If

!     De-allocations.
!     ---------------

      Call GetMem('Scr','Free','Real',ip_Scr,l_Scr)
      Call GetMem('VISqrt','Free','Real',ip_VISqrt,l_VISqrt)
      Call GetMem('VSqrt','Free','Real',ip_VSqrt,l_VSqrt)
      Call GetMem('V','Free','Real',ip_V,l_V)

      Return
      End
