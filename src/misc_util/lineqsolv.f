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
* Copyright (C) Thomas Bondo Pedersen                                  *
*               2017, Ignacio Fdez. Galvan                             *
************************************************************************
*  LinEqSolv
*
*> @brief
*>   Solve linear equations \f$ Ax=B \f$ or \f$ A^\text{T}x=B \f$ where
*>   \f$ A \f$ is a general nonsingular matrix
*> @author Thomas Bondo Pedersen
*> @modified_by Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> For \p TransA(1:1) = ``'N'`` or ``'n'``,
*> this routine solves the equations \f$ Ax = B \f$ for any number of
*> righthand sides stored as columns of the matrix \f$ B \f$. For
*> \p TransA(1:1) = ``'T'`` or ``'t'``, the equations \f$ A^\text{T}x = B \f$ are
*> solved.
*>
*> The solution vectors are stored as columns of \p B on exit. Note
*> that array \p A will be destroyed during the solution (the matrix
*> is replaced by its factors obtained by Gaussian elimination).
*>
*> Return codes are:
*>
*> - \p irc = ``-1``: input error(s) detected and nothing has been done.
*> - \p irc =  ``0``: successful completion.
*> - \p irc =  ``1``: \p A is estimated to be singular and no solution vectors have been computed.
*>
*> @note
*> 4 * \p nDim real*8 and 2 * \p nDim integer words of memory must be available
*> on entry.
*>
*> @param[out]    irc    Return code
*> @param[in]     TransA Transposition of \p A
*> @param[in,out] A      Array containing the nonsingular matrix \p A on entry and the factorized matrix (Gaussian elimination) on exit
*> @param[in]     ldA    Leading dimension of \p A
*> @param[in,out] B      Array containing righthand sides of equations on entry and solution vectors on exit
*> @param[in]     ldB    Leading dimension of \p B
*> @param[in]     nDim   Column dimension of \p A
*> @param[in]     nEq    Number of equations, i.e. column dimension of \p B
************************************************************************
      SubRoutine LinEqSolv(irc,TransA,A,ldA,B,ldB,nDim,nEq)
      Implicit None
      Character*(*) TransA
      Integer irc, ldA, ldB, nDim, nEq
      Real*8  A(ldA,nDim), B(ldB,nEq)
#include "WrkSpc.fh"

      Integer ip_iPivot, l_iPivot
      Integer ip_Scr, l_Scr
      Integer ip_iScr, l_iScr
      Integer iErr, lTransA
      Real*8  RC,AN
      Character*1 myTransA
      Real*8  DLANGE
      External DLANGE

C     Test input.
C     -----------

      irc = 0
      If (nEq.lt.1 .or. nDim.lt.1) Return
      If (ldA.lt.nDim .or. ldB.lt.nDim) Then
         irc = -1 ! dimension error
         Return
      End If

C     Translate TransA to integer Job to be used by DGESL.
C     ----------------------------------------------------

      lTransA = len(TransA)
      If (lTransA .gt. 0) Then
         myTransA = TransA(1:1)
      Else
         irc = -1 ! TransA error
         Return
      End If

      If (myTransA.ne.'N' .and. myTransA.ne.'n' .and.
     &    myTransA.ne.'T' .and. myTransA.ne.'t') Then
         irc = -1 ! TransA error
         Return
      End If

C     Allocate pivot array and scratch space
C     --------------------------------------

      l_iPivot = nDim
      l_Scr = 4*nDim
      l_iScr = nDim
      Call GetMem('LES_Pivot','Allo','Inte',ip_iPivot,l_iPivot)
      Call GetMem('LES_Scr','Allo','Real',ip_Scr,l_Scr)
      Call GetMem('LES_iScr','Allo','Inte',ip_iScr,l_iScr)

C     Factor A by Gaussian elimination and estimate reciprocal condition
C     number (RC). Check for singularity.
C     ------------------------------------------------------------------

      RC = 0.0d0
      AN = DLANGE('1',nDim,nDim,A,ldA,Work(ip_Scr))
      Call DGETRF_(nDim,nDim,A,ldA,iWork(ip_iPivot),iErr)
      Call DGECON('1',nDim,A,ldA,AN,RC,Work(ip_Scr),iWork(ip_iScr),iErr)
      If (((1.0d0+RC) .eq. 1.0d0) .or. (iErr .gt. 0)) Then
         irc = 1 ! error: A is (probably) singular
         Go To 1 ! exit after deallocations
      End If

C     Solve equations.
C     ----------------

      Call DGETRS_(myTransA,nDim,nEq,A,ldA,iWork(ip_iPivot),B,ldB,iErr)
      If (iErr .gt. 0) Then
         irc = 1 ! error: A is (probably) singular
         Go To 1 ! exit after deallocations
      End If

C     Deallocations.
C     --------------

    1 Call GetMem('LES_Pivot','Free','Inte',ip_iPivot,l_iPivot)
      Call GetMem('LES_Scr','Free','Real',ip_Scr,l_Scr)
      Call GetMem('LES_iScr','Free','Inte',ip_iScr,l_iScr)

      End
