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
************************************************************************
      SubRoutine LinEqSolv(irc,TransA,A,ldA,B,ldB,nDim,nEq)
************************************************************
*
*   <DOC>
*     <Name>LinEqSolv</Name>
*     <Syntax>Call LinEqSolv(irc,TransA,A,ldA,B,ldB,nDim,nEq)</Syntax>
*     <Arguments>
*       \Argument{irc}{Return code}{Integer}{out}
*       \Argument{TransA}{Transposition of A}{Character*(*)}{in}
*       \Argument{A}{Array containing the nonsingular matrix A on entry
*                    and the factorized matrix (Gaussian elimination) on
*                    exit}{Real*8(ldA,nDim)}{inout}
*       \Argument{ldA}{Leading dimension of A}{Integer}{in}
*       \Argument{B}{Array containing righthand sides of equations on
*                    entry and solution vectors on exit}
*                {Real*8(ldB,nEq)}{inout}
*       \Argument{ldB}{Leading dimension of B}{Integer}{in}
*       \Argument{nDim}{Column dimension of A}{Integer}{in}
*       \Argument{nEq}{Number of equations, i.e. column dimension of B}
*                {Integer}{in}
*     </Arguments>
*     <Purpose>Solve linear equations Ax=B or A(T)x=B where A is a
*              general nonsingular matrix
*     </Purpose>
*     <Dependencies>nDim real*8 and nDim integer words of memory must be
*     available on entry. Makes use of routines DGECO and DGESL.
*     </Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        For TransA(1:1)='N' or 'n',
*        this routine solves the equations Ax = B for any number of
*        righthand sides stored as columns of the matrix B. For
*        TransA(1:1)='T' or 't', the equations A(Transposed)x = B are
*        solved.
*        The solution vectors are stored as columns of B on exit. Note
*        that array A will be destroyed during the solution (the matrix
*        is replaced by its factors obtained by Gaussian elimination).
*        Return codes are:
*           irc=-1: input error(s) detected and nothing has been done,
*           irc=0: successful completion,
*           irc=1: A is estimated to be singular and no solution vectors
*                  have been computed.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Character*(*) TransA
      Integer irc, ldA, ldB, nDim, nEq
      Real*8  A(ldA,nDim), B(ldB,nEq)
#include "WrkSpc.fh"

      Integer ip_iPivot, l_iPivot
      Integer ip_Scr, l_Scr
      Integer iColB, Job, lTransA
      Real*8  RC
      Character*1 myTransA

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

      If (myTransA.eq.'N' .or. myTransA.eq.'n') Then
         Job = 0
      Else If (myTransA.eq.'T' .or. myTransA.eq.'t') Then
         Job = 1
      Else
         irc = -1 ! TransA error
         Return
      End If

C     Allocate pivot array and scratch space for dGeco.
C     -------------------------------------------------

      l_iPivot = nDim
      l_Scr = nDim
      Call GetMem('LES_Pivot','Allo','Inte',ip_iPivot,l_iPivot)
      Call GetMem('LES_Scr','Allo','Real',ip_Scr,l_Scr)

C     Factor A by Gaussian elimination and estimate reciprocal condition
C     number (RC). Check for singularity.
C     ------------------------------------------------------------------

      RC = 0.0d0
      Call DGECO(A,ldA,nDim,iWork(ip_iPivot),RC,Work(ip_Scr))
      If ((1.0d0+RC) .eq. 1.0d0) Then
         irc = 1 ! error: A is (probably) singular
         Go To 1 ! exit after deallocations
      End If

C     Solve equations.
C     ----------------

      Do iColB = 1,nEq
         Call DGESL(A,ldA,nDim,iWork(ip_iPivot),B(1,iColB),Job)
      End Do

C     Deallocations.
C     --------------

    1 Call GetMem('LES_Pivot','Free','Inte',ip_iPivot,l_iPivot)
      Call GetMem('LES_Scr','Free','Real',ip_Scr,l_Scr)

      End
