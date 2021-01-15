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
* Copyright (C) 2005, Per-Olof Widmark                                 *
************************************************************************
      Subroutine NIdiag(H,U,n,nv,iOpt)
************************************************************************
*                                                                      *
* This routine is a wrapper that calls appropriate routines to         *
* perform diagonalization of of symmetric matrices stored in lower     *
* triangular form.                                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written November 2005                                                *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
* n    - Dimension of matrix                                           *
* nv   - Length of eigenvectors nv>=n                                  *
* H    - Matrix to be diagonalized                                     *
* U    - Eigenvectors                                                  *
* iOpt - Option flag, for future improvements.                         *
*----------------------------------------------------------------------*
      External OrbPhase
      Integer n,nv,iOpt
      Real*8  H(*), U(nv,n)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer ierr, i
      Real*8  Tmp, OrbPhase
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      If (n.eq.0) Return
      Call Givens(H,U,n,nv)
      Call QLdiag(H,U,n,nv,ierr)
      If(ierr.eq.1) Then
c         Write (6,*) 'NIdiag: backup call to Jacob!'
         CALL Jacob(H,U,n,nv)
      End If
      Do i = 1, n
         Tmp = OrbPhase(U(1,i),nV)
      End Do
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
#ifdef _WARNING_WORKAROUND_
      If (.False.) Then
         Call Unused_integer(iOpt)
         Call Unused_real(Tmp)
      End If
#endif
      End
