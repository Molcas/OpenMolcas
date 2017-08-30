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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************
*  CD_InCore_p
*
*> @brief
*>   Cholesky decompose the symmetric positive (semi-)definite matrix \p X
*> @author Francesco Aquilante
*>
*> @details
*> The \p n &times; \p n matrix \p X is Cholesky decomposed and the resulting
*> \p NumCho Cholesky vectors are returned in array \p Vec. Note that \p X
*> is modified in this routine! The array \p iD(k) on exit contains
*> the index of the diagonal from which the \p k -th Cholesky
*> vector was generated. A non-zero return code signals
*> that an error has occured (\p X might e.g. be non-positive
*> (semi-)definite) and the output is ill-defined.
*>
*> @param[in,out] X      Matrix to be Cholesky decomposed
*> @param[in]     n      Linear dimension of \p X
*> @param[out]    Vec    Storage array for Cholesky vectors
*> @param[in]     MxVec  Max. number of vectors allowed
*> @param[in,out] iD     Index array for parent diagonals
*> @param[out]    NumCho Number of Cholesky vectors
*> @param[in]     Thr    Decomposition threshold (precision)
*> @param[out]    irc    Return code
************************************************************************
      SubRoutine CD_InCore_p(X,n,Vec,MxVec,iD,NumCho,Thr,irc)
      Implicit None
      Integer n, MxVec, NumCho, irc
      Real*8  X(n,n)
      Real*8  Vec(n,MxVec)
      Real*8  Thr
      Integer iD(MxVec)

      Character*11 SecNam
      Parameter (SecNam = 'CD_InCore_p')

      Real*8 DefThr
      Parameter (DefThr = 1.0d-6)
      Real*8 ThrNeg, ThrFail
      Parameter (ThrNeg = -1.0d-13, ThrFail = -1.0d-8)

      Call qEnter(SecNam)

      irc = 0
      NumCho = 0
      If (n .lt. 1) Go To 1 ! exit (nothing to do)
      If (Thr .lt. 0.0d0) Thr = DefThr

      If (MxVec .gt. 0) Then
         Call CD_InCore_1p(X,n,Vec,MxVec,NumCho,Thr,ThrNeg,ThrFail,
     &                     iD,irc)
      Else
         irc = -1
      End If

    1   Continue
       Call qExit(SecNam)
      End
