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
      SubRoutine CD_InCore(X,n,Vec,MxVec,NumCho,Thr,irc)
************************************************************
*
*   <DOC>
*     <Name>CD\_InCore</Name>
*     <Syntax>Call CD\_InCore(X,n,Vec,MxVec,NumCho,Thr,irc)</Syntax>
*     <Arguments>
*       \Argument{X}{Matrix to be Cholesky decomposed, dimension
*                    X(n,n)}{Real*8}{inout}
*       \Argument{n}{Linear dimension of X}{Integer}{in}
*       \Argument{Vec}{Storage array for Cholesky vectors, dimension
*                      Vec(n,MxVec)}{Real*8}{out}
*       \Argument{MxVec}{Max. number of vectors allowed}{Integer}{in}
*       \Argument{NumCho}{Number of Cholesky vectors}{Integer}{out}
*       \Argument{Thr}{Decomposition threshold (precision)}
*                {Real*8}{in}
*       \Argument{irc}{Return code}{Integer}{out}
*     </Arguments>
*     <Purpose>Cholesky decompose the symmetric positive (semi-)
*              definite matrix X</Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*        The n-by-n matrix X is Cholesky decomposed and the resulting
*        NumCho Cholesky vectors are returned in array Vec. Note that X
*        is modified in this routine! A non-zero return code signals
*        that an error has occured (X might f.ex. be non-positive
*        (semi-) definite) and the output is ill-defined.
*     </Description>
*    </DOC>
*
************************************************************

      Implicit None
      Integer n, MxVec, NumCho, irc
      Real*8  X(n,n)
      Real*8  Vec(n,MxVec)
      Real*8  Thr

      Character*9 SecNam
      Parameter (SecNam = 'CD_InCore')

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
         Call CD_InCore_1(X,n,Vec,MxVec,NumCho,Thr,ThrNeg,ThrFail,irc)
      Else
         irc = -1
      End If

    1   Continue
       Call qExit(SecNam)
      End
