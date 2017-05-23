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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Integer Function Cho_X_GetTol(iTolDef)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_GetTol</Name>
*     <Syntax>Cho\_X\_GetTol(iTolDef)</Syntax>
*     <Arguments>
*       \Argument{iTolDef}{Default tolerance}{Integer}{in}
*     </Arguments>
*     <Purpose>Set tolerance integer for use with Add\_Info (for
*              verification).
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by>Thomas Bondo Pedersen, October 2010: LDF support.
*     </Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>The tolerance in verification might depend on the
*                  threshold of the Cholesky decomposition. The integer
*                  used to specify tolerance in Add\_Info is computed
*                  here according to the formula:
*                     iTol = -log(Thr)
*                  where Thr is the threshold and "log" is the logarithm
*                  (base 10).
*                  If LDF (local DF) is used, Thr is the LDF target
*                  accuracy.
*                  If the integrals have not been Cholesky decomposed
*                  (or represented with DF or LDF), this function simply
*                  returns iTolDef.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit None
      Integer iTolDef
#include "cholesky.fh"
#include "choini.fh"

      Real*8   Get_LDFAccuracy
      External Get_LDFAccuracy

      Logical DidCholesky, DidLDF
      Real*8  ThrAbs, d
      Integer ChoIsIni

      Call DecideOnCholesky(DidCholesky)
      If (DidCholesky) Then
         Call DecideOnLocalDF(DidLDF)
         If (DidLDF) Then
            ThrAbs = Abs(Get_LDFAccuracy())
         Else
            Call Get_iScalar('ChoIni',ChoIsIni)
            If (ChoIsIni .ne. ChoIniCheck) Then ! not initialized
               Call Get_dScalar('Cholesky Threshold',ThrCom)
            End If
            ThrAbs = Abs(ThrCom)
         End If
         d = -log(ThrAbs)/log(1.0d1)
         Cho_X_GetTol = nint(d)
      Else
         Cho_X_GetTol = iTolDef
      End If

      End
      Real*8 Function Get_LDFAccuracy()
      Implicit None
#include "localdf.fh"
      Logical  LDF_X_IsSet
      External LDF_X_IsSet
      If (.not.LDF_X_IsSet()) Then
         Call Get_dScalar('LDF Accuracy',Thr_Accuracy)
      End If
      Get_LDFAccuracy=Thr_Accuracy
      End
