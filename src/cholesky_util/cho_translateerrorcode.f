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
      SubRoutine Cho_TranslateErrorCode(iErr,MolcasCode)
      Implicit None
      Integer iErr, MolcasCode
#include "warnings.fh"

      If (iErr .eq. 3) Then ! special code used in parallel
         MolcasCode = _RC_NOT_AVAILABLE_
      Else If (iErr .eq. 100) Then
         MolcasCode = _RC_CHO_DUM_
      Else If (iErr .eq. 101) Then
         MolcasCode = _RC_CHO_MEM_
      Else If (iErr .eq. 102) Then
         MolcasCode = _RC_CHO_INI_
      Else If (iErr .eq. 103) Then
         MolcasCode = _RC_CHO_LOG_
      Else If (iErr .eq. 104) Then
         MolcasCode = _RC_CHO_RUN_
      Else If (iErr .eq. 105) Then
         MolcasCode = _RC_CHO_INP_
      Else
         MolcasCode = _RC_GENERAL_ERROR_
      End If

      End
