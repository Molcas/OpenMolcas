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
      subroutine SlapAf(ireturn)
      use Slapaf_Parameters, only: CallLast, isFalcon
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
!     Is this usage permitted? - Yingjin
!#include "dmrginfo_mclr.fh"
      Character*8 ELOOP
************ columbus interface ****************************************
      Integer Columbus
#include "warnings.fh"
*                                                                      *
************************************************************************
*                                                                      *
      Lu=6
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************************************************************************
*                                                                      *
*     Let us go to work !
*
      iStop=-1
      Call RlxCtl(iStop)
      ! in case it is dmrg calculation
*      call read_dmrg_parameter_for_mclr()
*                                                                      *
************************************************************************
*                                                                      *
*
*     Epilogue
*
      If (nPrint(1).ge.6) Then
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put out the proper return code.

      If (iStop.eq.2) Then
************ columbus interface ****************************************
* prevent slapaf to startup alaska
* The last MRCI or MRAQCC energy is written to the RUNFILE by the
* Columbus modules, so no need to recompute it internally

        Call Get_iScalar('Columbus',Columbus)
        if (Columbus.eq.1) then
           iReturn=_RC_ALL_IS_WELL_
        else
*
*        Compute the last energy (module Last_Energy)
*
         if (CallLast) then
          Call Start_Last_Energy()
*
          iReturn=_RC_INVOKED_OTHER_MODULE_
         else
          iReturn=_RC_ALL_IS_WELL_
         endif
        endif
*

      Else If (iStop.eq.3) Then
*
*        Invoke Alaska
*
         Call Start_Alaska()
*
         iReturn=_RC_INVOKED_OTHER_MODULE_
*
      Else If (iStop.eq.6) Then
*
*        Saddle TS optimization
*
         iReturn=_RC_INVOKED_OTHER_MODULE_
*
      Else If (iStop.eq.8) Then
*
*        Predicted energy change too large
*
         iReturn=_RC_NOT_CONVERGED_
*
      Else If (iStop.eq.16) Then
*
*        Optimization didn't converge!
*
         iReturn=_RC_NOT_CONVERGED_
*
      Else If (iStop.eq.1) Then
*
*        Continue looping!
*
         Call GetEnvF('EMIL_InLoop',ELOOP)
         If (ELOOP.eq.' ') ELOOP='0'
         If (ELOOP(1:1).ne.'0') Then
            iReturn=_RC_CONTINUE_LOOP_
         Else
            iReturn=_RC_ALL_IS_WELL_
         End If
         If (isFalcon) Then
            iReturn=_RC_ALL_IS_WELL_
         End if
*
      Else If (iStop.eq.0) Then
*
*        ???
*
         iReturn=_RC_ALL_IS_WELL_
*
      Else
*
*        Bad bad boy!
*
         iReturn=_RC_INTERNAL_ERROR_
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
