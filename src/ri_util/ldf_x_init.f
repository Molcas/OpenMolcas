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
*  LDF_X_Init
*
*> @brief
*>   Initialize LDF (Local Density Fitting) environments
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Initialization of LDF environment. If \p Ini_Sew,
*> Seward is initialized as part of the LDF setup (if not, Seward
*> must be initialized before calling this routine). \p nDiff is only
*> used if Seward is to be initialized.
*>
*> \p BufFrac specifies the fraction of available memory (after all
*> other LDF data has been allocated) that can at most be used as
*> buffer for fitting coefficients. This buffer serves to minimize
*> the amount of I/O done in iterative procedures such as SCF.
*> A non-zero return code means that the initialization failed:
*>
*> Return codes:
*>
*> - \p irc = ``-2``: not LDF (according to Runfile)
*> - \p irc = ``-1``: illegal \p nDiff input parameter
*> - \p irc =  ``0``: initialization success
*> - \p irc =  ``1``: error in setup routine
*>
*> @note
*> Seward must have generated the LDF fitting coefficients.
*>
*> @param[in]  Ini_Sew Do initialization of Seward environment
*> @param[in]  nDiff   Initialize Seward to compute derivatives of order \p nDiff
*> @param[in]  BufFrac Memory fraction used for fitting coefficient buffer}
*> @param[out] irc     Return code
*>
*> @see ::LDF_X_Final
************************************************************************
      Subroutine LDF_X_Init(Ini_Sew,nDiff,BufFrac,irc)
      Implicit None
      Logical Ini_Sew
      Integer nDiff
      Real*8  BufFrac
      Integer irc
#include "localdf.fh"
#include "localdf_print.fh"

      Logical  LDF_X_IsSet, LDF_With2CF
      External LDF_X_IsSet, LDF_With2CF

      Integer  iPrintLevel
      External iPrintLevel

      Character*10 SecNam
      Parameter (SecNam='LDF_X_Init')

      Logical Verbose
      Parameter (Verbose=.False.)

      Logical isLocalDF
      Logical FirstCall
      Logical DoPairs
      Logical DoERI

      Integer LDF_Unset
      Integer iCLDF

      Save FirstCall
      Data FirstCall /.True./

      ! Init return code
      irc=0

      ! Define LDF_Unset
      LDF_Unset=LDF_Set-1

      ! Check that this is Local DF
      Call DecideOnLocalDF(isLocalDF)
      If (.not.isLocalDF) Then
         irc=-2
         Call WarningMessage(1,
     &                     'WARNING: '//SecNam//
     &                     ' called but this is not a Local DF run ?!?')
         Call LDF_SetStatusOnRunFile(LDF_Unset)
         Return
      End If

      ! Check if already initialized
      If (FirstCall) Then ! it cannot be already done
         FirstCall=.False.
      Else ! might be already done
         If (LDF_X_IsSet()) Then ! it is already done
            irc=0
            Return
         End If
      End If

      ! Define all variables in common blocks (dummy variables)
      Call LDF_SetInc()

      ! Set run mode to "external"
      LDF_Run_Mode=LDF_RUN_EXTERNAL

      ! Set LuPri in Cholesky (enables use of printing routines in
      ! Cholesky utility)
      Call LDF_Cho_SetLuPri(6)

      ! Initialize Seward if requested
      If (Ini_Sew) Then
         If (nDiff.lt.0) Then
            irc=-1
            Call WarningMessage(1,'WARNING: '//SecNam//': nDiff<0')
            Call LDF_SetStatusOnRunFile(LDF_Unset)
            Return
         Else
            DoERI=.True.
            Call IniSew(DoERI,nDiff)
         End If
      End If

      ! Get target accuracy from Runfile and set thresholds
      Call LDF_GetThrs()

      ! Get constraint type from Runfile
      Call Get_iScalar('LDF Constraint',iCLDF)
      Call LDF_AddConstraint(iCLDF) ! this also checks the value

      ! Initialize LDF, except pair information
      DoPairs=.False.
      Call LDF_Init(DoPairs,Verbose,irc)
      If (irc.ne.0) Then
         irc=1
         Call WarningMessage(1,
     &                       'WARNING: '//SecNam//': Error in LDF_Init')
         Call LDF_SetStatusOnRunFile(LDF_Unset)
         Return
      End If

      ! Read pair info from disk
      Call LDF_ReadAtomPairInfo(irc)
      If (irc.ne.0) Then
         irc=1
         Call WarningMessage(1,
     &           'WARNING: '//SecNam//': Error in LDF_ReadAtomPairInfo')
         Call LDF_SetStatusOnRunFile(LDF_Unset)
         Return
      End If

      ! Set flag LDF2 (.True. if two-center functions present)
      LDF2=LDF_With2CF()

      ! Initialize coefficient I/O:
      ! - open coefficient file
      ! - allocate and initialize coefficient buffer
      Call LDF_CIO_Init(BufFrac,irc)
      If (irc.ne.0) Then
         If (irc.eq.-1) Then
            Call WarningMessage(1,
     &             'WARNING: '//SecNam//': Coefficient file not found!')
            irc=0
         Else
            Call WarningMessage(1,
     &                   'WARNING: '//SecNam//': Error in LDF_CIO_Init')
            Call LDF_SetStatusOnRunFile(LDF_Unset)
         End If
      End If

      ! Set up info needed for constrained LDF
      ! TODO/FIXME: only really needed for debug options
      Call LDF_SetMltplCenters(max(LDF_Constraint,0))

      ! Get print level
      iPrint=iPrintLevel(-1)

      ! Set status on runfile
      Call LDF_SetStatusOnRunFile(LDF_Set)

      End
