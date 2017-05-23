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
      Subroutine LDF_X_Final(Term_Sew,irc)
************************************************************************
*
*   <DOC>
*     <Name>LDF\_X\_Final</Name>
*     <Syntax>Call LDF\_X\_Final(Term\_Sew,irc)</Syntax>
*     <Arguments>
*       \Argument{Term\_Sew}{Do termination of Seward environment}
*       {Logical}{in}
*       \Argument{irc}{Return code}{Integer}{out}
*     </Arguments>
*     <Purpose>Finalize LDF enviroment</Purpose>
*     <Dependencies>Call LDF\_X\_Init</Dependencies>
*     <Author>Thomas Bondo Pedersen</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description> Finalization of the LDF environment, including
*     proper deallocations of memory. If (Term\_Sew), terminate the
*     Seward environment as well. On exit, irc=0 signals successful
*     finalization.
*     </Description>
*    </DOC>
*
************************************************************************
      Implicit None
      Logical Term_Sew
      Integer irc
#include "localdf.fh"

      Character*11 SecNam
      Parameter (SecNam='LDF_X_Final')

      Integer LDF_Status

      ! Init return code
      irc=0

      ! Read init flag on runfile
      Call Get_iScalar('LDF Status',LDF_Status)

      ! Finalize if set
      If (LDF_Status.eq.LDF_Set) Then
         Call LDF_UnsetMltplCenters(max(LDF_Constraint,0))
         Call LDF_CIO_Final()
         Call LDF_Final(Term_Sew,irc)
         If (irc.ne.0) Then
            Write(6,'(A,A,I8)')
     &      SecNam,': LDF_Final returned code',irc
            irc=1
         End If
         If (Term_Sew) Then
            Call ClsSew()
         End If
         LDF_Status=LDF_Set-1
         Call LDF_SetStatusOnRunFile(LDF_Status)
      End If

      End
