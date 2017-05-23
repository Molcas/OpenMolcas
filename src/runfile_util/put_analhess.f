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
* Copyright (C) Roland Lindh                                           *
************************************************************************
      Subroutine Put_AnalHess(AnalHess,nAnalHess)
************************************************************
*
*   <DOC>
*     <Name>Put\_AnalHess</Name>
*     <Syntax>Call Put\_AnalHess(AnalHess,nAnalHess)</Syntax>
*     <Arguments>
*       \Argument{AnalHess}{Array with the symmetry blocked nuclear Hessian in cartesian coordinates.}{Real*8 (nAnalHess)}{in}
*       \Argument{nAnalHess}{Size of the array AnalHess}{Integer}{in}
*     </Arguments>
*     <Purpose>To write the symmetry blocked nuclear Hessian in cartesian coordinates on the run file.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>R. Lindh</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>The utility will write the symmetry blocked nuclear Hessian in cartesian coordinates on the run file.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 AnalHess(nAnalHess)
      Character Label*80
      Logical Found
      Integer iTmp, iSI1(7)

      Call Put_dArray('Analytic Hessian',AnalHess,nAnalHess)

*---- Add the iteration number corresponding to this Hessian
      iSI1(2)=0
      Call Qpg_iArray('Slapaf Info 1',Found,iTmp)
      If (Found) Call Get_iArray('Slapaf Info 1',iSI1,7)
      Call Getenvf("MOLCAS_ITER",Label)
      Read (Label,'(I80)') iTmp
      Call Getenvf("EMIL_InLoop",Label)
      Read (Label,*,IOStat=irderr) inLoop
      If (irderr.ne.0) inLoop=0
      If (inLoop.lt.1) iTmp=0
      If (iTmp.eq.0) Then
        Call Put_iScalar('HessIter',iTmp)
      Else
        Call Put_iScalar('HessIter',iSI1(2)+1)
      End If

      Return
      End
