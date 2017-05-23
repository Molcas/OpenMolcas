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
* Copyright (C) 2004, Giovanni Ghigo                                   *
*               2005, Francesco Aquilante                              *
************************************************************************
************************************************************************
* This Driver precedes the standard (TraCtl) and Cholesky (Cho_TraCtl) *
* AO/MO two-electrons transformation programs.                         *
* -------------------------------------------------------------------- *
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden                                    *
* Written:  August-September 2004                                      *
* Modified for Cholesky-MP2 May 2005                                   *
************************************************************************
      SUBROUTINE TraCtl_Drv(iType,DoExch2,iPart)
************************************************************
*
*   <DOC>
*     <Name>TraCtl\_Drv</Name>
*     <Syntax>Call TraCtl\_Drv(iType,DoExch2,iPart)</Syntax>
*     <Arguments>
*       \Argument{iType}{Caller program}{integer}{in}
*       \Argument{DoExch2}{Flag for the generation of Exch-2 integrals}
*       {logical}{in}
*       \Argument{iPart}{Partitioning of temp files}{integer}{in}
*     </Arguments>
*     <Purpose>
*      Driver that calls the Cholesky or the Conventional routine for the
*      generation of the Two-electrons integrals file in MO basis.
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>G. Ghigo</Author>
*     <Modified_by>F. Aquilante</Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      All programs that need the generation of the Two-electrons
*      integrals file in MO basis must tell to the Cholesky routine
*      who they are through iType:
*        \begin{itemize}
*          \item 1 MBPT2
*          \item 2 CASPT2
*          \item 3 MCLR
*          \item 4 CC
*        \end{itemize}
*      Programs must also tell to the Cholesky routine whether they need
*      exchange-2 integrals throught logical variable DoExch2. Both
*      values are not used by the Conventional routine.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer iType
#include "rasdim.fh"
#include "caspt2.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
      Logical DoCholesky, DoExch2
#include "chocaspt2.fh"
      Character*10 SECNAM

      SECNAM='TraCtl_Drv'

      Call DecideOnCholesky(DoCholesky)

      If (DoCholesky) then

        If (iType.EQ.1) then

          Call ChoMP2_TraCtl(LUINTM,Work(LCMO),NCMO)

        elseIf (iALGO.eq.0) Then

          Call Cho_TraCtl(iType,LUINTM,Work(LCMO),NCMO,DoExch2)

        elseIf (iALGO.eq.1) Then

* caspt2 with cholesky does no longer use call to tractl_drv/cho_caspt2_drv
*          Call Cho_caspt2_drv(Work(LCMO))

        else

          Call Cho_x_Quit(SecNam,' !!! Unknown algorithm !!! ',' ')

        EndIf

      Else

        Call TRACTL(iPart)

      EndIf


      RETURN
      END
