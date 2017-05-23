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
* Copyright (C) 2008, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine initialize the peek/poke utility.                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: May 2008                                                    *
*                                                                      *
************************************************************************
*
* <DOC>
*   <Name>Init\_ppu(Force)</Name>
*   <Syntax>Call Init\_ppu</Syntax>
*   <Arguments>
*     \Argument{Force}{Force initialization}{Logical}{in}
*   </Arguments>
*   <Purpose>To initialize the peek/poke utility.</Purpose>
*   <Dependencies></Dependencies>
*   <Author>Per-Olof Widmark</Author>
*   <Modified_by></Modified_by>
*   <Side_Effects></Side_Effects>
*   <Description>
*     This routine is used to initialize the peek/poke utility.
*   </Description>
* </DOC>
*
************************************************************************
      Subroutine Init_ppu(Force)
      Implicit None
#include "pp_ds_info.fh"
*----------------------------------------------------------------------*
* Dummy arguments.                                                     *
*----------------------------------------------------------------------*
      Logical Force
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Logical FirstTime
      Save    FirstTime
      Data    FirstTime/.true./
*----------------------------------------------------------------------*
* Initialize the peek/poke utility                                     *
*----------------------------------------------------------------------*
      If(FirstTime.or.Force) Then
         ds_no=0
      End If
      FirstTime=.false.
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
