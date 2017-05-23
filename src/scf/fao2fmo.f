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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
      Subroutine Fao2Fmo
************************************************************************
*                                                                      *
* This routine transforms a forck matrix in AO representation to An MO *
* representation.                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Called from:                                                         *
* Call to:                                                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden.                                    *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Logical Debug
*----------------------------------------------------------------------*
* Perform setup                                                        *
*----------------------------------------------------------------------*
#ifdef _DEBUG_
      Debug=.true.
#else
      Debug=.false.
#endif
      Call qEnter('Fao2Fmo')
      Write(6,'(a)') '+++ Entering Fao2Fmo'
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Write(6,'(a)') '+++ Exiting Fao2Fmo'
      Call qExit ('Fao2Fmo')
      Return
      End
