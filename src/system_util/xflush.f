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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine xflush(Lu)
**********************************************************************
*                                                                    *
*   purpose:                                                         *
*   Dump buffers                                                     *
*                                                                    *
*   calling arguments:                                               *
*   Lu      : integer                                                *
*             Logical unit number                                    *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
*   written by:                                                      *
*   M.P. Fuelscher                                                   *
*   University of Lund, Sweden, 1996                                 *
*                                                                    *
*--------------------------------------------------------------------*
*                                                                    *
*   history: none                                                    *
*                                                                    *
**********************************************************************

#if   defined(_CRAY_C90_)
      If (Lu.EQ.6) then
         Call Flush(101)
      else
         Call Flush(Lu)
      End If
#elif defined(__INTEL_COMPILER) || defined(NAGFOR)
       Return
       If (.False.) Call Unused_integer(Lu)
#elif defined(_SOLARIS_) || defined(_IRIX64_) || defined(_HP_UX_)
      Call Flush(Lu)
#elif defined(_PRIMEPOWER_) || defined(_LINUX_)
      Call Flush(Lu)
#endif
      Return
      End
