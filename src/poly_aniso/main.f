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
      program main
#ifdef _FPE_TRAP_
      Use, Intrinsic :: IEEE_Exceptions
#endif
      Integer iReturn
      Character(Len=20) Module_Name
      Parameter (Module_Name = 'poly_aniso')
#ifdef _FPE_TRAP_
      Call IEEE_Set_Halting_Mode(IEEE_Usual,.True._4)
#endif

      Call Start(Module_Name)
      Call poly_aniso_open(iReturn)
      Call QStat(' ')
      Call Finish(iReturn)
      End
