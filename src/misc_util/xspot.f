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
      Subroutine xSpot(Label)
      Character*(*) Label
*
      Write (6,*)
      Write (6,'(A)') Label
#if defined(_GA_)
      Call GetMem(Label,'Check','Real',iDum,iDum)
#else
      Call GetMem('Check','Check','Real',iDum,iDum)
#endif
*
      Return
      End
