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
       Subroutine xerror(Label,ix,ier,lvl)
       Character*(*) Label
       Integer ix,ier,lvl
*
       Write (6,*) 'Terminate in xerror!'
       Write (6,'(A)') Label
       Write (6,'(A,I5)') 'ix=',ix
       Write (6,'(A,I5)') 'ier=',ier
       Write (6,'(A,I5)') 'lvl=',lvl
       Call QTrace()
       Call Abend()
*
       Return
       End
