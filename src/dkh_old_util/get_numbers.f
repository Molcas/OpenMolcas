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
         Integer Function Get_sNumber(dkhunit3)
         Integer dkhunit3, snumber
         character*(43) scrtxt
#ifdef _MOLCAS_
         Call molcas_open(dkhunit3,'dkhops.13')
#endif
         Rewind(dkhunit3)
         Read (dkhunit3,'(A43)') scrtxt(1:43)
         Read (dkhunit3,'(A43)') scrtxt(1:43)
1050     Read (dkhunit3,'(A3)') scrtxt(1:3)
         If (scrtxt(1:3).ne.'+++') goto 1050
         Read (dkhunit3,'(I3)') snumber
#ifdef _MOLCAS_
         Close (dkhunit3)
#endif
*
         Get_sNumber = snumber
         Return
         End
         Integer Function Get_tNumber(dkhunit4)
         Integer dkhunit4, tnumber
         character*(3) scrtxt
*
#ifdef _MOLCAS_
         Call molcas_open(dkhunit4,'dkhops.14')
#endif
         Rewind(dkhunit4)
1051     Read (dkhunit4,'(A3)') scrtxt(1:3)
         If (scrtxt(1:3).ne.'+++') goto 1051
         Read (dkhunit4,'(I3)') tnumber
#ifdef _MOLCAS_
         Close (dkhunit4)
#endif
*
         Get_tNumber = tnumber
         Return
         End
         Integer Function Get_uNumber(dkhunit5)
         Integer dkhunit5, unumber
         character*(3) scrtxt
*
#ifdef _MOLCAS_
         Call molcas_open(dkhunit5,'dkhops.15')
#endif
         Rewind(dkhunit5)
1052     Read (dkhunit5,'(A3)') scrtxt(1:3)
         If (scrtxt(1:3).ne.'+++') goto 1052
         Read (dkhunit5,'(I3)') unumber
#ifdef _MOLCAS_
         Close (dkhunit5)
#endif
*
         Get_uNumber = unumber
         Return
         End
