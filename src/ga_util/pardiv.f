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
      INTEGER FUNCTION iParDiv(nMax,nMin)
*
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: nProcs, Is_Real_Par
      IF (Is_Real_Par()) THEN
        iParDiv = nMax/nProcs+1 + nMin
      ELSE
        iParDiv = nMax + nMin
      ENDIF
#else
      iParDiv = nMax + nMin
#endif
*
      Return
      End
