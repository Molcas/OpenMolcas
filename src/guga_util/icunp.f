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
      INTEGER FUNCTION ICUNP(ICSPCK,L)
      DIMENSION ICSPCK(*)

      INTW=ICSPCK((L+14)/15)
      IPOW=2**(28-2*MOD(L-1,15))
      ICUNP=MOD(INTW/IPOW,4)
      RETURN
      END
