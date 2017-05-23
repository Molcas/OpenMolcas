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
      INTEGER FUNCTION JSUNP(INTSYM,L)
      DIMENSION INTSYM(*)

c      INTW=INTSYM((L+9)/10)
c      IPOW=2**(27-3*MOD(L-1,10))
c      JSUNP=1+MOD(INTW/IPOW,8)
      JSUNP=1+MOD(INTSYM((L+9)/10)/(2**(27-3*MOD(L-1,10))),8)
      RETURN
      END
