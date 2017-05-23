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
* Copyright (C) 1986, Per E. M. Siegbahn                               *
************************************************************************
      SUBROUTINE ICPCK(ICSPCK,L,ICASE)
      DIMENSION ICSPCK(*)

c      INTW=ICSPCK((L+14)/15)
c      IPOW=2**(28-2*MOD(L-1,15))
c      INTW=INTW+ICASE*IPOW
c      ICSPCK((L+14)/15)=INTW
       MY=(L+14)/15
      IPOW=2**(28-2*MOD(L-1,15))
      ICSPCK(MY)=ICSPCK(MY)+ICASE*IPOW
      RETURN
      END
