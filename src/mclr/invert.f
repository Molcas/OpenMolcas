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
      Subroutine Invert(rmat,n)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Real*8 rmat(*),rdet(2),rcond
      Integer n,iptmp
      Call Getmem('Tmp','ALLO','REAL',iptmp,200*n)
      Call DGEICD(rmat,n,n,0,rcond,rdet,work(iptmp),200*n)
      Call Getmem('Tmp','FREE','REAL',iptmp,200*n)
      Return
      end
