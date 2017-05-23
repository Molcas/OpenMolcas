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
      subroutine getLIMIT(l1,l2,l3,l4,Lanf,Lend)
      implicit Integer (a-z)
cbs   get the minimum and maximum L-values
cbs   of the the coulomb-potential to interact
cbs   with l1-l4
      lower1=iabs(l1-l3)
      lower2=iabs(l2-l4)
      lupper1=l1+l3
      lupper2=l2+l4
      Lanf=max(lower1,lower2)
      Lend=min(lupper1,lupper2)
cbs     check for parity
      lsum=Lanf+l1+l3
      if (mod(lsum,2).eq.1) Lanf=Lanf+1
      lsum=Lend+l1+l3
      if (mod(lsum,2).eq.1) Lend=Lend-1
cbs   check the other parity
      lsum=Lanf+l2+l4
      if (mod(lsum,2).eq.1) then
      write(6,*) ' error in getLIMIT: '
      write(6,*) ' parity inconsistency for '
      write(6,*) 'l1,l2,l3,l4= ',l1,l2,l3,l4
      Call Abend()
      endif
      return
      end
