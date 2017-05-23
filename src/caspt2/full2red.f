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
      SUBROUTINE full2red(XLT,Xab)
      Implicit Real*8 (a-h,o-z)
      Integer  ISLT(8),cho_isao
      External cho_isao
      Dimension XLT(*)
      Dimension Xab(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

* Select table column for use with caspt2:
      iLoc=3
* jSym=1 always: Used for density matrices.
      jSym = 1
c Offsets to symmetry block in the LT matrix
      IS=0
      DO ISYM=1,NSYM
       ISLT(ISYM)=IS
       NB=NBAS(ISYM)
       IS=IS+(NB*(NB+1))/2
      END DO


      Do jRab=1,nnBstR(jSym,iLoc)
         kRab = iiBstr(jSym,iLoc) + jRab
         iRab = iWork(ip_IndRed-1+nnBstrT(1)*(iLoc-1)+kRab)
         idx = ip_iRS2F+2*(iRab-1)
         iag   = iWork(idx)
         ibg   = iWork(idx+1)
         iSyma = cho_isao(iag)
         ias   = iag - ibas(iSyma)
         ibs   = ibg - ibas(iSyma)
         if(ias.ge.ibs) then
          iab=(ias*(ias-1))/2+ibs
         else
          iab=(ibs*(ibs-1))/2+ias
         end if
         kfrom = isLT(iSyma) + iab
         Xab(jRab) = Xab(jRab)+XLT(kfrom)
      End Do

      Return
      End
      SUBROUTINE red2full(XLT,Xab)
      Implicit Real*8 (a-h,o-z)
      Integer  ISLT(8),cho_isao
      External cho_isao
      Dimension XLT(*)
      Dimension Xab(*)
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

* Select table column for use with caspt2:
      iLoc=3
* jSym=1 always: Used for density matrices.
      jSym = 1
c Offsets to symmetry block in the LT matrix
      IS=0
      DO ISYM=1,NSYM
       ISLT(ISYM)=IS
       NB=NBAS(ISYM)
       IS=IS+(NB*(NB+1))/2
      END DO

      Do jRab=1,nnBstR(jSym,iLoc)
         kRab = iiBstr(jSym,iLoc) + jRab
         iRab = iWork(ip_IndRed-1+nnBstrT(1)*(iLoc-1)+kRab)
         idx = ip_iRS2F+2*(iRab-1)
         iag   = iWork(idx)
         ibg   = iWork(idx+1)
         iSyma = cho_isao(iag)
         ias   = iag - ibas(iSyma)
         ibs   = ibg - ibas(iSyma)
         if(ias.ge.ibs) then
          iab=(ias*(ias-1))/2+ibs
         else
          iab=(ibs*(ibs-1))/2+ias
         end if
         kto = isLT(iSyma) + iab
         XLT(kto) = XLT(kto)+Xab(jRab)
      End Do

      Return
      End
