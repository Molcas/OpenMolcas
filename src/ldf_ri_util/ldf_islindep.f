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
      Logical Function LDF_isLinDep(i,iShell,iAtom,iAtomPair)
      Implicit None
      Integer i, iShell, iAtom, iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Integer iLD, nLD, ip_LD, ip

      Integer j, k
      Integer AP_1CLinDep
      AP_1CLinDep(j,k)=iWork(ip_AP_1ClinDep-1+2*(k-1)+j)

      LDF_isLinDep=.False.
      nLD=AP_1CLinDep(1,iAtomPair)
      If (nLD.lt.1) Return

      ip_LD=AP_1CLinDep(2,iAtomPair)
      iLD=0
      Do While (.not.LDF_isLinDep .and. iLD.lt.nLD)
         ip=ip_LD+3*iLD
         LDF_isLinDep=iWork(ip).eq.iAtom  .and.
     &                iWork(ip+1).eq.iShell .and.
     &                iWork(ip+2).eq.i
         iLD=iLD+1
      End Do

      End
