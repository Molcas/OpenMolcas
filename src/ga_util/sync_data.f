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
      Subroutine Sync_Data(Pren,Prem,nBtch,mBtch,kBtch)
      use definitions, only: iwp, wp
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: nProcs, Is_Real_Par
#endif
      Implicit None
      real(kind=wp), Intent(inout):: Pren,Prem
      integer(kind=iwp), intent(in):: nBtch,mBtch,kBtch
*
#ifdef _MOLCAS_MPP_

      Real(kind=wp) PrenV(2)
      Integer(kind=iwp) nBtchV(3)
*
      If (.Not. Is_Real_Par()) Return
      If (nProcs.eq.1) Return
      PrenV(1)=Pren
      PrenV(2)=Prem
      Call GADGOP(PrenV,2,'+')
      Pren=PrenV(1)
      Prem=PrenV(2)
*
      nBtchV(1)=nBtch
      nBtchV(2)=mBtch
      nBtchV(3)=kBtch
      Call GAIGOP(nBtchV,3,'+')
      nBtchV(1)=nBtch
      nBtchV(2)=mBtch
      nBtchV(3)=kBtch
#else
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real(Pren)
         Call Unused_real(Prem)
         Call Unused_integer(nBtch)
         Call Unused_integer(mBtch)
         Call Unused_integer(kBtch)
      End If
#endif
*
      End Subroutine Sync_Data
