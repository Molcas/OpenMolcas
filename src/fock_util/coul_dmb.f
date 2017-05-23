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
      Subroutine Coul_DMB(GetFM,nDM,Rep_EN,FM,DMA,DMB,lFDM)

      Implicit Real*8 (a-h,o-z)
      Logical GetFM
      Integer nDM, lFDM
      Real*8  FM(lFDM), DMA(lFDM), DMB(lFDM)
#include "WrkSpc.fh"
#include "cholesky.fh"
#include "choorb.fh"
      Character*16 NamRfil
*
      If (nDM.gt.2 .or. nDM.lt.1) Then
         write(6,*) ' In Coul_DMB: wrong value of nDM= ',nDM
         Call SysAbendMsg('Coul_DMB ',' nDM must be 1 or 2 ',' ')
      EndIf
*
      If (.not.GetFM) Go To 99
*
      Call Get_NameRun(NamRfil) ! save the old RUNFILE name
      Call NameRun('AUXRFIL')   ! switch RUNFILE name
*
      Call GetMem('DM','Allo','Real',ip_DM,lFDM)
      Call get_dArray('D1ao',Work(ip_DM),lFDM)
*
      Call FZero(FM,lFDM)
      Call CHO_FOCK_DFT_RED(irc,Work(ip_DM),FM)
      If (irc.ne.0) Then
         Call SysAbendMsg('Coul_DMB ',' non-zero rc ',' ')
      EndIf
*
      Call GetMem('DM','Free','Real',ip_DM,lFDM)
*
      Call NameRun(NamRfil)   ! switch back RUNFILE name
*
99    Continue
      Rep_EN = ddot_(lFDM,DMA,1,FM,1)
      If (nDM.eq.2) Then
         Rep_EN = Rep_EN + ddot_(lFDM,DMB,1,FM,1)
      EndIf
*
      Return
      End
