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
* Copyright (C) 1992, Roland Lindh                                     *
************************************************************************
      SubRoutine Get_Info_Static(ioffr)
************************************************************************
*                                                                      *
* Object: to read all input information stored in common blocks        *
*         INFO, RINFO, LINFO and CINFO.                                *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              OpnCom                                                  *
*              ClsCom                                                  *
*              RdCom                                                   *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "relae.fh"
#include "RelLight.fh"
      Integer iix(2)
      Real*8 rix(2)
*
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*
*     Prologue
*
*     Call qEnter('GetInfo')
*
*     Load the common INFO
*
      Len = iiLoc(ixEnd)-iiLoc(ixStrt)
      Len = (Len+nbyte_i)/nbyte_i
c      call izero(iwork(ixStrt),Len)
      Call Get_iArray('SewIInfo',ixStrt,Len) ! temporarely deactivated
      Call Get_iArray('nExp',nExp,Mx_Shll)
      Call Get_iArray('nBasis',nBasis,Mx_Shll)
      Call Get_iArray('nBasis_Cntrct',nBasis_Cntrct,Mx_Shll)
      Call Get_iArray('ipCff',ipCff,Mx_Shll)
      Call Get_iArray('ipCff_Cntrct',ipCff_Cntrct,Mx_Shll)
      Call Get_iArray('ipCff_Prim',ipCff_Prim,Mx_Shll)
      Call Get_iArray('ipExp',ipExp,Mx_Shll)
      Call Get_iArray('IndS',IndS,nShlls)
      Call Get_iArray('ip_Occ',ip_Occ,Mx_Shll)
      Call Get_iArray('ipAkl',ipAkl,Mx_Shll)
      Call Get_iArray('ipBk',ipBk,Mx_Shll)
      Call Get_iArray('nOpt',nOpt,nCnttp)
      Call Get_iArray('iCoSet',iCoSet,64*Mx_mdc)
      Call Get_iArray('iSOInf',iSOInf,3*4*MxAO)
      Call ICopy(MxUnq,1,0,IrrCmp,1)
      Call Get_iArray('IrrCmp',IrrCmp,Mx_Unq)
      Call Get_iArray('ipFockOp',ipFockOp,Mx_Shll)
*
*     And some in iAOtSO
*
      Call Allocate_iWork(ip_AS,8*Mx_AO)
      Call Get_iArray('iAOtSO',iWork(ip_AS),8*Mx_AO)
      jp_AS = ip_AS
      Do i = 1, Mx_AO
         Call ICopy(8,iWork(jp_AS),1,iAOtSO(i,0),MxAO)
         jp_AS = jp_AS + 8
      End Do
      Call Free_iWork(ip_AS)
*
      iRELAE=iRELAE_Info
*
*     Load the common LINFO
*
      Len = iiLoc(lxEnd)-iiLoc(lxStrt)
      Len = (Len+nbyte_i)/nbyte_i
      Call Get_iArray('SewLInfo',lxStrt,Len)
      Len = ilLoc(Prjct(Mx_Shll))-ilLoc(Prjct(1))
      Len = (Len+nByte_i)/nByte_i
      Call Get_lArray('Prjct',Prjct,Len)
      Call Get_lArray('Transf',Transf,Len)
      Call Get_lArray('AuxShell',AuxShell,Len)
      Call Get_lArray('FragShell',FragShell,Len)
*
*     Load the common RINFO
*
      Len = idLoc(rxEnd)-idLoc(rxStrt)
      Len = (Len+nByte_r)/nByte_r
      Call Get_dArray('SewRInfo',rxStrt,Len)
      Len = idLoc(RMax_Shll(Mx_Shll))-idLoc(RMax_Shll(1))
      Len = (Len+nByte_r)/nByte_r
      Call Get_dArray('RMax_Shll',RMax_Shll,Len)
      CLight=CLight_Info
*
*     Load the common CINFO
*
      Len = icLoc(cxEnd)-icLoc(cxStrt)
      Len = (Len+nByte_i)/nByte_i
      Call Get_iArray('SewCInfo',cxStrt,Len)
*
*     Epilogue, end
*
*     Call qExit('GetInfo')
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(ioffr)
      End
