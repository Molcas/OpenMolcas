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
      SubRoutine Get_Info_Static()
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
      use Symmetry_Info, only: Symmetry_Info_Get
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "relae.fh"
#include "RelLight.fh"
      Integer iix(2)
      Real*8 rix(2)
*
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*
      Call Get_Info_Static_Internal(cxStrt,ixStrt,lxStrt,rxStrt)
      Return
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine Get_Info_Static_Internal(cxStrt,ixStrt,lxStrt,rxStrt)
      Use Iso_C_Binding
      Integer, Target :: cxStrt,ixStrt,lxStrt
      Real*8, Target :: rxStrt
      Integer, Pointer :: p_cx(:),p_ix(:),p_lx(:)
      Real*8, Pointer :: p_rx(:)
*
*     Load the common INFO
*
      Len = iiLoc(ixEnd)-iiLoc(ixStrt)
      Len = (Len+nbyte_i)/nbyte_i
      Call C_F_Pointer(C_Loc(ixStrt),p_ix,[Len])
      Call Get_iArray('SewIInfo',p_ix,Len) ! temporarely deactivated
*
      iRELAE=iRELAE_Info
*
*     Load the common LINFO
*
      Len = iiLoc(lxEnd)-iiLoc(lxStrt)
      Len = (Len+nbyte_i)/nbyte_i
      Call C_F_Pointer(C_Loc(lxStrt),p_lx,[Len])
      Call Get_iArray('SewLInfo',p_lx,Len)
*
*     Load the common RINFO
*
      Len = idLoc(rxEnd)-idLoc(rxStrt)
      Len = (Len+nByte_r)/nByte_r
      Call C_F_Pointer(C_Loc(rxStrt),p_rx,[Len])
      Call Get_dArray('SewRInfo',p_rx,Len)
      CLightAU=CLight_Info
*
*     Load the common CINFO
*
      Len = icLoc(cxEnd)-icLoc(cxStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(cxStrt),p_cx,[Len])
      Call Get_iArray('SewCInfo',p_cx,Len)
*
      Nullify(p_ix,p_lx,p_rx,p_cx)
*
      Call Symmetry_Info_Get()
*
      Return
      End SubRoutine Get_Info_Static_Internal
*
      End
