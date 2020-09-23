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
      use Sizes_of_Seward, only: Size_Get
      use DKH_Info, only: DKH_Info_Get
      use Real_Info, only: Real_Info_Get
      use RICD_Info, only: RICD_Info_Get
      use Logical_Info, only: Logical_Info_Get
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
      Integer iix(2)
*
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
*
      Call Get_Info_Static_Internal(ixStrt)
      Call Symmetry_Info_Get()
      Call Size_Get()
      Call DKH_Info_Get()
      Call Real_Info_Get()
      Call RICD_Info_Get()
      Call Logical_Info_Get()
      Return
*
*     This is to allow type punning without an explicit interface
*
      Contains

      SubRoutine Get_Info_Static_Internal(ixStrt)
      Use Iso_C_Binding
      Integer, Target :: ixStrt
      Integer, Pointer :: p_ix(:)
*
*     Load the common INFO
*
      Len = iiLoc(ixEnd)-iiLoc(ixStrt)
      Len = (Len+nbyte_i)/nbyte_i
      Call C_F_Pointer(C_Loc(ixStrt),p_ix,[Len])
      Call Get_iArray('SewIInfo',p_ix,Len)
*
      Nullify(p_ix)
*
      Return
      End SubRoutine Get_Info_Static_Internal
*
      End
