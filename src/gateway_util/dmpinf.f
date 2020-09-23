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
*               Markus P. Fuelscher                                    *
************************************************************************
      SubRoutine DmpInf()
************************************************************************
*                                                                      *
* Object: to dump all input information on the file INFO.              *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    :                                                         *
*              Put_dArray                                              *
*              Get_dArray                                              *
*              Put_iArray                                              *
*              Get_iArray                                              *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
*                                                                      *
* modified by M.P. Fuelscher                                           *
* - changed to used communication file                                 *
************************************************************************
      use External_Centers
      use Basis_Info, only: Basis_Info_Dmp
      use Center_Info, only: Center_Info_Dmp
      use Symmetry_Info, only: Symmetry_Info_Dmp
      use SOAO_Info, only: SOAO_Info_Dmp
      use Sizes_of_Seward, only: Size_Dmp
      use DKH_Info, only: DKH_Info_Dmp
      use Real_Info, only: Real_Info_Dmp
      use RICD_Info, only: RICD_Info_Dmp
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
#include "rctfld.fh"
#include "nq_info.fh"
      Integer  iix(2)
      Real*8   rix(2)
      nbyte_i = iiloc(iix(2)) - iiloc(iix(1))
      nbyte_r = idloc(rix(2)) - idloc(rix(1))
*
      Call DmpInf_Internal(ixStrt,
     & cRFStrt,iRFStrt,lRFStrt,rRFStrt,cQStrt,iQStrt,rQStrt)
*
*     This is to allow type punning without an explicit interface
      Contains
      SubRoutine DmpInf_Internal(ixStrt,
     & cRFStrt,iRFStrt,lRFStrt,rRFStrt,cQStrt,iQStrt,rQStrt)
      Use Iso_C_Binding
      Integer, Target :: ixStrt,cRFStrt,iRFStrt,lRFStrt,
     &                   cQStrt,iQStrt
      Real*8, Target :: rRFStrt,rQStrt
      Integer, Pointer :: p_ix(:),p_cRF(:),p_iRF(:),
     &                    p_lRF(:),p_cQ(:),p_iQ(:)
      Real*8, Pointer :: p_rRF(:),p_rQ(:)
*
*     Save the common INFO
*
      Len = iiLoc(ixEnd)-iiLoc(ixStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(ixStrt),p_ix,[Len])
      Call Put_iArray('SewIInfo',p_ix,Len)
*                                                                      *
************************************************************************
*                                                                      *
      Call SOAO_Info_Dmp()
      Call Basis_Info_Dmp()
      Call Center_Info_Dmp()
      Call Symmetry_Info_Dmp()
      Call Size_Dmp()
      Call DKH_Info_Dmp()
      Call Real_Info_Dmp()
      Call RICD_Info_Dmp()
*                                                                      *
************************************************************************
*                                                                      *
*     Reaction field parameters
*
      Len = iiLoc(lRFEnd)-iiLoc(lRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(lRFStrt),p_lRF,[Len])
      Call Put_iArray('RFlInfo',p_lRF,Len)
*
      Len = idLoc(rRFEnd)-idLoc(rRFStrt)
      Len = (Len+nByte_r)/nByte_r
      Call C_F_Pointer(C_Loc(rRFStrt),p_rRF,[Len])
      Call Put_dArray('RFrInfo',p_rRF,Len)
*
      Len = iiLoc(iRFEnd)-iiLoc(iRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(iRFStrt),p_iRF,[Len])
      Call Put_iArray('RFiInfo',p_iRF,Len)
*
      Len = iiLoc(cRFEnd)-iiLoc(cRFStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(cRFStrt),p_cRF,[Len])
      Call Put_iArray('RFcInfo',p_cRF,Len)
*
      Nullify(p_lRF,p_rRF,p_iRF,p_cRF)
*                                                                      *
************************************************************************
*                                                                      *
*     Numerical integration information and parameters
*
      Len = idLoc(rQEnd)-idLoc(rQStrt)
      Len = (Len+nByte_r)/nByte_r
      Call C_F_Pointer(C_Loc(rQStrt),p_rQ,[Len])
      Call Put_dArray('Quad_r',p_rQ,Len)
*
      Len = iiLoc(iQEnd)-iiLoc(iQStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(iQStrt),p_iQ,[Len])
      Call Put_iArray('Quad_i',p_iQ,Len)
*
      Len = iiLoc(cQEnd)-iiLoc(cQStrt)
      Len = (Len+nByte_i)/nByte_i
      Call C_F_Pointer(C_Loc(cQStrt),p_cQ,[Len])
      Call Put_iArray('Quad_c',p_cQ,Len)
*
      Nullify(p_rQ,p_iQ,p_cQ)
*                                                                      *
************************************************************************
*                                                                      *
      Call External_Centers_Dmp()
*                                                                      *
************************************************************************
*                                                                      *
      Call DMP_EFP()
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End SubRoutine DmpInf_Internal
*
      End SubRoutine DmpInf
