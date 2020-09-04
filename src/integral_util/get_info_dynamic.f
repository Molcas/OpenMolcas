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
* Copyright (C) 1992,2020, Roland Lindh                                *
************************************************************************
      SubRoutine Get_Info_Dynamic()
************************************************************************
*                                                                      *
* Object: to read all input information on the file INFO.              *
*                                                                      *
* Called from: Alaska                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              GetMem                                                  *
*              OpnCom                                                  *
*              ClsCom                                                  *
*              RdCom                                                   *
*              SetUp_RW                                                *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             January 1992                                             *
************************************************************************
*                                                                      *
      Use Basis_Info
      Use Center_Info
*                                                                      *
************************************************************************
*                                                                      *
      Call Basis_Info_Get()
      Call Center_Info_Get()
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
