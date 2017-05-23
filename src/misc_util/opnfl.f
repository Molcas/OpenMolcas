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
* Copyright (C) 2000, Roland Lindh                                     *
************************************************************************
      SubRoutine OpnFl(Name,Lu,Exist)
************************************************************************
*                                                                      *
* Object:                                                              *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    : f_Inquire                                               *
*              isFreeUnit                                              *
*              Open                                                    *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemical Physics,                 *
*             University of Lund, SWEDEN                               *
*             February 2000                                            *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Character*(*) name
      Logical Exist
*
*---- Find an unused unit number
*
      lu_=lu
      lu=isFreeUnit(lu_)
      Exist = .False.
*
*---- Check that file exists
*
c      Open(unit=Lu,file=name,status='UNKNOWN',form='FORMATTED')
      Call F_Inquire(Name,Exist)
      Call Molcas_Open(Lu,Name)
*
      Return
      End
