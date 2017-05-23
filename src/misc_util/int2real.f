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
* Copyright (C) 1999, Thorstein Thorsteinsson                          *
************************************************************************
      Subroutine Int2Real(IntArg,RealArg)
************************************************************************
*                                                                      *
*    Purpose: Convert integer argument to Real and vice versa          *
*                                                                      *
*    Calling parameters:                                               *
*    IntArg : contains integer array to be converged                   *
*    RealArg: contains on return the corresponding real                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     T. Thorsteinsson                                                 *
*     University of Lund, Sweden, 1999                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit None
*
#ifndef _I8_
      Integer IntArg(2),iTemp(2)
#else
      Integer IntArg(1),iTemp(1)
#endif
      Real*8 RealArg,Temp
      Equivalence (iTemp,Temp)

      iTemp(1)=IntArg(1)
#ifndef _I8_
      iTemp(2)=IntArg(2)
#endif
      RealArg=Temp
      Return
      end

      Subroutine Real2Int(RealArg,IntArg)
*
#ifndef _I8_
      Integer IntArg(2),iTemp(2)
#else
      Integer IntArg(1),iTemp(1)
#endif
      Real*8 RealArg,Temp
      Equivalence (iTemp,Temp)

      Temp=RealArg
      IntArg(1)=iTemp(1)
#ifndef _I8_
      IntArg(2)=iTemp(2)
#endif
      Return
      End
