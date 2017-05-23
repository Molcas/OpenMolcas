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
* Copyright (C) 2017, Ignacio Fdez. Galvan                             *
************************************************************************
      Logical Function Is_First_Iter()
      Implicit None
      Character(Len=80) :: EnvStr
      Integer :: Iter, Iter_S
      Logical :: Found
      Integer, Dimension(7) :: Tmp

      ! If this is the first iteration in a "saddle" branch
      Call qpg_iScalar('Saddle Iter',Found)
      If (Found) Then
        Call Get_iScalar('Saddle Iter',Iter_S)
        If (Iter_S .eq. 0) Then
          Is_First_Iter = .True.
          Return
        Else
          Is_First_Iter = .False.
          Return
        End If
      End If

      ! If Slapaf information has been stripped out (e.g. IRC restart)
      Call qpg_iArray('Slapaf Info 1',Found,Tmp(1))
      If (Found) Then
        Call Get_iArray('Slapaf Info 1',Tmp,7)
        If (Tmp(1) .eq. -99) Then
          Is_First_Iter = .True.
          Return
        End If
      End If

      ! If MOLCAS_ITER <= 1
      Call Getenvf("MOLCAS_ITER",EnvStr)
      Read (EnvStr,*) Iter
      If (Iter .le. 1) Then
        Is_First_Iter = .True.
        Return
      End If

      Is_First_Iter = .False.

      End Function Is_First_Iter
