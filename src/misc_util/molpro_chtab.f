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
* Copyright (C) 2014, Naoki Nakatani                                   *
************************************************************************
      SubRoutine Molpro_ChTab(nIrrep,Label,iChMolpro)
************************************************************************
*                                                                      *
* Object: To convert MOLCAS character table to MOLPRO format           *
*         This is a part of integral_util/chtab.f                      *
*                                                                      *
*         Written by N.Nakatani Nov. 2014                              *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Integer iChMolpro(8)
      Integer iOper(nIrrep)
      Character*3 Label
      Logical Rot
************************************************************************
      Call Get_iArray('Symmetry operations',iOper,nIrrep)

      Do i=1,8
        iChMolpro(i)=0
      End Do

****** C1  symmetry ******
      If (nIrrep.eq.1) Then
        Label='c1 '
        iChMolpro(1)=1
      Else If (nIrrep.eq.2) Then
****** Ci  symmetry ******
        If (iOper(2).eq.7) Then
          Label='ci '
****** Cs  symmetry ******
        Else If (iOper(2).eq.1.or.iOper(2).eq.2.or.iOper(2).eq.4) Then
          Label='cs '
****** C2  symmetry ******
        Else
          Label='c2 '
        End If
        iChMolpro(1)=1
        iChMolpro(2)=2
      Else If (nIrrep.eq.4) Then
****** C2h symmetry ******
        If (iOper(2).eq.7.or.iOper(3).eq.7.or.iOper(4).eq.7) Then
          Label='c2h'
* MOLCAS::[ ag, bg, au, bu ] => MOLPRO::[ ag, au, bu, bg ]
          iChMolpro(1)=1
          iChMolpro(2)=4
          iChMolpro(3)=2
          iChMolpro(4)=3
        Else
          Rot = .True.
          Do i=1,nIrrep
            If (iOper(i).eq.1.or.iOper(i).eq.2.or.iOper(i).eq.4)
     &        Rot=.False.
          End Do
****** D2  symmetry ******
          If (Rot) Then
            Label='d2 '
* MOLCAS::[ a, b2, b1, b3 ] => MOLPRO::[ a, b3, b2, b1 ]
            iChMolpro(1)=1
            iChMolpro(2)=3
            iChMolpro(3)=4
            iChMolpro(4)=2
****** C2v symmetry ******
          Else
            Label='c2v'
* MOLCAS::[ a1, b1, a2, b2 ] => MOLPRO::[ a1, b1, b2, a2 ]
            iChMolpro(1)=1
            iChMolpro(2)=2
            iChMolpro(3)=4
            iChMolpro(4)=3
          End If
        End If
****** D2h symmetry ******
      Else If (nIrrep.eq.8) Then
        Label='d2h'
* MOLCAS::[ ag, b2g, b1g, b3g, au, b2u, b1u, b3u ] => MOLPRO::[ ag, b3u, b2u, b1g, b1u, b2g, b3g, au ]
        iChMolpro(1)=1
        iChMolpro(2)=6
        iChMolpro(3)=4
        iChMolpro(4)=7
        iChMolpro(5)=8
        iChMolpro(6)=3
        iChMolpro(7)=5
        iChMolpro(8)=2
      Else
         Call WarningMessage(2,'MOLPRO_ChTab: Illegal value of nIrrep')
         Write (6,*) 'nIrrep=',nIrrep
         Call Abend()
      End If

      Return
      End
