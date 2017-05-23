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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This procedure classify the atomic orbitals of atoms into 6 classes  *
* and return the count. Which shells that are returned are specified   *
* by the the option switch opt:                                        *
*                                                                      *
*  1 DeepCore                                                          *
*  2 Core                                                              *
*  4 SoftCore                                                          *
*  8 DeepValence                                                       *
* 16 Valence                                                           *
* 32 ExtraValence                                                      *
*                                                                      *
* The numbers are added up to get more than one shell reported.        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University                                             *
*          Sweden                                                      *
* Written: May 2003                                                    *
*                                                                      *
************************************************************************
      Subroutine OrbType(Z,List,opt)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Integer Z,List,opt
      Dimension List(4)
*----------------------------------------------------------------------*
* Type declare local variables                                         *
*----------------------------------------------------------------------*
      Integer DeepCore,Core,SoftCore
      Integer DeepValence,Valence,ExtraValence
      Integer i
*----------------------------------------------------------------------*
* Dimension local variables                                            *
*----------------------------------------------------------------------*
      Dimension DeepCore(4)
      Dimension Core(4)
      Dimension SoftCore(4)
      Dimension DeepValence(4)
      Dimension Valence(4)
      Dimension ExtraValence(4)
*----------------------------------------------------------------------*
* Is this a legal element?                                             *
*----------------------------------------------------------------------*
CNIKO If(Z.lt.1 .or. Z.gt.112) Then
CNIKO
      If(Z.lt.0 .or. Z.gt.112) Then
CNIKO
         Write(6,*) 'orbtype: do only know elements 1-112'
         Call Abend()
      End If
*----------------------------------------------------------------------*
* Initialize                                                           *
*----------------------------------------------------------------------*
      Do i=1,4
         DeepCore(i)=0
         Core(i)=0
         SoftCore(i)=0
         DeepValence(i)=0
         Valence(i)=0
         ExtraValence(i)=0
         List(i)=0
      End Do
*----------------------------------------------------------------------*
* How many shells are there for this atom                              *
*----------------------------------------------------------------------*
*
* Dummy
      If (Z.eq.0) Then
         Valence(1)=0 ! this is a redundant operation.
* H-He
      Else If(Z.le.2) Then
         Valence(1)=1
* Li-Be
      Else If(Z.le.4) Then
         Core(1)=1
         Valence(1)=1
         ExtraValence(2)=1
* B-Ne
      Else If(Z.le.10) Then
         Core(1)=1
         Valence(1)=1
         Valence(2)=1
* Na-Mg
      Else If(Z.le.12) Then
         DeepCore(1)=1
         SoftCore(1)=1
         SoftCore(2)=1
         Valence(1)=1
         ExtraValence(2)=1
* Al-Ar
      Else If(Z.le.18) Then
         DeepCore(1)=1
         Core(1)=1
         Core(2)=1
         Valence(1)=1
         Valence(2)=1
* K-Ca
      Else If(Z.le.20) Then
         DeepCore(1)=2
         DeepCore(2)=1
         SoftCore(1)=1
         SoftCore(2)=1
         Valence(1)=1
         ExtraValence(2)=1
* Sc-Zn
      Else If(Z.le.30) Then
         DeepCore(1)=2
         DeepCore(2)=1
         Core(1)=1
         Core(2)=1
         Valence(1)=1
         Valence(3)=1
         ExtraValence(2)=1
* Ga-Kr
      Else If(Z.le.36) Then
         DeepCore(1)=2
         DeepCore(2)=1
         Core(1)=1
         Core(2)=1
         Core(3)=1
         Valence(1)=1
         Valence(2)=1
* Rb-Sr
      Else If(Z.le.38) Then
         DeepCore(1)=3
         DeepCore(2)=2
         DeepCore(3)=1
         SoftCore(1)=1
         SoftCore(2)=1
         Valence(1)=1
         ExtraValence(2)=1
* Y-Cd
      Else If(Z.le.48) Then
         DeepCore(1)=3
         DeepCore(2)=2
         DeepCore(3)=1
         Core(1)=1
         Core(2)=1
         Valence(1)=1
         Valence(3)=1
         ExtraValence(2)=1
* In-Xe
      Else If(Z.le.54) Then
         DeepCore(1)=3
         DeepCore(2)=2
         DeepCore(3)=1
         Core(1)=1
         Core(2)=1
         Core(3)=1
         Valence(1)=1
         Valence(2)=1
* Cs-Ba
      Else If(Z.le.56) Then
         DeepCore(1)=4
         DeepCore(2)=3
         DeepCore(3)=2
         SoftCore(1)=1
         SoftCore(2)=1
         Valence(1)=1
         ExtraValence(1)=1
* La-Yb
      Else If(Z.le.70) Then
         DeepCore(1)=4
         DeepCore(2)=3
         DeepCore(3)=2
         Core(1)=1
         Core(2)=1
         Valence(1)=1
         Valence(4)=1
         ExtraValence(2)=1
* Lu-Hg
      Else If(Z.le.80) Then
         DeepCore(1)=4
         DeepCore(2)=3
         DeepCore(3)=2
         Core(1)=1
         Core(2)=1
         SoftCore(4)=1
         Valence(1)=1
         Valence(3)=1
         ExtraValence(2)=1
* Tl-Rn
      Else If(Z.le.86) Then
         DeepCore(1)=4
         DeepCore(2)=3
         DeepCore(3)=2
         Core(1)=1
         Core(2)=1
         Core(4)=1
         SoftCore(3)=1
         Valence(1)=1
         Valence(2)=1
* Fr-Ra
      Else If(Z.le.88) Then
         DeepCore(1)=5
         DeepCore(2)=4
         DeepCore(3)=3
         DeepCore(4)=1
         SoftCore(1)=1
         SoftCore(2)=1
         Valence(1)=1
         ExtraValence(2)=1
* Ac-No
      Else If(Z.le.102) Then
         DeepCore(1)=5
         DeepCore(2)=4
         DeepCore(3)=3
         DeepCore(4)=1
         Core(1)=1
         Core(2)=1
         Valence(1)=1
         Valence(4)=1
         ExtraValence(2)=1
* Lr-Cn
      Else If(Z.le.112) Then
         DeepCore(1)=5
         DeepCore(2)=4
         DeepCore(3)=3
         DeepCore(4)=1
         Core(1)=1
         Core(2)=1
         SoftCore(4)=1
         Valence(1)=1
         Valence(3)=1
         ExtraValence(2)=1
      Else
         Write(6,*) 'orbtype: element',Z,' not yet implemented'
         Call Abend()
      End If
*----------------------------------------------------------------------*
* Fill up the list to be returned                                      *
*----------------------------------------------------------------------*
      if(iAnd(opt,1).ne.0) Then
         Do i=1,4
            List(i)=List(i)+DeepCore(i)
         End Do
      End If
      if(iAnd(opt,2).ne.0) Then
         Do i=1,4
            List(i)=List(i)+Core(i)
         End Do
      End If
      if(iAnd(opt,4).ne.0) Then
         Do i=1,4
            List(i)=List(i)+SoftCore(i)
         End Do
      End If
      if(iAnd(opt,8).ne.0) Then
         Do i=1,4
            List(i)=List(i)+DeepValence(i)
         End Do
      End If
      if(iAnd(opt,16).ne.0) Then
         Do i=1,4
            List(i)=List(i)+Valence(i)
         End Do
      End If
      if(iAnd(opt,32).ne.0) Then
         Do i=1,4
            List(i)=List(i)+ExtraValence(i)
         End Do
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
