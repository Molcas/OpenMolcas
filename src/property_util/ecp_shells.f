************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine ECP_shells(iAtmNr,List)
#include "itmax.fh"
      Integer List(0:iTabMx)
*                                                                      *
************************************************************************
*                                                                      *
      Call ICopy(iTabMx+1,[0],0,List,1)
*                                                                      *
************************************************************************
*                                                                      *
* Dummy
      If (iAtmNr.eq.0) Then
         list(0)=0
         list(1)=0
         list(2)=0
         list(3)=0
* H-He
      Else If((iAtmNr.gt.2) .and. (iAtmNr.le.4)) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.10) Then
         list(0)=1
         list(1)=1
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.12) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.18) Then
         list(0)=1
         list(1)=1
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.20) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.30) Then
         list(0)=1
         list(1)=0
         list(2)=1
         list(3)=0
      Else If(iAtmNr.le.36) Then
         list(0)=1
         list(1)=1
         list(2)=1
         list(3)=0
      Else If(iAtmNr.le.38) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.48) Then
         list(0)=1
         list(1)=0
         list(2)=1
         list(3)=0
      Else If(iAtmNr.le.54) Then
         list(0)=1
         list(1)=1
         list(2)=1
         list(3)=0
      Else If(iAtmNr.le.56) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.70) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=1
      Else If(iAtmNr.le.80) Then
         list(0)=1
         list(1)=0
         list(2)=1
         list(3)=1
      Else If(iAtmNr.le.86) Then
         list(0)=1
         list(1)=1
         list(2)=1
         list(3)=1
      Else If(iAtmNr.le.88) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else If(iAtmNr.le.102) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=1
      Else If(iAtmNr.le.112) Then
         list(0)=1
         list(1)=0
         list(2)=1
         list(3)=1
      Else If(iAtmNr.le.118) Then
         list(0)=1
         list(1)=1
         list(2)=1
         list(3)=1
      Else If(iAtmNr.le.120) Then
         list(0)=1
         list(1)=0
         list(2)=0
         list(3)=0
      Else
          Write (6,*) 'ECP_shells can not handle atom numbers beyond'
     &               //' 112.'
          Call Abend()
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
