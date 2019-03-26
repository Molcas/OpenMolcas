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
      Integer Function nProp_Int(Do_Index,Index,nIndex)
      Implicit real*8 (a-h,o-z)
#include "real.fh"
      Character*8 Label
      Logical Do_Index
      Integer Index(4,nIndex)
      Dimension nInt(1)
*                                                                      *
************************************************************************
*                                                                      *
      nProp_Int=0
*                                                                      *
************************************************************************
*                                                                      *
c     Scan the ONEINT file for multipole moment operators
c
C     Write (*,*) ' Starting scan of ONEINT for multipole moments'
      do iMltpl=1,99
         nComp=(iMltpl+1)*(iMltpl+2)/2
c
         Write (Label,'(a,i2)') 'MLTPL ',iMltpl
         irc=-1
         iopt=1
         iComp=1
         Call iRdOne (irc,iopt,Label,iComp,nInt,iSmLbl)
         If (irc.ne.0) Go To 110
         If (Do_Index) Then
            Do iComp=1,nComp
               nProp_Int = nProp_Int + 1
               Index(1,nProp_Int) = 1
               Index(2,nProp_Int) = iMltpl
               Index(3,nProp_Int) = iComp
               Index(4,nProp_Int) = 0
            End Do
         Else
            nProp_Int=nProp_Int + nComp
         End If
      End Do
110   Continue
*                                                                      *
************************************************************************
*                                                                      *
c     Scan 'ONEINT' for electric field integrals
c
C     Write (*,*) ' Starting scan of ONEINT for various elec. field integrals'
c
      Do iEF=0,2
         nComp=(iEF+1)*(iEF+2)/2
c
         maxCen=9999
         Do iCent=1,maxCen
            write (Label,'(a,i1,i5)') 'EF',iEF,iCent
            irc=-1
            iopt=1
            iComp=1
            Call iRdOne (irc,iopt,Label,iComp,nInt,iSmLbl)
            If (irc.ne.0) Go To 201
            If (Do_Index) Then
               Do iComp=1,nComp
                  nProp_Int = nProp_Int + 1
                  Index(1,nProp_Int) = 2
                  Index(2,nProp_Int) = iEF
                  Index(3,nProp_Int) = iComp
                  Index(4,nProp_Int) = iCent
               End Do
            Else
               nProp_Int=nProp_Int + nComp
            End If
         End Do
201      Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
c     Scan 'ONEINT' for contact term integrals
c
C     Write (*,*) ' Starting scan of ONEINT for various contact term integrals'
c
      nComp = 1
c
      maxCen=9999
      Do iCent=1,maxCen
         write (Label,'(a,i5)') 'Cnt',iCent
         irc=-1
         iopt=1
         iComp=1
         Call iRdOne (irc,iopt,Label,iComp,nInt,iSmLbl)
         If (irc.ne.0) Go To 301
         If (Do_Index) Then
            Do iComp=1,nComp
               nProp_Int = nProp_Int + 1
               Index(1,nProp_Int) = 3
               Index(2,nProp_Int) = 1
               Index(3,nProp_Int) = 1
               Index(4,nProp_Int) = iCent
            End Do
         Else
            nProp_Int=nProp_Int + nComp
         End If
      End Do
301   Continue
*                                                                      *
************************************************************************
*                                                                      *
C     If (Do_Index) Call iVcPrt('Index',' ',Index,4*nProp_Int)
      Return
      End
