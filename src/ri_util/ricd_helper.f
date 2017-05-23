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
      Subroutine RICD_Helper(Do_nacCD_Basis,nTest,iAngMin_,iAngMax_,
     &                       jAngMin_,jAngMax_,nBS,iAng,jAng,list,
     &                       nBS_Max)
      Implicit Real*8 (a-h,o-z)
      Logical Do_nacCD_Basis
      Integer iAngMin_(0:nBS_Max-1),
     &        iAngMax_(0:nBS_Max-1)
      Parameter (iTabMx=15)
      Integer jAngMin_(0:nBS_Max-1,0:nBS_Max-1),
     &        jAngMax_(0:nBS_Max-1,0:nBS_Max-1)
      Integer list(2,0:((nTest+1)*(nTest+2))/2,0:nTest*2)
      Integer list2(0:nTest**2)
*                                                                      *
************************************************************************
*                                                                      *
      if(.Not.Do_nacCD_Basis) Then
*                                                                      *
************************************************************************
*                                                                      *
         nBS=(nTest+2)/2
         Do iBS=0, nBS-1
            iAngMin_(iBS)=iBS
            iAngMax_(iBS)=nTest-iBS
            Do iAng=0, iAngMax_(iBS)
               jAngMax_(iBS,iAng)=Min(iAng,iAngMin_(iBS))
               If (iAng.eq.iAngMax_(iBS))
     &             jAngMax_(iBS,iAng)=iAngMax_(iBS)
               If (iAng.lt.iAngMin_(iBS)) jAngMax_(iBS,iAng)=0
               jAngMin_(iBS,iAng)=iAngMin_(iBS)
               If (iAng.le.iAngMin_(iBS)) jAngMin_(iBS,iAng)=0
               Do jAng=jAngMin_(iBS,iAng), jAngMax_(iBS,iAng)
                  list(1,0,iAng)=iAng
                  list(2,0,iAng)=jAng
               End Do
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
            nBS=1
            iPair=0
            Do iBS=0, nBS-1
            iAngMax_(iBS)=nTest*2
            Do iAng=iAngMin_(iBS), iAngMax_(iBS)
               jAngMax_(iBS,iAng)=0
               jAngMin_(iBS,iAng)=0
               Do jAng=jAngMin_(iBS,iAng), jAngMax_(iBS,iAng)
                  list2(iAng)=0
                  Do k=0, nTest
                     Do l=0, k
                        Do m=iAng,0,-2
                           n=k-l
                           If ((n.eq.m).and.((k+l).ge.iAng)) Then
                              iPair=list2(iAng)
                              list(1,iPair,iAng)=l
                              list(2,iPair,iAng)=k
                              list2(iAng)=list2(iAng)+1
                           End If
                        End Do ! m
                     End Do    ! l
                  End Do       ! k
              End Do           ! jAng
            End Do             ! iAng
         End Do                ! iBS
      End if
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
