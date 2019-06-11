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
      Subroutine WR_VibRot_Info1(Lu,iOpt,iDisk,
     &                           ntit1,J1A,J2A,lambda,n0,
     &                           nvib1,Redm,Umax,Umin,ngrid,
     &                           isn1,isn2,Req,xMass1,xMass2)
      Implicit Real*8 (a-h,o-z)
      Dimension dum(1),idum(1)
*
      idum(1)=ntit1
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      ntit1=idum(1)
      idum(1)=J1A
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      J1A=idum(1)
      idum(1)=J2A
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      J2A=idum(1)
      idum(1)=lambda
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      lambda=idum(1)
      idum(1)=n0
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      n0=idum(1)
      idum(1)=nvib1
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      nvib1=idum(1)
      dum(1)=Redm
      Call dDaFile(Lu,iOpt,dum,1,iDisk)
      Redm=dum(1)
      dum(1)=Umax
      Call dDaFile(Lu,iOpt,dum,1,iDisk)
      Umax=dum(1)
      dum(1)=Umin
      Call dDaFile(Lu,iOpt,dum,1,iDisk)
      Umin=dum(1)
      idum(1)=ngrid
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      ngrid=idum(1)
      idum(1)=isn1
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      isn1=idum(1)
      idum(1)=isn2
      Call iDaFile(Lu,iOpt,idum,1,iDisk)
      isn2=idum(1)
      dum(1)=Req
      Call dDaFile(Lu,iOpt,dum,1,iDisk)
      Req=dum(1)
      dum(1)=xMass1
      Call dDaFile(Lu,iOpt,dum,1,iDisk)
      xMass1=dum(1)
      dum(1)=xMass2
      Call dDaFile(Lu,iOpt,dum,1,iDisk)
      xMass2=dum(1)
*
      Return
      End
