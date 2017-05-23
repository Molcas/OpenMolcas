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
*
      Call iDaFile(Lu,iOpt,ntit1, 1,iDisk)
      Call iDaFile(Lu,iOpt,J1A,   1,iDisk)
      Call iDaFile(Lu,iOpt,J2A,   1,iDisk)
      Call iDaFile(Lu,iOpt,lambda,1,iDisk)
      Call iDaFile(Lu,iOpt,n0,    1,iDisk)
      Call iDaFile(Lu,iOpt,nvib1, 1,iDisk)
      Call dDaFile(Lu,iOpt,Redm,  1,iDisk)
      Call dDaFile(Lu,iOpt,Umax,  1,iDisk)
      Call dDaFile(Lu,iOpt,Umin,  1,iDisk)
      Call iDaFile(Lu,iOpt,ngrid, 1,iDisk)
      Call iDaFile(Lu,iOpt,isn1,  1,iDisk)
      Call iDaFile(Lu,iOpt,isn2,  1,iDisk)
      Call dDaFile(Lu,iOpt,Req,   1,iDisk)
      Call dDaFile(Lu,iOpt,xMass1,1,iDisk)
      Call dDaFile(Lu,iOpt,xMass2,1,iDisk)
*
      Return
      End
