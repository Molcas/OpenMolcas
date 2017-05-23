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
*----------------------------------------------------------------------*
* A routine inspired by the wr_motra_info utility.                     *
*----------------------------------------------------------------------*
      Subroutine WrRdSim(iLu,iOpt,iDisk,iTcSim,nTcSim,Etot
     &,Radie,nPart,Gamma,Gam,Esav)
      Implicit Real*8 (a-h,o-z)
      Dimension iTcSim(nTcSim)

      Call iDaFile(iLu,iOpt,iTcSim,nTcSim,iDisk)
      Call dDaFile(iLu,iOpt,Etot,1,iDisk)
      Call dDaFile(iLu,iOpt,Radie,1,iDisk)
      Call iDaFile(iLu,iOpt,nPart,1,iDisk)
      Call dDaFile(iLu,iOpt,Gamma,1,iDisk)
      Call dDaFile(iLu,iOpt,Gam,1,iDisk)
      Call dDaFile(iLu,iOpt,Esav,1,iDisk)

      Return
      End
