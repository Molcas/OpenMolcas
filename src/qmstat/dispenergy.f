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
      Subroutine DispEnergy(EEDisp,BoMaH,BoMaO,dAtO1,dAtH1,dAtH2
     &                     ,Rab13i,Rab23i,Rab33i,indQAt)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      Dimension BoMaH(MxAt),BoMaO(MxAt)
      Dimension DampBarn(3)

*
*--- If Damping, do it, otherwise easy stuff
*
      If(DispDamp) then
*
*--- Get the damping, for the relevant QM-atom with all solvent
*    atoms.
*
        kFac=1
        DampBarn(1)=1.0d0
        Do 661, k=1,6
          kFac=kFac*k
          DampBarn(1)=DampBarn(1)+(BoMaH(indQAt)*dAtH1)**k/dble(kFac)
661     Continue
        DampBarn(1)=1.0d0-DampBarn(1)*exp(-BoMaH(indQAt)*dAtH1)

        kFac=1
        DampBarn(2)=1.0d0
        Do 662, k=1,6
          kFac=kFac*k
          DampBarn(2)=DampBarn(2)+(BoMaH(indQAt)*dAtH2)**k/dble(kFac)
662     Continue
        DampBarn(2)=1.0d0-DampBarn(2)*exp(-BoMaH(indQAt)*dAtH2)

        kFac=1
        DampBarn(3)=1.0d0
        Do 663, k=1,6
          kFac=kFac*k
          DampBarn(3)=DampBarn(3)+(BoMaO(indQAt)*dAtO1)**k/dble(kFac)
663     Continue
        DampBarn(3)=1.0d0-DampBarn(3)*exp(-BoMaO(indQAt)*dAtO1)

*
*--- If not damping, set factors to 1.0d0
*
      Else
        DampBarn(1)=1.0d0
        DampBarn(2)=1.0d0
        DampBarn(3)=1.0d0
      Endif

*
*--- Now evaluate the Dispersion energy.
*
      EfromO1=Rab13i**2*DampBarn(3)*uDisp(indQAt,1)
      EfromH1=Rab23i**2*DampBarn(1)*uDisp(indQAt,2)
      EfromH2=Rab33i**2*DampBarn(2)*uDisp(indQAt,2)

      EEDisp=EEdisp+EfromO1+EfromH1+EfromH2

      Return
      End
