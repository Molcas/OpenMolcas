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
* Copyright (C) Anders Ohrn                                            *
************************************************************************
*  HelState
*
*> @brief
*>   Couple the electrostatic part of the solvent with the QM-region.
*>   Only include the static part, no polarization at this moment
*> @author A. Ohrn
*>
*> @note
*> The quadrupoles are put in 'Buckingham-style'.
*>
*> @details
*> Rather easy to follow. This subroutine is a slightly modified
*> copy of hel.f. The interesting quantities are collected in
*> \p Vmat and are later to be added to the 'RASSI-matrix'.
*>
*> @param[in]  Eint    Field from static part of solvent on the Qm-molecule centers
*> @param[in]  nrstate Number of states in RASSI
*> @param[in]  ici     Number of MME-centers
*> @param[in]  Cha     Charges
*> @param[in]  Dip     Dipoles
*> @param[in]  Qua     Quadrupoles
*> @param[out] Vmat    The electrostatic part of the solute-solvent interaction matrix
*> @param[in]  iPrint  Print-level
************************************************************************
      Subroutine HelState(Eint,nrstate,ici,Cha,Dip,Qua,Vmat,iPrint)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension Eint(MxQCen,10),Vmat(MxOt)
      Dimension Cha(MxStOT,MxQCen),Dip(MxStOT,3,MxQCen)
      Dimension Qua(MxStOT,6,MxQCen)

      kaunt=0
      Do 9, i=1,nrState
        Do 99, j=1,i
          kaunt=kaunt+1
          Vmat(kaunt)=0.0d0
99      Continue
9     Continue

      kaunt=0  !The interaction between the distributed multipoles
               !and the generalized field from the solvent.
      Do 10, i=1,nrState
        Do 11, j=1,i
          kaunt=kaunt+1
          Do 12, k=1,ici
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,1)*Cha(kaunt,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,2)*Dip(kaunt,1,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,3)*Dip(kaunt,2,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,4)*Dip(kaunt,3,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,5)*Qua(kaunt,1,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,7)*Qua(kaunt,3,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,10)*Qua(kaunt,6,k)
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,6)*Qua(kaunt,2,k)*2
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,8)*Qua(kaunt,4,k)*2
            Vmat(kaunt)=Vmat(kaunt)+Eint(k,9)*Qua(kaunt,5,k)*2
12        Continue
11      Continue
10    Continue

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iPrint)
      End
