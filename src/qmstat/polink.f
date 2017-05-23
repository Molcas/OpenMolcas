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
*----------------------------------------------------------------------*
      Subroutine Polink(Energy,iCall,iAtom2,iCi,iFil,VpolMat,fil,polfac
     &,poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ,xyzMyI,xyzMyP
     &,RoMat,xyzQuQ,CT)
************************************************************
*
*   <DOC>
*     <Name>Polink</Name>
*     <Syntax>Call Polink(Energy,iCall,iAtom2,iCi,iFil,VpolMat,fil,polfac,poli,iCstart,iTri,iAtom)</Syntax>
*     <Arguments>
*       \Argument{Energy}{The energy of the electrostatic interaction}{}{out}
*       \Argument{iCall}{An integer that tells if this is the first call in the iteration. Necessary for the copy of the one-particle hamiltonian}{}{inout}
*       \Argument{iAtom2}{Number of particles in the solvent, times number of polarizabilities per solvent molecule}{}{in}
*       \Argument{iCi}{Number of centers in QM-molecule}{}{in}
*       \Argument{iFil}{Pointer to the static field from the solvent}{}{in}
*       \Argument{VpolMat}{The matrix due to polarization}{}{out}
*       \Argument{fil}{The field from the induced dipoles in the solvent}{}{inout}
*       \Argument{polfac}{A factor for the computation of the image}{}{in}
*       \Argument{poli}{The solvent polarized field on QM-region}{}{out}
*       \Argument{iCstart}{Number to keep track of solvent molecules}{}{in}
*       \Argument{iTri}{iOrb(1)*(iOrb(1)+1)/2}{}{in}
*       \Argument{iAtom}{Number of atoms in QM-molecule}{}{in}
*     </Arguments>
*     <Purpose>
*    Add the field from the QM-region onto the solvent. Include the
*    field from the polarizabilities in the solvent onto the QM-region.
*    (The effect of the static field is taken care of in hel.f.)
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>A.Ohrn</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*    To begin with we obtain the charge distribution of the QM-region
*    as it exists due to the pressent density matrix (recall that it
*    is the changes in the density matrix that causes the QM-region to
*    be polarized). The field from these new multipoles are added on to
*    solvent. We also include the reaction field from the QM-region.
*    Then, with the new field from the QM-region included, we compute
*    the field from the polarizabiolities in the solvent onto the QM-region,
*    which is done just like in hel.f.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qm1.fh"
#include "WrkSpc.fh"

      Dimension Fil(npart*npol,3),Qm(MxQCen),Dm(MxQCen,3),QQm(MxQCen,6)
     &,Poli(MxQCen,10),Gunnar(10),Eil(MxPut*MxPol,3),xyzMyC(3),CofC(3)
     &,VpolMat(MxOt),ChaNuc(MxAt),xyzMyQ(3),xyzMyI(3),xyzMyP(3)
     &,RoMat(MxOT)
      Dimension xyzQuQ(6),qQ(6),qD(6),qK(6),CT(3)
      Dimension iFil(MxQCen,10)

*----------------------------------------------------------------------*
* Begin with some zeros.                                               *
*----------------------------------------------------------------------*
      iCnum=iCStart/Ncent
      Do 2, i=1,iCi
        Qm(i)=0.0d0
        If(i.le.iQ_Atoms) Qm(i)=-ChaNuc(i) !Here is nuclear contribution
        Do 3, j=1,3                        !added to atoms.
          Dm(i,j)=0.0d0
          QQm(i,j)=0.0d0
          QQm(i,j+3)=0.0d0
3       Continue
2     Continue
      Do 6644, i=1,MxPut*MxPol
        Do 6645, j=1,3
          Eil(i,j)=0.0d0
6645    Continue
6644  Continue
*----------------------------------------------------------------------*
* Below we compute how the MME of the QM-molecule changes with the new *
* density matrix Romat. What we actually do is a HF-SCF procedure with *
* a MME-expanded density. A change in the density has the effect that  *
* the set of multipoles in the MME are slightly perturbed.             *
*----------------------------------------------------------------------*
      Do 4, i=1,iTri
        Do 41, j=1,iCi
          Qm(j)=Cha(i,j)*Romat(i)+Qm(j)
          Dm(j,1)=Dm(j,1)+DipMy(i,1,j)*Romat(i)
          Dm(j,2)=Dm(j,2)+DipMy(i,2,j)*Romat(i)
          Dm(j,3)=Dm(j,3)+DipMy(i,3,j)*Romat(i)
          QQm(j,1)=QQm(j,1)+Quad(i,1,j)*Romat(i)
          QQm(j,3)=QQm(j,3)+Quad(i,3,j)*Romat(i)
          QQm(j,6)=QQm(j,6)+Quad(i,6,j)*Romat(i)
          QQm(j,2)=QQm(j,2)+Quad(i,2,j)*Romat(i)
          QQm(j,4)=QQm(j,4)+Quad(i,4,j)*Romat(i)
          QQm(j,5)=QQm(j,5)+Quad(i,5,j)*Romat(i)
41      Continue
4     Continue
      Do 775, kk=1,3
        xyzMyQ(kk)=0
        xyzMyC(kk)=0
        CofC(kk)=0
775   Continue
      Do 776, kk=1,6
        xyzQuQ(kk)=0
        qQ(kk)=0
        qD(kk)=0
        qK(kk)=0
776   Continue
      !Observe one trixy thing about xyzmyq: the electric multipoles
      !we use above are actually of opposite sign, so how can xyzmyq be
      !the dipole in the qm-region unless we change sign (which we does
      !not)? The reason is that the density matrix elements will also
      !have opposite sign, which in turn has not physical meaning.
      !We also compute the quadupole moment - a mezzy formula.
      Do 866,i=1,iCi
        xyzMyQ(1)=xyzMyQ(1)+Dm(i,1)+Qm(i)*outxyz(i,1)
        xyzMyQ(2)=xyzMyQ(2)+Dm(i,2)+Qm(i)*outxyz(i,2)
        xyzMyQ(3)=xyzMyQ(3)+Dm(i,3)+Qm(i)*outxyz(i,3)
        qQ(1)=qQ(1)+Qm(i)*(outxyz(i,1)-CT(1))*(outxyz(i,1)-CT(1))
        qQ(2)=qQ(2)+Qm(i)*(outxyz(i,1)-CT(1))*(outxyz(i,2)-CT(2))
        qQ(3)=qQ(3)+Qm(i)*(outxyz(i,1)-CT(1))*(outxyz(i,3)-CT(3))
        qQ(4)=qQ(4)+Qm(i)*(outxyz(i,2)-CT(2))*(outxyz(i,2)-CT(2))
        qQ(5)=qQ(5)+Qm(i)*(outxyz(i,2)-CT(2))*(outxyz(i,3)-CT(3))
        qQ(6)=qQ(6)+Qm(i)*(outxyz(i,3)-CT(3))*(outxyz(i,3)-CT(3))
        qD(1)=qD(1)+2*Dm(i,1)*(outxyz(i,1)-CT(1))
        qD(2)=qD(2)+Dm(i,1)*(outxyz(i,2)-CT(2))
     &             +Dm(i,2)*(outxyz(i,1)-CT(1))
        qD(3)=qD(3)+Dm(i,1)*(outxyz(i,3)-CT(3))
     &             +Dm(i,3)*(outxyz(i,1)-CT(1))
        qD(4)=qD(4)+2*Dm(i,2)*(outxyz(i,2)-CT(2))
        qD(5)=qD(5)+Dm(i,2)*(outxyz(i,3)-CT(3))
     &             +Dm(i,3)*(outxyz(i,2)-CT(2))
        qD(6)=qD(6)+2*Dm(i,3)*(outxyz(i,3)-CT(3))
        qK(1)=qK(1)+QQm(i,1)
        qK(2)=qK(2)+QQm(i,2)
        qK(3)=qK(3)+QQm(i,4)
        qK(4)=qK(4)+QQm(i,3)
        qK(5)=qK(5)+QQm(i,5)
        qK(6)=qK(6)+QQm(i,6)
866   Continue
      Trace1=qQ(1)+qQ(4)+qQ(6)
      Trace2=qD(1)+qD(4)+qD(6)
      Trace1=Trace1/3
      Trace2=Trace2/3
      xyzQuQ(1)=1.5*(qQ(1)+qD(1)-Trace1-Trace2)+qK(1)
      xyzQuQ(2)=1.5*(qQ(2)+qD(2))+qK(2)
      xyzQuQ(3)=1.5*(qQ(3)+qD(3))+qK(3)
      xyzQuQ(4)=1.5*(qQ(4)+qD(4)-Trace1-Trace2)+qK(4)
      xyzQuQ(5)=1.5*(qQ(5)+qD(5))+qK(5)
      xyzQuQ(6)=1.5*(qQ(6)+qD(6)-Trace1-Trace2)+qK(6)
      If(ChargedQM) then  !If charged system, then do...
        qs=0
        Do 721, i=1,iCi
          CofC(1)=CofC(1)+abs(qm(i))*outxyz(i,1) !Center of charge
          CofC(2)=CofC(2)+abs(qm(i))*outxyz(i,2)
          CofC(3)=CofC(3)+abs(qm(i))*outxyz(i,3)
          qs=qs+abs(qm(i))
721     Continue
        CofC(1)=CofC(1)/qs
        CofC(2)=CofC(2)/qs
        CofC(3)=CofC(3)/qs
        Gx=CofC(1)-outxyz(1,1)+Cordst(1,1) !Where C-of-C is globally
        Gy=CofC(2)-outxyz(1,2)+Cordst(1,2)
        Gz=CofC(3)-outxyz(1,3)+Cordst(1,3)
        xyzMyC(1)=xyzMyC(1)+qtot*Gx !Dipole
        xyzMyC(2)=xyzMyC(2)+qtot*Gy
        xyzMyC(3)=xyzMyC(3)+qtot*Gz
      Endif
      Do 9977, i=1,3
          !Change sign on both the dipoles, which in effect gives
          !no sign change, all in order with Boettcher, p.145.
        Energy=Energy+Polfac*xyzMyQ(i)*(xyzMyQ(i)+xyzMyi(i))
9977  Continue
*----------------------------------------------------------------------*
* The multipoles of the QM-region, modified due to the polarization,   *
* now interacts with each polarizability in the solvent.               *
*----------------------------------------------------------------------*
      Do 5, i=1,iCi
        Do 6, j=1+(nPol*iCnum),iAtom2
          Do 7, k=1,3
            Eil(j,k)=Eil(j,k)+Work(iFil(i,1)-1+j+(k-1)*nPart*nPol)*Qm(i)
            Do 8, l=1,3
              Eil(j,k)=Eil(j,k)+Work(iFil(i,l+1)-1+j+(k-1)*nPart*nPol)
     &                *Dm(i,l)
8           Continue
         Eil(j,k)=Eil(j,k)+Work(iFil(i,5)-1+j+(k-1)*nPart*nPol)*QQm(i,1)
         Eil(j,k)=Eil(j,k)+Work(iFil(i,7)-1+j+(k-1)*nPart*nPol)*QQm(i,3)
        Eil(j,k)=Eil(j,k)+Work(iFil(i,10)-1+j+(k-1)*nPart*nPol)*QQm(i,6)
       Eil(j,k)=Eil(j,k)+Work(iFil(i,6)-1+j+(k-1)*nPart*nPol)*QQm(i,2)*2
       Eil(j,k)=Eil(j,k)+Work(iFil(i,8)-1+j+(k-1)*nPart*nPol)*QQm(i,4)*2
       Eil(j,k)=Eil(j,k)+Work(iFil(i,9)-1+j+(k-1)*nPart*nPol)*QQm(i,5)*2
7         Continue
6       Continue
5     Continue
C...THIS IS LEBENSGEFAHRLICH (original comment says it all!)
*----------------------------------------------------------------------*
* We add up the field from the QM-region to the field on all the       *
* solvent polarizabilities.                                            *
*----------------------------------------------------------------------*
      Do 10, i=1+(nPol*iCNum),iAtom2
        Iu=i-((i-1)/nPol)*nPol
        Do 11, j=1,3  !Here we add the QM-molecule image to the solvent
                      !polarizabilities. Good old classical dielectric
                      !cavity model!
          Fil(i,j)=Fil(i,j)+PolFac*xyzMyQ(j)+Eil(i,j)
          Energy=Energy+Fil(i,j)*Eil(i,j)*Pol(iu) !How much the induced
                            !dipoles in solvent interacts with the field
                            !from the QM-region.
11      Continue
10    Continue
*----------------------------------------------------------------------*
* Now we wish to make the induced field from the solvent interact with *
* the QM-region. The static field has already interacted in helstate.f.*
* The reaction field of the QM-region in the dielectric cavity is      *
* also included, excluding the quadrupoles and higher; they are        *
* small anyway, so this is not a major restriction.                    *
*----------------------------------------------------------------------*
      Do 1801, i=1,10
        Gunnar(i)=0
1801  Continue
      Gunnar(2)=PolFac*(xyzMyP(1)+xyzMyQ(1)+xyzMyI(1)+xyzMyC(1))
      Gunnar(3)=PolFac*(xyzMyP(2)+xyzMyQ(2)+xyzMyI(2)+xyzMyC(2))
      Gunnar(4)=PolFac*(xyzMyP(3)+xyzMyQ(3)+xyzMyI(3)+xyzMyC(3))
      Do 1304, l=1,iCi
        Gunnar(1)=Gunnar(2)*outxyz(l,1)+Gunnar(3)*outxyz(l,2)
     &           +Gunnar(4)*outxyz(l,3) !Potential from the apparent
        Do 1305, i=1,10          !surface charge, see Boettcher (4.22).
          Poli(l,i)=Gunnar(i)
          Do 1306, j=1+(nPol*iCnum),iAtom2
            Iu=j-((j-1)/nPol)*nPol
            Do 1307, k=1,3            !Compute the generalized field
                              !from induced dipoles in solvent on the
                              !QM-region cites.
              Poli(l,i)=Poli(l,i)-Fil(j,k)*Pol(iu)
     &                 *Work(iFil(l,i)-1+j+(k-1)*nPart*nPol)
1307        Continue
1306      Continue
1305    Continue
1304  Continue
      Do 201, i=1,iTri
        VpolMat(i)=0
201   Continue
      Do 300, i=1,iTri
        Do 301, j=1,iCi
          Vpolmat(i)=Vpolmat(i)+Poli(j,1)*Cha(i,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,2)*DipMy(i,1,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,3)*DipMy(i,2,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,4)*DipMy(i,3,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,5)*Quad(i,1,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,7)*Quad(i,3,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,10)*Quad(i,6,j)
          Vpolmat(i)=Vpolmat(i)+Poli(j,6)*Quad(i,2,j)*2
          Vpolmat(i)=Vpolmat(i)+Poli(j,8)*Quad(i,4,j)*2
          Vpolmat(i)=Vpolmat(i)+Poli(j,9)*Quad(i,5,j)*2
301     Continue
300   Continue
      Do 400, i=1,iQ_Atoms  !This is how the nuclei interact with the
                      !induced field (in equil2 exists a corresponding
                      !term for the interaction with the static field).
                      !This way interaction between a charged molecule
                      !and the induced/permanent potential is included.
        Energy=Energy-2*Poli(i,1)*ChaNuc(i)
400   Continue

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iCall)
      End
