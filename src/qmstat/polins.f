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
      Subroutine Polins(Energy,iCall,iAtom2,iCi,iFil,VpolMat,fil,polfac
     &     ,poli,xyzmyq,xyzmyi,xyzmyp,iCstart,iQ_Atoms,qtot,ChaNuc
     &     ,RoMatSt,xyzQuQ,CT)
************************************************************
*
*   <DOC>
*     <Name>Polins</Name>
*     <Syntax>Call Polins(Energy,iCall,iAtom2,iCi,iFil,VpolMat,fil,polfac,poli,xyzmyq,xyzmyi,xyzmyp,qtot,iCstart,iAtom)</Syntax>
*     <Arguments>
*       \Argument{Energy}{The energy of the electrostatic interaction}{}{out}
*       \Argument{iCall}{An integer that tells if this is the first call in the iteration. Necessary for the copy of the one-particle hamiltonian}{}{inout}
*       \Argument{iAtom2}{Number of particles in the solvent, times number of polarizabilities per solvent molecule}{}{in}
*       \Argument{iCi}{Number of centers in QM-molecule}{}{in}
*       \Argument{iFil}{Pointer to the static field from the solvent}{}{in}
*       \Argument{VpolMat}{The polarization matrix}{}{out}
*       \Argument{fil}{The field from the induced dipoles in the solvent}{}{inout}
*       \Argument{polfac}{A factor for the computation of the image}{}{in}
*       \Argument{poli}{The solvent polarized field on QM}{}{out}
*       \Argument{xyzmyq}{Total dipole of QM-region}{}{in}
*       \Argument{xyzmyi}{Total induced dipole of solvent}{}{in}
*       \Argument{xyzmyp}{Total permanent dipole of solvent}{}{in}
*       \Argument{qtot}{Total charge of QM-region}{}{in}
*       \Argument{iCstart}{Number to keep track of solvent molecules}{}{in}
*       \Argument{iAtom}{Number of atoms in QM-molecule}{}{in}
*     </Arguments>
*     <Purpose>
*    Add the field from the QM-region onto the solvent. Include the
*    field from the polarizabilities in the solvent onto the QM-region.
*    (The effect of the static field is taken care of in helstate.f.)
*     </Purpose>
*     <Dependencies></Dependencies>
*     <Author>A.Ohrn</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*    First we compute the field from the QM-region onto the solvent. The
*    central quantity is the density matrix, which gives us how the
*    QM-region polarizes. The reaction field due to the QM-region is also
*    accounted for. Then we allow the polarized field from the solvent
*    to interact with the QM-region. The static field from the solvent,
*    in other word that from the charges, is already coupled in
*    helstate.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qm2.fh"
#include "WrkSpc.fh"

      Dimension Fil(npart*npol,3),Qm(MxQCen),Dm(MxQCen,3),QQm(MxQCen,6)
     &,Poli(MxQCen,10),Gunnar(10),ChaNuc(MxAt)
     &,Eil(MxPut*MxPol,3),xyzMyC(3),CofC(3),xyzmyq(3),xyzmyi(3)
     &,xyzmyp(3),VpolMat(MxStOt),RoMatSt(MxStOT)
      Dimension xyzQuQ(6),qQ(6),qD(6),qK(6),CT(3)
      Dimension iFil(MxQCen,10)

*----------------------------------------------------------------------*
* Begin with some zeros.                                               *
*----------------------------------------------------------------------*
      iCnum=iCStart/Ncent
      Do 2, i=1,iCi
        Qm(i)=0
        If(i.le.iQ_Atoms) Qm(i)=-ChaNuc(i)
        Do 3, j=1,3
          Dm(i,j)=0
          QQm(i,j)=0
          QQm(i,j+3)=0
3       Continue
2     Continue
      Do 6644, i=1,MxPut*MxPol
        Do 6645, j=1,3
          Eil(i,j)=0
6645    Continue
6644  Continue
*----------------------------------------------------------------------*
* Below we compute how the MME of the QM-molecule changes with the new *
* density matrix Romat. In this step we connect a state density with   *
* various multipoles, which will go on and interact with the solvent.  *
* This is one of the pivotal steps in coupling the QM-region with the  *
* classical region.                                                    *
*----------------------------------------------------------------------*
      kaunt=0
      Do 14, i=1,nState
        Do 15, j=1,i
          kaunt=kaunt+1
          Do 41, k=1,iCi
            Qm(k)=Qm(k)+RasCha(kaunt,k)*RomatSt(kaunt)
            Dm(k,1)=Dm(k,1)+RasDip(kaunt,1,k)*RomatSt(kaunt)
            Dm(k,2)=Dm(k,2)+RasDip(kaunt,2,k)*RomatSt(kaunt)
            Dm(k,3)=Dm(k,3)+RasDip(kaunt,3,k)*RomatSt(kaunt)
            QQm(k,1)=QQm(k,1)+RasQua(kaunt,1,k)*RomatSt(kaunt)
            QQm(k,3)=QQm(k,3)+RasQua(kaunt,3,k)*RomatSt(kaunt)
            QQm(k,6)=QQm(k,6)+RasQua(kaunt,6,k)*RomatSt(kaunt)
            QQm(k,2)=QQm(k,2)+RasQua(kaunt,2,k)*RomatSt(kaunt)
            QQm(k,4)=QQm(k,4)+RasQua(kaunt,4,k)*RomatSt(kaunt)
            QQm(k,5)=QQm(k,5)+RasQua(kaunt,5,k)*RomatSt(kaunt)
41        Continue
15      Continue
14    Continue
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
      Do 866,i=1,iCi
        Do 8661, kk=1,3
          xyzMyQ(kk)=xyzMyQ(kk)+Dm(i,kk)+Qm(i)*outxyzRAS(i,kk)
8661    Continue
        qQ(1)=qQ(1)+Qm(i)*(outxyzRAS(i,1)-CT(1))*(outxyzRAS(i,1)-CT(1))
        qQ(2)=qQ(2)+Qm(i)*(outxyzRAS(i,1)-CT(1))*(outxyzRAS(i,2)-CT(2))
        qQ(3)=qQ(3)+Qm(i)*(outxyzRAS(i,1)-CT(1))*(outxyzRAS(i,3)-CT(3))
        qQ(4)=qQ(4)+Qm(i)*(outxyzRAS(i,2)-CT(2))*(outxyzRAS(i,2)-CT(2))
        qQ(5)=qQ(5)+Qm(i)*(outxyzRAS(i,2)-CT(2))*(outxyzRAS(i,3)-CT(3))
        qQ(6)=qQ(6)+Qm(i)*(outxyzRAS(i,3)-CT(3))*(outxyzRAS(i,3)-CT(3))
        qD(1)=qD(1)+2*Dm(i,1)*(outxyzRAS(i,1)-CT(1))
        qD(2)=qD(2)+Dm(i,1)*(outxyzRAS(i,2)-CT(2))
     &             +Dm(i,2)*(outxyzRAS(i,1)-CT(1))
        qD(3)=qD(3)+Dm(i,1)*(outxyzRAS(i,3)-CT(3))
     &             +Dm(i,3)*(outxyzRAS(i,1)-CT(1))
        qD(4)=qD(4)+2*Dm(i,2)*(outxyzRAS(i,2)-CT(2))
        qD(5)=qD(5)+Dm(i,2)*(outxyzRAS(i,3)-CT(3))
     &             +Dm(i,3)*(outxyzRAS(i,2)-CT(2))
        qD(6)=qD(6)+2*Dm(i,3)*(outxyzRAS(i,3)-CT(3))
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
          CofC(1)=CofC(1)+abs(qm(i))*outxyzRAS(i,1) !Center of charge
          CofC(2)=CofC(2)+abs(qm(i))*outxyzRAS(i,2)
          CofC(3)=CofC(3)+abs(qm(i))*outxyzRAS(i,3)
          qs=qs+abs(qm(i))
721     Continue
        CofC(1)=CofC(1)/qs
        CofC(2)=CofC(2)/qs
        CofC(3)=CofC(3)/qs
        Gx=CofC(1)-outxyzRAS(1,1)+Cordst(1,1) !Where C-of-C is globally
        Gy=CofC(2)-outxyzRAS(1,2)+Cordst(1,2)
        Gz=CofC(3)-outxyzRAS(1,3)+Cordst(1,3)
        xyzMyC(1)=xyzMyC(1)+qtot*Gx !Dipole
        xyzMyC(2)=xyzMyC(2)+qtot*Gy
        xyzMyC(3)=xyzMyC(3)+qtot*Gz
      Endif
      Do 9977, i=1,3  !The energy of the induced dipole in its reaction
                      !field. It is ok since polfac*(xyzMyQ+xyzMyi) is
                      !the field from the induced dipole according to
                      !the image approximation. And the sought energy
                      !is -0.5*my_perm*R_ind, which is the thing below,
                      !see Boethcer p. 145.
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
                      !cavity model! Fil now contains the field on each
                      !solvent polarizability from all different
                      !sources.
          Fil(i,j)=Fil(i,j)+PolFac*xyzMyQ(j)+Eil(i,j)
          Energy=Energy+Fil(i,j)*Eil(i,j)*Pol(iu)
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
        Gunnar(1)=Gunnar(2)*outxyzRAS(l,1)+Gunnar(3)*outxyzRAS(l,2)
     &           +Gunnar(4)*outxyzRAS(l,3) !The potential from the
                                           !dipole (the 1/r*r*r is
                                           !in Work(iFil...).
        Do 1305, i=1,10
          Poli(l,i)=Gunnar(i)
          Do 1306, j=1+(nPol*iCnum),iAtom2
            Iu=j-((j-1)/nPol)*nPol
            Do 1307, k=1,3
       !Poli is the polarized field of the solvent. Remember that
       !Fil() is the total field on the polarizabilities.
              Poli(l,i)=Poli(l,i)-Fil(j,k)*Pol(iu)
     &                 *Work(iFil(l,i)-1+j+(k-1)*nPart*nPol)
1307        Continue
1306      Continue
1305    Continue
1304  Continue
      Do 299, i=1,nState*(nState+1)/2
        VpolMat(i)=0
299   Continue
      kaunt=0
      Do 300, iS=1,nState  !Attention! The reason we use RasCha etc. and
        Do 301, jS=1,iS   !not the computed Qm, Dm etc. from above is
          kaunt=kaunt+1  !that the density we want to describe is the
          Do 302, j=1,iCi !density of the basis-functions. Compare with
                          !ordinary <psi_i|V_el|psi_j>.
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,1)*RasCha(kaunt,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,2)*RasDip(kaunt,1,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,3)*RasDip(kaunt,2,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,4)*RasDip(kaunt,3,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,5)*RasQua(kaunt,1,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,7)*RasQua(kaunt,3,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,10)*RasQua(kaunt,6,j)
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,6)*RasQua(kaunt,2,j)*2
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,8)*RasQua(kaunt,4,j)*2
            Vpolmat(kaunt)=Vpolmat(kaunt)+Poli(j,9)*RasQua(kaunt,5,j)*2
302       Continue
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
