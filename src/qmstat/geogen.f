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
*  Geogen
*
*> @brief
*>   Generate the new configuration in typical random manner
*> @author A. Ohrn
*>
*> @details
*> The creator of random geometries in typical Monte-Carlo fashion.
*> The changes made are:
*>
*> 1. The dielectric radius is modified
*> 2. Each coordinate *except* the \c iSta-1 -th molecule is changed
*> 3. All molecules (with the above exception) are rotated around the oxygen (which approximately equals the CM)
*> 4. Every molecule except the fixed ones are rotated slightly around one of the global \f$ x \f$-, \f$ y \f$- or \f$ z \f$-axes;
*>    the purpose of this is to emulate a rotation of the central molecule and therefore make the system more dynamic.
*>
*> @param[in,out] Ract  The dielectric radius on input and the slightly perturbed radius on output
*> @param[out]    Rold  Stores the input dielectric radius
*> @param[in]     iCNum How many solvent places that are taken up by the QM-molecule
************************************************************************
      Subroutine Geogen(Ract,Rold,iCNum,iQ_Atoms)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qmcom.fh"
#include "qminp.fh"

      Dimension Dq(3)
      External Ranf
*------------------------------------------------------------------------*
* Store old configuration.                                               *
*------------------------------------------------------------------------*
      Do 12,i=1,nPart*nCent
        Do 13,j=1,3
          OldGeo(i,j)=Cordst(i,j)
13      Continue
12    Continue
*------------------------------------------------------------------------*
* Query which type of simulation this is, and if quantum then change     *
* the quantum molecule.                                                  *
*------------------------------------------------------------------------*
      If(Qmeq.or.QmProd) then !Which coordinates to keep fixed.
        iSta=iCNum+1  !This sees to that the QM-molecule is excluded
                    !from the moves below.
        Dq(1)=delX*(ranf(iseed)-0.5)
        Dq(2)=delX*(ranf(iseed)-0.5)
        Dq(3)=delX*(ranf(iseed)-0.5)
        Do 21,iAt=1,iQ_Atoms
          Do 22,ii=1,3
            Cordst(iAt,ii)=Cordst(iAt,ii)+Dq(ii) !Move QM-mol.
22        Continue
21      Continue
      Endif
*-----------------------------------------------------------------------*
* Obtain the random-stuff and make small geometry change.               *
*-----------------------------------------------------------------------*
      Rold=Ract
      Ract=Ract+(ranf(iseed)-0.5)*DelR !Change in cavity radius
      Do 100, i=iSta,nPart !Which molecules to give new coordinates.
        ij=(i-1)*nCent
        Do 101, j=1,3
          Dx=DelX*(ranf(iseed)-0.5)
          Do 102, k=1,nCent
            ii=ij+k
            Cordst(ii,j)=Cordst(ii,j)+Dx  !Make translation
102       Continue
101     Continue
        Cx=Cordst(ij+1,1)  !The oxygen, around which we rotate
        Cy=Cordst(ij+1,2)
        Cz=Cordst(ij+1,3)
        B=(ranf(iseed)-0.5)*DelFi
        CB=Cos(B)
        SB=Sin(B)
        Do 111, k=2,nCent !Rotate around the oxygen in yz-plane, i.e.
          y=Cordst(ij+k,2)-Cy   !around x-axis.
          z=Cordst(ij+k,3)-Cz
          yNy=y*CB+z*SB  !This is a rotation matrix
          zNy=z*CB-y*SB
          Cordst(ij+k,2)=yNy+Cy
          Cordst(ij+k,3)=zNy+Cz
111     Continue
        B=(ranf(iseed)-0.5)*DelFi
        CB=Cos(B)
        SB=Sin(B)
        Do 112, k=2,nCent !And now rotate in xz-plane
          x=Cordst(ij+k,1)-Cx
          z=Cordst(ij+k,3)-Cz
          xNy=x*CB+z*SB
          zNy=z*CB-x*SB
          Cordst(ij+k,1)=xNy+Cx
          Cordst(ij+k,3)=zNy+Cz
112     Continue
        B=(ranf(iseed)-0.5)*DelFi
        CB=Cos(B)
        SB=Sin(B)
        Do 113, k=2,nCent  !To your surprise, here we rotate in the
          x=Cordst(ij+k,1)-Cx !xy-plane.
          y=Cordst(ij+k,2)-Cy
          xNy=x*CB+y*SB
          yNy=y*CB-x*SB
          Cordst(ij+k,1)=xNy+Cx
          Cordst(ij+k,2)=yNy+Cy
113     Continue
100   Continue
*Here all other water molecules rotate around one of the three axes,
*except the ones we fix, which in a qm-simulation is the
*quantum particle.
      A=ranf(iseed)
      B=(ranf(iseed)-0.5)*DelFi*0.1
      CB=Cos(B)
      SB=Sin(B)
      If(A.le.0.33333333) then  !make it random whether we rotate around
                                !x, y or z.
        Do 201, i=iSta,nPart
          ij=(i-1)*nCent
          Do 202, k=1,nCent
            ii=ij+k
            Cy=Cordst(ii,2)
            Cz=Cordst(ii,3)
            Cordst(ii,2)=Cy*CB+Cz*SB
            Cordst(ii,3)=Cz*CB-Cy*SB
202       Continue
201     Continue
      Else
        If(A.le.0.66666667) then
          Do 211, i=iSta,nPart
            ij=(i-1)*nCent
            Do 212, k=1,nCent
              ii=ij+k
              Cx=Cordst(ii,1)
              Cz=Cordst(ii,3)
              Cordst(ii,1)=Cx*CB+Cz*SB
              Cordst(ii,3)=CB*Cz-SB*Cx
212         Continue
211       Continue
        Else
          Do 221, i=iSta,nPart
            ij=(i-1)*nCent
            Do 222, k=1,nCent
              ii=ij+k
              Cx=Cordst(ii,1)
              Cy=Cordst(ii,2)
              Cordst(ii,1)=Cx*CB+Cy*SB
              Cordst(ii,2)=CB*Cy-SB*Cx
222         Continue
221       Continue
        Endif
      Endif
*---------------------------------------------------------------------*
* Generate the image points that correspond with the new coordinates. *
* We follow Friedman. Since no image is created here for the qm-      *
* molecule, start with query.                                         *
*---------------------------------------------------------------------*
      Ind=0    !The SM-defaults.
      iImage=1
      If(Qmeq.or.QmProd) then
        Ind=iCNum*nCent  !Makes sure that the first slots in CordIm
                         !are empty
        iImage=iSta
      Endif
      iQsta=nCent-nCha+1
      A2=Ract**2
      DiFac=-(DiEl-1.0)/(DiEl+1.0)
      Do 301, i=1,3
        xyzMyp(i)=0
301   Continue
      Do 302, i=iImage,nPart
        Do 303, j=1,nCent
          Ind=Ind+1
          S2=0
          Do 304, k=1,3
            S2=S2+Cordst(Ind,k)**2
304       Continue
          S2=A2/S2
          Sqrts2=Sqrt(S2)
          Sqrs(Ind)=Sqrts2
          If(j.le.nPol) then
            QImp(Ind)=0
            Do 305, k=1,3
              Dim(Ind,k)=0
305         Continue
          Endif
          If(j.ge.iQsta) then
            qq=Qsta(j-nCent+nCha)
            q=DiFac*Sqrts2*qq
            Qim(Ind)=q
          Else
            qq=0
            Qim(Ind)=0
          Endif
          Do 306, k=1,3
            xyzMyp(k)=xyzMyp(k)-qq*Cordst(Ind,k) !Total dipole of
                              !the cavity; used in polink.f.
            CordIm(Ind,k)=Cordst(Ind,k)*S2
306       Continue
303     Continue
302   Continue
*-----------------------------------------------------------------------*
* Exit.                                                                 *
*-----------------------------------------------------------------------*
      Return
      End
