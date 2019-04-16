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
*  MultiNew
*
*> @brief
*>   Perform the MME in contracted AO-basis
*> @author A. Ohrn
*>
*> @details
*> (i) Read in the multipole integrals from Seward. (ii) Construct
*> some data to simplify accessing the computed data. (iii) Make
*> the actual MME.
*>
*> @note
*> Requires numbers taken from ::qfread. We also need some integrals
*> that supposedly have been computed by Seward.
*>
*> @param[in]  nAt      Number of atoms in QM-molecule
*> @param[in]  nBas     Number of contracted basis functions
*> @param[in]  nOcc     Number of basis functions of the \f$ i \f$ -th atom-type
*> @param[in]  natyp    Number of atoms of the \f$ i \f$ -th atom-type
*> @param[in]  nntyp    Number of atom-types
*> @param[out] iMME     Pointer to the multicenter multipole expanded densities of unique pairs of contracted basis functions
*> @param[out] iCenTri  Set of indices that tells to which center the \f$ i \f$ -th unique pair of basis functions in a lower triangulary stored matrix belongs
*> @param[out] iCenTriT Just like \p iCenTri, but in square shape
*> @param[out] nMlt     Highest multipole in MME
*> @param[out] outxyz   Expansion centers in molecule
************************************************************************
      Subroutine MultiNew(nAt,nBas,nOcc,natyp,nntyp,iMME,iCenTri
     &,iCenTriT,nMlt,outxyz,SlExpQ,lSlater)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension xyz(MxAt,MxAt,3),CordMul(MxMltp,3),outxyz(MxQCen,3)
      Dimension nOcc(MxAt),natyp(MxAt),nBasAt(MxBas)
      Dimension iCenTri(*),iCenTriT(*)
      Dimension iX(6),iY(6),iMult(MxMltp,MxComp)
      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6)
      Dimension SlExpQ(MxMltp+1,MxQCen)
*Jose.No Nuclear charges in Salter ,SlPQ(MxQCen)
      Character MemLab*20,MMElab*20,ChCo*2,ChCo2*2
      Character*9 Integrals(3)
      Logical Lika, Changed1, Changed2, lSlater
      Data iX/1,1,1,2,2,3/
      Data iY/1,2,3,2,3,3/
      Data Integrals/'MLTPL  0','MLTPL  1','MLTPL  2'/
      Dimension iDum(1)
*----------------------------------------------------------------------*
* Read the multipole integrals in contracted AO-basis.                 *
*----------------------------------------------------------------------*
      irc=-1
      Lu_One=49
      Lu_One=IsFreeUnit(Lu_One)
      Call OpnOne(irc,0,'ONEINT',Lu_One)
      If(irc.ne.0) then
        Write(6,*)
        Write(6,*)'ERROR! Could not open one-electron integral file.'
        Call Quit(_RC_IO_ERROR_READ_)
      Endif

*
*-- This loop will terminate when no more multipole integrals are
*   available, hence there is not a problem that we apparently loop
*   over MxMltpl, which is a fixed number.
*
      Do 100, iMlt=1,MxMltp
        nComp=iMlt*(iMlt+1)/2
        Do 101, iComp=1,nComp
          irc=-1
          iOpt=1
          iSmLbl=1
          Call iRdOne(irc,iOpt,integrals(iMlt),iComp,iDum,iSmLbl)
          nSize=iDum(1)
          If(irc.ne.0) then
            If(iComp.ne.1) then
              Write(6,*)
              Write(6,*)'ERROR! Failed to read number of one-electron i'
     &//'ntegrals.'
              Call Quit(_RC_IO_ERROR_READ_)
            Else  !Normal exit here.
              nMlt=iMlt-1
              Go To 199
            Endif
          Endif
          If(nSize.ne.0) then
            Write(ChCo,'(I2.2)')iMlt
            Write(ChCo2,'(I2.2)')iComp
            Write(MemLab,*)'MEM'//ChCo//ChCo2
            Call GetMem(MemLab,'Allo','Real',iMult(iMlt,iComp),nSize+4)
            irc=-1
            iOpt=0
            iSmLbl=0
            Call RdOne(irc,iOpt,integrals(iMlt),iComp
     &               ,Work(iMult(iMlt,iComp)),iSmLbl) !Collect integrals
          Else
            Write(6,*)
            Write(6,*)'ERROR! Problem reading ',integrals(iMlt)
            Call Quit(_RC_IO_ERROR_READ_)
          Endif
101     Continue
        Do 110, i=1,3
          CordMul(iMlt,i)=Work(iMult(iMlt,1)+nSize+i-1)
110     Continue
        nMlt=MxMltp
100   Continue
199   Continue

*----------------------------------------------------------------------*
* Collect centers from preceeding MpProp calculation. Compute two      *
* index vectors. First one gives index of atom on which the i:th basis *
* function is centered. The other (iCenTri) gives to which center the  *
* i:th unique basis function product belong.                           *
*----------------------------------------------------------------------*
      Call Get_Centers(nAt,xyz)
      kaunt=0
      Do 2001, i=1,nAt
        kaunt=kaunt+1
        outxyz(kaunt,1)=xyz(i,i,1)
        outxyz(kaunt,2)=xyz(i,i,2)
        outxyz(kaunt,3)=xyz(i,i,3)
2001  Continue
      kaunt=nAt
      Do 2002, i=1,nAt
        Do 2003, j=1,i-1
          kaunt=kaunt+1
          outxyz(kaunt,1)=xyz(i,j,1)
          outxyz(kaunt,2)=xyz(i,j,2)
          outxyz(kaunt,3)=xyz(i,j,3)
2003    Continue
2002  Continue
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
*Jose. Collect data of the Slater representation of the Quantum System *
* Prefactors, Exponents, PointNuclearCharges.                          *
*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*
      If(lSlater) then
        Call Get_Slater(SlExpQ,LMltSlQ,outxyz,nAt)

        If(LMltSlQ+1.ne.nMlt) then
        Write(6,*)'ERROR! Multipole order',LMltSlQ,' in DiffPr file is'
     &//' different from order',nMlt-1,' in One-electron file. Check'
     &//' your files.'
          Call Quit(_RC_GENERAL_ERROR_)
        Endif
      Endif

      kaunter=0
      iAt=0
      Do 201, i=1,nntyp
        nBasA=nOcc(i)/natyp(i)
        Do 202, j=1,natyp(i)
          iAt=iAt+1
          Do 203,k=1,nBasA
            kaunter=kaunter+1
            nBasAt(kaunter)=iAt
203       Continue
202     Continue
201   Continue

      kaunter=0
      Indie=nAt
      IndiePrev=1
      nB1Prev=1
      nB2Prev=1
      Do 204, iB1=1,nBas  !Count over unique pairs of bas.func.
        Do 205, iB2=1,iB1
          kaunter=kaunter+1
          Lika=nBasAt(iB1).eq.nBasAt(iB2)
          If(Lika) then  !If equal indeces, then take that number.
            iCenTri(kaunter)=nBasAt(iB1)
            nB1Prev=nBasAt(iB1)
            nB2Prev=nBasAt(iB2)
          Else
            Changed1=nB1Prev.ne.nBasAt(iB1) !Check if changed atom.
            Changed2=nB2Prev.ne.nBasAt(iB2)
            If(Changed1.and..not.Changed2) then !Case when from center
              Indie=Indie+1                     !1 to nAt+1.
              nB1Prev=nBasAt(iB1)
              IndiePrev=Indie
            Elseif(Changed2.and..not.Changed1) then !Moving to new
              If(iB2.eq.1) then   !atom horizontally in lower triangular
                Indie=IndiePrev   !matrix. If it is a jump back to left
              Else                !corner, do not increase index, but
                Indie=Indie+1     !get old one.
              Endif
              nB2Prev=nBasAt(iB2)
            Elseif(Changed1.and.Changed2) then  !Changing both atoms.
              Indie=Indie+1
              nB1Prev=nBasAt(iB1)
              nB2Prev=nBasAt(iB2)
              IndiePrev=Indie
            Endif
            iCenTri(kaunter)=Indie
          Endif
205     Continue
204   Continue
      Ind=0
      Do 206, i=0,nBas-1  !Lets be square.
        Do 207, j=0,i
          Ind=Ind+1
          iCenTriT(1+i+j*nBas)=iCenTri(ind)
          iCenTriT(1+j+i*nBas)=iCenTri(ind)
207     Continue
206   Continue

*----------------------------------------------------------------------*
* Start the MME. To get a MME-dipole, we want <psi_i|x-x_o|psi_j> but  *
* we have <psi_i|x-x_M|psi_j> where x_o is the chosen MME-center, while*
* x_M is the center that Molcas uses. We transform in this manner:     *
* <psi_i|x-x_o|psi_j>=<psi_i|x-x_M|psi_j>+(x_M-x_o)*<psi_i|psi_j>.     *
* The quadrupole contains a further complication: not only must we     *
* include more terms, but the dipole correction may not be the MME-    *
* dipole due to that Molcas may not have used the same center for      *
* dipoles and quadrupoles. Let have a look: <psi_i|(x-x_o)(y-y_o)|psi_j>=*
* <psi_i|(x-x_M)(y-y_M)|psi_j>+(x_M-x_o)*<psi_i|y-y_M|psi_j>+...       *
* But the last dipole term may need to be transformed further if y_M   *
* for the quadrupoles are not the same as y_M for the dipoles. This is *
* the explanation for the somewhat "sliskiga" expression below for the *
* quadrupoles.                                                         *
*----------------------------------------------------------------------*
      If(nMlt.gt.3) then  !This number is connected to for how high
                          !order of multipole we have implemented below.
        Write(6,*)
        Write(6,*)'Too high order of multipole in MME.'
        Call Quit(_RC_INTERNAL_ERROR_)
      Endif
      nMul=0
      Do 210, i=1,nMlt
        nMul=nMul+i*(i+1)/2
210   Continue
      Do 2101, iMlt=1,nMul
        Write(ChCo,'(I2.2)')iMlt
        Write(MMElab,*)'MME'//ChCo
        Call GetMem(MMElab,'Allo','Real',iMME(iMlt),nSize)
2101  Continue

*
*-- The MME.
*
      kaunt=0
      Do 211, iB1=1,nBas
        Do 212, iB2=1,iB1
*
*-- The charge. No translation.
*
          Work(iMME(1)+kaunt)=Work(iMult(1,1)+kaunt)
*
*-- The dipole. Translation gives rise to charge.
*
          Do 221, i=1,3
            Corr=(CordMul(2,i)-xyz(nBasAt(iB1),nBasAt(iB2),i))
     &          *Work(iMult(1,1)+kaunt)
            Work(iMME(i+1)+kaunt)=Work(iMult(2,i)+kaunt)+Corr
221       Continue
*
*-- The quadrupole. Translation gives rise to dipoles and charges.
*   Also we have to keep track of the centers for the integrals computed
*   by Seward.
*
          Do 222, i=1,6
          CorrDip1=(CordMul(3,iX(i))-xyz(nBasAt(iB1),nBasAt(iB2),iX(i)))
     &            *(Work(iMult(2,iY(i))+kaunt)
     &            +(CordMul(2,iY(i))-CordMul(3,iY(i)))
     &            *Work(iMult(1,1)+kaunt))
          CorrDip2=(CordMul(3,iY(i))-xyz(nBasAt(iB1),nBasAt(iB2),iY(i)))
     &            *(Work(iMult(2,iX(i))+kaunt)
     &            +(CordMul(2,iX(i))-CordMul(3,iX(i)))
     &            *Work(iMult(1,1)+kaunt))
          CorrOvl=(CordMul(3,iX(i))-xyz(nBasAt(iB1),nBasAt(iB2),iX(i)))
     &           *(CordMul(3,iY(i))-xyz(nBasAt(iB1),nBasAt(iB2),iY(i)))
     &           *Work(iMult(1,1)+kaunt)
          Work(iMME(i+4)+kaunt)=Work(iMult(3,i)+kaunt)
     &                         +CorrDip1+CorrDip2+CorrOvl
222       Continue
          kaunt=kaunt+1
212     Continue
211   Continue

*
*-- Deallocations.
*
      Do 301,iMlt=1,nMlt
        nComp=iMlt*(iMlt+1)/2
        Do 302, iComp=1,nComp
          Write(ChCo,'(I2.2)')iMlt
          Write(ChCo2,'(I2.2)')iComp
          Write(MemLab,*)'MEM'//ChCo//ChCo2
          Call GetMem(MemLab,'Free','Real',iMult(iMlt,iComp),nSize+4)
302     Continue
301   Continue
      Call ClsOne(irc,Lu_One)

      Return
      End
