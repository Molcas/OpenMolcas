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
*  OverLq
*
*> @brief
*>   Compute overlap between primitve bases, which for bases of other type than s,
*>   will mean that several overlaps between basis-functions are computed. One
*>   function is on the solvent, the other in the QM-region. Observe that we do not
*>   care about overlaps within the QM-region or among the solvent molecules.
*> @author A. Ohrn
*>
*> @details
*> Uses the formulas in \cite Tak1966-JPSJ-21-2313. It is hard to give any
*> easy explanation, so if you want to understand exactly what is
*> going on below, see the article, especially equations (2.4) and
*> (2.12); then the source-code comments will provide you with
*> sufficient information.
*>
*> @param[in]  Bori   Center for the QM-region contracted basis-function
*> @param[in]  Cori   Like Bori, but for the solvent basis-function
*> @param[in]  Alfa   Exponents for the primitive basis-functions that build this contracted function
*> @param[in]  Beta   Like \p alfa, but for solvent
*> @param[in]  iQ1    = ``1`` if s-type, = ``2`` if p-type, etc. for the function in the QM-region
*> @param[in]  iQ2    Like \p iQ1, but for solvent function
*> @param[in]  nExp1  How many primitives there are in this contracted function
*> @param[in]  nExp2  Like \p nExp1, but for (surprise) the solvent
*> @param[out] iPSint Pointer to the matrix of overlaps
*> @param[in]  Trans  Transition matrix between Cartesian and spherical basis functions
************************************************************************
      Subroutine OverLq(Bori,Cori,Alfa,Beta,iQ1,iQ2,nExp1,nExp2,iPSint
     &                 ,Trans)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

* MaxAngqNr=4 means f-function is top. There is no limit in the
* algorithm though, so if higher is needed, change this number.
      Parameter(MaxAr=MxAngqNr*(MxAngqNr+1)/2)
      Dimension Bori(3),Cori(3),Alfa(MxCont),Beta(MxCont)
      Dimension Trans((3*MxAngqNr**2-2*MxAngqNr-10+8*MxAngqNr**3
     &                +3*MxAngqNr**4)/12)
      Dimension PAxyz(3),PBxyz(3),TheCent(3)
      Dimension FactorX(2*MxAngqNr+1),FactorY(2*MxAngqNr+1)
      Dimension FactorZ(2*MxAngqNr+1)
      Dimension nCartxQ(MaxAr),nCartyQ(MaxAr),nCartzQ(MaxAr)
      Dimension nCartxC(MaxAr),nCartyC(MaxAr),nCartzC(MaxAr)
*----------------------------------------------------------------------*
* Prepare some numbers for later.                                      *
*----------------------------------------------------------------------*
      Call Qenter('OverLQ')
      ind=0
      nSpecific1=iQ1*(iQ1+1)/2 !Remember that each base consist of
      nSpecific2=iQ2*(iQ2+1)/2 !many functions - this is how many.
      icompo=nSpecific1
      Do 14, ix=0,iQ1-1  !Then follows the NEW loops that compute how
        Do 15, iy=0,iQ1-1-ix !the cartesian components are ordered in
          iz=iQ1-1-iy-ix    !various basis functions.
          ncartxQ(icompo)=ix
          ncartyQ(icompo)=iy
          ncartzQ(icompo)=iz
          icompo=icompo-1
15       Continue
14     Continue
      icompo=nSpecific2
      Do 16, ix=0,iQ2-1       !And for the solvent orbitals.
        Do 17, iy=0,iQ2-1-ix
          iz=iQ2-1-iy-ix
          ncartxC(icompo)=ix
          ncartyC(icompo)=iy
          ncartzC(icompo)=iz
          icompo=icompo-1
17       Continue
16     Continue
      nSph1=2*iQ1-1
      nSph2=2*iQ2-1
      nSizeCart=nSpecific1*nSpecific2
      nSizeSph=nSph1*nSph2
      nBigP=nExp1*nExp2*nSph1*nSph2
      Call GetMem('PrimCar','Allo','Real',iPpS,nSizeCart)
      Call GetMem('PrimSph','Allo','Real',iPsphS,nSizeSph)
      Call GetMem('AllPrims','Allo','Real',iPSint,nBigP)
      Do 8, i=1,nBigP
        Work(iPSint+i-1)=0
8     Continue
      Separation=((Bori(1)-Cori(1))**2+(Bori(2)-Cori(2))**2
     &          +(Bori(3)-Cori(3))**2)
*----------------------------------------------------------------------*
* Start loop over primitives.                                          *
*----------------------------------------------------------------------*
      Kaunt=0
      Do 101, iP1=1,nExp1
        Do 102, iP2=1,nExp2
          Kaunt=Kaunt+1
          TheCent(1)=(Alfa(iP1)*Bori(1)+Beta(iP2)*Cori(1))
          TheCent(2)=(Alfa(iP1)*Bori(2)+Beta(iP2)*Cori(2))
          TheCent(3)=(Alfa(iP1)*Bori(3)+Beta(iP2)*Cori(3))
          Divide=1/(Alfa(iP1)+Beta(iP2))  !gamma in article
          TheCent(1)=TheCent(1)*Divide  !The new center, P
          TheCent(2)=TheCent(2)*Divide
          TheCent(3)=TheCent(3)*Divide
          Piconst=Pi*Divide
          SqPiconst=sqrt(Piconst)
          Piconst=Piconst*SqPiconst  !That constant to the power of 3/2
          Expo=Alfa(iP1)*Beta(iP2)*Separation*Divide
          TheFirstFac=Piconst*exp(-Expo) !This is the exponential factor
*Now we should get those difficult f-functions.
          PAxyz(1)=TheCent(1)-Bori(1)
          PAxyz(2)=TheCent(2)-Bori(2)
          PAxyz(3)=TheCent(3)-Bori(3)
          PBxyz(1)=TheCent(1)-Cori(1)
          PBxyz(2)=TheCent(2)-Cori(2)
          PBxyz(3)=TheCent(3)-Cori(3)
          kaunter=0
          Do 103, iSp1=1,nSpecific1
            Do 104, iSp2=1,nSpecific2
              kaunter=kaunter+1
              loneX=ncartxQ(iSp1)
              ltwoX=ncartxC(iSp2)
              lsumX=loneX+ltwoX
              loneY=ncartyQ(iSp1)
              ltwoY=ncartyC(iSp2)
              lsumY=loneY+ltwoY
              loneZ=ncartzQ(iSp1)
              ltwoZ=ncartzC(iSp2)
              lsumZ=loneZ+ltwoZ
              Call fFactor(loneX,ltwoX,lsumX,loneY,ltwoY,lsumY,loneZ
     &                    ,ltwoZ,lsumZ,PAxyz,PBxyz,FactorX,FactorY
     &                    ,FactorZ)
*Now we have the f-factors for this specific angular type of this
*specific primitive basis-function. Now put things together.
              iUpX=lsumX/2 !Yes, it should be like this, even when lsumX
              iUpY=lsumY/2 !is odd.
              iUpZ=lsumZ/2
              SummaX=0
              SummaY=0
              SummaZ=0
              Do 131, ixxx=0,iUpX
                Extra=iDubFac(2*ixxx-1)*(0.5*Divide)**ixxx !This is
                SummaX=SummaX+FactorX(2*ixxx+1)*Extra  !just a matter of
131           Continue                               !putting things
              Do 132, iyyy=0,iUpY             !together according to the
                Extra=iDubFac(2*iyyy-1)*(0.5*Divide)**iyyy  !formula
                SummaY=SummaY+FactorY(2*iyyy+1)*Extra
132           Continue
              Do 133, izzz=0,iUpZ
                Extra=iDubFac(2*izzz-1)*(0.5*Divide)**izzz
                SummaZ=SummaZ+FactorZ(2*izzz+1)*Extra
133           Continue
              Primequals=TheFirstFac*SummaX*SummaY*SummaZ
              Work(iPpS+kaunter-1)=Primequals
104         Continue
103       Continue
*----------------------------------------------------------------------*
* This was the overlap for the primitives in terms of cartesian        *
* functions, but in the new qmstat we use spherical functions, so we   *
* need to transform if any d-function or higher is involved. In the    *
* matrix Trans the numbers for how spherical functions are expressed in*
* cartesian functions are stored, including the extra normalization so *
* that all d-functions (and higher) have the same combined contraction *
* and normalization coefficient. The rest is just a matter of getting  *
* the matrix multiplications right. The convention I use is this: The  *
* matrix with the overlaps contain elements such as <psi_QM|psi_Solv>  *
* in other words, the QM-orbitals count over the rows and the solvent  *
* orbitals over the columns; observe however that this is NOT the way  *
* the matrix enters from above, since there the fastest couting index  *
* is over solvent orbitals (iSp2), so given this and the knowledge of  *
* how Fortran stores multidimensional matrices, we can figure out when *
* to transpose. All this means that if it is the QM-orbitals that are  *
* to be transformed, the transformation matrix is multiplied from left,*
* while it is multiplied from the right -- transposed of course -- if  *
* it is the solvent orbitals that are to be transformed. In the case   *
* that both orbitals are to be transformed we simply apply the trans-  *
* formation matrix from both directions.                               *
*----------------------------------------------------------------------*
          If(iQ1.ge.3.or.iQ2.ge.3) then !Check if any transformations
                                        !are necessary.
            If(iQ2.lt.3) then  !If only the base of the QM-region needs
                               !to be transformed.
              ind=1+(iQ1-3)*(3*iQ1**3+5*iQ1**2+12*iQ1+40)/12
              Call Dgemm_('N','T',nSph1,nSpecific2,nSpecific1,ONE
     &                  ,Trans(ind),nSph1,Work(iPps),nSpecific2,ZERO
     &                  ,Work(iPsphS),nSph1)
            Elseif(iQ1.lt.3) then !If only solvent base needs to be
                                  !transformed.
              ind=1+(iQ2-3)*(3*iQ2**3+5*iQ2**2+12*iQ2+40)/12
              Call Dgemm_('T','T',nSph1,nSph2,nSpecific2,ONE
     &                  ,Work(iPps),nSpecific2,Trans(ind),nSph2,ZERO
     &                  ,Work(iPsphS),nSph1)
            Else !Both QM-region and Solvent need to be transformed.
              ind1=1+(iQ1-3)*(3*iQ1**3+5*iQ1**2+12*iQ1+40)/12
              ind2=1+(iQ2-3)*(3*iQ2**3+5*iQ2**2+12*iQ2+40)/12
              Call GetMem('Intmd','Allo','Real',iPInte,nSph1*nSpecific2)
              Call Dgemm_('N','T',nSph1,nSpecific2,nSpecific1,ONE
     &                  ,Trans(ind1),nSph1,Work(iPps),nSpecific2,ZERO
     &                  ,Work(iPInte),nSph1)
              Call Dgemm_('N','T',nSph1,nSph2,nSpecific2,ONE
     &                  ,Work(iPInte),nSph1,Trans(ind2),nSph2,ZERO
     &                  ,Work(iPsphS),nSph1)
              Call GetMem('Intmd','Free','Real',iPInte,nSph1*nSpecific2)
            Endif
          Else  !Here we only transpose to get integrals in right
            krakna=0  !order.
            Do 136, i=0,nSph1-1
              Do 137, j=0,nSph2-1
                Work(iPsphS+krakna)=Work(iPps+i+j*nSph1)
                krakna=krakna+1
137           Continue
136         Continue
          Endif
*----------------------------------------------------------------------*
* Put this thing in the slowly growing overlap matrix for the          *
* primitive basis functions. The reason the index is so nasty is that  *
* we compute small blocks of the matrix and now have to fit it in the  *
* right place in the growing, much larger, matrix. Nasty!              *
*----------------------------------------------------------------------*
          krakna=0
          Do 141, j=0,nSph2-1
            Do 142, i=0,nSph1-1
              jndex=i+j*nExp1*nSph1+(iP1-1)*nSph1
     &              +(iP2-1)*nExp1*nSph1*nSph2
              Work(iPSint+jndex)=Work(iPsphS+krakna)
              krakna=krakna+1
142         Continue
141       Continue
102     Continue
101   Continue
*----------------------------------------------------------------------*
* Deallocate and ta'ta!                                                *
*----------------------------------------------------------------------*
      Call GetMem('PrimCar','Free','Real',iPpS,nSizeCart)
      Call GetMem('PrimSph','Free','Real',iPsphS,nSizeSph)
      Call Qexit('OverLQ')

      Return
      End

*----------------------------------------------------------------------*
* A function that returns the binomial coefficient. The coefficients   *
* are stored since N and P will not under normal circumstances be      *
* so large.                                                            *
*----------------------------------------------------------------------*
      Integer Function NoverP_Q(N,P)
      Integer N,P,Bino(22)
      Data (Bino(i),i=1,21)/1,1,1,1,2,1,1,3,3,1
     &        ,1,4,6,4,1,1,5,10,10,5,1/
      NoverP_Q=1
      If(N.ge.6) then
        Write(6,*)'Must extend NoverP_Q!'
        Call Quit(_RC_INTERNAL_ERROR_)
      Else
        ind=(N+1)*(N+2)/2-(N-P)
        NoverP_Q=Bino(ind)
      Endif
      Return
      End
*----------------------------------------------------------------------*
* A function that will return the double factorial. We do not expect   *
* big numbers, so we do it brute-force. Observe that N must be odd, but*
* to skip the if-sentence, we assume that the one who calls this       *
* function has seen to that.                                           *
*----------------------------------------------------------------------*
      Integer Function iDubFac(N)
      Integer N
      iDubFac=1
      Do 1101, k=3,N,2
        iDubFac=iDubFac*k
1101  Continue
      Return
      End
*----------------------------------------------------------------------*
* A subroutine that computes those darn f-factors. They are definied   *
* in equation (2.4) in the article cited above. As can be seen from    *
* that equation, the computation of the f-factors is actually a matter *
* of using the binomial theorem. This is what we do below and to make  *
* the computation efficient the expression (2.4) is written as a       *
* succint double sum.                                                  *
*----------------------------------------------------------------------*
      Subroutine fFactor(loneX,ltwoX,lsumX,loneY,ltwoY,lsumY,loneZ
     &                 ,ltwoZ,lsumZ,PAxyz,PBxyz,FactorX,FactorY,FactorZ)
      Implicit Real*8 (a-h,o-z)

      Parameter(MaxAngqNr=6,MaxAr=MaxAngqNr*(MaxAngqNr+1)/2)
      Dimension FactorX(2*MaxAngqNr+1),FactorY(2*MaxAngqNr+1)
      Dimension FactorZ(2*MaxAngqNr+1)
      Dimension PAxyz(3),PBxyz(3)

      Do 105, ia=0,lsumX !We use unrolled loops with regard to x,y and z
        fff2=0          !therefore, here we start with the x-factors.
        iLowB=max(0,ia-ltwoX)  !These lower and upper bounds have to do
        iUpB=min(ia,loneX)   !with the allowed numbers in the binomial
        Do 106, i=iLowB,iUpB  !coefficients.
          fff1=NoverP_Q(loneX,i)*NoverP_Q(ltwoX,ia-i)
          If(i.ne.0) then  !This is needed for some compilers (NAG_64)
            PAraise=PAxyz(1)**i
          Else
            PAraise=1.0d0
          Endif
          If(ia-i.ne.0) then
            PBraise=PBxyz(1)**(ia-i)
          Else
            PBraise=1.0d0
          Endif
          fff2=fff2+fff1*PAraise*PBraise
106     Continue
        FactorX(lsumX-ia+1)=fff2
105   Continue
      Do 115, ia=0,lsumY  !y-factors.
        fff2=0
        iLowB=max(0,ia-ltwoY)
        iUpB=min(ia,loneY)
        Do 116, i=iLowB,iUpB
          fff1=NoverP_Q(loneY,i)*NoverP_Q(ltwoY,ia-i)
          If(i.ne.0) then  !This is needed for some compilers (NAG_64)
            PAraise=PAxyz(2)**i
          Else
            PAraise=1.0d0
          Endif
          If(ia-i.ne.0) then
            PBraise=PBxyz(2)**(ia-i)
          Else
            PBraise=1.0
          Endif
          fff2=fff2+fff1*PAraise*PBraise
116     Continue
        FactorY(lsumY-ia+1)=fff2
115   Continue
      Do 125, ia=0,lsumZ  !z-factorz.
        fff2=0
        iLowB=max(0,ia-ltwoZ)
        iUpB=min(ia,loneZ)
        Do 126, i=iLowB,iUpB
          fff1=NoverP_Q(loneZ,i)*NoverP_Q(ltwoZ,ia-i)
          If(i.ne.0) then  !This is needed for some compilers (NAG_64)
            PAraise=PAxyz(3)**i
          Else
            PAraise=1.0d0
          Endif
          If(ia-i.ne.0) then
            PBraise=PBxyz(3)**(ia-i)
          Else
            PBraise=1.0d0
          Endif
          fff2=fff2+fff1*PAraise*PBraise
126     Continue
        Factorz(lsumZ-ia+1)=fff2
125   Continue
      Return
      End
