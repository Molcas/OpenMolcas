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
*  ContractOvl
*
*> @brief
*>   Compute the overlaps between solvent and solute in contracted basis-functions
*> @author A. Ohrn
*>
*> @details
*> Here the overlap between the QM-region contracted AO-basis
*> functions and the present solvent molecule contracted AO-basis
*> functions are computed. In order to use the fact that we use
*> contracted functions to the maximum, we compute the overlaps with
*> primitive functions only once, then we transform this matrix to
*> all relevant contracted overlaps. After that, the old primitive
*> integrals are discarded and a new set of primitive are computed.
*> This is very nice since ::OverLq is rather slow. The problems we
*> get are that we must use rather elaborate schemes to get right
*> digit in right place.
*>
*> @param[out] Sint     The contracted basis function overlaps
*> @param[out] SintPar  The contracted basis function overlaps with extra atom--atom weights *if* this has been requested by user, otherwise unchanged
*> @param[in]  nBaseQ   Number of AO-basis functions in QM-region
*> @param[in]  nBaseC   Like \p nBaseQ but for solvent
*> @param[in]  N        Which solvent molecule this is
*> @param[in]  nCent    How many centers the solvent molecule has
*> @param[in]  iEl      Number of elements in QM-region
*> @param[in]  nAtomsCC How many solvent atoms
*> @param[in]  iPrint   Print level
************************************************************************
      Subroutine ContractOvl(Sint,SintPar,nBaseQ,nBaseC
     &,N,nCent,iEl,iQ_Atoms,nAtomsCC,iPrint,Inside)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "integral.fh"
#include "WrkSpc.fh"

      Parameter(MxSphAng=2*MxAngqNr-1)
      Dimension Sint(MxBas,MxBasC),SintPar(MxBas,MxBasC)
      Dimension Bori(3),Cori(3),Alf(MxCont),Bet(MxCont)
      Dimension Conkort(MxCont),Donkort(MxCont),ContrI(MxSphAng**2)
      Dimension Inside(MxAt,3)
      Logical Inside

      nSph1=0
      nSph2=0
      iQcontBSAV=0
      iCcontBSAV=0
      iQcontBSAV1=0
      iCcontBSAV1=0
      iQcontBSAV2=0
      iCcontBSAV2=0
      Do 11, iA1=1,iQ_Atoms  !The atoms
        Do 12, iA2=1,nAtomsCC
          If(.not.Inside(iA1,iA2)) then  !when atom-pair too far from
                                       !each other, do this then skip.
            iCcontBSAV=iCcontBSAV+nBonA_C(iA2)
            iCcontBSAV1=iCcontBSAV
            iCcontBSAV2=iCcontBSAV
            If(iA2.eq.nAtomsCC) then
              iQcontBSAV=iQcontBSAV+nBonA_Q(iA1)
              iQcontBSAV1=iQcontBSAV
            Endif
            Go To 12
          Endif
          Do 13, iB1=1,nBA_Q(iA1) !The basis functions on this specific
            Do 14, iB2=1,nBA_C(iA2) !atom
              iQcontB=iQcontBSAV
              Do 15, iNcB1=1,nCBoA_Q(iA1,iB1) !The basis of angular
                iQcontB=iQcontB+1             !type.
                iCcontB=iCcontBSAV
                Bori(1)=BasOri(1,iQcontB) !Suck-out proper coord
                Bori(2)=BasOri(2,iQcontB) !for QM.
                Bori(3)=BasOri(3,iQcontB)
                iqqqQ=iQang(iQcontB)  !Various integers, see qfread
                nExp1=nPrimus(iQcontB) !to understand their meaning.
                nSph1=2*iqqqQ-1
                Do 5411, i=1,nPrimus(iQcontB) !Suck-out the proper
                  Alf(i)=alfa(iQcontB,i)  !exponents for QM-region
                  Conkort(i)=cont(iQcontB,i)
5411            Continue
                Do 16, iNcB2=1,nCBoA_C(iA2,iB2)
                  iCcontB=iCcontB+1
                  Cori(1)=CasOri(1,iCcontB)  !Coord. of the atoms of
                  Cori(2)=CasOri(2,iCcontB)  !this solvent mol.
                  Cori(3)=CasOri(3,iCcontB)
                  iqqqC=iQn(iCcontB)
                  nExp2=mPrimus(iCcontB)
                  nSph2=2*iqqqC-1
                  Do 5412, j=1,mPrimus(iCcontB) !Exponents and stuff.
                    Bet(j)=beta(iCcontB,j)
                    Donkort(j)=dont(iCcontB,j)
5412              Continue
*------ Now call on the routine that computes a block of primitive
*       integrals. So if we are integrating the np-mp overlap we
*       compute ALL primitive p-p integrals, in the first call, then
*       they are merely contracted. This is an economical procedure
*       for both general and ordinary contracted basis sets since all
*       primitve overlaps are needed at some point in the contracted
*       overlaps, the difference between general and ordinary is that
*       in the former primitve overlaps are needed at all instances,
*       while in the latter primitve overlaps are needed only once.
                  If(iNcB1.eq.1.and.iNcB2.eq.1) then
                    Call OverLq(Bori,Cori,Alf,Bet,iqqqQ,iqqqC,nExp1
     &                         ,nExp2,iPSint,Trans)
                  Endif
                  kaunter=0
                  Do 501, i=1,nSph2 !contract
                    Do 502, j=1,nSph1
                      kaunter=kaunter+1
                      DaNumber=0
                      Do 503, iCC=1,nExp2
                        Do 504, iCQ=1,nExp1
                          iindex=(i-1)*nSph1*nExp1+j-1+nSph1*(iCQ-1)
     &                          +nSph1*nExp1*nSph2*(iCC-1)
                          DaNumber=DaNumber+Conkort(iCQ)*Donkort(iCC)
     &                            *Work(iPSint+iindex)
504                     Continue
503                   Continue
                      ContrI(kaunter)=DaNumber
502                 Continue
501               Continue
                  If(iPrint.ge.30) then
                    Write(6,*)'Basis',iQcontB,iCcontB
                    Write(6,*)'Coord.',Bori(1),Bori(2),Bori(3)
                    Write(6,*)'Coord.',Cori(1),Cori(2),Cori(3)
                    Write(6,*)'Alfa',(Alf(i),i=1,nPrimus(iQcontB))
                    Write(6,*)'Beta',(Bet(i),i=1,mPrimus(iCcontB))
                    Write(6,*)'ConQ',(Conkort(i),i=1,nPrimus(iQcontB))
                    Write(6,*)'ConC',(Donkort(i),i=1,mPrimus(iCcontB))
                    Write(6,*)'Angular',iqqqQ,iqqqC
                    Write(6,*)'#primitive',nExp1,nExp2
                    Write(6,*)(ContrI(k),k=1,nSph1*nSph2)
                  Endif
                  kreichner=0
                  Do 5421, iC=1,nSph2
                    Do 5422, iQ=1,nSph1
                    kreichner=kreichner+1
                    Sint(iWoGehenQ(iQcontB,iQ),iWoGehenC(iCcontB,iC))=
     &                                ContrI(kreichner)
5422                Continue
5421              Continue
16              Continue
15            Continue
              iCcontBSAV=iCcontB
         !This vector is allocated in OverLq.
              nSize=nExp1*nExp2*nSph1*nSph2
              Call GetMem('AllPrims','Free','Real',iPSint,nSize)
14          Continue
            iQcontBSAV=iQcontB   !OH NO!, these things have to do with
            iCcontBSAV1=iCcontBSAV  !getting the right number in right
            iCcontBSAV=iCcontBSAV2  !place. This should probably be
13        Continue                !changed sometime to a less cumbersome
          iQcontBSAV1=iQcontBSAV  !method.
          iQcontBSAV=iQcontBSAV2
          iCcontBSAV2=iCcontBSAV1
          iCcontBSAV=iCcontBSAV1
12      Continue
        iQcontBSAV2=iQcontBSAV1
        iQcontBSAV=iQcontBSAV1
        iCcontBSAV=0
        iCcontBSAV1=0
        iCcontBSAV2=0
11    Continue
      If(iPrint.ge.30) then  !Optional print-out.
        Write(6,*)
        Write(6,*)'OVERLAP BETWEEN QM-SYSTEM AND SOLVENT MOLECULE'
     &,N/nCent
        Write(6,*)'QM-AO  SOLV-AO  OVERLAP'
        Do 545, i=1,nBaseQ
          Do 546, j=1,nBaseC
            Write(6,8888)i,j,Sint(i,j)
546       Continue
545     Continue
      Endif
8888  Format(I3,'    ',I3,'       ',F12.10)

      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(SintPar)
        Call Unused_integer(iEl)
      End If
      End
