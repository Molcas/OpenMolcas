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
*  RassiHandM
*
*> @brief
*>   Make a multicenter multipole expansion of the various densities in the RASSI-state Hamiltonian to be.
*>   We also construct the gas-phase RASSI-state Hamiltonian
*> @author A. Ohrn
*>
*> @details
*> First construct the unperturbed RASSI-state Hamiltonian. The \c RasEne
*> are given in input (could be changed later). Then we obtain the
*> MME of the densities of each unique pair of AO-basis functions,
*> which we with the transition density matrices transform to their
*> RASSI-state counterparts---a process that requires some knowledge
*> about that matrix. For example, it should be noted that the matrix
*> is triangularily stored *with* corrections made for the difference
*> between diagonal and non-diagonal elements; therefore we do not
*> need to treat them differently, like we have to in the subroutine
*> ::scfhandm. If requested, we compute total charges and dipoles of
*> every state and print.
*>
*> @note
*> Requires Qfread and of course a RASSI-computation and also MPPROP.
*>
*> @param[in] nBas  Number of contracted AO-basis functions
*> @param[in] nOcc  Number of basis functions on the \f$ i \f$ th atom-type.
*> @param[in] natyp Number of atoms of the \f$ i \f$ th atom-type
*> @param[in] nntyp Number of atom-types
************************************************************************
      Subroutine RassiHandM(nBas,iQ_Atoms,nOcc,natyp,nntyp)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "files_qmstat.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension nBas(MxSym),natyp(MxAt),nOcc(MxBas)
      Dimension iMME(MxMltp*(MxMltp+1)*(MxMltp+2)/6),iCent(MxBas*MxBas)
      Character MMElab*20,ChCo*2

*
*-- A modest entrance.
*

*
*-- Zeros.
*
      kaunt=0
      iCi=iQ_Atoms*(iQ_Atoms+1)/2
      Do 1, i=1,nState
        Do 2, j=1,i
          kaunt=kaunt+1
          Do 3, k=1,iCi
            RasCha(kaunt,k)=0
            RasDip(kaunt,1,k)=0
            RasDip(kaunt,2,k)=0
            RasDip(kaunt,3,k)=0
            RasQua(kaunt,1,k)=0
            RasQua(kaunt,2,k)=0
            RasQua(kaunt,3,k)=0
            RasQua(kaunt,4,k)=0
            RasQua(kaunt,5,k)=0
            RasQua(kaunt,6,k)=0
3         Continue
2       Continue
1     Continue

*
*-- Construct H_0 with external perturbation if requested. Construct
*   a copy also.
*
      Write(6,*)
      If(.not.AddExt) then
        Write(6,*)'     Constructs H_0.'
      Else
        Write(6,*)'     Constructs H_0 with external perturbation.'
      Endif
      Call Chk_OneHam(nBas)
      Call RasH0(nBas(1))
      kaunter=0
      Do 11, i=1,nState
        Do 12, j=1,i
          kaunter=kaunter+1
          HMatSOld(kaunter)=HmatState(kaunter)
12      Continue
11    Continue
      Write(6,*)'     ...Done!'
      Write(6,*)

*
*-- Here the MME is performed. First the expansion is performed in the
*   AO-basis. Then depending on whether we use a reduced MO-basis or
*   not, the proper approach is chosen.
*
*
*-- Say what we do.
*
      Write(6,*)
      Write(6,*)'     Multicenter multipole expanding the charge'
     &//' density expressed in RASSI eigenstates.'

*
*-- Frist obtain MME in AO-basis.
*
      Call GetMem('Dummy','Allo','Inte',iDum,nBas(1)**2)
      Call MultiNew(iQ_Atoms,nBas(1),nOcc,natyp,nntyp,iMME,iCent
     &             ,iWork(iDum),nMlt,outxyzRAS
     &             ,SlExpQ,lSlater)
      Call GetMem('Dummy','Free','Inte',iDum,nBas(1)**2)

*
*-- Set nTyp, which is number of unique multipole components.
*
      nTyp=0
      Do 100, i=1,nMlt
        nTyp=nTyp+i*(i+1)/2
100   Continue

*
*-- Transform to State-basis. The logical flag MoAveRed decides which
*   path to go in this subroutine.
*
      Call StateMME(MoAveRed,nBas(1),nRedMO,nState,nTyp,iCi,iBigT,iMME
     &             ,iCent,ipAvRed,RasCha,RasDip,RasQua)

*
*-- Deallocate the MME in AO-basis.
*
      Do 106,i=1,nTyp
        Write(ChCo,'(I2.2)')i
        Write(MMElab,*)'MME'//ChCo
        Call GetMem(MMElab,'Free','Real',iMME(i),nSize)
106   Continue

*
*-- Buckinghamification of the quadrupoles.
*
        kaunter=0
        Do 991, i=1,nState
          Do 992, j=1,i
            kaunter=kaunter+1
            Do 993, k=1,ici
              Do 994, l=1,6
                RasQua(kaunter,l,k)=RasQua(kaunter,l,k)*1.5
994           Continue
              Tra=RasQua(kaunter,1,k)+RasQua(kaunter,3,k)
     &           +RasQua(kaunter,6,k)
              Tra=Tra/3
              RasQua(kaunter,1,k)=RasQua(kaunter,1,k)-Tra
              RasQua(kaunter,3,k)=RasQua(kaunter,3,k)-Tra
              RasQua(kaunter,6,k)=RasQua(kaunter,6,k)-Tra
993         Continue
992       Continue
991     Continue
*
*--- And do some printing if asked for.
*
      If(iPrint.ge.10) then
        Write(6,*)
        Write(6,*)'    Distributed multipoles for each state'
        Do 31, i=1,nState
          k=i*(i+1)/2
          Write(6,*)'     State ',i
          Do 32, j=1,iCi
            Write(6,*)'        Center ',j
            Q=-RasCha(k,j)
            If(j.le.iQ_Atoms) Q=Q+ChaNuc(j)
            Write(6,*)'          ',Q
            D1=-RasDip(k,1,j)
            D2=-RasDip(k,2,j)
            D3=-RasDip(k,3,j)
            Write(6,*)'          ',D1,D2,D3
32        Continue
31      Continue
      Endif
      If(iPrint.ge.5) then
        Write(6,*)
        Write(6,*)'    Summed multipoles for each state (not state-over'
     &//'laps)'
        Write(6,*)'              Charge  Dipole(x)   Dipole(y)   Dipole'
     &//'(z)   Quadrup(xx) Quadrup(xy) Quadrup(xz) Quadrup(yy) Quadrup('
     &//'yz) Quadrup(zz)'
      Endif
      Do 26, i=1,nState  !Total charge and dipole.
        k=i*(i+1)/2
        qEl=0
        dipx=0
        dipy=0
        dipz=0
        dipx0=0
        dipy0=0
        dipz0=0
        quaxx=0
        quaxy=0
        quaxz=0
        quayy=0
        quayz=0
        quazz=0
        quaDxx=0
        quaDxy=0
        quaDxz=0
        quaDyx=0
        quaDyy=0
        quaDyz=0
        quaDzx=0
        quaDzy=0
        quaDzz=0
        quaQxx=0
        quaQxy=0
        quaQxz=0
        quaQyy=0
        quaQyz=0
        quaQzz=0
        qtot=0
        dTox=0
        dToy=0
        dToz=0
        dQxx=0
        dQxy=0
        dQxz=0
        dQyy=0
        dQyz=0
        dQzz=0
        Trace1=0
        Trace2=0
        Trace3=0
        Do 27, j=1,iCi
          qEl=qEl+RasCha(k,j)
          dipx=dipx+RasDip(k,1,j)
          dipy=dipy+RasDip(k,2,j)
          dipz=dipz+RasDip(k,3,j)
          dipx0=dipx0+RasCha(k,j)*outxyzRAS(j,1)
          dipy0=dipy0+RasCha(k,j)*outxyzRAS(j,2)
          dipz0=dipz0+RasCha(k,j)*outxyzRAS(j,3)
          quaxx=quaxx+RasQua(k,1,j)
          quaxy=quaxy+RasQua(k,2,j)
          quaxz=quaxz+RasQua(k,4,j)
          quayy=quayy+RasQua(k,3,j)
          quayz=quayz+RasQua(k,5,j)
          quazz=quazz+RasQua(k,6,j)
          quaDxx=quaDxx+RasDip(k,1,j)*(outxyzRAS(j,1)-CT(1))
          quaDxy=quaDxy+RasDip(k,1,j)*(outxyzRAS(j,2)-CT(2))
          quaDxz=quaDxz+RasDip(k,1,j)*(outxyzRAS(j,3)-CT(3))
          quaDyx=quaDyx+RasDip(k,2,j)*(outxyzRAS(j,1)-CT(1))
          quaDyy=quaDyy+RasDip(k,2,j)*(outxyzRAS(j,2)-CT(2))
          quaDyz=quaDyz+RasDip(k,2,j)*(outxyzRAS(j,3)-CT(3))
          quaDzx=quaDzx+RasDip(k,3,j)*(outxyzRAS(j,1)-CT(1))
          quaDzy=quaDzy+RasDip(k,3,j)*(outxyzRAS(j,2)-CT(2))
          quaDzz=quaDzz+RasDip(k,3,j)*(outxyzRAS(j,3)-CT(3))
          quaQxx=quaQxx+RasCha(k,j)*(outxyzRAS(j,1)-CT(1))
     &                             *(outxyzRAS(j,1)-CT(1))
          quaQxy=quaQxy+RasCha(k,j)*(outxyzRAS(j,1)-CT(1))
     &                             *(outxyzRAS(j,2)-CT(2))
          quaQxz=quaQxz+RasCha(k,j)*(outxyzRAS(j,1)-CT(1))
     &                             *(outxyzRAS(j,3)-CT(3))
          quaQyy=quaQyy+RasCha(k,j)*(outxyzRAS(j,2)-CT(2))
     &                             *(outxyzRAS(j,2)-CT(2))
          quaQyz=quaQyz+RasCha(k,j)*(outxyzRAS(j,2)-CT(2))
     &                             *(outxyzRAS(j,3)-CT(3))
          quaQzz=quaQzz+RasCha(k,j)*(outxyzRAS(j,3)-CT(3))
     &                             *(outxyzRAS(j,3)-CT(3))
27      Continue
        Do 28, kk=1,iQ_Atoms
          qtot=qtot+ChaNuc(kk)
          dTox=dTox+ChaNuc(kk)*outxyzRAS(kk,1)
          dToy=dToy+ChaNuc(kk)*outxyzRAS(kk,2)
          dToz=dToz+ChaNuc(kk)*outxyzRAS(kk,3)
          dQxx=dQxx+ChaNuc(kk)*(outxyzRAS(kk,1)-CT(1))
     &                        *(outxyzRAS(kk,1)-CT(1))
          dQxy=dQxy+ChaNuc(kk)*(outxyzRAS(kk,1)-CT(1))
     &                        *(outxyzRAS(kk,2)-CT(2))
          dQxz=dQxz+ChaNuc(kk)*(outxyzRAS(kk,1)-CT(1))
     &                        *(outxyzRAS(kk,3)-CT(3))
          dQyy=dQyy+ChaNuc(kk)*(outxyzRAS(kk,2)-CT(2))
     &                        *(outxyzRAS(kk,2)-CT(2))
          dQyz=dQyz+ChaNuc(kk)*(outxyzRAS(kk,2)-CT(2))
     &                        *(outxyzRAS(kk,3)-CT(3))
          dQzz=dQzz+ChaNuc(kk)*(outxyzRAS(kk,3)-CT(3))
     &                        *(outxyzRAS(kk,3)-CT(3))
28      Continue
*--- Observe! qTot is not just a check. It is used later as a sign
*             to see if QM-region is charged.
        qtot=qtot-qEl
        dTox=dTox-dipx-dipx0
        dToy=dToy-dipy-dipy0
        dToz=dToz-dipz-dipz0
        Trace1=dQxx+dQyy+dQzz
        Trace2=-quaDxx-quaDyy-quaDzz
        Trace3=-quaQxx-quaQyy-quaQzz
        Trace1=Trace1/3
        Trace2=2*Trace2/3
        Trace3=Trace3/3
        dQxx=1.5*(dQxx-2*quaDxx-quaQxx-Trace1-Trace2-Trace3)-quaxx
        dQyy=1.5*(dQyy-2*quaDyy-quaQyy-Trace1-Trace2-Trace3)-quayy
        dQzz=1.5*(dQzz-2*quaDzz-quaQzz-Trace1-Trace2-Trace3)-quazz
        dQxy=1.5*(dQxy-quaDxy-quaDyx-quaQxy)-quaxy
        dQxz=1.5*(dQxz-quaDxz-quaDzx-quaQxz)-quaxz
        dQyz=1.5*(dQyz-quaDyz-quaDzy-quaQyz)-quayz
        If(iPrint.ge.5) then
          Write(6,9001)'      State ',i,'  ',qtot,dTox,dToy,dToz,dQxx
     &                                      ,dQxy,dQxz,dQyy,dQyz,dQzz
        Endif
        If(.not.abs(qtot).le.0.0001) ChargedQM=.true.
26    Continue
*Jose This format has problems to print anions
*9001  Format(A,i2,A,F5.3,9(F12.8))
9001  Format(A,i2,A,F5.1,9(F12.8))
      Write(6,*)'     ...Done!'

*----------------------------------------------------------------------*
* The end has come.                                                    *
*----------------------------------------------------------------------*
      Return
      End
