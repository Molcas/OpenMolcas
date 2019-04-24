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
      Subroutine Qmstat(ireturn)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension Coord(MxAt*3)
      Dimension nBas(1),nBas_C(1),nCnC_C(MxBasC)
*******JoseMEP New variables to perform the MEP calculation
      Dimension Eint(MxQCen,10),Poli(MxQCen,10)
      Dimension SumElcPot(MxQCen,10)
      Dimension PertElcInt(MxBas*(MxBas+1)/2)
      Dimension nOcc(MxBas),natyp(MxAt)
      Character Labjhr*4
****************
      Logical  Is_Real_Par
      External Is_Real_Par

*
*-- The journey begins. Set non-zero return code.
*
      Call Qenter('QMSTAT')
      ireturn=99

*
*-- If parallel compilation, terminate, but gracefully.
*
#ifdef _MOLCAS_MPP_
      If (Is_Real_Par()) Then
        Write(6,*)
        Write(6,*)' QMStat does not run in parallel!'
        Write(6,*)
        Call Quit(_RC_NOT_AVAILABLE_)
      End If
#endif

*
*-- Set defaults, zeros and initial values.
*
      Call Qmstat_init

*
*-- Make discrete banner.
*

*
*-- Read ('process') input. To do that we need the number of atoms in
*   the QM-region, therefore that initial call to the RUNFILE.
*
      Call Get_iScalar('Unique atoms',iQ_Atoms)
      Call Get_Qmstat_Input(iQ_Atoms)

*
*-- If only the startfile is to be edited or the sampfile analyzed,
*   go here, then terminate.
*
      If(EdSt) then
        Call EditStart
        Go To 666
      Endif
      If(Anal) then
        Call Analyze_Q(iQ_Atoms)
        Go To 666
      Endif

*
*-- Read in orbitals, basis functions, integrals and bla bla bla. This
*   is the centre for communicating with the rest of Molcas.
*
*******JoseMEP*** Qfread is called with more variables to the MEP calculation
      Call Qfread(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C
     &,nOcc,natyp,nntyp)
********

*
*-- The turning point for the single-point calculation.
*
      nCalls=0
3333  Continue

*
*-- If user request a set of single point calcualtions, then go in to
*   a separate routine and do a 'reintrepreation' of the input.
*
      If(SingPoint) Call SingP(nCalls,iQ_Atoms,ipStoreCoo,nPart2)

*
*-- Decide which type of calculation we are to run. At this stage only
*   QM equilibration and QM production is available. Should probably
*   be so in the future, hence all-classical should be removed
*   at some stage. Also, ordinary sampfile analysis is performed
*   seperately.
*
      If(Smeq.or.Smprod) then
        Write(6,*)
        Write(6,*)'All classical simulation is currently not available.'
        Call Quit(_RC_GENERAL_ERROR_)
      Elseif(Qmeq.or.Qmprod) then  !Qmeq=.true. is default option.
        If(QmType(1:3).eq.'SCF') then
          Call EqScf(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C
     &              ,iSupDeAll,iV1DeAll)
        Elseif(QmType(1:4).eq.'RASS') then
          Call EqRas(iQ_Atoms,nAtomsCC,Coord,nBas,nBas_C,nCnC_C
     &              ,iBigDeAll,nSizeBig,ip,nRM)
        Endif
      Endif
*
*********JoseMEP
* If MEP option true. Calculate the Interaction energy coming from the
* Interaction with Mean Elec. Potential, Field and Field-Gradient
* E=q*Pot-dipol*field+quad*field-grad
* Also is included the non_Electr perturbation
* For a quantum AO i nonE=d*Sum_k <Xi|Chi_k>*Sum_k <Xi|Chi_k>*E_k
* These Energies are added to the SEWARD One-electron file
* as a perturbation to the One-e Hamiltoniam
*
      If(lExtr(8)) then
        Labjhr='Pert'
        Call AverMEP(Labjhr,Eint,Poli,iCi,SumElcPot
     &               ,NCountField,PertElcInt
     &               ,iQ_Atoms,nBas(1),nOcc(1),natyp(1),nntyp)
      Endif
*********
*
*-- If we are doing single point calculations, then jump back up.
*
      If(SingPoint) then
        If(nCalls.lt.nPart2) then
          Write(6,*)
          Write(6,*)
          Write(6,*)
          Write(6,*)'     *********** << NEW CALCULATION >> ***********'
          Write(6,'(A,I3,A)')'Single point-calculation nr.',nCalls+1,' '
     &//'follows below.'
          Write(6,*)
          Write(6,*)
          Go To 3333
        Else
          Call GetMem('Store','Free','Real',ipStoreCoo,nPart2*nCent*3)
          Write(6,*)
          Write(6,*)
          Write(6,*)' Single-point calculation ended.'
          Write(6,*)
        Endif
      Endif

*
*-- Deallocations of stuff from qfread to simulations. These are unique
*   for the QM-method.
*
      If(QmType(1:4).eq.'RASS') then
        Call GetMem('ALLES','Free','Real',iBigDeAll,nSizeBig)
        If(MoAveRed) then
          Call GetMem('UncleMoe','Free','Real',ip,nBas(1)*nRM)
        Endif
      Elseif(QmType(1:3).eq.'SCF') then
        nS=iOrb(1)*(iOrb(1)+1)/2
        Call GetMem('SUPER','Free','Real',iSupDeAll,nS**2)
        Call GetMem('OrbCoeffQ','Free','Real',iV1DeAll,iOrb(1)*nBas(1))
      Endif

*
*-- Exit
*
666   Continue

      ireturn=0
      Call QExit('QMSTAT')

      Return
      End
