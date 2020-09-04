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
      SubRoutine Prepare(nGrdPt,ipGrid,ipB,ipGrdI)
      use Basis_Info
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
*
*     Some stuff for preparing the gradient integral computation
*
#include "espf.fh"
*
#include "disp.fh"
      Logical TstFnc, DoRys
      Character*1 xyz(0:2)
      Data xyz/'x','y','z'/
*
      Call qEnter('prepare')
*
      LuWr=6
      DoRys = .True.
      nDiff = 3
      Call IniSew(DoRys,nDiff)
*
*     Copy the grid coordinates and weights in ONE array
*     This is the only solution I found to pass info trough oneel_g !
*
      Do iPnt = 1, nGrdPt
         Work(ipGrdI+(iPnt-1)*4  ) = Work(ipGrid+(iPnt-1)*3  )
         Work(ipGrdI+(iPnt-1)*4+1) = Work(ipGrid+(iPnt-1)*3+1)
         Work(ipGrdI+(iPnt-1)*4+2) = Work(ipGrid+(iPnt-1)*3+2)
         Work(ipGrdI+(iPnt-1)*4+3) = Work(ipB+iPnt-1)
      EndDo
*
      nCnttp_Valence=0
      Do iCnttp = 1, nCnttp
         If (dbsc(iCnttp)%Aux) Exit
         nCnttp_Valence = nCnttp_Valence+1
      End Do
*
*---- Compute number of centers and displacements. Ignore pseudo centers.
*     If any pseudo centers disable use of translational and rotational
*     invariance.
*
      mDisp = 0
      mdc = 0
      Do 10 iCnttp = 1, nCnttp_Valence
         If (dbsc(iCnttp)%pChrg) Then
             mdc = mdc + dbsc(iCnttp)%nCntr
             Go To 10
         End If
         Do 20 iCnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
            mDisp = mDisp + 3*(nIrrep/nStab(mdc))
 20      Continue
 10   Continue
*
*
*     Initialize the Direct array. Why? I don't know.
*
      Do i = 1, 3*mxdc
         Direct(i) = .True.
      EndDo
*
*     Generate symmetry adapted cartesian displacements
*
      Call ICopy(mxdc*8,[0],0,IndDsp,1)
      Call ICopy(mxdc*3,[0],0,InxDsp,1)
      call dcopy_(3*MxSym*mxdc,[One],0,Disp_Fac,1)
      Call ICopy(3*mxdc,[1],0,mult_Disp,1)
      nDisp = 0
      Do iIrrep = 0, nIrrep-1
         lDisp(iIrrep) = 0
*        Loop over basis function definitions
         mdc = 0
         mc = 1
         Do iCnttp = 1, nCnttp_Valence
*           Loop over unique centers associated with this basis set.
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
               IndDsp(mdc,iIrrep) = nDisp
*              Loop over the cartesian components
               Do iCar = 0, 2
                  iComp = 2**iCar
                  If ( TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                nIrrep/nStab(mdc),iChTbl,iIrrep,
     &                iComp,nStab(mdc)) .and.
     &                .Not.dbsc(iCnttp)%pChrg ) Then
                      nDisp = nDisp + 1
                      If (iIrrep.eq.0) InxDsp(mdc,iCar+1) = nDisp
                      lDisp(iIrrep) = lDisp(iIrrep) + 1
                      mult_Disp(nDisp)=nIrrep/nStab(mdc)
                      If (iIrrep.eq.0) Then
                         Do jOper = 0, nIrrep-1
                            Disp_Fac(iCar+1,jOper,mdc)=
     &                        iPrmt( jOper ,iComp) *iChTbl(iIrrep,jOper)
                         End Do
                      End If
                      Write (ChDisp(nDisp),'(A,1X,A1)')
     &                       dc(mdc)%LblCnt,xyz(iCar)

                  End If
               End Do
               mc = mc + nIrrep/nStab(mdc)
            End Do
         End Do
*
      End Do
*
      If (nDisp.ne.mDisp) Then
         Call WarningMessage(2,'Error in espf/prepare')
         Write (LuWr,*)
     &      ' Wrong number of symmetry adapted displacements',
     &       nDisp,'=/=',mDisp
         Call Abend()
      End If
*
      Call qExit('prepare')
      Return
      End
