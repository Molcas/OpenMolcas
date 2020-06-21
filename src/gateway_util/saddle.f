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
* Copyright (C) 2009, Roland Lindh                                     *
*               2010, Mickael G. Delcey                                *
************************************************************************
      SubRoutine Saddle(DInf,nDinf)
************************************************************************
*                                                                      *
* Object: to set up for a TS optimization with the SADDLE approach.    *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : qEnter                                                  *
*              OpnCom                                                  *
*                                                                      *
*     Author: Roland Lindh, Dep. of Theoretical Chemistry,             *
*             University of Lund, SWEDEN                               *
*             January 2009                                             *
************************************************************************
      use external_centers
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Character*1 Mode
      Character*16 StdIn
      Logical Not_First_Iter, FindTS,Found,quadratic,Invar
#include "angstr.fh"
#include "warnings.fh"
      Character*2, Dimension(:), Allocatable :: Elm
      Real*8, Dimension(:), Allocatable :: TanVec, TmpA, W
      Real*8, Dimension(:,:), Allocatable :: Vec, MEP
      Integer, Dimension(:), Allocatable :: iStab
      Integer ipX2, ipX3
      Integer ipRef ipOpt
      Real*8, Dimension(:,:), Allocatable :: XYZ
      Real*8 DInf(nDInf)
#include "periodic_table.fh"
************************************************************************
*                                                                      *
*                            Prologue                                  *
*                                                                      *
************************************************************************
*
      nSaddle_Max=100
      Delta_Max=0.10D0
      ratio=Zero
      quadratic=.false.
*
**    Avoid warnings
*
      R11=Zero
      R22=Zero
      R1R2=Zero
      iX0=-1
      iX1=-1
*
**    If lRP true in Info and Saddle block active but set to zero then
**    this is after the Saddle procedure is terminated and lRP should
**    be ignored.
*
      If (lRP) Then
         Call qpg_dArray('Saddle',Not_First_Iter,nData)
         If (Not_First_Iter.and.nData.eq.0) lRP=.False.
      End If
      lRP_Post=lRP
*                                                                      *
************************************************************************
*                                                                      *
*     Get informations about Saddle in RunFile
*
      If (lRP) Then
*
**       Find out whether the energy is assumed invariant to trans. & rot.
*
         Call Get_iScalar('System BitSwitch',iSBS)
         Invar=(iAnd(iSBS,2**7).eq.0).and.(iAnd(iSBS,2**8).eq.0)
*
         nAt = nRP / 3
         Call qpg_dArray('Saddle',Not_First_Iter,nData)
         nSaddle=2*nRP+5
         Call mma_allocate(TmpA,nSaddle,label='TmpA')
         If (Not_First_Iter) Then
            Call Get_dArray('Saddle',TmpA,nSaddle)
            Update=TmpA(2*nRP+5)
*
            If (Update.ne.One) Then
               HSR=TmpA(2*nRP+4)
               Call mma_deallocate(TmpA)
               lRP_Post=.False.
               Go To 100
            End If
            Update=Zero
*
            call dcopy_(3*nAt,TmpA(1      ),1,RP_Centers(1,1,1),1)
            call dcopy_(3*nAt,TmpA(1+3*nAt),1,RP_Centers(1,1,2),1)
            E1    =TmpA(6*nAt+1)
            E2    =TmpA(6*nAt+2)
            HSR0  =TmpA(6*nAt+3)
         Else
            dHSR=Zero ! Dummy initialize
            Update=Zero
            HSR0  = Zero ! Dummy initialize
            iOff_Iter=0
            Call Put_iScalar('iOff_Iter',iOff_Iter)
*                                                                      *
************************************************************************
*                                                                      *
*           Retrieve the weights even if the structures are not
*           going to be aligned explicitly
*
            Call mma_Allocate(XYZ,3*nAt*8,2,label='XYZ')
            iReac=1
            iProd=2
            Call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iReac),mAt,
     &                       nIrrep,iOper)
            Call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iProd),mAt,
     &                       nIrrep,iOper)
            Call Qpg_dArray('Weights',Found,nData)
            If (Found.And.(nData.ge.mAt)) Then
              Call mma_allocate(W,nData,label='W')
              Call Get_dArray('Weights',W,nData)
            Else
              Call SysAbendMsg('Saddle',
     &             'No or wrong weights were found in the RUNFILE.','')
            End If
*
*           Get the symmetry stabilizers for each center
*
            Call mma_allocate(iStab,nAt,label='iStab')
            iAt=1
            nsc=0
            Do i=1,nCnttp
               Do iCnt=1,nCntr(i)
                  nsc=nsc+1
                  If (.Not.(pChrg(i).Or.FragCnttp(i).Or.AuxCnttp(i)))
     &                Then
                     iStab(iAt)=jStab(1,nsc)
                     iAt=iAt+1
                  End If
               End Do
            End Do
*
*           Align the reactant and product structures the first time.
*           Only if energy is invariant
*
            If (Do_Align .and. Invar) Then
* Note: this might break symmetry
               Call Superpose_w(XYZ(1,iReac),XYZ(1,iProd),W,mAt,
     &                          RMS,RMSMax)
               Call Fix_Symmetry(XYZ(1,iReac),nAt,iStab)
               Call Add_Info('RMSD',[RMS],1,6)
               Call Add_Info('RMSMax',[RMSMax],1,6)
               call dcopy_(3*nAt,XYZ(1,iReac),1,RP_Centers(1,1,1),1)
               call dcopy_(3*nAt,XYZ(1,iProd),1,RP_Centers(1,1,2),1)
*
               If (Align_Only) Then
*
**           Get the atomic symbols with symmetry unfolded
*
                  Call mma_allocate(Elm,nAt)
                  iAt=0
                  iAtSym=nAt
                  ndc=0
                  Do iCnttp=1,nCnttp
                    Do iCnt=1,nCntr(iCnttp)
                      ndc=ndc+1
                      If (.Not.(pChrg(iCnttp).Or.
     &                          FragCnttp(iCnttp).Or.
     &                          AuxCnttp(iCnttp))) Then
                        iAt=iAt+1
                        Elm(iAt)=PTab(iAtmNr(iCnttp))
                        Do i=1,nIrrep/nStab(ndc)-1
                          iAtSym=iAtSym+1
                          Elm(iAtSym)=Elm(iAt)
                        End Do
                      End If
                    End Do
                  End Do
*
                  Write (6,*)
                  Write (6,*) 'Aligned Reactants and Products'
                  Write (6,*) '=============================='
                  Write (6,*)
                  Write (6,*)
                  Write (6,*) ' Reactants / Angstrom'
                  Write (6,*) '====================='
                  Write (6,*)
                  Do iAt = 1, mAt
                     Write (6,'(A,1X,3F15.8)') Elm(iAt),
     &                     (XYZ((iAt-1)*3+jAt,iReac)*Angstr,jAt=1,3)
                  End Do
                  Write (6,*)
                  Write (6,*)
                  Write (6,*) ' Products / Angstrom'
                  Write (6,*) '===================='
                  Write (6,*)
                  Do iAt = 1, mAt
                     Write (6,'(A,1X,3F15.8)') Elm(iAt),
     &                     (XYZ((iAt-1)*3+jAt,iProd)*Angstr,jAt=1,3)
                  End Do
                  Write (6,*)
                  Write (6,*)
                  Call WarningMessage(2,'Molecular alignment completed')
                  Write (6,*)
                  iReturn=_RC_ALL_IS_WELL_
                  Call mma_deallocate(XYZ)
                  Call mma_deallocate(iStab)
                  Call mma_deallocate(W)
                  Call mma_deallocate(Elm)
                  Call ClsSew()
                  Call xQuit(iReturn)
               End If
            End If
*
            Call mma_deallocate(XYZ)
            Call mma_deallocate(iStab)
            Call mma_deallocate(W)
         End If
*                                                                      *
************************************************************************
*                                                                      *
*        Determine the distance between the two structures
*
         Call mma_allocate(iStab,nAt,label='iStab')
         iAt=1
         nsc=0
         Do i=1,nCnttp
            Do iCnt=1,nCntr(i)
               nsc=nsc+1
               If (.Not.(pChrg(i).Or.FragCnttp(i).Or.AuxCnttp(i))) Then
                  iStab(iAt)=jStab(1,nsc)
                  iAt=iAt+1
               End If
            End Do
         End Do
*
         Call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
         iRA1=1
         iRA2=2
         Call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iRA1),mAt,nIrrep,
     &                    iOper)
         Call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iRA2),mAt,nIrrep,
     &                    iOper)
         Call Qpg_dArray('Weights',Found,nData)
         If (Found.And.(nData.ge.mAt)) Then
           Call mma_allocate(W,nData,label='W')
           Call Get_dArray('Weights',W,nData)
         Else
           Call SysAbendMsg('Saddle',
     &          'No or wrong weights were found in the RUNFILE.','')
         End If
         If (.Not.Invar) Then
*
**         If the energy is not trans/rot invariant, compute the weighted
**         RMS with the current structures
*
           HSR = Zero
           wTot = Zero
           iOff = 1
           Do i = 1, mAt
             Do ixyz = 1, 3
               diff = XYZ(iOff,iRA2)-XYZ(iOff,iRA1)
               HSR = HSR + W(i)*diff**2
               iOff = iOff+1
             End Do
             wTot = wTot + W(i)
           End Do
           HSR = Sqrt( HSR/wTot )
         Else
           Call Get_RMSD_w(XYZ(1,iRA2),XYZ(1,iRA1),W,mAt,HSR)
         End If
         Call mma_deallocate(XYZ)
*
         If (E1.le.E2) Then
            Mode='R'
         Else
            Mode='P'
         End If
*
************************************************************************
*                                                                      *
*                   Determine the desired HSR                          *
*                                                                      *
************************************************************************
         If (Not_First_Iter) Then
*
**   FindTS
*
            FindTS=HSR.le.(1.5d0*SadStep)
            If (FindTS) Then
               Write (6,*) '**************************'
               Write (6,*) '* Enable TS optimization *'
               Write (6,*) '**************************'
               Update=2.0d0
               Delta=0.5D0 - 6.25D0 * (E2-E1)
               Delta=Min(Delta,0.75D0)
               Delta=Max(Delta,0.25D0)
               If (Mode.eq.'R') Then
                  HSR=(One-Delta)*HSR
               Else
                  HSR=Delta*HSR
               End If
*
**        Calculate tangent vector
*
               Call mma_allocate(TanVec,3*nAt,label='TanVec')
               If (Mode.eq.'R') Then
                 Call Calc_LSTvec(3*nAt,RP_Centers(1,1,2),
     &                                  RP_Centers(1,1,1),TanVec,Invar)
               Else
                 Call Calc_LSTvec(3*nAt,RP_Centers(1,1,1),
     &                                  RP_Centers(1,1,2),TanVec,Invar)
               End If
               Call Put_dArray('TanVec',TanVec,3*nAt)
               Call mma_deallocate(TanVec)
            Else
*
**         Try to be smart:
**         be slow the first iterations and the last ones
**         and also if large reorganisation the previous iteration
*
               Call Qpg_iScalar('nMEP',Found)
               If (Found) Then
                  quadratic=.true.
                  Call Get_iScalar('nMEP',iSaddle)
*
**        Compute some distances, used for quadratic interpolation
*
                  Call mma_Allocate(MEP,nRP,nSaddle_Max,label='MEP')
                  Call Get_dArray('MEP-Coor    ',MEP,nRP*nSaddle_Max)
                  call mma_allocate(Vec,nRP,2,label='Vec')
                  iX0=iSaddle-2
                  iX1=iSaddle-1
                  if (Mode.eq.'R') Then
                     ipX2 = 1
                     ipX3 = 2
                  Else
                     ipX2 = 2
                     ipX3 = 1
                  End If
                  If (iSaddle.lt.3) iX0=iX1
*
**        Align everything with the current structure (X2)
*
                  Call mma_Allocate(XYZ,3*nAt*8,4)
                  iXA0=1
                  iXA1=2
                  iXA2=3
                  iXA3=4
                  Call Expand_Coor(MEP(1,iX0),nAt,XYZ(1,iXA0),
     &                             mAt,nIrrep,iOper)
                  Call Expand_Coor(MEP(1,iX1),nAt,XYZ(1,iXA1),
     &                             mAt,nIrrep,iOper)
                  Call Expand_Coor(RP_Centers(1,1,ipX2),nAt,XYZ(1,iXA2),
     &                             mAt,nIrrep,iOper)
                  Call Expand_Coor(RP_Centers(1,1,ipX3),nAt,XYZ(1,iXA3),
     &                             mAt,nIrrep,iOper)
                  If (Invar) Then
                    Call Superpose_w(XYZ(1,iXA0),XYZ(1,iXA2),W,
     &                               mAt,RMSD,RMax)
                    Call Fix_Symmetry(XYZ(1,iXA0),nAt,iStab)
                    Call Superpose_w(XYZ(1,iXA1),XYZ(1,iXA2),W,
     &                               mAt,RMSD,RMax)
                    Call Fix_Symmetry(XYZ(1,iXA1),nAt,iStab)
                    Call Superpose_w(XYZ(1,iXA3),XYZ(1,iXA2),W,
     &                               mAt,RMSD,RMax)
                    Call Fix_Symmetry(XYZ(1,iXA3),nAt,iStab)
                  End If
                  iX0=iXA0
                  iX1=iXA1
                  iX3=iXA3
*
**       Compute deviation=(X1-X0)*(X2-X1).
**       If it is lower than 0.8, slow down and do linear interpolation
**       (the direction is probably broken)
*
                  If (iSaddle.gt.2) Then
                   call dcopy_(nRP,XYZ(1,iX1),1,Vec(1,1),1)
                   Call daxpy_(nRP,-One,XYZ(1,iX0),1,Vec(1,1),1)
                   call dcopy_(nRP,RP_Centers(1,1,ipX2),1,Vec(1,2),1)
                   Call daxpy_(nRP,-One,XYZ(1,iX1),1,Vec(1,2),1)
                   deviation=dmwdot(nAt,mAt,Vec(1,1),Vec(1,2))
                   R11=dmwdot(nAt,mAt,Vec(1,2),Vec(1,2))
                   R22=dmwdot(nAt,mAt,Vec(1,1),Vec(1,1))
                   deviation=deviation/Sqrt(R11*R22)
                   If (deviation.lt.0.85d0) Then
                     quadratic=.false.
                     dHSR=SadStep*0.8d0
                     Go To 35
                   EndIf
                  EndIf
*
**       Compute R11=(X3-X1)**2 and R22=(X2-X1)**2
*
                  call dcopy_(nRP,XYZ(1,iX3),1,Vec(1,1),1)
                  Call daxpy_(nRP,-One,XYZ(1,iX1),1,Vec(1,1),1)
                  R11 = dmwdot(nAt,mAt,Vec(1,1),Vec(1,1))
                  call dcopy_(nRP,RP_Centers(1,1,ipX2),1,Vec(1,2),1)
                  Call daxpy_(nRP,-One,XYZ(1,iX1),1,Vec(1,2),1)
                  R22 = dmwdot(nAt,mAt,Vec(1,2),Vec(1,2))
                  R1R2= dmwdot(nAt,mAt,Vec(1,2),Vec(1,1))
*
**       The direction of the previous iteration is far from the R-P direction
*
                  tmp=R1R2/(Sqrt(R11)*Sqrt(R22))
                  If (tmp.lt.0.3D0) Then
                     quadratic=.false.
                     dHSR=SadStep*0.8d0
                  ElseIf (tmp.lt.Zero) Then
                     quadratic=.false.
                     dHSR=SadStep*0.6
                  Else
                     dHSR=SadStep*1.0d0
                     If (isaddle.gt.1) Then
                        dHSR=SadStep*1.3d0
                     EndIf
                  EndIf
  35              Continue
                  Call mma_deallocate(XYZ)
                  Call mma_deallocate(MEP)
                  If (.not.quadratic) Call mma_deallocate(Vec)
*
**      Slow down close to the TS or in the first iteration
*
                  If (HSR.le.2.5d0*SadStep) Then
                     dHSR=SadStep*0.55d0
                  ElseIf (quadratic.and.(HSR.le.4.0d0*SadStep)) Then
                     dHSR=SadStep*0.75d0
                  EndIf
               Else
                  dHSR=SadStep*0.7d0
               EndIf
               Delta=dHSR/HSR
               If (Mode.eq.'P') Delta=One-Delta
               HSR=HSR-dHSR
            End If
            TmpA(6*nAt+5)=Update
            TmpA(6*nAt+4)=HSR
         Else
*                                                                      *
************************************************************************
*                                                                      *
*           This is written here only the first time around
*
            FindTS=.False.
            HSR0=HSR
*
**          Already findTS!
*
            FindTS=HSR.le.(1.5d0*SadStep)
            If (FindTS) Then
               Write (6,*) '**************************'
               Write (6,*) '* Enable TS optimization *'
               Write (6,*) '**************************'
               Update=2.0d0
               Delta=0.5D0 - 6.25D0 * (E2-E1)
               Delta=Min(Delta,0.75D0)
               Delta=Max(Delta,0.25D0)
               If (Mode.eq.'R') Then
                  HSR=(One-Delta)*HSR
               Else
                  HSR=Delta*HSR
               End If
               ratio=1.0d0
*
**        Calculate tangent vector
*
               Call mma_allocate(TanVec,3*nAt,label='TanVec')
               If (Mode.eq.'R') Then
                 Call Calc_LSTvec(3*nAt,RP_Centers(1,1,2),
     &                                  RP_Centers(1,1,1),TanVec,Invar)
               Else
                 Call Calc_LSTvec(3*nAt,RP_Centers(1,1,1),
     &                                  RP_Centers(1,1,2),TanVec,Invar)
               End If
               Call Put_dArray('TanVec',TanVec,3*nAt)
               Call mma_deallocate(TanVec)
            Else
*
**          Do not go too fast the first time
*
               dHSR=SadStep*0.7d0
               Delta=dHSR/HSR
               HSR=HSR-dHSR
               If (Mode.eq.'P') Delta=One-Delta
            End If
            call dcopy_(3*nAt,RP_Centers(1,1,1),1,TmpA(1      ),1)
            call dcopy_(3*nAt,RP_Centers(1,1,2),1,TmpA(1+3*nAt),1)
            TmpA(6*nAt+1)=E1
            TmpA(6*nAt+2)=E2
            TmpA(6*nAt+3)=HSR0
            TmpA(6*nAt+4)=HSR
            TmpA(6*nAt+5)=Update
         End If
*
         Call Put_dArray('Saddle',TmpA,nSaddle)
         Call mma_deallocate(TmpA)
*
************************************************************************
*                                                                      *
*                      Some paper work                                 *
*                                                                      *
************************************************************************
*        Set the point with the highest energy as the reference
*        structure. Put the reference structure on the runfile.
*
         Write (6,*)
         Write (6,'(A)') ' -- TS optimization a la the Saddle approach'
         If (FindTS) Then
            Write (6,'(A)') '   Last Macro iteration'
            Call Merge_Lists(Mode,nAt)
         End If
         If (Mode.eq.'R') Then
             ipRef=2
             ipOpt=1
            Write (6,'(A)') '     Reference structure: product side'
            Write (6,'(A,F15.8)') '       Associated Energy: ',E2
            Write (6,'(A)') '     Optimized structure: reactant side'
            Write (6,'(A,F15.8)') '       Associated Energy: ',E1
         Else
             ipRef=1
             ipOpt=2
            Write (6,'(A)') '     Reference structure: reactant side'
            Write (6,'(A,F15.8)') '       Associated Energy: ',E1
            Write (6,'(A)') '     Optimized structure: product side'
            Write (6,'(A,F15.8)') '       Associated Energy: ',E2
         End If
         Write (6,*)
*
*        Align the reference structure with the current structure
*
         Call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
         iRefAlign=1
         iOptExp  =2
         Call Expand_Coor(RP_Centers(1,1,ipRef),nAt,XYZ(1,iRefAlign),
     &                    mAt,nIrrep,iOper)
         Call Expand_Coor(RP_Centers(1,1,ipOpt),nAt,XYZ(1,iOptExp),
     &                    mAt,nIrrep,iOper)
         If (Invar) Then
           Call Superpose_w(XYZ(1,iRefAlign),XYZ(1,iOptExp),W,
     &                      mAt,RMSD,RMax)
           Call Fix_Symmetry(XYZ(1,iRefAlign),nAt,iStab)
         End If
         Call Put_dArray('Ref_Geom',XYZ(1,iRefAlign),3*nAt)
         Call mma_deallocate(XYZ)
*                                                                      *
************************************************************************
*                                                                      *
*        Fix so that two copies of runfile, vector files, jobiph, etc.
*        are maintained correctly. This is only done the very first
*        iteration.
*
         If (.Not.Not_First_Iter) Then
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_Open(LuInput,StdIn)
*
*           Make two copies of the virgin runfile
*
            Write (LuInput,'(A)')
     &            '> CLONE $Project.RunFile $Project.Reac.RunFile'
            Write (LuInput,'(A)')
     &            '> CLONE $Project.RunFile $Project.Prod.RunFile'
            Write (LuInput,'(A)') '> RM $Project.RunFile'
*
*           Make the appropriate links for the other files.
*
            If (Mode.eq.'R') Then
               Write (LuInput,'(A)') '> EXPORT SubProject=.Reac'
            Else
               Write (LuInput,'(A)') '> EXPORT SubProject=.Prod'
            End If
            Write (LuInput,'(A)')
     &            '> CLONE $Project.GssOrb $Project$SubProject.GssOrb'
            Write (LuInput,'(A)') '> RM -FORCE $Project.GssOrb'
*
*           Signal that this is a saddle run, and the first iteration in a branch
*
            Write (LuInput,'(A)') '> EXPORT MOLCAS_SADDLE=1'
            Write (LuInput,'(A)') '> EXPORT SADDLE_FIRST=1'
            Close(LuInput)
         Else
            lRP_Post=.False.
         End If
*
************************************************************************
*                                                                      *
*       Compute the starting structure and put it into the runfile     *
*                                                                      *
************************************************************************
 30      Continue
         If (quadratic) Then
*
**       Quadratic interpolation:
**       Find the solutions of the system
*
            R1_2=R11-Two*R1R2+R22
            Delta=HSR**2*(R11*R22*R1_2+HSR**2*(R1R2**2-R11*R22))
            If (Delta.lt.Zero) Then
               Write(6,*) 'Delta is negative!!!'
               quadratic=.false.
               Go To 30
            End If
            C=(HSR**2*(Two*R1R2**2-(R11+R1R2)*R22)+
     &                  (-Two*R1R2+R22)*Sqrt(Delta))/(R11*R22*R1_2)
            D=(HSR**2*(R22-R1R2)+Sqrt(Delta))/(R22*R1_2)
*
**      Compute the next structure
*
            Call daxpy_(nRP,(C-One),Vec(1,1),1,Vec(1,1),1)
            Call daxpy_(nRP,D,Vec(1,2),1,Vec(1,1),1)
            Call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
            iRA1=1
            iRA2=2
            Call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iRA1),mAt,
     &                       nIrrep,iOper)
            Call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iRA2),mAt,
     &                       nIrrep,iOper)
            If (Mode.eq.'R') Then
              If (Invar) Then
                Call Superpose_w(XYZ(1,iRA2),XYZ(1,iRA1),W,mAt,
     &                           RMSD,RMax)
                Call Fix_Symmetry(XYZ(1,iRA2),nAt,iStab)
              End If
              Call daxpy_(nRP,One,XYZ(1,iRA2),1,Vec(1,1),1)
            Else
              If (Invar) Then
                Call Superpose_w(XYZ(1,iRA1),XYZ(1,iRA2),W,mAt,
     &                           RMSD,RMax)
                Call Fix_Symmetry(XYZ(1,iRA1),nAt,iStab)
              End If
              Call daxpy_(nRP,One,XYZ(1,iRA1),1,Vec(1,1),1)
            End If
            Call mma_deallocate(XYZ)
            j=1
            Do iCnttp=1,nCnttp
              If (.Not.pChrg(iCnttp).and..Not.FragCnttp(iCnttp) .and.
     &            .Not.AuxCnttp(iCnttp)) Then
                ixyz=ipCntr(iCnttp)
                Do iCnt=1,nCntr(iCnttp)
                  Do i=0,2
                    DInf(ixyz)=Vec(j,1)
                    ixyz=ixyz+1
                    j=j+1
                  End Do
                End Do
              End If
            End Do
            If (Not_First_Iter) Call Put_Coord_New(Vec(1,1),nAt)
            Call mma_deallocate(Vec)
         Else
*                                                                      *
************************************************************************
*                                                                      *
*        Linear interpolation
*
            Call mma_allocate(TmpA,nRP,label='TmpA')
            Call DZero(TmpA,nRP)
            Call mma_allocate(XYZ,3*nAt*8,2,label='XYZ')
            iRA1=1
            iRA2=2
            Call Expand_Coor(RP_Centers(1,1,1),nAt,XYZ(1,iRA1),mAt,
     &                       nIrrep,iOper)
            Call Expand_Coor(RP_Centers(1,1,2),nAt,XYZ(1,iRA2),mAt,
     &                       nIrrep,iOper)
            If (Invar) Then
              If (Mode.eq.'R') Then
                Call Superpose_w(XYZ(1,iRA2),XYZ(1,iRA1),W,mAt,
     &                           RMSD,RMax)
                Call Fix_Symmetry(XYZ(1,iRA2),nAt,iStab)
              Else
                Call Superpose_w(XYZ(1,iRA1),XYZ(1,iRA2),W,mAt,
     &                           RMSD,RMax)
                Call Fix_Symmetry(XYZ(1,iRA1),nAt,iStab)
              End If
            End If
            Call daxpy_(nRP,(One-Delta),XYZ(1,iRA1),1,TmpA,1)
            Call daxpy_(nRP,(    Delta),XYZ(1,iRA2),1,TmpA,1)
            Call mma_deallocate(XYZ)
            j=1
            Do iCnttp=1,nCnttp
              If (.Not.pChrg(iCnttp).and..Not.FragCnttp(iCnttp) .and.
     &            .Not.AuxCnttp(iCnttp)) Then
                ixyz=ipCntr(iCnttp)
                Do iCnt=1,nCntr(iCnttp)
                  Do i=0,2
                    DInf(ixyz)=TmpA(j)
                    ixyz=ixyz+1
                    j=j+1
                  End Do
                End Do
              End If
            End Do
            If (Not_First_Iter) Call Put_Coord_New(TmpA,nAt)
            Call mma_deallocate(TmpA)
         End If
*
         Call mma_deallocate(iStab)
         Call mma_deallocate(W)
*                                                                      *
************************************************************************
*                                                                      *
*        Write constraint on the runfile for slapaf to pick up later.
*
 100     Continue
         Lu_UDC=97
         Lu_UDC = IsFreeUnit(Lu_UDC)
         Call Molcas_Open(Lu_UDC,'UDC.Saddle')
         Write (Lu_UDC,*) 'R = Sphere'
         Write (Lu_UDC,*) 'Value'
         Write (Lu_UDC,*) 'R = ', HSR, ' soft'
         Write (Lu_UDC,*) 'END'
         Close(Lu_UDC)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
************************************************************************
************************************************************************
************************************************************************
      Real*8 function dmwdot(nAt,mAt,A,B)
      Implicit Real*8 (a-h,o-z)
      Integer nAt,mAt
      Real*8 A(3,nAt),B(3,nAt)
      Logical Found
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
      Real*8, Dimension(:), Allocatable :: W
************************************************************************
*                                                                      *
*      Object : compute the weighted dot product of two vectors        *
*                                                                      *
************************************************************************
      tmp=Zero
      TMass=Zero
      Call Qpg_dArray('Weights',Found,nData)
      If (Found.And.(nData.ge.mAt)) Then
        Call mma_allocate(W,nData,label='W')
        Call Get_dArray('Weights',W,nData)
      Else
        Call SysAbendMsg('dmwdot',
     &       'No or wrong weights were found in the RUNFILE.','')
      End If
      iAt=0
      Do iCnttp = 1, nCnttp
         If (.Not.pChrg(iCnttp).and..Not.FragCnttp(iCnttp) .and.
     &       .Not.AuxCnttp(iCnttp)) Then
             Do iCnt = 1, nCntr(iCnttp)
                iAt = iAt + 1
                Fact=DBLE(iDeg(A(1,iAt),iOper,nIrrep))
                xMass=Fact*W(iAt)
                TMass=TMass+xMass
                Do i = 1, 3
                   tmp=tmp + xMass*A(i,iAt)*B(i,iAt)
                End Do
             End Do
         End If
      End Do
      Call mma_deallocate(W)
      dmwdot=tmp/TMass
      return
      end
************************************************************************
************************************************************************
************************************************************************
      subroutine calc_LSTvec(mynRP,Reac,Prod,TanVec,Invar)
      Implicit Real*8 (a-h,o-z)
      Real*8 Reac(mynRP),Prod(mynRP),TanVec(mynRP),norm
      Logical Found,Invar
#include "itmax.fh"
#include "info.fh"
#include "real.fh"
#include "stdalloc.fh"
      Integer, Dimension(:), Allocatable :: iStab
      Real*8, Dimension(:), Allocatable :: W
      Real*8, Dimension(:,:), Allocatable :: XYZ
************************************************************************
*                                                                      *
*      Object : compute the LST vector, that is Reac-Prod              *
*      Of course it is not that simple, because of metric problems     *
*                                                                      *
*      M.G. Delcey     June 2010                                       *
*      Lund University                                                 *
*                                                                      *
************************************************************************
*
*     Get the symmetry stabilizers for each center
*
      nAt=mynRP/3
      Call mma_Allocate(iStab,nAt,label='iStab')
      iAt=1
      nsc=0
      Do i=1,nCnttp
         Do iCnt=1,nCntr(i)
            nsc=nsc+1
            If (.Not.(pChrg(i).Or.FragCnttp(i).Or.AuxCnttp(i))) Then
               iStab(iAt)=jStab(1,nsc)
               iAt=iAt+1
            End If
         End Do
      End Do
*
**    Align the structures and compute the vector (un-weighted)
*
      Call mma_allocate(XYZ,3*nAt*8,2)
      iReacA=1
      iProdA=2
      Call Expand_Coor(Reac,nAt,XYZ(1,iReacA),mAt,nIrrep,iOper)
      Call Expand_Coor(Prod,nAt,XYZ(1,iProdA),mAt,nIrrep,iOper)
      Call Qpg_dArray('Weights',Found,nData)
      If (Found.And.(nData.ge.mAt)) Then
        Call mma_allocate(W,nData,label='W')
        Call Get_dArray('Weights',W,nData)
      Else
        Call SysAbendMsg('calc_LSTvec',
     &       'No or wrong weights were found in the RUNFILE.','')
      End If
      If (Invar) Then
        Call Superpose_w(XYZ(1,iReacA),XYZ(1,iProdA),W,mAt,RMSD,RMax)
        Call Fix_Symmetry(XYZ(1,iReacA),nAt,iStab)
      End If
      call dcopy_(mynRP,XYZ(1,iReacA),1,TanVec,1)
      Call daxpy_(mynRP,-One,XYZ(1,iProdA),1,TanVec,1)
      Call mma_deallocate(XYZ)
      Call mma_deallocate(iStab)
      Call mma_deallocate(W)
*
**    And normalize it
*
      norm=DDot_(mynRP,TanVec,1,TanVec,1)
      Call dScal_(mynRP,One/Sqrt(norm),TanVec,1)
*     Call RecPrt('TanVec',' ',TanVec,3,nAt)
      end
