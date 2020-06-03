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
      Subroutine Convrg(iter,kIter, nInter, qInt, Shift, Grad,
     &                  Lbl,GNrm,Energy,Stat,MaxItr,Stop,iStop,ThrCons,
     &                  ThrEne, ThrGrd, MxItr, UpMeth, HUpMet, mIntEff,
     &                  Baker, Cx,Gx,nAtom,mTtAtm,iOper,nSym,ed,iNeg,
     &                  GoOn,Step_Trunc,GrdMax,StpMax,GrdLbl,StpLbl,
     &                  Analytic_hessian,rMEP,MEP,nMEP,Numerical,
     &                  Just_Frequencies,FindTS,ipCoor,eMEPTest,nLambda,
     &                  TSReg)
      Use Chkpnt
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "weighting.fh"
#include "nadc.fh"
#include "print.fh"
#include "warnings.fh"
      Real*8 Shift(nInter,iter),Grad(nInter,iter),Cx(3*nAtom,iter+1),
     &       Gx(3*nAtom,iter+1),GNrm(iter),Energy(iter+1),
     &       qInt(nInter,iter+1),Maxed,MaxErr
      Character Lbl(nInter)*8, Stat(0:MaxItr)*128, GrdLbl*8, StpLbl*8
      Character*6 UpMeth, HUpMet, ConLbl(5)*5
      Character*1 Step_Trunc
      Character*16 StdIn
      Character*80 Point_Desc
      Integer   iOper(0:nSym-1), iNeg(2)
      Logical Stop, Conv1, Baker, GoOn,Analytic_hessian, MEP,
     &        Found, Terminate, Numerical, Last_Energy, rMEP,
     &        Just_Frequencies, Saddle, FindTS, eMEPTest, eTest,
     &        IRCRestart, Conv2, ConvTmp, TSReg, BadConstraint
      Character*8 Temp
*
      Lu=6
      nSaddle_Max=100
      iRout=116
      iPrint=nPrint(116)
      If (iPrint.ge.99) Then
         Call RecPrt('Convrg: Energy',' ',Energy,1,iter)
         Call RecPrt('Convrg: Grad',' ',Grad,nInter,iter)
         Call RecPrt('Convrg: Shift',' ',Shift,nInter,iter)
         Call RecPrt('Convrg: qInt',' ',qInt,nInter,iter+1)
         Call RecPrt('Convrg: Cx',' ',Cx,3*nAtom,iter+1)
         Call RecPrt('Convrg: Gx',' ',Gx,3*nAtom,iter+1)
      End If
      Call QEnter('Convrg')
*
      Call Get_iScalar('Saddle Iter',iter_S)
      If (iter_S.eq.0) Then
         iter_S=1
         Call Put_iScalar('Saddle Iter',iter_S)
         Call f_Inquire('RUNFILE2',Found)
         If (Found) Then
            Call NameRun('RUNFILE2')
            Call Put_iScalar('Saddle Iter',iter_S)
            Call NameRun('RUNFILE')
         End If
      End If
      If (Analytic_Hessian) Then
         If (HUPMET.eq.'  No  '.or.HUPMET.eq.' None ') Then
*           Temp='Analytic'
            Temp='Computed'
         Else
            Temp(1:6)= HUPMET(1:6)
            Temp(7:8)='  '
         End if
      Else
         Temp(1:6)= HUPMET(1:6)
         Temp(7:8)='  '
      End If
*
      Call GetMem('x','Allo','Real',ipx,3*mTtAtm)
      Call GetMem('y','Allo','Real',ipy,3*mTtAtm)
      Call AtmLst(Cx(1,iter  ),nAtom,Work(ipx),iOper,nSym,mTtAtm)
      Call AtmLst(Cx(1,iter+1),nAtom,Work(ipy),iOper,nSym,mTtAtm)
      Call OptRMS_Slapaf(Work(ipx),Work(ipy),mTtAtm,RMS,RMSMax)
      Call GetMem('y','Free','Real',ipy,3*mTtAtm)
      Call GetMem('x','Free','Real',ipx,3*mTtAtm)
*
      If (kIter.ne.iter .and. kIter.eq.1) Then
         Fabs   = GNrm(kiter)
         E = Energy(kiter)
      Else
         Fabs   = GNrm(iter)
         E = Energy(iter)
      End If
      Fabs = Max(Zero,Fabs)
      E0 = E + ed
      Energy(iter+1)=E0
      Call FZero(Gx(1,iter+1),3*nAtom)
      If (kiter.eq.1) Then
         eChng=Zero
      Else
         If (kIter.ne.iter .and. kIter.eq.2) Then
            eChng=Energy(iter)-Energy(1)
         Else
            eChng=Energy(iter)-Energy(iter-1)
         End If
      End If
*
      If (MEP.or.rMEP) Then
         Saddle=.False.
         iMEP=0
         Call Qpg_iScalar('nMEP',Found)
         If (Found) Call Get_iScalar('nMEP',iMEP)
         If (iMEP.eq.0) Then
            iOff_Iter=0
            Call Put_iScalar('iOff_Iter',iOff_Iter)
         Else
            Call Get_iScalar('iOff_Iter',iOff_Iter)
         End If
      Else
         iOff_Iter=0
         iSaddle=0
         Call qpg_dArray('Saddle',Saddle,nSaddle)
         Saddle=Saddle.and..NOT.Just_Frequencies
         If (Saddle) Then
            Call Qpg_iScalar('nMEP',Found)
            If (Found) Call Get_iScalar('nMEP',iSaddle)
            Call Get_iScalar('iOff_Iter',iOff_Iter)
         End If
      End If
*
*----- Convergence criteria
*
*      Too many iterations
*
*      or
*
* 1)   a la Baker
*      Abs(GrdMax).lt.ThrGrd
*      and
*      Abs(eChng).lt.ThrEne or RMSMax.lt.ThrGrd
*
* 2)   a la Gaussian
*      Abs(Fabs/Sqrt(mIntEff)).lt.ThrGrd
*      and
*      Abs(GrdMax).lt.ThrGrd*1.5
*      and
*      ((RMS.lt.ThrGrd*4
*        and
*        RMSMax.lt.ThrGrd*6)
*       or
*       Abs(eChng).lt.ThrEne)
*
      If (Baker) Then
         Val1= Abs(eChng)
         Thr1= ThrEne
         If (kIter.le.1) Then
            ConLbl(1)=' --- '
         Else If (Val1.lt.Thr1) Then
            ConLbl(1)=' Yes '
         Else
            ConLbl(1)=' No  '
         End If
         Val2= RMSMax
         Thr2= ThrGrd
         If (Val2.lt.Thr2) Then
            If (Step_Trunc.eq.' ') Then
               ConLbl(2)=' Yes '
            Else
               ConLbl(2)=' No *'
            End If
         Else
            ConLbl(2)=' No  '
         End If
         Val3= Abs(GrdMax)
         Thr3= ThrGrd
         If (Val3.lt.Thr3) Then
            ConLbl(3)=' Yes '
         Else
            ConLbl(3)=' No  '
         End If
         Conv1= Val1.lt.Thr1.and.kIter.gt.1
         Conv1= Conv1.or. (Val2.lt.Thr2 .and. Step_Trunc.eq.' ')
         Conv1= Conv1.and. Val3.lt.Thr3
      Else
         Val2=Abs(Fabs/Sqrt(DBLE(mIntEff)))
         Thr2=ThrGrd
         Conv1=Val2.lt.Thr2
         If (Conv1) Then
            ConLbl(2)=' Yes '
         Else
            ConLbl(2)=' No  '
         End If
         Conv1= Conv1.and.Abs(GrdMax).lt.ThrGrd*1.5D0
         Val4=Abs(GrdMax)
         Thr4=ThrGrd*1.5D0
         ConvTmp=Val4.lt.Thr4
         Conv1=Conv1.and.ConvTmp
         If (ConvTmp) Then
            ConLbl(4)=' Yes '
         Else
            ConLbl(4)=' No  '
         End If
         Conv2= RMS.lt.ThrGrd*4.D0 .and. Step_Trunc.eq.' '
         Val1=RMS
         Thr1=ThrGrd*4.0D0
         ConvTmp=Val1.lt.Thr1
         Conv2=ConvTmp .and. Step_Trunc.eq.' '
         If (ConvTmp) Then
            If (Step_Trunc.eq.' ') Then
               ConLbl(1)=' Yes '
            Else
               ConLbl(1)=' No *'
            End If
         Else
            ConLbl(1)=' No  '
         End If
         Val3=RMSMax
         Thr3=ThrGrd*6.0D0
         ConvTmp=Val3.lt.Thr3
         Conv2=Conv2.and.ConvTmp
         If (ConvTmp) Then
            If (Step_Trunc.eq.' ') Then
               ConLbl(3)=' Yes '
            Else
               ConLbl(3)=' No *'
            End If
         Else
            ConLbl(3)=' No  '
         End If
         Val5=Abs(eChng)
         Thr5=ThrEne
         ConvTmp=Val5.lt.Thr5 .and. kIter.gt.1
         If (ConvTmp) Then
            ConLbl(5)=' Yes '
         Else
            If (kIter.gt.1) Then
               ConLbl(5)=' No  '
            Else
               ConLbl(5)=' --- '
            End If
         End If
         Conv2=Conv2.or.ConvTmp
         Conv1=Conv1.and.Conv2
      End If
*
      Stop = Conv1 .or. (kIter-iOff_Iter).ge.MxItr ! CGG
      iStop=1
      If (kIter-iOff_Iter.ge.MxItr) iStop=16       ! CGG
      If (Conv1.or.Just_Frequencies)  iStop= 0
*
      If (GoOn) Then
         Stop=.False.
         iStop=1
      Else
         If (Just_Frequencies) Stop=.True.
         nPrint(52)=nPrint(52)+1
         nPrint(54)=nPrint(54)+1
         iPrint    =iPrint+1
         nPrint(53)=nPrint(53)+1
      End If
      If (.Not.Just_Frequencies) Then
         Call Status(kIter-iOff_Iter,E,Fabs,GrdMax,GrdLbl,StpMax,StpLbl,
     &               E0,Stat,MaxItr-1,eChng,iNeg,UpMeth,Temp,Step_Trunc,
     &               .NOT.Numerical)
      End If
*
      If (Baker) Then
         If (iPrint.ge.5) Then
            Write (Lu,'(A)')
     &        '                +----------------------------------+'
            Write (Lu,'(A)')
     &        '                +  Value      Threshold Converged? +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
            Write (Lu,4)
     &        '+ Max. gradient +',Val3,' ',Thr3,'    ',ConLbl(3),'  +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
            Write (Lu,4)
     &        '+ Max. disp.    +',Val2,' ',Thr2,'    ',ConLbl(2),'  +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
            Write (Lu,4)
     &        '+ Energy diff.  +',Val1,' ',Thr1,'    ',ConLbl(1),'  +'
            Write (Lu,'(A)')
     &        '+---------------+----------------------------------+'
 4          Format(A,2(ES11.4,A),A,A)
         End If
      Else
         If (iPrint.ge.5) Then
            Write (Lu,'(A)')
     &        '       +----------------------------------+'
     &      //'----------------------------------+'
            Write (Lu,'(A)')
     &        '       +    Cartesian Displacements       +'
     &      //'    Gradient in internals         +'
            Write (Lu,'(A)')
     &        '       +  Value      Threshold Converged? +'
     &      //'  Value      Threshold Converged? +'
            Write (Lu,'(A)')
     &        ' +-----+----------------------------------+'
     &      //'----------------------------------+'
            Write (Lu,5)
     &        ' + RMS +',Val1,   ' ',Thr1,'    ',ConLbl(1),  '  +',
     &                Val2,   ' ',Thr2,'    ',ConLbl(2),  '  +'
            Write (Lu,'(A)')
     &       ' +-----+----------------------------------+'
     &          //'----------------------------------+'
            Write (Lu,5)
     &        ' + Max +',Val3,   ' ',Thr3,'    ',ConLbl(3),  '  +',
     &                Val4,   ' ',Thr4,'    ',ConLbl(4),  '  +'
            Write (Lu,'(A)')
     &       ' +-----+----------------------------------+'
     &      //'----------------------------------+'
            If (ThrEne.gt.Zero) Then
               Write (Lu,5)
     &           ' + dE  +',Val5,   ' ',Thr5,'    ',ConLbl(5),  '  +'
               Write (Lu,'(A)')
     &          ' +-----+----------------------------------+'
            End If
            Write (Lu,*)
 5          Format(A,2(2(ES11.4,A),A,A))
         End If
      End If
*
      If (Stop.and.Conv1) Then
         Call Qpg_dScalar('Max error',Found)
         If (Found) Call Get_dScalar('Max error',MaxErr)
         If (MaxErr.gt.ThrCons) Then
            iStop=1
            Conv1=.False.
            Stop=.False.
            Write(Lu,'(A,ES11.4)') 'Maximum constraint error: ',MaxErr
            Write(Lu,*)
         End If
      End If
*
      nConst=0
      If (iNeg(1).eq.0) Then
         If (EDiffZero) Then
            If (NADC) Then
               nConst=2
            Else
               nConst=1
            End If
            Point_Desc='Minimum Energy Crossing Point Structure'
         Else
            Point_Desc='Minimum Structure'
         End If
      Else If (iNeg(1).eq.1) Then
         Point_Desc='Transition State Structure'
      Else
         Point_Desc='Higher Order Saddle Point Structure'
      End If
      If (nLambda.gt.nConst) Point_Desc='Constrained '//Trim(Point_Desc)
      If (iPrint.ge.5) Then
         If (Stop) Then
            If (Conv1) Then
               Write (Lu,'(A,I3,A)') ' Geometry is converged in ',
     &            kIter-iOff_iter,' iterations to a '//Trim(Point_Desc)
            Else
               Write (Lu,'(A)') ' No convergence after max iterations'
               If (Lu.ne.6) Write (6,'(/A)')
     &                          ' No convergence after max iterations'
            End If
         Else
            Write (Lu,'(A)') ' Convergence not reached yet!'
         End If
      End If
      If (FindTS.and.Stop.and.Conv1) Then
         If (.Not.TSReg) Then
            If (iPrint.ge.5) Then
               Write (Lu,*)
               Write (Lu,'(A)')
     &' FindTS was requested, but the TS regime was not reached.'
               Write (Lu,'(A)')
     &' The converged structure is probably not the desired TS.'
            End If
            iStop=16
         End If
      End If
c      If (iPrint.eq.7) Then
c         Write (Lu,*)
c         Write (Lu,'(A)') '*********************************'//
c     &      ' Geometry Statistics for Geometry Optimization '//
c     &                    '*********************************'
c         nPrint(118) = 7
c         Call List(' Internal coordinates ',Lbl,qInt,nInter,iter+1)
c         Call List(' Internal forces    ',Lbl,Grad,nInter,iter)
c      End If
*
*     The energy change should not be too large
      Maxed=1.0d2
      If (Abs(ed).gt.Maxed) Then
         Write (6,*) 'The predicted energy change is too large: ',ed
         Write (6,'(A)') ' This can''t be right!'
         Write (6,'(A)') ' This job will be terminated.'
         iStop=8
         Stop=.True.
      End If
      If (iPrint.ge.5) Then
         Write (Lu,*)
         Write (Lu,'(A)') '*********************'//
     &      '*******************************************************'
     &      //'*************************************'
         Write (Lu,'(A)') '*********************'//
     &      '*******************************************************'
     &      //'*************************************'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Terminate=.False.
      IRCRestart=.False.
*                                                                      *
************************************************************************
*                                                                      *
*-----Write summary of conical intersection characterization data
*
      If (Conv1.and.NADC.and.EDiffZero) Then
         If (iPrint.ge.5) Then
            If (.Not.ApproxNADC) Call CI_Summary(Lu)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Book keeping for Saddle optimization for a TS. Note that this will
*     have to be done on both of the runfiles!
*
      If (Conv1.and.Saddle) Then
*
*        Here if a macro iteration in the Saddle TS optimization is
*        completed.
*
         ENew=Energy(iter)+ed
         Call Allocate_Work(ipTmp,nSaddle)
*
*        Store the info for later generation of MOLDEN formated files
*
         Call Get_dArray('Saddle',Work(ipTmp),nSaddle)
         E_Reac=Work(ipTmp+6*nAtom  )
         E_Prod=Work(ipTmp+6*nAtom+1)
*
         Call Allocate_Work(ipE,nSaddle_Max)
         Call Allocate_Work(ipC,3*nAtom*nSaddle_Max)
         Call Allocate_Work(ipG,3*nAtom*nSaddle_Max)
         If (iSaddle.eq.0) Then
*
*           Initiate with data from the starting points
*
            Call FZero(Work(ipE),nSaddle_Max)
            Call FZero(Work(ipC),3*nAtom*nSaddle_Max)
            Call FZero(Work(ipG),3*nAtom*nSaddle_Max)
            iSaddle=1
            If (E_Reac.le.E_Prod) Then
               Work(ipE+(iSaddle-1))=E_Reac
               call dcopy_(3*nAtom,Work(ipTmp        ),1,
     &                    Work(ipC+(iSaddle-1)*3*nAtom),1)
            Else
               Work(ipE+(iSaddle-1))=E_Prod
               call dcopy_(3*nAtom,Work(ipTmp+3*nAtom),1,
     &                    Work(ipC+(iSaddle-1)*3*nAtom),1)
            End If
            call dcopy_(3*nAtom,[Zero],0,
     &                 Work(ipG+(iSaddle-1)*3*nAtom),1)
*
         Else
*
            Call Get_dArray('MEP-Energies',Work(ipE),nSaddle_Max)
            Call Get_dArray('MEP-Coor',Work(ipC),3*nAtom*nSaddle_Max)
            Call Get_dArray('MEP-Grad',Work(ipG),3*nAtom*nSaddle_Max)
*
         End If
*
*        Add the new data
*
         iSaddle=iSaddle+1
         Work(ipE+(iSaddle-1))=Energy(iter)
         call dcopy_(3*nAtom,Cx(1,iter),1,
     &              Work(ipC+(iSaddle-1)*3*nAtom),1)
         call dcopy_(3*nAtom,Gx(1,iter),1,
     &              Work(ipG+(iSaddle-1)*3*nAtom),1)
*
*        Put data on RUNFILE
*
         Call Put_dArray('MEP-Energies',Work(ipE),nSaddle_Max)
         Call Put_dArray('MEP-Coor',Work(ipC),3*nAtom*nSaddle_Max)
         Call Put_dArray('MEP-Grad',Work(ipG),3*nAtom*nSaddle_Max)
         Call Put_iScalar('nMEP',iSaddle)
*
         Call Free_Work(ipG)
         Call Free_Work(ipC)
         Call Free_Work(ipE)
*                                                                      *
************************************************************************
*                                                                      *
*        Now update the "Saddle" field on both runfiles.
*
         Do iFile = 1, 2
*
         If (iFile.eq.1) Then
            Call NameRun('RUNREAC')
         Else
            Call NameRun('RUNPROD')
         End If
*
*        Update info on the runfile.
*
         Call Get_dArray('Saddle',Work(ipTmp),nSaddle)
         E1=Work(ipTmp+6*nAtom  )
         E2=Work(ipTmp+6*nAtom+1)
C        Write (6,*) 'ENew=',ENew
C        Write (6,*) 'E1,E2=',E1,E2
         If (E1.le.E2) Then
C           Write (6,*) 'Update reactant'
            Work(ipTmp+6*nAtom  )=Energy(iter)
            E1=Energy(iter)
            call dcopy_(3*nAtom,Cx(1,iter),1,Work(ipTmp        ),1)
         Else
C           Write (6,*) 'Update product'
            Work(ipTmp+6*nAtom+1)=Energy(iter)
            E2=Energy(iter)
            call dcopy_(3*nAtom,Cx(1,iter),1,Work(ipTmp+3*nAtom),1)
         End If
*        Set flag that seward should process the info! This should not
*        be done for the final macro iteration.
         If (.Not.FindTS) Work(ipTmp+6*nAtom+4)=One
         Call Put_dArray('Saddle',Work(ipTmp),nSaddle)
*
         End Do
*
*                                                                      *
************************************************************************
*                                                                      *
*        Converged or not, create the saddle.molden file
*        after each macro iteration
*
         Call NameRun('RUNREAC')
         Call Allocate_Work(ipE_r,nSaddle_Max)
         Call Allocate_Work(ipC_r,3*nAtom*nSaddle_Max)
         Call Allocate_Work(ipG_r,3*nAtom*nSaddle_Max)
         Call Qpg_iScalar('nMEP',Found)
         if(Found) Then
            Call Get_dArray('MEP-Energies',Work(ipE_r),nSaddle_Max)
            Call Get_dArray('MEP-Coor',Work(ipC_r),3*nAtom*nSaddle_Max)
            Call Get_dArray('MEP-Grad',Work(ipG_r),3*nAtom*nSaddle_Max)
            Call Get_iScalar('nMEP',iSaddle_r)
         Else
            Work(ipE_r)=E_Reac
            call dcopy_(3*nAtom,Work(ipTmp),1,Work(ipC_r),1)
            Call FZero(Work(ipG_r),3*nAtom)
            iSaddle_r=1
         End If
*
         Call NameRun('RUNPROD')
         Call Allocate_Work(ipE_p,nSaddle_Max)
         Call Allocate_Work(ipC_p,3*nAtom*nSaddle_Max)
         Call Allocate_Work(ipG_p,3*nAtom*nSaddle_Max)
         Call Qpg_iScalar('nMEP',Found)
         if(Found) Then
            Call Get_dArray('MEP-Energies',Work(ipE_p),nSaddle_Max)
            Call Get_dArray('MEP-Coor',Work(ipC_p),3*nAtom*nSaddle_Max)
            Call Get_dArray('MEP-Grad',Work(ipG_p),3*nAtom*nSaddle_Max)
            Call Get_iScalar('nMEP',iSaddle_p)
         Else
            Work(ipE_p)=E_Prod
            call dcopy_(3*nAtom,Work(ipTmp+3*nAtom),1,Work(ipC_p),1)
            Call FZero(Work(ipG_p),3*nAtom)
            iSaddle_p=1
         End If
*
*        Merge the two lists
*
         jSaddle=iSaddle_r
         Do iSaddle=iSaddle_p, 1, -1
            jSaddle=jSaddle+1
            Work(ipE_r+(jSaddle-1))=Work(ipE_p+(iSaddle-1))
            call dcopy_(3*nAtom,Work(ipC_p+(iSaddle-1)*3*nAtom),1,
     &                         Work(ipC_r+(jSaddle-1)*3*nAtom),1)
            call dcopy_(3*nAtom,Work(ipG_p+(iSaddle-1)*3*nAtom),1,
     &                         Work(ipG_r+(jSaddle-1)*3*nAtom),1)
         End Do
         Call Free_Work(ipG_p)
         Call Free_Work(ipC_p)
         Call Free_Work(ipE_p)
*
*        Align the structures sequentially, only for visualization
*        (gradients are not changed, though)
* TODO   Rotate the gradients too
*
         jindex=ipC_r
         Do iSaddle=1,iSaddle_r+iSaddle_p-1
            iindex=jindex+3*nAtom
            Call Align(Work(iindex),Work(jindex),nAtom)
            jindex=iindex
         End Do
*
         Call Intergeo('MD_SADDLE',Work(ipE_r),Work(ipC_r),
     &                 Work(ipG_r),nAtom,iSaddle_r+iSaddle_p)
*
         Call Free_Work(ipG_r)
         Call Free_Work(ipC_r)
         Call Free_Work(ipE_r)
*
*                                                                      *
************************************************************************
*                                                                      *
*        If the Saddle TS optimization is not yet completed set up
*        data for the next macro iteration.
*
         If (.Not.FindTS) Then
            Call Free_Work(ipTmp)
C           Write (*,*) 'Reset $SubProject'
*
*           Reset $SubProject for the next macro iteration.
*
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_Open(LuInput,StdIn)
            If (E1.le.E2) Then
               Write (LuInput,'(A)') '> EXPORT SubProject=.Reac'
C              Write (6,*) 'SubProject=.Reac'
               Call NameRun('RUNREAC')
            Else
               Write (LuInput,'(A)') '> EXPORT SubProject=.Prod'
C              Write (6,*) 'SubProject=.Prod'
               Call NameRun('RUNPROD')
            End If
*
*           Signal whether next iteration will be the first in the branch
*
            Call Qpg_iScalar('nMEP',Found)
            If (.Not.Found) Then
               Write (LuInput,'(A)') '> EXPORT SADDLE_FIRST=1'
            End If
            Write (LuInput,'(A,I3)') '> EXIT ',_RC_CONTINUE_LOOP_
            Close(LuInput)
*
*           Set flags to request yet another macro iteration.
*
            Terminate=.False.
            iStop=6
            Stop=.False.
         Else
            Call NameRun('RUNFILE')
            Call Free_Work(ipTmp)
            nSaddle=0
            Call Put_dArray('Saddle',[Zero],nSaddle)
            Call Put_iScalar('nMEP',nSaddle)
         End If
*
*        Update the active runfile wrt the total number of micro
*        iterations done in all macro iterations of this branch.
*
         Call NameRun('RUNFILE')
         Call Put_iScalar('iOff_Iter',iter)
      End If
*
*     Disable first iteration signal right after the first iteration
*     (in each branch)
*
      If ((.Not.Conv1).and.Saddle) Then
         If ((iter.Eq.1).and.(iStop.eq.1)) Then
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_Open(LuInput,StdIn)
            Write (LuInput,'(A)') '> EXPORT SADDLE_FIRST=0'
            Write (LuInput,'(A,I3)') '> EXIT ',_RC_CONTINUE_LOOP_
            Close(LuInput)
            iStop=6
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Book keeping for minimum energy path search
*
      Call Qpg_iScalar('IRC',Found)
      If (Found) Then
         Call Get_iScalar('IRC',IRC)
      Else
         IRC=0
      End If
*
      If (Conv1.and.(MEP.or.rMEP)) Then
*
*        Is this the first iteration or not?
*
         iMEP=iMEP+1
*
*        Save information for the current step
*
         Call Allocate_Work(ipE,nMEP+1)
         Call Allocate_Work(ipC,3*nAtom*(nMEP+1))
         Call Allocate_Work(ipG,3*nAtom*(nMEP+1))
         If (iMEP.gt.1) Then
            Call Get_dArray('MEP-Energies',Work(ipE),nMEP+1)
            Call Get_dArray('MEP-Coor',Work(ipC),3*nAtom*(nMEP+1))
            Call Get_dArray('MEP-Grad',Work(ipG),3*nAtom*(nMEP+1))
         Else
            Call FZero(Work(ipE),nMEP+1)
            Call FZero(Work(ipC),3*nAtom*(nMEP+1))
            Call FZero(Work(ipG),3*nAtom*(nMEP+1))
            Work(ipE+(iMEP-1))=Energy(iOff_iter+1)
            call dcopy_(3*nAtom,Cx(1,iOff_iter+1),1,
     &                 Work(ipC+(iMEP-1)*3*nAtom),1)
            call dcopy_(3*nAtom,Gx(1,iOff_iter+1),1,
     &                 Work(ipG+(iMEP-1)*3*nAtom),1)
         End If
*
         Work(ipE+iMEP)=Energy(iter)
         call dcopy_(3*nAtom,Cx(1,iter),1,Work(ipC+(iMEP)*3*nAtom),1)
         call dcopy_(3*nAtom,Gx(1,iter),1,Work(ipG+(iMEP)*3*nAtom),1)
         Call Put_dArray('MEP-Energies',Work(ipE),nMEP+1)
         Call Put_dArray('MEP-Coor',Work(ipC),3*nAtom*(nMEP+1))
         Call Put_dArray('MEP-Grad',Work(ipG),3*nAtom*(nMEP+1))
         Call Put_iScalar('nMEP',iMEP)
*
*        Save the path so far (energies, coordinates and forces)
*
         Call Intergeo('MD_MEP',Work(ipE),Work(ipC),Work(ipG),nAtom,
     &                 iMEP+1)
*
*        Should we terminate or not? Not done on the first iteration.
*
         If (iMEP.gt.0) Then
*
*           Test if the energy increase (optionally disabled).
            If (eMEPTest) Then
               eTest= Work(ipE+(iMEP  )).gt.Work(ipE+(iMEP-1))
            Else
               eTest=.FALSE.
            End If
            If ( MEP .and. eTest) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  If (IRC.eq.0) Then
                     Write (6,'(A)')
     &                     ' MEP-search terminated'//
     &                     ' due to energy increase!'
                  Else If (IRC.eq.1) Then
                     Write (6,'(A)')
     &                     ' IRC(forward)-search terminated'//
     &                     ' due to energy increase!'
                  Else
                     Write (6,'(A)')
     &                     ' IRC(backward)-search terminated'//
     &                     ' due to energy increase!'
                  End If
                  Write (6,*)
               End If
            End If
*
*           Test if the energy decrease (optionally disabled).
            If (eMEPTest) Then
               eTest= Work(ipE+(iMEP  )).lt.Work(ipE+(iMEP-1))
            Else
               eTest=.FALSE.
            End If
            If (rMEP .and. eTest) Then
               Terminate=.True.
               If (iPrint.ge.5) Then
                  Write (6,*)
                  Write (6,'(A)')
     &                  ' rMEP-search terminated'//
     &                  ' due to energy decrease!'
                  Write (6,*)
               End If
            End If
         End If
*
*        Test on max number of points.
         If ((iMEP.ge.nMEP).and.(.not.Terminate)) Then
            Terminate=.True.
            If (iPrint.ge.5) Then
               Write (6,*)
               If (MEP) Then
                  If (IRC.eq.0) Then
                     Write (6,'(A)') ' MEP-search '//
     &                   'terminated due to max number of path points!'
                  Else If (IRC.eq.1) Then
                     Write (6,'(A)') ' IRC(forward)-search '//
     &                   'terminated due to max number of path points!'
                  Else
                     Write (6,'(A)') ' IRC(backward)-search '//
     &                   'terminated due to max number of path points!'
                  End If
               Else If (rMEP) Then
                  Write (6,'(A)') ' rMEP-search '//
     &                'terminated due to max number of path points!'
               End If
               Write (6,*)
            End If
         End If
*
*        If IRC reset for backward IRC search.
*
         If (Terminate) Then
            If (IRC.ne.0) Then
               If (IRC.eq.1) Then
                  IRCRestart=.True.
               End If
            End If
         End If
*
         Call Chkpnt_update_MEP(IRCRestart)
*
         Call Free_Work(ipE)
         Call Free_Work(ipG)
         Call Free_Work(ipC)
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----List internal coordinates and gradients
*
      kkIter=iter+1
      If (iPrint.ge.8) Then
         Call List(' Internal coordinates ',Lbl,qInt,nInter,kkIter)
         Call List(' Internal forces    ',Lbl,Grad,nInter,iter)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Put out the new reference structure and the new starting
*     structure to be used for the next MEP point.
*     For rMEP keep the reference structure!
*     Note that this is done in weighted Cartesian coordinates!
*
      If ((Conv1.or.(iter.eq.1)).and.(MEP.or.rMEP)) Then
         If ((iMEP.ge.1).and.(iPrint.ge.5)) Then
            Write (6,*)
            Call CollapseOutput(1,'IRC/Minimum Energy Path Information')
         End If
*
         Call MEP_Dir(Cx,Gx,nAtom,iMEP,iOff_iter,iPrint,IRCRestart,
     &                BadConstraint)
         Call Put_iScalar('iOff_Iter',iter)
*
*        Test on constraint misbehavior
         If (BadConstraint.and.(.not.Terminate)) Then
            Terminate=.True.
            If (iPrint.ge.5) Then
               Write (6,*)
               If (MEP) Then
                  If (IRC.eq.0) Then
                     Write (6,'(A)') ' MEP-search '//
     &                   'terminated due to problematic constraint!'
                  Else If (IRC.eq.1) Then
                     Write (6,'(A)') ' IRC(forward)-search '//
     &                   'terminated due to problematic constraint!'
                  Else
                     Write (6,'(A)') ' IRC(backward)-search '//
     &                   'terminated due to problematic constraint!'
                  End If
               Else If (rMEP) Then
                  Write (6,'(A)') ' rMEP-search '//
     &                'terminated due to problematic constraint!'
               End If
               Write (6,*)
            End If
         End If
*
         If (Conv1.and.Terminate) Then
            If (IRC.ne.0) Then
               Call Allocate_Work(ipE,nMEP+1)
               Call Allocate_Work(ipC,3*nAtom*(nMEP+1))
               Call Allocate_Work(ipG,3*nAtom*(nMEP+1))
               Call Get_dArray('MEP-Energies',Work(ipE),nMEP+1)
               Call Get_dArray('MEP-Coor',Work(ipC),3*nAtom*(nMEP+1))
               Call Get_dArray('MEP-Grad',Work(ipG),3*nAtom*(nMEP+1))
               If (IRC.eq.1) Then
                  IRCRestart=.True.
                  IRC=-1
                  Call Put_iScalar('IRC',IRC)
*
*                 Store away data for IRC molden file. Forward part.
*
                  Call Put_dArray('IRC-Energies',Work(ipE),iMEP+1)
                  Call Put_dArray('IRC-Coor',Work(ipC),3*nAtom*(iMEP+1))
                  Call Put_dArray('IRC-Grad',Work(ipG),3*nAtom*(iMEP+1))
                  Call Put_dArray('Ref_Geom',Cx,3*nAtom)
*
*                 Write a temporary file
*                 (will be overwritten when the backward part is done)
*
                  Call Intergeo('MD_IRC',Work(ipE),Work(ipC),Work(ipG),
     &                 nAtom,iMEP+1)
*
                  Terminate=.False.
*
               Else If (IRC.eq.-1) Then
*
*                 Assemble molden file for IRC
*
                  nBackward=iMEP+1
                  Call qpg_dArray('IRC-Energies',Found,nForward)
                  nIRC=nForward+nBackward-1
                  Call GetMem('IRCE','Allo','Real',ipE_IRC,nIRC)
                  Call GetMem('IRCC','Allo','Real',ipC_IRC,3*nAtom*nIRC)
                  Call GetMem('IRCG','Allo','Real',ipG_IRC,3*nAtom*nIRC)
*
                  j=0
                  Do i = nBackward, 1, -1
                     j = j+1
                     Work(ipE_IRC+(j-1))=Work(ipE+(i-1))
                     call dcopy_(3*nAtom,
     &                          Work(ipC+(i-1)*3*nAtom),1,
     &                          Work(ipC_IRC+(j-1)*3*nAtom),1)
                     call dcopy_(3*nAtom,
     &                          Work(ipG+(i-1)*3*nAtom),1,
     &                          Work(ipG_IRC+(j-1)*3*nAtom),1)
                  End Do
*
                  Call Get_dArray('IRC-Energies',
     &                            Work(ipE_IRC+(nBackward-1)),nForward)
                  Call Get_dArray('IRC-Coor',
     &                            Work(ipC_IRC+(nBackward-1)*3*nAtom),
     &                            nForward*3*nAtom)
                  Call Get_dArray('IRC-Grad',
     &                            Work(ipG_IRC+(nBackward-1)*3*nAtom),
     &                            nForward*3*nAtom)
*
                  Call Intergeo('MD_IRC',Work(ipE_IRC),Work(ipC_IRC),
     &                          Work(ipG_IRC),nAtom,nIRC)
*
                  Call Free_Work(ipG_IRC)
                  Call Free_Work(ipC_IRC)
                  Call Free_Work(ipE_IRC)
               End If
               Call Free_Work(ipG)
               Call Free_Work(ipC)
               Call Free_Work(ipE)
            End If
         End If
*
         If (.Not.Terminate) Then
           iStop=1
           Stop=.False.
         End If
*
*        Print out the path so far
*
         If ((iMEP.ge.1).and.(iPrint.ge.5)) Then
            Call Allocate_Work(ipE,nMEP+1)
            Call Allocate_Work(ipC,3*nAtom*(nMEP+1))
            Call Allocate_Work(ipLen,nMEP+1)
            Call Allocate_Work(ipCur,nMEP+1)
            Call Get_dArray('MEP-Energies',Work(ipE),nMEP+1)
            Call Get_dArray('MEP-Coor',Work(ipC),3*nAtom*(nMEP+1))
            Call Get_dArray('MEP-Lengths',Work(ipLen),nMEP+1)
            Call Get_dArray('MEP-Curvatures',Work(ipCur),nMEP+1)
            Write (6,*)
            CumLen=Zero
            If (Work(ipCur+iMEP).ge.Zero) Then
               Write(6,*) '         Cumul.'
               Write(6,*) 'Point  Length (bohr)       Energy  Curvature'
               Write(6,*) '--------------------------------------------'
               Do i = 0, iMEP
                  CumLen=CumLen+Work(ipLen+i)
                  Write (6,200) i,CumLen,Work(ipE+i),Work(ipCur+i)
               End Do
            Else
               Write(6,*) '         Cumul.'
               Write(6,*) 'Point  Length (bohr)       Energy'
               Write(6,*) '---------------------------------'
               Do i = 0, iMEP
                  CumLen=CumLen+Work(ipLen+i)
                  Write (6,200) i,CumLen,Work(ipE+i)
               End Do
            End If
200         Format (1X,I5,1X,F10.6,1X,F16.8,1X,F10.6)
            If (iPrint.gt.6) Then
               Write (6,*)
               Do i = 0, iMEP
                  Call RecPrt(' Coordinates',' ',Work(ipC+i*3*nAtom),
     &                        3,nAtom)
               End Do
            End If
            Call CollapseOutput(0,'IRC/Minimum Energy Path Information')
            Write(6,*)
            Call Free_Work(ipE)
            Call Free_Work(ipC)
            Call Free_Work(ipLen)
            Call Free_Work(ipCur)
         End If
*
      End If
*
      If (IRCRestart) Then
*
*         Prepare the runfile to start from Scratch
*
          iMEP=0
          Call Put_iScalar('nMEP',iMEP)
          Call GetMem(' iter','Allo','Inte',ipItr,7)
          iWork(ipItr)=-99     ! Deactivate the record
          Call Put_iArray('Slapaf Info 1',iWork(ipItr),7)
          Call GetMem(' iter','Free','Inte',ipItr,7)
          iOff_Iter=0
          Call Put_iScalar('iOff_Iter',iOff_Iter)
*
*         Restore data
*
          Call Put_dScalar('Last Energy',Energy(1))
          Call Put_Grad(Gx,3*nAtom)
          Call Put_dArray('Unique Coordinates',Cx,3*nAtom)
          Call Put_Coord_New(Cx,nAtom)
          call dcopy_(3*nAtom,Cx,1,Work(ipCoor),1)
          call dcopy_(3*nAtom,Cx,1,Cx(1,iter+1),1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Figure out if the last energy should be computed!
*
      Last_Energy = Stop .and. iStop.ne.16 .and. iStop.ne.8
      Last_Energy = Last_Energy .and. .Not.MEP .and. .Not.rMEP
      Last_Energy = Last_Energy .and.
     &             .Not. (Numerical .and. kIter.eq.1)
      If (Last_Energy) iStop = 2
*                                                                      *
************************************************************************
*                                                                      *
*     Write geometry file for MOLDEN. Note that this file should not be
*     generated if Slapaf is running new geometries for any numerical
*     procedure!
*
      If (
     &     .NOT. Just_Frequencies  .AND.
     &     .NOT. Numerical
     &   ) Then
         Call Write_QMMM(Cx,nAtom,iter)
         Call Intergeo('MD_GEO',Energy,Cx,Gx,nAtom,iter+1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('Convrg')
      Return
      End
