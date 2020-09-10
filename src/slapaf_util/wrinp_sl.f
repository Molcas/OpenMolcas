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
      Subroutine WrInp_sl(iRow)
      use kriging_mod
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "info_slapaf.fh"
#include "print.fh"
#include "constants.fh"
*
      iRout=3
      iPrint=nPrint(iRout)
*
      Lu=6
*
      Call QEnter('WrInp')
*
      If (lNmHss) Then
         lOld = .False.
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.5) Then
*                                                                      *
************************************************************************
*                                                                      *
      Write (Lu,*)
      Write (Lu,*)
      Call CollapseOutput(1,'      Slapaf input parameters:')
      Write (Lu,'(3X,A)')   '      ------------------------'
      Write (Lu,*)
*                                                                      *
************************************************************************
*                                                                      *
      Write (Lu,'(A,I5)')  ' Max iterations:                           '
     &      ,MxItr
      If (Baker) Then
         Write (Lu,'(A)')  ' Convergence test a la Baker.'
      Else
         Write (Lu,'(A)')  ' Convergence test a la Schlegel.'
      End If
      Write (Lu,'(A,ES8.1)')
     &    ' Convergence criterion on gradient/para.<=:',ThrGrd
      Write (Lu,'(A,ES8.1)')
     &    ' Convergence criterion on step/parameter<=:',ThrGrd
      Write (Lu,'(A,ES8.1)')
     &    ' Convergence criterion on energy change <=:',ThrEne
      Write (Lu,'(A)')
     &    ' Parameters for step-restricted optimization'
      If (.NOT.Kriging) Then
      Write (Lu,'(A,ES9.2)')
     &    ' Max step length (initial seed):          ',Beta
      Else
      Write (Lu,'(A,ES9.2)')
     &    ' Max step length (micro iterations):      ',Beta
      End If
      Write (Lu,*)
*                                                                      *
************************************************************************
*                                                                      *
      If (Kriging) Then
       Write (Lu,*) '-RVO activated with parameters:'
!      Write (Lu,'(A,I6)')
!    &    '   GEK starts at iteration:                   ',nspAI
       Write (Lu,'(A,I6)')
     &    '   Maximum number of data points used in GEK: ',nWndw/2
!      Write (Lu,'(A,I6)')
!   &     '   Parameter of diff. for Matern (p):         ',pAI
       Write (Lu,'(A,I6)')
     &    '   Maximum number of micro iterations:        ',miAI
       If (set_l) Then
          Call Get_dScalar('Value_l',Value_l)
          Write (Lu,*) '  Global characteristic length scale, l:     ',
     &              Value_l
       Else
          Write (Lu,*) '  Individual characteristic length scales set '
     &          //'to reproduce HMF Hessian.'
       End If
*
       If (blaAI) then
          write (6,'(A,F10.5,A)')
     &          '   Baseline is highest energy plus: ',blavAI,' a.u'
       Else
          if (mblAI) then
             write (6,*) '  Baseline set to maximum value of the energy'
          else if (blAI) then
             write (6,'(A,F9.5,A,/,A,F9.5,A)')
     &              '  Baseline (trend function) changed to value:',
     &              blvAI, 'a.u.',
     &              '                                             ',
     &              blvAI * CONV_AU_TO_KJ_PER_MOLE_,
     &              ' kJ/mol'
          endif
       Endif
       write (6,'(A,F10.5,A)')
     &       '   Maximum dispersion accepted:     ',Beta_disp,
     &       ' * abs(g.max.comp)'
      Else
       Write (Lu,*) '-RFO activated with parameters:'
       Write (Lu,'(A,I6)')
     &    '   Maximum number of data points used in RFO: ',nWndw
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (Line_Search) Then
         Write (Lu,'(A)') ' Line search is performed'
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iAnd(iOptC,256).eq.256) Then
*
         Write (Lu,'(1X,A)') '-Constrained optimization.'
*
         If (MEP) Then
            If (IRC.eq.0) Then
               Write (Lu,'(1X,A)') ' Minimum Energy Path (MEP) search'
               Write (Lu,'(1X,A,I5)') ' Max number of points:',nMEP
            Else If (IRC.eq.1) Then
               Write (Lu,'(1X,A)') ' IRC forward search'
               Write (Lu,'(1X,A,I5)') ' Max number of points:',nMEP
            Else
               Write (Lu,'(1X,A)') ' IRC backward search'
               Write (Lu,'(1X,A,I5)') ' Max number of points:',nMEP
            End If
            If (eMEPtest)
     &         Write (Lu,'(1X,A)') ' Stop when energy increases'
            If (MEP_Algo.eq.'GS') Then
               Write (Lu,'(1X,A)') ' MEP optimization algorithm:'
     &                           //' Gonzalez-Schlegel'
            Else If (MEP_Algo.eq.'MB') Then
               Write (Lu,'(1X,A)') ' MEP optimization algorithm:'
     &                           //' Mueller-Brown'
            End If
            If (MEP_Type.eq.'SPHERE    ') Then
               Write (Lu,'(1X,A)') ' Type of constraint: Hypersphere'
            Else If (MEP_Type.eq.'TRANSVERSE') Then
               Write (Lu,'(1X,A)') ' Type of constraint: Hyperplane'
            End If
         End If
*
         If (rMEP) Then
            Write (Lu,'(1X,A)') ' Reverse Minimum Energy Path '
     &                        //'(rMEP) search'
            Write (Lu,'(1X,A,I3)') ' Max number of points:',nMEP
            If (eMEPtest)
     &         Write (Lu,'(1X,A)') ' Stop when energy decreases'
            If (MEP_Type.eq.'SPHERE    ') Then
               Write (Lu,'(1X,A)') ' Type of constraint: Hypersphere'
            Else If (MEP_Type.eq.'TRANSVERSE') Then
               Write (Lu,'(1X,A)') ' Type of constraint: Hyperplane'
            End If
         End If
*
         If (FindTS) Then
            Write (Lu,'(1X,A)')
     &       '-The optimization will home in on a transition state if:'
            Write (Lu,'(A)')
     &       '  a) Negative curvature is encountered, and'
            Write (Lu,'(A,F10.4)')
     &       '  b) the norm of the gradient is below:',GNrm_Threshold
            If (iAnd(iOptC,512).eq.512) Then
               Write (Lu,'(A)') '  TS-search by RS-I-RFO.'
               !If (Kriging) Then
               !   Write (Lu,'(A)') '  TS-search by RV-I-RFO.'
               !Else
               !   Write (Lu,'(A)') '  TS-search by RS-I-RFO.'
               !End If
            Else
               Write (Lu,'(A)') '  TS-search by RS-P-RFO.'
               !If (Kriging) Then
               !   Write (Lu,'(A)') '  TS-search by RV-P-RFO.'
               !Else
               !   Write (Lu,'(A)') '  TS-search by RS-P-RFO.'
               !End If
            End If
         End If
*
      End If
      Write (Lu,*)
*
      If (iAnd(iOptC,128).eq.128) Then
         Write (Lu,'(1X,A)') '-Optimization for minimum.'
         If (iAnd(iOptC,1).eq.1) Then
            Write (Lu,'(A)') '  Optimization method: quasi-NR.'
         Else If (iAnd(iOptC,2).eq.2) Then
            Write (Lu,'(A)') '  Optimization method: C1-DIIS.'
         Else If (iAnd(iOptC,4).eq.4) Then
            Write (Lu,'(A)') '  Optimization method: C2-DIIS.'

         Else If (iAnd(iOptC,8).eq.8) Then
            If (Kriging) Then
               Write (Lu,'(A)') '  Optimization method: RVO.'
            Else
               Write (Lu,'(A)') '  Optimization method: RS-RFO.'
            End If
         Else
            Call WarningMessage(2,' WrInp: Wrong iOptC setting!')
            Write (Lu,*) ' iOptC=',iOptC
            Call Abend()
         End If
      Else
         Write (Lu,'(1X,A)') '-Optimization for transition state.'
         If (iAnd(iOptC,512).eq.512) Then
            Write (Lu,'(A)') '  Optimization method: RS-I-RFO'
            !If (Kriging) Then
            !   Write (Lu,'(A)') '  Optimization method: RV-I-RFO'
            !Else
            !   Write (Lu,'(A)') '  Optimization method: RS-I-RFO'
            !End If
         Else
            Write (Lu,'(A)') '  Optimization method: RS-P-RFO'
            !If (Kriging) Then
            !   Write (Lu,'(A)') '  Optimization method: RV-P-RFO'
            !Else
            !   Write (Lu,'(A)') '  Optimization method: RS-P-RFO'
            !End If
         End If
         If (Mode.gt.0) Then
            Write (Lu,'(A,I2)') '  Original mode to follow:',Mode
         Else
            Write (Lu,'(A)') '  No mode to follow is specified!'
            Write (Lu,'(A)') '  Optimization will follow mode with'
     &                     //' the lowest eigenvalue.'
         End If
      End If
      Write (Lu,*)
*                                                                      *
************************************************************************
*                                                                      *
*.....Special options on DIIS
*
      If (iAnd(iOptC,2).eq.2.or.iAnd(iOptC,4).eq.4) Then
         If (iAnd(iOptC,16).eq.16) Then
            Write (Lu,'(1X,A)') '-DIIS based on <dx|dx>.'
         Else If (iAnd(iOptC,32).eq.32) Then
            Write (Lu,'(1X,A)') '-DIIS based on <g|dx>.'
         Else If (iAnd(iOptC,64).eq.64) Then
            Write (Lu,'(1X,A)') '-DIIS based on <g|g>.'
         Else
            Call WarningMessage(2,' WrInp: Wrong iOptC setting!')
            Write (Lu,*) ' iOptC=',iOptC
            Call Abend()
         End If
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*.....Hessian information
*
      If (Analytic_Hessian) Then
         Write (Lu,'(1X,A)') '-The Hessian is analytic.'
         Write (Lu,'(1X,A)') ' Hessian from either input or runfile.'
      Else
         If (.Not.lOld.and..Not.lNmHss) Then
            If (DDV_Schlegel) Then
              Write (Lu,'(1X,A)') '-Initial Hessian guessed a'
     &                //' la Schlegel.'
            Else
              If (Kriging) Then
                 Write (Lu,'(1X,A)') '-Hessian guessed by'
     &                   //' Kriging surrogate surface.'
              Else
                 Write (Lu,'(1X,A)') '-Initial Hessian guessed by'
     &                   //' Hessian Model Function (HMF).'
                 If (iAnd(iOptC,1024).eq.1024) Then
                   Write (Lu,'(A)') '  HMF augmented with'
     &                   //' weak interactions.'
                 End If
              End If
            End If
         Else If (lOld.and..Not.lNmHss) Then
            Write (Lu,'(1X,A)')
     &       '-Initial Hessian guess was read from a RUNFILE file.'
         Else
            Write (Lu,'(1X,A,/,A,E9.2)')
     &         '-Initial Hessian guess is estimated with finite'
     &       //' differences.',
     &         '    Two point symmetric formula, Delta=',Delta
            If (Cubic) Then
               Write (Lu,'(1X,A)')
     &          '-Cubic force constants evaluated numerically.'
            End If
         End If
      End If
      Write (Lu,*)
*                                                                      *
************************************************************************
*                                                                      *
*.....Hessian update method
*
      If (.NOT.Kriging) Then
      If (iAnd(iOptH,1).eq.1) Then
         Write (Lu,'(1X,A)')
     &       '-Hessian update method: Fletcher-Meyer'
      Else If (iAnd(iOptH,2).eq.2) Then
         Write (Lu,'(1X,A)')
     &       '-Hessian update method: Broyden-Powell'
      Else If (iAnd(iOptH,4).eq.4) Then
         Write (Lu,'(1X,A)')
     &       '-Hessian update method: Broyden-Fletcher-Goldfarb-Shanno'
      Else If (iAnd(iOptH,8).eq.8) Then
         Write (Lu,'(1X,A)') '-Hessian update method: none'
      Else If (iAnd(iOptH,16).eq.16) Then
         Write (Lu,'(1X,A)')
     &       '-Hessian update method: Murtagh-Sargent-Powell'
      Else If (iAnd(iOptH,64).eq.64) Then
         Write (Lu,'(1X,A)')
     &       '-Hessian update method: EU update by Bofill'
      Else If (iAnd(iOptH,128).eq.128) Then
         Write (Lu,'(1X,A)')
     &       '-Hessian update method: TS-BFGS update by Bofill'
      Else
         Call WarningMessage(2,' WrInp: Wrong iOptH setting!')
         Write (Lu,*) ' Nonrecognizable iOptH setting:',iOptH
         Call Abend()
      End If
      If (.Not.(iAnd(iOptH,8).eq.8)) Then
         Write (Lu,'(A,I3)')
     &         '  Max number of points in Hessian update:',
     &                        nWndw
      End If
      If (iAnd(iOptH,32).eq.32) Then
         Write (Lu,'(A)')
     &       '  Hessian update order according to Schlegel'
      End If
      Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (rHidden.ge.Two) Then
         Write(Lu,'(1X,A,/,1X,A,F6.2,A)')
     &    '-Improved QM/MM Hessian.',
     &    ' Hidden atoms until ',rHidden,' bohrs are included.'
         Write (Lu,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iRow.gt.0) Then
         If (Redundant) Then
            Write (Lu,'(1X,A)') '-Relaxation will be done in user'
     &             //' supplied redundant internal coordinates.'
         Else
            Write (Lu,'(1X,A)') '-Relaxation will be done in user'
     &             //' supplied non-redundant internal coordinates.'
         End If
      Else
         If (CurviLinear) Then
            If (HWRS) Then
               If (Redundant) Then
                  Write (Lu,'(1X,A)') '-Relaxation will be done on'
     &                   //' redundant internal coordinates,'
     &                   //' based on'
                  Write (Lu,*) ' force constant weighted redundant'
     &                   //' internal coordinates.'
               Else
                  Write (Lu,'(1X,A)') '-Relaxation will be done on'
     &                   //' non-redundant internal coordinates,'
     &                   //' based on'
                  Write (Lu,*) ' force constant weighted redundant'
     &                   //' internal coordinates.'
               End If
            Else
               If (Redundant) Then
                  Write (Lu,'(1X,A)')
     &               '-Relaxation will be done in redundant'
     &                   //' delocalized internal coordinates.'
               Else
                  Write (Lu,'(1X,A)')
     &               '-Relaxation will be done in non-redundant'
     &                   //' delocalized internal coordinates.'
               End If
            End If
         Else
            If (Redundant) Then
               Write (Lu,'(1X,A)')
     &            '-Relaxation will be done in redundant'
     &                //' Cartesian coordinates.'
            Else
               Write (Lu,'(1X,A)')
     &            '-Relaxation will be done in approximate'
     &                //' non-redundant'
     &                //' Cartesian normal mode coordinates.'
            End If
         End If
      End If
      Write (Lu,*)
*                                                                      *
************************************************************************
*                                                                      *
      If (Ref_Geom) Then
         Write (Lu,'(1X,A)')
     &         '-The origin of the hyper sphere is defined implicitly.'
      End If
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.6) Then
         Write (Lu,*)
         Write (Lu,'(A)') ' Header from ONEINT:'
         Call Banner(Header,2,Len(Header(1))+12)
         Write (Lu,*)
         Call PrList('Symmetry Distinct Nuclear Coordinates / bohr',
     &                AtomLbl,nsAtom,Work(ipCoor),3,nsAtom)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.5)
     &   Call CollapseOutput(0,'      Slapaf input parameters:')
*
      Call QExit('WrInp')
      Return
      End
