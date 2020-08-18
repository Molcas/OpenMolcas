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
      Subroutine PrInp_MCLR(iPL)
************************************************************************
*                                                                      *
*     Echo input                                                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "disp_mclr.fh"
#include "sa.fh"
      Character*8   Fmt1,Fmt2
      Character XYZ(3)
      Character*100  Line,BlLine,StLine
      Data XYZ / 'X','Y','Z' /
*----------------------------------------------------------------------*
*     Start and define the paper width                                 *
*----------------------------------------------------------------------*
      lPaper=110
*     lPaper=80
*----------------------------------------------------------------------*
*     Initialize blank and header lines                                *
*----------------------------------------------------------------------*

      lLine=Len(Line)
      Do i=1,lLine
         BlLine(i:i)=' '
         StLine(i:i)='*'
      End Do
      lPaper=132
*     left=(lPaper-lLine)/2
      left=5
      Write(Fmt1,'(A,I3.3,A)') '(',left,'X,A)'
      Write(Fmt2,'(A,I3.3,A)') '(',left,'X,'
*----------------------------------------------------------------------*
*     Print the project title                                          *
*----------------------------------------------------------------------*
*     If ( mTit.gt.0 ) then
         If (iPL.ge.3) Write(6,*)
         nLine=mTit+5
         Do i=1,nLine
            Line=BlLine
            If ( i.eq.1 .or. i.eq.nLine ) Line=StLine
            If ( i.eq.3 ) Line='Project:'
            If ( i.ge.4 .and. i.le.nLine-2 )
     &         Write(Line,'(18A4)')(TitleIN((i-4)*18+j),j=1,18)
            If (iPL.ge.3) Then
               Call Center(Line)
               Write(6,Fmt1) '*'//Line//'*'
            End If
         End Do
         If (iPL.ge.3) Write(6,*)
*     End If
*----------------------------------------------------------------------*
*     Print file identifiers                                           *
*----------------------------------------------------------------------*
      If (iPL.ge.3) Then
         Write(6,*)
         Write(6,Fmt1) 'Header of the ONEINT file:'
         Write(6,Fmt1) '--------------------------'
         Write(Line,Fmt1)  Header1I(1)
         Write(6,'(A)') Line(:mylen(Line))
         Write(Line,Fmt1)  Header1I(2)
         Write(6,'(A)') Line(:mylen(Line))
         Write(6,*)
*----------------------------------------------------------------------*
*     Print cartesian coordinates of the system                        *
*----------------------------------------------------------------------*
         Write(6,*)
         Write(6,Fmt1)'Cartesian coordinates:'
         Write(6,Fmt1)'----------------------'
         Write(6,*)
       Write(6,Fmt1)'----------------------------------------------'
       Write(6,Fmt1)' No.    Label       X         Y         Z     '
       Write(6,Fmt1)'----------------------------------------------'
         Do 10 iAt=1,nAtoms
            Write(6,Fmt2//'I3,5X,A6,2X,3F10.5)')
     &      iAt,AtLbl(iAt),Coor(1,iAt),Coor(2,iAt),Coor(3,iAt)
10       Continue
       Write(6,Fmt1)'----------------------------------------------'
       Write(6,Fmt2//'A,F20.10)')'Nuclear repulsion energy =',PotNuc
      End If
*----------------------------------------------------------------------*
*     Print orbital and wavefunction specifications                    *
*----------------------------------------------------------------------*
      If (iMethod.eq.iCASSCF) Then
         ntIsh=0
         ntAsh=0
         ntSsh=0
         ntBas=0
         Do 20 iSym=1,nSym
            ntIsh=ntIsh+nIsh(iSym)
            ntAsh=ntAsh+nAsh(iSym)
            ntBas=ntBas+nBas(iSym)
            ntSsh=ntSsh+nBas(iSym)
     &           -nFro(iSym)-nDel(iSym)-nIsh(iSym)-nAsh(iSym)
20       Continue
         If (iPL.ge.2) Then
            Write(6,*)
            Line=''
            Write(Line(left-2:),'(A)') 'Wave function specifications:'
            Call CollapseOutput(1,Line)
            Write(6,Fmt1)              '-----------------------------'
            Write(6,*)
            Write(6,Fmt2//'A,T47,I6)')
     &                  'Number of closed shell electrons',
     &                              2*ntIsh
            Write(6,Fmt2//'A,T47,I6)')
     &                  'Number of electrons in active shells',
     &                              nActEl
            Write(6,Fmt2//'A,T47,I6)')
     &               'Max number of holes in RAS1 space',
     &                           nHole1
            Write(6,Fmt2//'A,T47,I6)')
     &                  'Max number of electrons in RAS3 space',
     &                              nElec3
            Write(6,Fmt2//'A,T47,I6)')
     &                 'Number of inactive orbitals',
     &                              ntIsh
            Write(6,Fmt2//'A,T47,I6)')
     &                 'Number of active orbitals',
     &                              ntAsh
            Write(6,Fmt2//'A,T47,I6)')
     &                 'Number of secondary orbitals',
     &                              ntSsh
            Write(6,Fmt2//'A,T47,F6.1)')
     &                'Spin quantum number',
     &                              (dble(iSpin)-1.)/2.
            Write(6,Fmt2//'A,T47,I6)')
     &               'State symmetry',
     &                             State_Sym
            Write(6,Fmt2//'A,T47,I6)') 'Number of roots',nroots
            Write(6,Fmt2//'A,(T47,10I6))')
     &       'States considered ',(iroot(i),i=1,nroots)
            Write(6,Fmt2//'A,(T47,10F6.3))') 'Weights ',
     &           (weight(i),i=1,nroots)
            Write(6,*)
            Write(6,Fmt2//'A,T47,8I6)')
     &           'Symmetry species',
     &                            (i,i=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)')
     &           'Skiped sym. species',
     &            (nSkip(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)')
     &            'Frozen orbitals',
     &            (nFro(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)')
     &            'Inactive orbitals',
     &                               (nIsh(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)')
     &             'Active orbitals',
     &                               (nAsh(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)') 'RAS1 orbitals',
     &                              (nRs1(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)') 'RAS2 orbitals',
     &                              (nRs2(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)') 'RAS3 orbitals',
     &                              (nRs3(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)') 'Deleted orbitals',
     &                              (nDel(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)')
     &               'Number of basis functions',
     &                              (nBas(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I6)')
     &               'Number of Orbitals',
     &                              (nOrb(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T47,8I8)')
     &            'Number of configurations',
     &                             (ncsf(isym),isym=1,nsym)

            Write(6,Fmt2//'A,T47,8I6)')
     &            'Number of combinations',
     &                             (nint(xispsm(isym,1)),isym=1,nsym)
*
            If ( iPt2.eq.0 ) then
               Write(6,Fmt1)
     &              'Natural orbitals are used in the last CI'
            Else
               Write(6,Fmt1)
     &          'Pseudo canonical orbitals are used in the last CI'
            End If
*
            Write(6,Fmt2//'A,T33,F20.10)')
     &           'RASSCF state energy = ',ERASSCF(istate)
            Write(6,Fmt2//'A,T47,I6)')
     &          'Size of explicit Hamiltonian in PCG: ',nExp_Max
            Call CollapseOutput(0,'Wave function specifications:')
         End If
      Else
         If (iPL.ge.2) Then
            Write(6,*)
            Line=''
            Write(Line(left-2:),'(A)') 'Wave function specifications:'
            Call CollapseOutput(1,Line)
            Write(6,Fmt1)              '-----------------------------'
            Write(6,Fmt2//'A,T50,A)')
     &                'Wavefunction type:','SCF'
            Write(6,Fmt2//'A,T45,I6)')
     &                'Number of irreducible symmetry groups',
     &                         nSym
            Write(6,Fmt2//'A,T45,8I6)')
     &                'Number of basis functions',
     &                         (nBas(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T45,8I6)')
     &                'Number of occupied orbitals',
     &                         (nish(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T31,F20.10)')
     &             'SCF energy = ',ESCF
            Call CollapseOutput(0,'Wave function specifications:')
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPL.ge.2) Then
*                                                                      *
************************************************************************
*                                                                      *
         Write(6,*)
         Write(6,Fmt2//'A,T42,ES11.4)')
     &      'Convergence threshold= ',Epsilon
         Write(6,Fmt2//'A,T45,I8)')
     &      'Max number of iterations in PCG: ',nIter
*
      If (SPINPOL) Then
         Write(6,Fmt1) 'CALCULATING SPIN POLARIZATION'
      Else If (PT2) Then
         Write(6,Fmt2//'A,A)') 'CALCULATING LAGRANGIAN MULTIPLIER',
     &                      ' FOR CASPT2'
      Else If (SA.or.iMCPD) Then
         If (isNAC) Then
            Write(6,Fmt2//'A,I3,"/",I3)')'Lagrangian multipliers '//
     &                            'are calculated for states no. ',
     &                            NACstates(1),NACstates(2)
            If ((NSSA(1).ne.NACstates(1)).or.
     &          (NSSA(2).ne.NACstates(2))) Then
               Write(6,Fmt2//'39X,A,I3,"/",I3,A)')'(SA roots no. ',
     &                               NSSA(1),NSSA(2),')'
            End If
         Else
            Write(6,Fmt2//'A,I3)')'Lagrangian multipliers are '//
     &                            'calculated for state no. ',irlxroot
            If (istate.ne.irlxroot) Then
               Write(6,Fmt2//'39X,A,I3,A)')'(SA root no. ',istate,')'
            End If
         End If
         If(TwoStep) Then
            If(StepType.eq.'RUN1') Write(6,Fmt1)
     &                      'TwoStep activated. Run 1 (preparation).'
            If(StepType.eq.'RUN2') Write(6,Fmt1)
     &                      'TwoStep activated. Run 2 (final run).'
         End If
      Else
         If (ndisp.ne.0) Then
            Write(6,*)
            Line=''
            Write(Line(left-2:),'(A)') 'Perturbation specifications:'
            Call CollapseOutput(1,Line)
            Write(6,Fmt1)              '----------------------------'
            Write(6,*)
            Write(6,Fmt2//'A,T47,8I4)')
     &             'Number of perturbations in each symmetry',
     &                           (ldisp(iSym),iSym=1,nSym)
            Write(6,Fmt2//'A,T50,A)') 'Type of perturbation:',
     &                            Perturbation
            Call CollapseOutput(0,'Perturbation specifications:')
            Write(6,*)
            Line=''
            Write(Line(left-2:),'(A)') 'Perturbations:'
            Call CollapseOutput(1,Line)
            Write(6,Fmt1)              '--------------'
            Write(6,*)
            Write(6,Fmt1)
     &             '-------------------------------------'
            Write(6,Fmt1)
     &          ' No.    Symmetry    Center Direction '
            Write(6,Fmt1)
     &             '-------------------------------------'
            jDisp=0
            Do iSym=1,nSym
               Do iDisp=1,lDisp(iSym)
                  jDisp=jDisp+1
                  If (iand(ntpert(jdisp),16).eq.16) Then
                     Write(6,Fmt2//'I3,T16,A3,T29,A)')
     &                  jDisp,chIrr(isym),ChDisp(dspvec(jDisp))
                  Else
                     Write(6,Fmt2//'I3,T16,A3,T29,A8,A,A)')
     &                  jDisp,chIrr(isym),swlbl(jDisp),' ',
     &                  XYZ(dspvec(jDisp))
                  End If
               End Do
            End Do
            Write(6,Fmt1) '-------------------------------------'
            Call CollapseOutput(0,'Perturbations:')
         End If
      End If
      Write(6,*)
*----------------------------------------------------------------------*
*     Print reference state information                                *
*----------------------------------------------------------------------*
      If (.Not. SA) Then
         Write(6,*)
         If (iMethod.eq.iCASSCF) Then
             Write(6,Fmt2//'A,I3)')
     &          'Linear response function is computed '//
     &                        'for root no. = ',irlxroot
         Else
            Write(6,Fmt2//'A,I3)')
     &         'Linear response function is computed '//
     &                       'for Restricted Hartree-Fock wavefunction'
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Write(6,*)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
