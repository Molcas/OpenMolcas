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
* Copyright (C) 2006, Roland Lindh                                     *
************************************************************************
      SubRoutine Output1_Seward(lOPTO)
************************************************************************
*                                                                      *
*     Object: to write the output of seward            .               *
*                                                                      *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
*             September '06                                            *
************************************************************************
      use Basis_Info
      use Center_Info
      use Period
      use GeoList
      use MpmC
      use EFP_Module
      use External_centers
      use Temporary_Parameters
      use DKH_Info
      use Sizes_of_Seward, only: S
      use Real_Info, only: ThrInt, CutInt, RPQMin
      use RICD_Info, only: iRI_Type, LDF, Do_RI, Cholesky,
     &                     Do_acCD_Basis, Skip_High_AC, Cho_OneCenter
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "rinfo.fh"
#include "real.fh"
#include "rmat.fh"
#include "rctfld.fh"
#include "relmp.fh"
#include "relae.fh"
#include "print.fh"
#include "gateway.fh"
#include "localdf.fh"
      Logical l_aCD_Thr, lOPTO, Found
      Logical lNoPair, lPam2, lECP, lPP
      Character(LEN=80) Title(10)
#include "angstr.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint=nPrint(iRout)
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
      lNoPair = .False.
      lPam2   = .False.
      lECP    = .False.
      lPP     = .False.
      Do i = 1, nCnttp
         lNoPair = lNoPair .or. dbsc(i)%NoPair
         lPam2   = lPam2   .or. dbsc(i)%lPam2
         lECP    = lECP    .or. dbsc(i)%ECP
         lPP     = lPP     .or. dbsc(i)%nPP.ne.0
         lECP    = lECP    .or. dbsc(i)%ECP
         lPP     = lECP    .or. dbsc(i)%ECP
      End Do

*                                                                      *
************************************************************************
*                                                                      *
*     Start of output
*
      If (Test) Then
         Write (LuWr,*)
         Write (LuWr,'(15X,88A)') ('*',i=1,45)
         Write (LuWr,'(15X,A)')
     &   '* TEST: SEWARD will only process the input! *'
         Write (LuWr,'(15X,88A)') ('*',i=1,45)
         Go To 99
      End If
*
      iDKH_X_Order = iRELAE/10000
      iParam = iRELAE-1000 - iDKH_X_Order*10000
      iDKH_H_Order = iParam/10
      iParam = iParam - iDKH_H_Order*10
*
      Write (LuWr,'(15X,A)') 'SEWARD will generate:'
      Write (LuWr,'(15X,A,I2)')
     &        '   Multipole Moment integrals up to order ',S%nMltpl
      If (.Not.Prprt) Then
         Write (LuWr,'(15X,A)')    '   Kinetic Energy integrals'
         If (Nuclear_Model.eq.Gaussian_Type) Then
            Write (LuWr,'(15X,A)')  '   Nuclear Attraction integrals'
     &                              //' (finite nuclei - Gaussian type)'
         Else If (Nuclear_Model.eq.Point_Charge) Then
            Write (LuWr,'(15X,A)')  '   Nuclear Attraction integrals'
     &                              //' (point charge)'
         Else
            Write (LuWr,'(15X,A)')  '   Nuclear Attraction integrals'
     &                    //' (finite nuclei -  Modified Gaussian type)'
         End If
         If (lECP) Then
            If (lNoPair) Then
               If (IRELMP.EQ.0) Then
                  Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &               //' modified with ECP and No-Pair contributions'
                Else If (IRELMP.EQ.1) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &         //' modified with ECP and No-Pair (DK1) contributions'
                Else If (IRELMP.EQ.2) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &         //' modified with ECP and No-Pair (DK2) contributions'
                Else If (IRELMP.EQ.3) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &         //' modified with ECP and No-Pair (DK3) contributions'
                Else If (IRELMP.EQ.11) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &               //' modified with ECP and RESC contributions'
                Else If (IRELMP.EQ.21) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &               //' modified with ECP and ZORA contributions'
                Else If (IRELMP.EQ.22) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &               //' modified with ECP and ZORA-FP contributions'
                Else If (IRELMP.EQ.23) Then
                   Write (LuWr,'(15X,A)')
     &                 '   One-Electron Hamiltonian integrals'
     &               //' modified with ECP and IORA contributions'
                End If
            Else
               Write (LuWr,'(15X,A)')
     &                '   One-Electron Hamiltonian integrals'
     &              //' modified with ECP contributions'
            End If
         Else
            Write (LuWr,'(15X,A)')
     &             '   One-Electron Hamiltonian integrals'
         End If
         If(FNMC) Then
            Write (LuWr,'(15X,A)')
     &             '   Finite nuclear mass correction added'
         End If
         If (lPAM2) Then
               Write (LuWr,'(15X,A)')
     &                '   Include potentials for DMFT calculation'
         End If
         If (lRel) Then
            Write (LuWr,'(15X,A)')    '   Mass-Velocity integrals'
            Write (LuWr,'(15X,A)')
     &             '   Darwin One-Electron Contact Term integrals'
         End If
         If (Vlct) Write (LuWr,'(15X,A)') '   Velocity integrals'
         If (DKroll) Then
            If (iRELAE.lt.1000) Then
               If (iRELAE.EQ.0) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic Douglas-Kroll integrals'
               Else If (IRELAE.EQ.1) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic Douglas-Kroll (DK1) integrals'
               Else If (IRELAE.EQ.2) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic Douglas-Kroll (DK2) integrals'
               Else If (IRELAE.EQ.3) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic Douglas-Kroll (DK3) integrals'
               Else If (IRELAE.EQ.4) Then
                  Write (LuWr,'(15X,A)')
     &          '   full Relativistic Douglas-Kroll (DK3) integrals'
               Else If (IRELAE.EQ.11) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic RESC integrals'
               Else If (IRELAE.EQ.21) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic ZORA integrals'
               Else If (IRELAE.EQ.22) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic ZORA-FP integrals'
               Else If (IRELAE.EQ.23) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic IORA integrals'
               Else If (IRELAE.EQ.101) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic X2C integrals'
               Else If (IRELAE.EQ.102) Then
                  Write (LuWr,'(15X,A)')
     &               '   Relativistic BSS integrals'
               Else If (BSS) Then
                  Write  (LuWr,'(15X,A)')
     &                '   Relativistic Barysz-Sadlej-Snijders integrals'
               End If
            Else
               If (LDKroll) Then
                  Write (LuWr,'(17X,A)')
     &               ' Relativistic Local-Douglas-Kroll-Hess integrals:'
                  If (nCtrLD.eq.0) Then
                     If (radiLD.eq.0.0d0) Then
                        Write(LuWr,'(17X,A)')
     &                  '   - Atomic approximation'
                     Else
                        Write(LuWr,'(17X,A)')
     &                  '   - Full local approximation'
                     End If
                  Else
                     Write(LuWr,'(17X,A)')
     &               '   - Partial local approximation:'
                     Write(LuWr,'(17X,A,10(A4))')
     &          '     - Centers: ', (dc(iCtrLD(i))%LblCnt,i=1,nCtrLD)
                     Write(LuWr,'(17X,A,F6.2,A)')
     &          '     - Cutoff radius: ', radiLD, ' Bohr'
                 End If
               Else
                  Write (LuWr,'(17X,A)')
     &               ' Relativistic Douglas-Kroll-Hess integrals:'
               Endif
               If (iParam.eq.1) Then
                  Write (LuWr,'(17X,A)')
     &                  '   - Parametrization         : OPT'
               Else If (iParam.eq.2) Then
                  Write (LuWr,'(17X,A)')
     &                  '   - Parametrization         : EXP'
               Else If (iParam.eq.3) Then
                  Write (LuWr,'(17X,A)')
     &                  '   - Parametrization         : SQR'
               Else If (iParam.eq.4) Then
                  Write (LuWr,'(17X,A)')
     &                  '   - Parametrization         : MCW'
               Else If (iParam.eq.5) Then
                  Write (LuWr,'(17X,A)')
     &                  '   - Parametrization         : CAY'
               End If
               Write (LuWr,'(17X,A,I2)')
     &               '   - DKH order of Hamiltonian:',iDKH_H_order
               Write (LuWr,'(17X,A,I2)')
     &               '   - DKH order of Properties :',iDKH_X_order
               Write (LuWr,'(17X,A)')
     &               '        - multipole moment operators'
               Write (LuWr,'(17X,A)')
     &               '        - electric potential operators'
               Write (LuWr,'(17X,A)')
     &               '        - contact operators'
            End If
         End If
         If (lRF) Then
            If (PCM) Then
               Write (LuWr,'(15X,A)')
     &               '   Reaction Field integrals (PCM)'
            Else If (lLangevin) Then
               Write (LuWr,'(15X,A)')
     &               '   Reaction Field integrals (Langevin)'
            Else
               Write (LuWr,'(15X,A)')
     &               '   Reaction Field integrals (KirkWood-Onsager)'
            End If
         End If
         If (lPP)Write (LuWr,'(15X,A)')
     &                '   Pseudo Potential integrals'
      End If
      If (Allocated(XF)) Write (LuWr,'(15X,A,I6,A)')
     &                       '   External field from',
     &                           nXF, ' point(s) added to the'
     &                           //' one-electron Hamiltonian'
      If (nEF.gt.0 .and. nOrdEF.ge.0) Write (LuWr,'(15X,A,I6,A)')
     &                       '   Electric potential for',
     &                           nEF, ' points'
      If (nEF.gt.0 .and. nOrdEF.ge.1) Write (LuWr,'(15X,A,I6,A)')
     &                       '   Electric field integrals for',
     &                           nEF, ' points'
      If (nEF.gt.0 .and. nOrdEF.ge.2) Write (LuWr,'(15X,A,I6,A)')
     &            '   Electric field gradient integrals for',
     &                           nEF, ' points'
      If (nEF.gt.0 .and. nOrdEF.ge.2) Write (LuWr,'(15X,A,I6,A)')
     &            '   Contact term integrals for',
     &                           nEF, ' points'
      If (Allocated(DMS_Centers)) Write (LuWr,'(15X,A,I6,A)')
     &            '   Diamagnetic shielding integrals for',
     &                           nDMS, ' points'
      If (Allocated(OAM_Center)) Write (LuWr,'(15X,A,3(F7.4,1X),A)')
     &                       '   Orbital angular momentum around (',
     &   (OAM_Center(i),i=1,3),')'
      If (Allocated(OMQ_Center)) Write (LuWr,'(15X,A,3(F7.4,1X),A)')
     &                       '   Orbital magnetic quadrupole around (',
     &   (OMQ_Center(i),i=1,3),')'
      If (Vlct.and.(S%nMltpl.ge.2)) Write (LuWr,'(15X,A,3(F7.4,1X),A)')
     &                       '   Velocity quadrupole around (',
     &   (Coor_MPM(i,3),i=1,3),')'
      If (Allocated(AMP_Center)) Write (LuWr,'(15X,A,3(F7.4,1X),A)')
     & '   Products of Orbital angular momentum operators around (',
     &   (AMP_Center(i),i=1,3),')'
      If (nWel.ne.0) Write (LuWr,'(15X,A,I4,A)')
     &             '   Spherical well for', nWel,
     &             ' exponent(s) added to the'
     &           //' one-electron Hamiltonian'
      If (lAMFI) Write (LuWr,'(15X,A)') '   Atomic mean-field integrals'
      If (lPSOI) Write (LuWr,'(15X,A)')
     & '   (PSO) Paramagnetic Spin-Orbit integrals'
     &     //' calculated from Gen1Int F90 library'
      If (DoFMM) Then
         Write (LuWr,'(15X,A)')
     &     '   Integral environment set up for FMM option'
         Write (LuWr,'(15X,A,F10.5)')
     &     '    - RPQMin: ',RPQMin
      End If
      If (lEFP) Then
#ifdef _EFP_
         Write (LuWr,'(15X,A)')
     &     '   Effective Fragment potentials added       '
         Write (LuWr,'(15X,A,I4)')
     &     '    - # of fragments: ',nEFP_fragments
#else
         Write (LuWr,'(15X,A)')
     &     '   EFP input specified but code not enabled for the option'
         Call Abend()
#endif
      End If
      If (.Not.Onenly) Then
         If (Cholesky) Then
            Write (LuWr,'(15X,A)')
     &        '   Cholesky decomposed two-electron'
     &           //' repulsion integrals'
           If (Cho_OneCenter) Then
            Write (LuWr,'(17X,A,G10.2)')
     &                 '  - 1C-CD Threshold: ',Thrshld_CD
           Else
            Write (LuWr,'(17X,A,G10.2)')
     &                 '  - CD Threshold: ',Thrshld_CD
           EndIf
         Else If (Do_RI) Then
            If (LocalDF) Then
               If (LDF_Constraint.eq.-1) Then
                  Write(LuWr,'(15X,A)')
     &                  '   Local Density Fitting coefficients'
               Else
                  Write(LuWr,'(15X,A)')
     &               '   Constrained Local Density Fitting coefficients'
                  If (LDF_Constraint.eq.0) Then
                     Write (LuWr,'(17X,A)')
     &                 '  - constraint type: charge'
                  Else
                     Call WarningMessage(2,'Unknown constraint!')
                     Write(6,'(A,I10)') 'LDF_Constraint=',LDF_Constraint
                     Call LDF_Quit(-1)
                  End If
               End If
               If (LDF2) Then
                  Write (LuWr,'(17X,A,G10.2)')
     &                 '  - two-center auxiliary functions included'
     &                 //' (when needed); target accuracy: ',
     &                 Thr_Accuracy
               Else
                  Write (LuWr,'(17X,A)')
     &                 '  - two-center auxiliary functions not included'
               End If
            Else If (LDF) Then
               Write (LuWr,'(15X,A)')
     &               '   LDF decomposed two-electron'
     &               //' repulsion integrals stored Cholesky style'
               Write (LuWr,'(15X,A)') '    Concept demonstration only!'
            Else
               Write (LuWr,'(15X,A)')
     &               '   RI decomposed two-electron'
     &               //' repulsion integrals stored Cholesky style'
            End If
            If (iRI_Type.eq.1) Then
               Write (LuWr,'(17X,A)')
     &                 '  - RIJ auxiliary basis'
            Else If (iRI_Type.eq.2) Then
               Write (LuWr,'(17X,A)')
     &                 '  - RIJK auxiliary basis'
            Else If (iRI_Type.eq.3) Then
               Write (LuWr,'(17X,A)')
     &                 '  - RIC auxiliary basis'
            Else If (iRI_Type.eq.5) Then
               Write (LuWr,'(17X,A)')
     &                 '  - External RICD auxiliary basis'
            Else
               If (Do_nacCD_Basis) Then
                  Write (LuWr,'(17X,A)')
     &                    '  - nacCD auxiliary basis'
               Else
                  If (Do_acCD_Basis) Then
                     Write (LuWr,'(17X,A)')
     &                    '  - acCD auxiliary basis'
                  Else
                     Write (LuWr,'(17X,A)')
     &                       '  - aCD auxiliary basis'
                  End If
               End If
               Write (LuWr,'(17X,A,G10.2)')
     &                 '  - CD Threshold: ',Thrshld_CD
               l_aCD_Thr=.False.
               Do iCnttp = 1, nCnttp
                  l_aCD_Thr=l_aCD_Thr .or. dbsc(iCnttp)%aCD_Thr.ne.One
               End Do
               If (l_aCD_Thr) Then
                  Write (LuWr,'(17X,A)')
     &                    '     Note that the threshold for individual'
     &                  //' basis sets might be modified!'
               End If
               If (Skip_High_AC) Then
                  Write (LuWr,'(17X,A)')
     &                    '  - Skip high angular momentum combinations'
               End If
            End If
         Else
            Write (LuWr,'(15X,A)')
     &        '   Two-Electron Repulsion integrals'
         End If
      End If
      If (RMat_On) Then
         Write (LuWr,*)
         Write (LuWr,'(15X,A,A)')
     &        '   OBSERVE that some integrals are modified to',
     &        ' enable variational R-matrix calculations!'
      End If
      If (GIAO) Then
         Write (LuWr,'(15X,A)')
     &        '   GIAO integrals differentiated with respect to B'
         Write (LuWr,'(15X,A)')
     &        '     dS/dB                                        '
         Write (LuWr,'(15X,A)')
     &        '     dT/dB                                        '
         Write (LuWr,'(15X,A)')
     &        '     dV/dB                                        '
      End If
*
*     Transition moment integrals for oscillator strengths of
*     electronic transitions.
*
      If (EMFR) Then
         Write (LuWr,'(15X,A)')
     &        '   Transition moment intergrals'
         Write (LuWr,'(15X,A,3(F7.4,1X),A)')
     &                       '   The wavevector k: (',
     &   (kVector(i),i=1,3),')'
         temp=Sqrt(KVector(1)**2+KVector(2)**2+kVector(3)**2)
         temp = (Two*Pi)/temp
         Write (LuWr,'(15X,A,(F10.4,1X),A)')
     &                       '   Wavelength:        ',
     &   Temp,'a.u.'
         Write (LuWr,'(15X,A,(F10.4,1X),A)')
     &                       '                      ',
     &   Temp*Angstr,'Angstrom'
         Write (LuWr,'(15X,A,(F10.4,1X),A)')
     &                       '                      ',
     &   Temp*Angstr/Ten,'nm'
      End If
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*
      Call Qpg_cArray('SewardXTitle',Found,nTtl)
      If (Found) Then
         nTtl=nTtl/80
         Call Get_cArray('SewardXTitle',Title(1),nTtl*80)
         If (iPrint.ge.6) Then
            Write (LuWr,*)
            Write (LuWr,'(15X,88A)') ('*',i=1,88)
            Write (LuWr,'(15X,88A)') '*', (' ',i=1,86), '*'
            Do iTtl = 1, nTtl
               Write (LuWr,'(15X,A,A,A)') '*   ',Title(iTtl),'   *'
            End Do
            Write (LuWr,'(15X,88A)') '*', (' ',i=1,86), '*'
            Write (LuWr,'(15X,88A)') ('*',i=1,88)
         Else
            Write (LuWr,*)
            Write (LuWr,'(A)') ' Title:'
            Do iTtl = 1, nTtl
               Write (LuWr,'(8X,A)') Title(iTtl)
            End Do
            Write (LuWr,*)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Write (LuWr,*)
      Write (LuWr,'(19X,A,E9.2)')
     &      'Integrals are discarded if absolute value <:',ThrInt
      Write (LuWr,'(19X,A,E9.2)')
     &      'Integral cutoff threshold is set to       <:',CutInt
*                                                                      *
************************************************************************
*                                                                      *
      If (Run_Mode.eq.GS_Mode) Call Print_Symmetry()
*                                                                      *
************************************************************************
*                                                                      *
      If (nIrrep.gt.1) Then
         If (MolWgh.eq.0) Then
            Write (LuWr,*)
            Write (LuWr,'(19X,A,A)') ' Symmetry adaptation a la',
     &                          ' DCR.'
            Write (LuWr,*)
         End If
         If (MolWgh.eq.1) Then
            Write (LuWr,*)
            Write (LuWr,'(19X,A,A)') ' Symmetry adaptation a la',
     &                          ' MOLECULE.'
            Write (LuWr,*)
         End If
         If (MolWgh.eq.2) Then
            Write (LuWr,*)
            Write (LuWr,'(19X,A)') ' Unitary symmetry adaptation'
            Write (LuWr,*)
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Print Cell unit information
*
      If (Cell_l) Then
          Write(LuWr,'(6X,30(''-''))')
          Write(LuWr,'(6X,A)') '* - the centers of the Unit Cell'
          Write(LuWr,'(A)') ' '
          Write(LuWr,'(A,3I3)') 'Spread of the unit cell:',
     &      (ispread(i),i=1,3)
          Write(LuWr,'(A)') ' '
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write out basis set information
*
      If (Run_Mode.eq.GS_Mode) Then
         Call Print_Basis(lOPTO)
*                                                                      *
************************************************************************
*                                                                      *
*     Write out coordinates, bond, angles and torsional angles
*
         If (lOPTO) then
            Call Print_Geometry(1)
         else
            Call Print_Geometry(0)
         EndIf
         Call Print_Isotopes()
*                                                                      *
************************************************************************
*                                                                      *
*     Rigid Rotor analysis etc.
*
         Call RigRot(Centr,Mass,S%kCentr)
*                                                                      *
************************************************************************
*                                                                      *
         Call Print_Basis2()
*                                                                      *
************************************************************************
*                                                                      *
         Call Print_OpInfo()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
