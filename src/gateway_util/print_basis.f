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
      SubRoutine Print_Basis(lOPTO)
************************************************************************
*                                                                      *
*     Object: to print the basis set                                   *
*                                                                      *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
*             September 2006                                           *
************************************************************************
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "relmp.fh"
#include "real.fh"
#include "print.fh"
      Character DBas*4
      Character ChCo*1, ChCa*1, ChSph*1
      Logical Output, Type(0:7), lOPTO
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint = nPrint(iRout)
      If (iPrint.eq.0) Return
      LuWr = 6
*                                                                      *
************************************************************************
*                                                                      *
      If (Show) Then
         Write (LuWr,*)
         Call CollapseOutput(1,'   Basis set information:')
         Write (LuWr,'(3X,A)') '   ----------------------'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write out basis set information
*                                                                      *
************************************************************************
*                                                                      *
*     Valence basis set
*
      Do iCnttp = 1, nCnttp
         mdc = dbsc(iCnttp)%mdci
         lSh= 0
         output=Show
         If (dbsc(iCnttp)%Aux.or.dbsc(iCnttp)%Frag)
     &     output = output .and. iPrint.ge.10
     &                     .and. iCnttp.ne.iCnttp_Dummy
         If (output) Then
            Write (LuWr,*)
            Write (LuWr,*)
            Write (LuWr,'(6X,A,1X,A)') 'Basis set label:',
     *        Trim(Bsl(iCnttp))
            If (lOPTO) GoTo 100
            Write (LuWr,*)
            dbas=LblCnt(mdc+1)(1:4)
            Call Upcase(dbas)
            If (dbas.eq.'DBAS') Then
               Write (LuWr,'(6X,A)') 'Diffuse basis set for R-matrix:'
               Write (LuWr,'(6X,A)') '==============================='
               If (dbsc(iCnttp)%nCntr.ne.1) Then
                  Call WarningMessage(2,
     &                        'Too many centers, should only be one!')
                  Call Quit_OnUserError()
               End If
            Else
               If (dbsc(iCnttp)%Aux) Then
                  Write (LuWr,'(6X,A)') 'Auxiliary basis set:'
                  Write (LuWr,'(6X,A)') '=================='
                  If (aCD_Thr(iCnttp).ne.One) Then
                     Write (LuWr,'(6X,A,G9.2)') 'Threshold in the '
     &                    //'auxiliary basis set generation is '
     &                    //'modified to ',
     &                    aCD_Thr(iCnttp)*Thrshld_CD
                  End If
               Else If (dbsc(iCnttp)%Frag) Then
                  Write (LuWr,'(6X,A)') 'Fragment basis set:'
                  Write (LuWr,'(6X,A)') '=================='
               Else
                  If (fmass(iCnttp).eq.One) Then
                     Write (LuWr,'(6X,A)')
     &                      'Electronic valence basis set:'
                     Write (LuWr,'(6X,A)') '------------------'
                  Else
                     Write (LuWr,'(6X,A)') 'Muonic valence basis set:'
                     Write (LuWr,'(6X,A)') '------------------'
                  End If
               End If
            End If
            If (Fixed(iCnttp)) Write (LuWr,'(6X,A)')
     &          'Centers of this basis set are frozen!'
            If (dbsc(iCnttp)%IsMM.eq.1) Then
               Write (LuWr,'(6X,A)') 'This is a MM atom: no basis set'
            Else
               If (pChrg(iCnttp)) Then
                  Write (LuWr,'(6X,A,F10.6,A)')
     &                'Associated Effective Charge ',
     &               Charge(iCnttp), ' au (this is a pseudo charge)'
               Else
                  Write (LuWr,'(6X,A,F10.6,A)')
     &                'Associated Effective Charge ',
     &               Charge(iCnttp), ' au'
               End If
               Write (LuWr,'(6X,A,F10.6,A)')
     &               'Associated Actual Charge    ',
     &               Max(Zero,DBLE(iAtmNr(iCnttp))), ' au'
*
               If (Nuclear_Model.eq.Point_Charge) Then
                  Write (LuWr,'(6X,A)') 'Nuclear Model: Point charge'
               Else If (Nuclear_Model.eq.Gaussian_type) Then
                  Write (LuWr,'(6X,A)')
     &                  'Nuclear Model: Finite nucleus - '
     &                             //'Gaussian distribution'
                  Write (LuWr,'(6X,A,E12.5)')
     &               '  Gaussian exponent, Xi/bohr**(-2): ',
     &               ExpNuc(iCnttp)
               Else If (Nuclear_Model.eq.mGaussian_type) Then
                  Write (LuWr,'(6X,A)')
     &                'Nuclear Model: Finite nucleus - '
     &                             //'Modified Gaussian distribution'
                  Write (LuWr,'(6X,A,E12.5,A,E12.5)')
     &               '  Parameters, Xi/bohr**(-2), w/bohr**(-2): ',
     &               ExpNuc(iCnttp),', ',w_mGauss(iCnttp)
               Else
                  Call WarningMessage(2,'Illegal Nuclear Model!')
                  Call Abend()
               End If
            End If
            Write (LuWr,*)
  100    Continue
*
         End If
         kShStr = dbsc(iCnttp)%iVal
         kShEnd = dbsc(iCnttp)%iVal+dbsc(iCnttp)%nVal-1
         Type(0) = .False.
         Do kSh = kShStr, kShEnd
            nExpk=Shells(kSh)%nExp
            Type(0) = Type(0) .or. nExpk*Shells(kSh)%nBasis.ne.0
         End Do
         If (output.and.Type(0) .AND..NOT.lOPTO) Then
            Write (LuWr,'(6X,A)')
     &          'Shell  nPrim  nBasis  Cartesian Spherical Contaminant'
         End If
         Do kSh = kShStr, kShEnd
            nExpk=Shells(kSh)%nExp
            ChCa=' '
            ChSph='X'
            ChCo=' '
            If (.not.Shells(kSh)%Transf) Then
               ChCa='X'
               ChSph=' '
            End if
            If (Shells(kSh)%Transf .and.
     &    .not. Shells(kSh)%Prjct ) ChCo='X'
            If (output.and.nExpk*Shells(kSh)%nBasis.ne.0.AND..NOT.lOPTO)
     &         Write (LuWr,'(9X,A,5X,I3,5X,I3,8X,A,8X,A,8X,A)')
     &            AngTp(lSh),nExpk,Shells(kSh)%nBasis,ChCa,ChSph,ChCo
            If (Shells(kSh)%Prjct) Then
               kComp = 2*lSh + 1
            Else
               kComp = (lSh+1)*(lSh+2)/2
            End If
            lSh = lSh + 1
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*        Process PP part, if any.
*
         kShStr=ipPP(iCnttp)
         kShEnd = kShStr + nPP_Shells(iCnttp)-1
         If (output.and.nPP_Shells(iCnttp).ne.0 .AND..NOT.lOPTO) Then
            Write (LuWr,*)
            Write (LuWr,'(6X,A)')
     &             'Pseudo Potential specification:'
            Write (LuWr,'(6X,A)')
     &             '======================================='
            Type(0)=.False.
            Do kSh = kShStr, kShEnd
               nExpk=Shells(kSh)%nExp/3
               Type(0)=Type(0).or.nExpk.ne.0
            End Do
            If (Type(0)) Then
               Write (LuWr,*)
               Write (LuWr,'(6X,A)') 'Potential  nTerms    '
            End If
         End If
         lSh= 0
         Do kSh = kShStr, kShEnd
            nExpk=Shells(kSh)%nExp/3
C           Write (*,*) 'kSh,lSh=',kSh,lSh
            If (output.and.nExpk.ne.0 .AND..NOT.lOPTO) Then
               If (lSh.eq.0) Then
                  Write (LuWr,'(9X,A,6X,I2)')
     &               '  H',nExpk
               Else
                  Write (LuWr,'(9X,A,6X,I2)')
     &               AngTp(lSh-1)//'-H',nExpk
               End If
            End If
            lSh = lSh + 1
         End Do
*                                                                      *
************************************************************************
*                                                                      *
*        Process ECP part, if any.
*
         kShStr = dbsc(iCnttp)%iPrj
         kShEnd = kShStr + dbsc(iCnttp)%nPrj-1
         If (output.and.dbsc(iCnttp)%ECP .AND..NOT.lOPTO) Then
            Write (LuWr,*)
            Write (LuWr,'(6X,A)')
     &             'Effective Core Potential specification:'
            Write (LuWr,'(6X,A)')
     &             '======================================='
            If (dbsc(iCnttp)%nM1.gt.0) Then
               Write (LuWr,*)
               Write (LuWr,'(6X,A,I5)')
     &               ' Number of M1 terms:',dbsc(iCnttp)%nM1
            End If
            If (dbsc(iCnttp)%nM2.gt.0) Then
               Write (LuWr,*)
               Write (LuWr,'(6X,A,I5)')
     &               ' Number of M2 terms:',dbsc(iCnttp)%nM2
            End If
            Type(0)=.False.
            Do kSh = kShStr, kShEnd
               nExpk=Shells(kSh)%nExp
               nBasisk=Shells(kSh)%nBasis
               Type(0)=Type(0).or.nExpk*nBasisk.ne.0
            End Do
            If (Type(0)) Then
               Write (LuWr,*)
               Write (LuWr,'(6X,A)') 'Projection basis set '
               Write (LuWr,'(6X,A)') 'Shell  nPrim  nBasis '
            End If
         End If
         lSh= 0
         Do kSh = kShStr, kShEnd
            nExpk=Shells(kSh)%nExp
            nBasisk=Shells(kSh)%nBasis
            If (output.and.nExpk*nBasisk.ne.0 .AND..NOT.lOPTO)
     &         Write (LuWr,'(9X,A,6X,I2,6X,I2)')
     &            AngTp(lSh),nExpk,nBasisk
            kComp = 2*lSh + 1
            lSh = lSh + 1
         End Do
*--------Spectral resolution basis set
         kShStr = dbsc(iCnttp)%iSRO
         kShEnd = kShStr + dbsc(iCnttp)%nSRO-1
         If (output.and.dbsc(iCnttp)%ECP .AND..NOT.lOPTO) Then
            Type(0)=.False.
            Do kSh = kShStr, kShEnd
               nExpk=Shells(kSh)%nExp
               Type(0)=Type(0).or.nExpk.ne.0
            End Do
            If (Type(0)) Then
               If (dbsc(iCnttp)%nOpt.ne.0) Then
                  Write (LuWr,*)
                  Write (LuWr,'(6X,A)') 'Spectral Resolvent Operators :'
                  If (iAnd(2**0,dbsc(iCnttp)%nOpt).ne.0)
     &               Write (LuWr,'(8X,A)') ' Exchange'
                  If (iAnd(2**1,dbsc(iCnttp)%nOpt).ne.0)
     &               Write (LuWr,'(8X,A)') ' Mass-Velocity'
                  If (iAnd(2**2,dbsc(iCnttp)%nOpt).ne.0)
     &               Write (LuWr,'(8X,A)') ' Darwin 1-electron'
     &                      //' contact term'
                  If (iAnd(2**3,dbsc(iCnttp)%nOpt).ne.0) Then
                     If (IRELMP.EQ.0) Then
                        Write (LuWr,'(8X,A)') ' No-Pair approximation'
                     Else If (IRELMP.EQ.1) Then
                        Write (LuWr,'(8X,A)')
     &                        ' No-Pair approximation (DK1)'
                     Else If (IRELMP.EQ.2) Then
                        Write (LuWr,'(8X,A)')
     &                        ' No-Pair approximation (DK2)'
                     Else If (IRELMP.EQ.3) Then
                        Write (LuWr,'(8X,A)')
     &                        ' No-Pair approximation (DK3)'
                     Else If (IRELMP.EQ.4) Then
                        Write (LuWr,'(8X,A)')
     &                        ' No-Pair approximation (DK3)'
                     Else If (IRELMP.EQ.11) Then
                        Write (LuWr,'(8X,A)')
     &                        ' RESC approximation'
                     Else If (IRELMP.EQ.21) Then
                        Write (LuWr,'(8X,A)')
     &                        ' ZORA approximation'
                     Else If (IRELMP.EQ.22) Then
                        Write (LuWr,'(8X,A)')
     &                        ' ZORA-FP approximation'
                     Else If (IRELMP.EQ.23) Then
                        Write (LuWr,'(8X,A)')
     &                        ' IORA approximation'
                     End If
                  End If
               End If
               Write (LuWr,*)
               Write (LuWr,'(6X,A)') 'Spectral Resolvent basis set '
               Write (LuWr,'(6X,A)') 'Shell  nPrim '
            End If
         End If
         lSh= 0
         Do kSh = kShStr, kShEnd
            nExpk=Shells(kSh)%nExp
            If (output.and.nExpk.ne.0 .AND..NOT.lOPTO)
     &         Write (LuWr,'(9X,A,6X,I2)')
     &            AngTp(lSh),nExpk
            kComp = 2*lSh + 1
            lSh = lSh + 1
         End Do
*
*--------Auxilliary SO core
*
         kShStr = ipSOC(iCnttp)
         kShEnd = kShStr + nSOC_Shells(iCnttp)-1
         If (output.and.dbsc(iCnttp)%ECP .AND..NOT.lOPTO) Then
            Type(0)=.False.
            Do kSh = kShStr, kShEnd
               nExpk=Shells(kSh)%nExp
               Type(0)=Type(0).or.nExpk.ne.0
            End Do
            If (Type(0)) Then
               If (dbsc(iCnttp)%nOpt.ne.0) Then
                  Write (LuWr,*)
                  Write (LuWr,'(6X,A)') 'Auxilliary core basis'
               End If
               Write (LuWr,*)
               Write (LuWr,'(6X,A)') 'SOC basis set '
               Write (LuWr,'(6X,A)') 'Shell  nPrim '
            End If
         End If
         lSh= 0
         Do kSh = kShStr, kShEnd
            nExpk=Shells(kSh)%nExp
            If (output.and.nExpk.ne.0)
     &         Write (LuWr,'(9X,A,6X,I2)')
     &            AngTp(lSh),nExpk
            kComp = 2*lSh + 1
            lSh = lSh + 1
         End Do
*
         If (output.and.iPrint.ge.6) Then
            Write (LuWr,*)
            Write (LuWr,'(6X,A)')
     &              ' Label   Cartesian Coordinates / Bohr'
            Write (LuWr,*)
            Do iCnt = 1, dbsc(iCnttp)%nCntr
               Call Write_LblCnt(LuWr,LblCnt(mdc+iCnt),
     &                           dbsc(iCnttp)%Coor(1,iCnt))
            End Do
         End If
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (Show) Then
         Call CollapseOutput(0,'   Basis set information:')
         Write (LuWr,*)
      End If
*
      Return
      End
