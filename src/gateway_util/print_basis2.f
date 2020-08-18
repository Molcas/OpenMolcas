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
      Subroutine Print_Basis2()
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "angtp.fh"
#include "info.fh"
#include "print.fh"
      Logical output
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint=nPrint(iRout)
*
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.6) Then
         Write (LuWr,*)
         Call CollapseOutput(1,'   Primitive basis info:')
         Write (LuWr,'(3X,A)') '   ---------------------'
         Write (LuWr,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Generate list of primitive basis functions
*
      If (iPrint.ge.6) Then
         Write (LuWr,*)
         Write (LuWr,'(19X,A)')
     &         ' *****************************************************'
         Write (LuWr,'(19X,A)')
     &         ' ******** Primitive Basis Functions (Valence) ********'
         Write (LuWr,'(19X,A)')
     &         ' *****************************************************'
      End If
*     Loop over distinct shell types
      jExp  =0
      iPrim = 0
      iPrim_Aux = -1
      iPrim_Frag = 0
      iBas  = 0
      iBas_Aux  = -1
      iBas_Frag = 0
      iShell=0
*     Loop over basis sets
      Do iCnttp = 1, nCnttp
         mdc = mdciCnttp(iCnttp)
         output=iPrint.ge.6
         If (AuxCnttp(iCnttp).or.FragCnttp(iCnttp))
     &     output = output .and. iPrint.ge.10
     &                     .and. iCnttp.ne.iCnttp_Dummy
         If(output) then
            Write (LuWr,*)
            Write (LuWr,*)
            Write (LuWr,'(A,A)') ' Basis set:',Bsl(iCnttp)
         End If
         iShSrt = ipVal(iCnttp)
*        Loop over distinct centers
         Do icnt = 1, dbsc(iCnttp)%nCntr
            mdc = mdc + 1
            if (mdc.gt.mxdc) then
               Call WarningMessage(2,'mxdc too small')
               write(LuWr,*) 'mxdc=',mxdc
               write(LuWr,*) 'Increase mxdc in info.fh and',
     &                       ' recompile the code!'
               Call Abend()
            end if
*           Loop over shells associated with this center
*           Start with s type shells
            jSh = iShSrt
            Do iAng = 0, nVal_Shells(iCnttp)-1
               iShell = iShell + 1
               nExpj=Shells(jSh)%nExp
               nBasisj=Shells(jSh)%nBasis
               If (MaxPrm(iAng).gt.0 .and. nExpj.gt.0 .and.
     &             nBasisj.gt.0 .and. output .and.
     &             iCnt.eq.1) Then
                  Write (LuWr,*)
                  Write (LuWr,*) '                 Type         '
                  Write (LuWr,'(19X,A)')
     &                                         AngTp(iAng)
                  Write (LuWr,*) '          No.      Exponent   ',
     &                        ' Contraction Coefficients'
               End If
*
               If (nBasisj.gt.0 .and. output) Then
                  Do kExp = 1, nExpj
                     jExp  = jExp  + 1
                     If (iCnt.eq.1)
     &               Write (LuWr,'( 9X,I4,1X,D16.9,10(1X,F10.6),'//
     &                        '1X,3(/,30X,10(1X,F10.6)))')
     &                     jExp , Shells(jSh)%Exp(kExp),
     &                     ( Shells(jSh)%Cff_c(kExp,ib,2),
     &                     ib=1,nBasisj)
                  End Do
               End If
               If (iShell.gt.MxShll) Then
                  Call WarningMessage(2,'iShell.gt.MxShll')
                  Write (LuWr,*) ' Change MxShll in info.fh and re'
     &                      //'compile the code!'
                  Call Abend()
               End If
               kCmp=(iAng+1)*(iAng+2)/2
               If (Shells(jSh)%Prjct) kCmp=2*iAng+1
               If (nBasisj.ne.0 ) Then
                  If (Shells(jSh)%Aux) Then
                     iPrim_Aux = iPrim_Aux + nExpj   * kCmp
     &                         * nIrrep/nStab(mdc)
                     iBas_Aux  = iBas_Aux  + nBasisj * kCmp
     &                         * nIrrep/nStab(mdc)
                  Else If (Shells(jSh)%Frag) Then
                     iPrim_Frag = iPrim_Frag + nExpj   * kCmp
     &                          * nIrrep/nStab(mdc)
                     iBas_Frag = iBas_Frag  + nBasisj * kCmp
     &                         * nIrrep/nStab(mdc)
                  Else
                     iPrim = iPrim + nExpj   * kCmp
     &                      * nIrrep/nStab(mdc)
                     iBas  = iBas  + nBasisj * kCmp
     &                     * nIrrep/nStab(mdc)
                  End If
               End If
               jSh = jSh + 1
            End Do
         End Do
*
      End Do
      If (iBas.ge.2*MaxBfn) Then
         Call WarningMessage(2,'MaxBfn too small')
         Write (LuWr,*) 'Input: Increase 2*MaxBfn to ', iBas
         Call Abend()
      End If
      If (iPrint.ge.6) Then
         Write(LuWr,*)
         Write(LuWr,*) ' Number of primitives                ', iPrim
         Write(LuWr,*) ' Number of basis functions           ', iBas
         If (lAux .and. iPrint.ge.10) Then
             Write(LuWr,*) ' Number of primitive aux. functions  ',
     &                      iPrim_Aux
             Write(LuWr,*) ' Number of auxiliary basis functions ',
     &                       iBas_Aux
         End If
         If (lFAIEMP .and. iPrint.ge.10) Then
             Write(LuWr,*) ' Number of primitive frag. functions ',
     &                      iPrim_Frag
             Write(LuWr,*) ' Number of fragment basis functions  ',
     &                       iBas_Frag
         End If
         Write (LuWr,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (lPAM2) Then
         Write (LuWr,*)
         Write (LuWr,'(19X,A)')
     &         ' *************************************************'
         Write (LuWr,'(19X,A)')
     &         ' ******** Primitive Basis Functions (PAM) ********'
         Write (LuWr,'(19X,A)')
     &         ' *************************************************'
         Do iCnttp=1,nCnttp
         If (PAM2(iCnttp)) Then
*            If (iPrint.ge.10) Then
            Write (LuWr,*)
            Write (LuWr,*)
            Write (LuWr,'(A,A)') ' Basis set:',Bsl(iCnttp)
            If (dbsc(iCnttp)%nPAM2.ne.-1) Then
               Write (LuWr,*)
               Write (LuWr,*) 'Angular momentum of PAM operator: ',
     &                         AngTp(dbsc(iCnttp)%nPAM2)
               iAddr=1
*
               Do iAngl=0,dbsc(iCnttp)%nPAM2
                  iPrimm = INT(dbsc(iCnttp)%PAM2(iAddr)  )
                  iBass =  INT(dbsc(iCnttp)%PAM2(iAddr+1))
                  Write(LuWr,'(14H Ang. moment: ,3x,a1)') AngTp(iAngl)
                  Write(LuWr,'(22H Number of  primitive:,i4,'
     &                      //'22H Number of contracted:,i4)')
     &                         iPrimm,iBass
                  Write (LuWr,*)
                  Write(LuWr,'(21H  N.        Exponents,'
     &                      //'23H            Coefficent:)')
*
                  Do ir=0,iPrimm - 1
                     Write (LuWr,'(i4,1x,f18.12,1x,12(1x,f12.8))')
     &               ir+1,dbsc(iCnttp)%PAM2(iAddr+2+ir),
     &                   (dbsc(iCnttp)%PAM2(iAddr+2+iPrimm+ic),
     &               ic=ir,(iPrimm)*iBass-1,iPrimm)
                  End Do
*
               iAddr = iAddr + 2 + iPrimm + iPrimm*iBass
               End Do
            End If
*            End If
         End If
         End Do
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.10 .and. (lECP.or.lPP)) Then
         Write (LuWr,*)
         Write (LuWr,'(19X,A)')
     &         ' *************************************************'
         Write (LuWr,'(19X,A)')
     &         ' ******** Primitive Basis Functions (ECP) ********'
         Write (LuWr,'(19X,A)')
     &         ' *************************************************'
C     Else
Cstart Molcas
C        If (iPrint.lt.6) Write (LuWr,*)
C    &                  'To display basis set information use the key',
C    &                  ' "BSSHOW" in the input.'
C        If ((lECP.or.lPP).and.iPrint.lt.10) Write (LuWr,*)
C    &                  'To display ECP information use the key',
C    &                  ' "ECPSHOW" in the input.'
C        If (lAUX .and. iPrint.lt.10 ) Write (LuWr,*)
C    &      'To display auxiliary basis information use the key',
C    &      ' "AUXSHOW" in the input.'
Cend
      End If
*
      Do iCnttp = 1, nCnttp
*
*------- Pseudo potential type ECP
*
         If (nPP_Shells(iCnttp).ne.0 .and. iPrint.ge.10) Then
            Write (LuWr,*)
            Write (LuWr,*)
            Write (LuWr,'(A,A)') ' Basis set:',Bsl(iCnttp)
            Write (LuWr,*)
            Write (LuWr,'(A)') ' Pseudo Potential'
            Write (LuWr,*)
            kShStr=ipPP(iCnttp)
            kShEnd = kShStr + nPP_Shells(iCnttp)-1
            lSh= 0
            Do kSh = kShStr, kShEnd
               nExpk=Shells(kSh)%nExp
               If (nExpk.ne.0.and.iCnttp.le.21) Then
                  If (lSh.eq.0) Then
                     Write (LuWr,'(4X,A)') '  H Potential'
                  Else
                     Write (LuWr,'(4X,A)') AngTp(lSh-1)//
     &                      '-H Potential'
                  End If
               End If
               lSh = lSh + 1
               Write (LuWr,'(A)') '  n     Exponent      Coefficient'
               iOff = 1
               Do iExp = 1, nExpk
                  ncr=Int(Shells(kSh)%Exp(iOff  ))
                  zcr=    Shells(kSh)%Exp(iOff+1)
                  ccr=    Shells(kSh)%Exp(iOff+2)
                  Write (LuWr,'(2x,I1,3X,2F15.10)') ncr,zcr,ccr
                  iOff = iOff + 3
               End Do
               Write (LuWr,*)
*
            End Do
         End If
*
*------- Huzinaga type ECP
*
         If (ECP(iCnttp)) Then
            If (iPrint.ge.10) Then
               Write (LuWr,*)
               Write (LuWr,*)
               Write (LuWr,'(A,A)') ' Basis set:',Bsl(iCnttp)
*
               If (dbsc(iCnttp)%nM1.ne.0) Then
                  Write (LuWr,*)
                  Write (LuWr,*) ' M1 operator       Exponent   ',
     &                           ' Contraction Coefficients'
                  Do irow = 1, dbsc(iCnttp)%nM1
                     Write (LuWr,'(14X,D16.9,1X,D19.9)')
     &                  dbsc(iCnttp)%M1xp(irow),
     &                  dbsc(iCnttp)%M1cf(irow)
                  End Do
               End If ! If (dbsc(iCnttp)%nM1.ne.0) Then
*
               If (dbsc(iCnttp)%nM2.ne.0) Then
                  Write (LuWr,*)
                  Write (LuWr,*) ' M2 operator       Exponent   ',
     &                           ' Contraction Coefficients'
                  Do irow = 1, dbsc(iCnttp)%nM2
                     Write (LuWr,'(14X,D16.9,1X,D19.9)')
     &                  dbsc(iCnttp)%M2xp(irow),
     &                  dbsc(iCnttp)%M2cf(irow)
                  End Do
               End If ! If (dbsc(iCnttp)%nM2.ne.0) Then
            End If ! If (iPrint.ge.10) Then
*
*--------------Projection Basis Set
*
            iSh = ipPrj(iCnttp)
            nSumB = 0
            jSh = iSh
            Do iAng = 0, nPrj_Shells(iCnttp)-1
               nSumB = nSumB + Shells(jSh)%nBasis
               jSh = jSh + 1
            End Do
            If (nSumB.ne.0.and.iPrint.ge.10) Then
               Write (LuWr,*)
               Write (LuWr,*)
               Write (LuWr,*) ' Proj. Operator'
            End If
            Do iAng = 0, nPrj_Shells(iCnttp)-1
               If (Shells(iSh)%nBk.ne.0) Then
                  If (iPrint.ge.10) Then
                     Write (LuWr,*)
                     Write (LuWr,'(19X,A,A)')
     &                     '        Angular Type: ', AngTp(iAng)
                     Write (LuWr,*) '                   Exponent   ',
     &                              ' Contraction Coefficients'
                     Write (LuWr,*)
                     Write (LuWr,
     &                  '(A,18X,8(G12.5),/,5(32X,8(G12.5),/))')
     &                  '     Bk-values',
     &                  (Shells(iSh)%Bk(i),i=1,Shells(iSh)%nBk)
                     Write (LuWr,
     &                  '(A,18X,8(G12.5),/,5(32X,8(G12.5),/))')
     &                  '     Frac.Occ.',
     &                  (Shells(iSh)%Occ(i),i=1,Shells(iSh)%nBk)
                  End If
*
                  Do i=1,Shells(iSh)%nBk
                     Shells(iSh)%Bk(i)=Shells(iSh)%Bk(i)
     &                                *Shells(iSh)%Occ(i)
                  End Do
*
                  If (iPrint.ge.10) Then
                     Do kExp = 1, Shells(iSh)%nExp
                        jExp  = jExp  + 1
                        Write (LuWr,'(14X,D16.9,8(G12.5),'//
     &                        '3(/,30X,8(G12.5)))')
     &                          Shells(ish)%Exp(kExp),
     &                        ( Shells(ish)%Cff_c(kExp,ib,2),
     &                     ib=1,Shells(iSh)%nBk)
                     End Do
                  End If ! If (iPrint.ge.10) Then
               End If ! If (Shells(iSh)%nBk.ne.0) Then
               iSh = iSh + 1
            End Do    ! iAng

*
*--------------Auxilliary core basis
*
            If (iPrint.ge.10) Then
               iSh = ipSOC(iCnttp)
               nSumB = 0
               jSh = iSh
               Do iAng = 0, nSOC_Shells(iCnttp)-1
                  nSumB = nSumB + Shells(jSh)%nBasis
                  jSh = jSh + 1
               End Do
               If (nSumB.ne.0.and.iPrint.ge.10) Then
                  Write (LuWr,*)
                  Write (LuWr,*)
                  Write (LuWr,*) ' SOC Basis'
               End If
               Do iAng = 0, nSOC_Shells(iCnttp)-1
                  If (Shells(iSh)%nBasis.ne.0) Then
                     Write (LuWr,*)
                     Write (LuWr,'(19X,A,A)')
     &                     '        Angular Type: ', AngTp(iAng)
                     Write (LuWr,*) '                   Exponent   ',
     &                              ' Contraction Coefficients'
                     Write (LuWr,*)
                     Do kExp = 1, Shells(iSh)%nExp
                        jExp  = jExp  + 1
                        Write (LuWr,'(14X,D16.9,10(1X,F10.6),'//
     &                           '3(/,30X,10(1X,F10.6)))')
     &                        Shells(iSh)%Exp(kExp),
     &                      ( Shells(iSh)%Cff_c(kExp,ib,1),
     &                        ib=1,Shells(iSh)%nBasis)
                     End Do
                  End If
                  iSh = iSh + 1
               End Do
            End If ! If (iPrint.ge.10) Then
*
*-----------Spectral Resolution Basis Set
*
            If (iPrint.ge.10) Then
               iSh = ipSRO(iCnttp)
               nSumA = 0
               jSh = iSh
               Do iAng = 0, nSRO_Shells(iCnttp)-1
                  nSumA = nSumA + Shells(jSh)%nExp
                  jSh = jSh + 1
               End Do
               If (nSumA.ne.0) Then
                  Write (LuWr,*)
                  Write (LuWr,*)
                  Write (LuWr,*) ' Spectral Resolution Basis Set'
               End If
               Do iAng = 0, nSRO_Shells(iCnttp)-1
                  nExpi=Shells(iSh)%nExp
                  If (nExpi.ne.0) Then
                     Write (LuWr,*)
                     Write (LuWr,'(19X,A,A)')
     &                     '        Angular Type: ', AngTp(iAng)
                     Call RecPrt(' Exponents',' ',
     &                           Shells(iSh)%Exp,nExpi,1)
                     If (iPrint.ge.11) Then
                        Call RecPrt(' The Akl matrix','(5D20.13)',
     &                              Shells(iSh)%Akl(1,1,1),nExpi,
     &                                                     nExpi)
                        Call RecPrt(' The Adl matrix','(5D20.13)',
     &                              Shells(iSh)%Akl(1,1,2),nExpi,
     &                                                     nExpi)
                     End If
                  End If
                  iSh = iSh + 1
               End Do ! iAng
            End If  ! If (iPrint.ge.10) Then
         End If ! If (ECP(iCnttp)) Then
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.6) Then
         Call CollapseOutput(0,'   Primitive basis info:')
         Write (LuWr,*)
      End If
*
      Return
      End
