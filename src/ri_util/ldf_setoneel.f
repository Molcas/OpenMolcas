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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetOneEl(Label)
C
C     Thomas Bondo Pedersen, February 2011
C
C     Purpose: set up data needed for evaluation of one-electron
C              integrals for operator OperatorLabel.
C              To unset, call LDF_UnsetOneEl().
C
      Implicit None
      Character*8 Label
#include "ldf_oneel.fh"

      Character*12 SecNam
      Parameter (SecNam='LDF_SetOneEl')

      If (OperatorLabel.eq.'IS_UNSET') Then
         OperatorLabel=Label
         If (OperatorLabel(1:6).eq.'Mltpl ') Then
            Call LDF_SetOneEl_Mltpl()
         Else
            Call WarningMessage(2,SecNam//': Unknown operator label')
            Write(6,'(A,A)') 'Label=',Label
            Call LDF_Quit(1)
         End If
      Else
         Call WarningMessage(2,SecNam//': info exists!')
         Write(6,'(A,A)') 'OperatorLabel=',OperatorLabel
         Write(6,'(A,A)') 'Label=',Label
         Call LDF_Quit(1)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_SetOneEl_Mltpl()
C
C     Thomas Bondo Pedersen, February 2011.
C     - based OneEl() and OneEl_() by R. Lindh.
C
C     Purpose: set up info for calculating multipole integrals.
C
      use MpmC
      use Symmetry_Info, only: iChBas
      use Sizes_of_Seward, only: S
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "rmat_option.fh"
#include "ldf_oneel.fh"

      Character*18 SecNam
      Parameter (SecNam='LDF_SetOneEl_Mltpl')

      External IrrFnc, MltLbl

      Integer iTwoj(0:7)
      Data iTwoj/1,2,4,8,16,32,64,128/

      CCoor(i,j)=Coor_MPM(i,j+1)

      ! Check that active operator is a multipole operator
      ! (OperatorLabel stored in ldf_oneel.fh)
      If (OperatorLabel(1:6).ne.'Mltpl ') Then
         Call WarningMessage(2,SecNam//': not multipole operator!')
         Write(6,'(A,A)') 'Operator=',OperatorLabel
         Call LDF_Quit(1)
      End If

      ! Get multipole order from operator label
      Read(OperatorLabel(7:8),'(I2)') iMltpl

      ! Set info (transcript of drv1el() in Seward)
      rHrmt=One
      nComp=(iMltpl+1)*(iMltpl+2)/2
      l_lOper=nComp
      Call GetMem('lOper','Allo','Inte',ip_lOper,l_lOper)
      l_kOper=nComp
      Call GetMem('kOper','Allo','Inte',ip_kOper,l_kOper)
      l_CCoor=3*nComp
      Call GetMem('CCoor','Allo','Real',ip_CCoor,l_CCoor)
      l_xZeta=S%m2Max
      Call GetMem('Zeta','Allo','Real',ip_xZeta,l_xZeta)
      l_xZI=S%m2Max
      Call GetMem('ZI','Allo','Real',ip_xZI,l_xZI)
      l_xKappa=S%m2Max
      Call GetMem('Kappa','Allo','Real',ip_xKappa,l_xKappa)
      l_xPCoor=3*S%m2Max
      Call GetMem('PCoor','Allo','Real',ip_xPCoor,l_xPCoor)
      iComp=0
      Do ix=iMltpl,0,-1
         If (Mod(ix,2).eq.0) Then
            iSymX=1
         Else
            ixyz=1
            iSymX=2**IrrFnc(ixyz)
            If (CCoor(1,iMltpl).ne.Zero) Then
               iSymX=iOr(iSymX,1)
            End If
         End If
         Do iy=iMltpl-ix,0,-1
            If (Mod(iy,2).eq.0) Then
               iSymY=1
            Else
               ixyz=2
               iSymY=2**IrrFnc(ixyz)
               If (CCoor(2,iMltpl).ne.Zero) Then
                  iSymY=iOr(iSymY,1)
               End If
            End If
            iz = iMltpl-ix-iy
            If (Mod(iz,2).eq.0) Then
               iSymZ=1
            Else
               ixyz=4
               iSymZ=2**IrrFnc(ixyz)
               If (CCoor(3,iMltpl).ne.Zero) Then
                  iSymZ = iOr(iSymZ,1)
               End If
            End If
            iChO=Mod(ix,2)*iChBas(2)
     &          +Mod(iy,2)*iChBas(3)
     &          +Mod(iz,2)*iChBas(4)
            iWork(ip_lOper+iComp)=MltLbl(iSymX,
     &                                   MltLbl(iSymY,iSymZ,nIrrep),
     &                                   nIrrep)
            iWork(ip_kOper+iComp)=iChO
            Call dCopy_(3,Coor_MPM(1,iMltpl+1),1,
     &                   Work(ip_CCoor+3*iComp),1)
            iComp=iComp+1
         End Do
      End Do
      nIC=0
      llOper=0
      Do iComp=1,nComp
         llOper=iOr(llOper,iWork(ip_lOper-1+iComp))
         Do iIrrep=0,nIrrep-1
            If (iAnd(iWork(ip_lOper-1+iComp),iTwoj(iIrrep)).ne.0) Then
               nIC=nIC+1
            End If
         End Do
      End Do
      Call SOS(iStabO,nStabO,llOper)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_UnsetOneEl(Label)
C
C     Thomas Bondo Pedersen, February 2011.
C
C     Purpose: unset data needed for evaluation of one-electron
C              integrals.
C
      Implicit None
      Character*8 Label
#include "ldf_oneel.fh"
      Character*14 SecNam
      Parameter (SecNam='LDF_UnsetOneEl')
      Character*8 ucLabel

      If (OperatorLabel.eq.'IS_UNSET') Then
#if defined (_DEBUG_)
         Call WarningMessage(0,SecNam//': nothing to do, returning')
         Call xFlush(6)
#endif
         Return
      Else
         ucLabel=Label
         Call UpCase(ucLabel)
         If (ucLabel.eq.'DUMMYLBL' .or. Label.eq.OperatorLabel) Then
            If (l_xPCoor.gt.0) Then
               Call GetMem('PCoor','Free','Real',ip_xPCoor,l_xPCoor)
               ip_xPCoor=0
               l_xPCoor=0
            End If
            If (l_xKappa.gt.0) Then
               Call GetMem('Kappa','Free','Real',ip_xKappa,l_xKappa)
               ip_xKappa=0
               l_xKappa=0
            End If
            If (l_xZI.gt.0) Then
               Call GetMem('ZI','Free','Real',ip_xZI,l_xZI)
               ip_xZI=0
               l_xZI=0
            End If
            If (l_xZeta.gt.0) Then
               Call GetMem('Zeta','Free','Real',ip_xZeta,l_xZeta)
               ip_xZeta=0
               l_xZeta=0
            End If
            If (l_CCoor.gt.0) Then
               Call GetMem('CCoor','Free','Real',ip_CCoor,l_CCoor)
               ip_CCoor=0
               l_CCoor=0
            End If
            If (l_kOper.gt.0) Then
               Call GetMem('kOper','Free','Inte',ip_kOper,l_kOper)
               ip_kOper=0
               l_kOper=0
            End If
            If (l_lOper.gt.0) Then
               Call GetMem('lOper','Free','Inte',ip_lOper,l_lOper)
               ip_lOper=0
               l_lOper=0
            End If
            Call iZero(iStabO,8)
            nStabO=0
            nIC=0
            nComp=0
            rHrmt=-9.9d9
            OperatorLabel='IS_UNSET'
         Else
            Call WarningMessage(2,SecNam//': Label mismatch!')
            Write(6,'(A,A,/,A,A)')
     &      'Label=',Label,'OperatorLabel=',OperatorLabel
            Call LDF_Quit(1)
         End If
      End If

      End
