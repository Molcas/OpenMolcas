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
      Subroutine Misc_Seward(iBas,iBas_Aux,iBas_Frag,DInf,nDInf)
      use Basis_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
      Real*8 DInf(nDInf)
*                                                                      *
************************************************************************
*                                                                      *
*---- Statement Function
*
      IndSOff(iCnttp,iCnt)=(iCnttp-1)*Max_Cnt+iCnt
*                                                                      *
************************************************************************
*                                                                      *
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over distinct shell types
*
      iBas  = 0
      iBas_Aux  = 0
      iBas_Frag = 0
*
      IndShl=0
      iShell=0
      mc = 1
*     Loop over basis sets
      iCnttp = 0
      Do jCnttp = 1, nCnttp
*
*        Make sure that we process the dummy shell last
*
         If (jCnttp.eq.iCnttp_Dummy .and. jCnttp.ne.nCnttp) Then
            iCnttp = iCnttp + 2
         Else If (jCnttp.eq.nCnttp .and. iCnttp.eq.jCnttp) Then
            iCnttp = iCnttp_Dummy
         Else
            iCnttp = iCnttp + 1
         End If
*
*        Loop over distinct centers
*
         Do icnt = 1, dbsc(iCnttp)%nCntr
            If (IndSOff(iCnttp,iCnt).gt.MxShll) Then
               Call WarningMessage(2,'MxShll too small:')
               write(LuWr,*) 'MxShll=',MxShll
               write(LuWr,*) 'Increase MxShll in info.fh and',
     &                       ' recompile the code!'
            End If
            Ind_Shell(IndSOff(iCnttp,iCnt)) = iShell
            mdc = iCnt + mdciCnttp(iCnttp)
            if(mdc.gt.mxdc) then
               Call WarningMessage(2,'mxdc too small:')
               write(LuWr,*) 'mxdc=',mxdc
               write(LuWr,*) 'Increase mxdc in info.fh and',
     &                       ' recompile the code!'
               Call Abend()
            end if
*           Loop over shells associated with this center
*           Start with s type shells
            jSh = ipVal(iCnttp)
            Do iAng = 0, nVal_Shells(iCnttp)-1
               iShell = iShell + 1
*              Pointer to the untouched contraction matrix as after input.
               iCff = ipCff(jSh)+nExp(jSh)*nBasis(jSh)
*
               If (nBasis_Cntrct(jSh).gt.0 )
     &            Call RdMx(RadMax,Shells(jSh)%Exp,nExp(jSh),
     &                      DInf(ipCff_Cntrct(jSh)),nBasis_Cntrct(jSh),
     &                      cdMax,EtMax)
               If (iShell.gt.MxShll) Then
                  Call WarningMessage(2,'iShell.gt.MxShll;'
     &                    //' Change MxShll in info.fh and re'
     &                    //'compile the code!')
                  Call Abend()
               End If
               IndS(iShell) = IndShl
               kCmp=(iAng+1)*(iAng+2)/2
               If (Prjct(jSh)) kCmp=2*iAng+1
               IndShl = IndShl + kCmp
*
               If (nBasis(jSh).ne.0 ) Then
                  If (AuxShell(jSh)) Then
                     iBas_Aux  = iBas_Aux  + nBasis(jSh) * kCmp
     &                         * nIrrep/nStab(mdc)
                  Else If (FragShell(jSh)) Then
                     iBas_Frag = iBas_Frag  + nBasis(jSh) * kCmp
     &                         * nIrrep/nStab(mdc)
                  Else
                     iBas  = iBas  + nBasis(jSh) * kCmp
     &                     * nIrrep/nStab(mdc)
                  End If
               End If
               jSh = jSh + 1
            End Do
            mc = mc + nIrrep/nStab(mdc)
         End Do
         nShlls = iShell
*
      End Do
      If (iBas.ge.2*MaxBfn) Then
         Call WarningMessage(2,'MaxBfn too small')
         Write (LuWr,*) 'Increase 2*MaxBfn to ', iBas
         Call Abend()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
