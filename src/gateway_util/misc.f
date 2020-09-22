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
      Subroutine Misc_Seward(iBas,iBas_Aux,iBas_Frag)
      use Basis_Info
      use Center_Info
      use Sizes_of_Seward, only: S
      use Real_Info, only: RadMax, cdMax, EtMax
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
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
      iShell=0
      mc = 1
      kdc = 0
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
         Do iCnt = 1, dbsc(iCnttp)%nCntr
            kdc = kdc + 1
            mdc = iCnt + dbsc(iCnttp)%mdci
            if(Max(mdc,kdc).gt.MxAtom) then
               Call WarningMessage(2,'MxAtom too small:')
               write(LuWr,*) 'MxAtom=',MxAtom
               write(LuWr,*) 'Increase mxAtom in Molcas.fh and',
     &                       ' recompile the code!'
               Call Abend()
            end if
*           Loop over shells associated with this center
*           Start with s type shells
            jSh = dbsc(iCnttp)%iVal
            Do iAng = 0, dbsc(iCnttp)%nVal-1
               iShell = iShell + 1
*
               If (Shells(jSh)%nBasis_C.gt.0 )
     &            Call RdMx(RadMax,Shells(jSh)%Exp,
     &                             Shells(jSh)%nExp,
     &                      Shells(jSh)%Cff_c(1,1,1),
     &                      Shells(jSh)%nBasis_C,
     &                      cdMax,EtMax)
               If (iShell.gt.MxShll) Then
                  Call WarningMessage(2,'iShell.gt.MxShll;'
     &                    //' Change MxShll in info.fh and re'
     &                    //'compile the code!')
                  Call Abend()
               End If
               kCmp=(iAng+1)*(iAng+2)/2
               If (Shells(jSh)%Prjct) kCmp=2*iAng+1
*
               If (Shells(jSh)%nBasis.ne.0 ) Then
                  If (Shells(jSh)%Aux) Then
                     iBas_Aux  = iBas_Aux  + Shells(jSh)%nBasis * kCmp
     &                         * nIrrep/dc(mdc)%nStab
                  Else If (Shells(jSh)%Frag) Then
                     iBas_Frag = iBas_Frag  + Shells(jSh)%nBasis * kCmp
     &                         * nIrrep/dc(mdc)%nStab
                  Else
                     iBas  = iBas  + Shells(jSh)%nBasis * kCmp
     &                     * nIrrep/dc(mdc)%nStab
                  End If
               End If
               jSh = jSh + 1
            End Do
            mc = mc + nIrrep/dc(mdc)%nStab
         End Do
*
      End Do
      S%nShlls = iShell
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
