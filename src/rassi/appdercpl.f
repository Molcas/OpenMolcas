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
      Subroutine AppDerCpl(natom,nST,ChgNuc,Prop,DerCpl,HAM)
      Implicit Real*8(A-H,O-Z)
*
*     Approximate derivative couplings:         <\Psi_I|\nabla H|\Psi_J>
*                                        f_IJ =  ----------------------
*                                                     E_J - E_I
*
*     If the wfn are real-valued: f_II = 0 ; f_JI = - f_IJ -> lower triangular storage
*
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='AppDerCpl')
#include "rassi.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "real.fh"
      Dimension ChgNuc(natom),Prop(nState,nState,NProp),
     &          DerCpl(nST,3,natom),Ham(Nstate,Nstate)
      Character*3 Label
      Save Label
      Data Label/'EF1'/
*
*
      nST = nState*(nState+1)/2
      Call FZero(DerCpl,3*natom*nST)
      Do iSta = 1, nState-1
         Ei = Ham(iSta,iSta)
         Do jSta = iSta + 1, nState
            Ej = Ham(jSta,jSta)
            iST = iSta * (jSta -1) / 2 + jSta
            Write(6,1000) iSta, jSta, Ej-Ei
            Do kProp = 1, nProp
               If (PName(kProp)(1:3) .eq. Label) Then
                  Read(PName(kProp)(5:8),'(i4)') lAt
                  DerCpl(iST,IComp(kProp),lAt) =
     &                         Prop(iSta,jSta,kProp)*ChgNuc(lAt)/(Ej-Ei)
               End If
            End Do
            SumX = Zero
            SumY = Zero
            SumZ = Zero
            Do lAt = 1, natom
               Write(6,1100) lAt,(DerCpl(iST,kXYZ,lAt),kXYZ=1,3)
               SumX = SumX + DerCpl(iST,1,lAt)
               SumY = SumY + DerCpl(iST,2,lAt)
               SumZ = SumZ + DerCpl(iST,3,lAt)
            End Do
            If(IPGLOB .ge. DEBUG) Write(6,1200) SumX,SumY,SumZ
         End Do
      End Do
1000  Format(/,' Approximate derivative couplings for states ',2i3,/,
     &       ' Energy difference = ',F15.8,/,
     &       '   Atom          X              Y              Z')
1100  Format(i7,3f15.8)
1200  Format('   Sum:',3f15.8)
*
      Return
      End
