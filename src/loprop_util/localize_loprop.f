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
      Subroutine Localize_LoProp(Ttot,Ttot_Inv,nBas,SMatrix,iCenter,
     &                           iType)
*                                                                      *
************************************************************************
*                                                                      *
* H2 molecule with 2s basis functions on each H: A and B are the two
* centers
* Original S matrix  (nz=non zero; 10-8 should be zero but it is not for
* numerical imprecision
*       1     10-8  nz   nz
* S0=   10-8  1     nz   nz
*       nz    nz    1    10-8
*       nz    nz    10-8 1
*
* Step 1. LO S0 in order to generate S1 (orthogonalization of all the
* blocks for the same center
*       1     0     nz   nz
* S1=   0     1     nz   nz
*       nz    nz    1    0
*       nz    nz    0    1
*
* Step 2. LO S1 -> S2 (orthogonalization of the OAOB and VAVB blocks
*       1     0     0    nz
* S2=   0     1     nz   0
*       0     nz    1    0
*       nz    0     0    1
*
* Step 3. GS S2 -> S3 (orthogonalization of the VBOA and VAOB blocks
*       1     0     0    0
* S3=   0     1     0    0
*       0     0     1    0
*       0     0     0    1
*                                                                      *
************************************************************************
*                                                                      *
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
      Real*8 Ttot(nBas,nBas), Ttot_Inv(nBas,nBas), SMatrix(nBas,nBas)
      Integer iCenter(nBas), iType(nBas)
      Integer IndType(7), Occ, Vir
      Parameter(Occ=1,Vir=0)
      Character OrbName*128, Note*80, Filename*6
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate temporary memory
*
      Call Allocate_Work(ip_T1,nBas**2)
      Call Allocate_Work(ip_T2,nBas**2)
      Call Allocate_Work(ip_T3,nBas**2)
      Call Allocate_Work(ip_T4,nBas**2)
      Call Allocate_Work(ip_tmp,nBas**2)
      Call Allocate_Work(ip_S,nBas**2)
      Call Allocate_Work(ip_Save,nBas**2)
*
*     Save S because GS will destroy it!
*
      call dcopy_(nBas**2,SMatrix,1,Work(ip_S),1)
*                                                                      *
************************************************************************
*                                                                      *
*     Step 1. GS S0 ->S1
*
      call dcopy_(nBas**2,[Zero],0,Work(ip_T1),1)
      call dcopy_(nBas,[One],0,Work(ip_T1),nBas+1)
*
      Call Step1(iCenter,Work(ip_S),nBas,Work(ip_T1),iType,
     &           SMatrix,Work(ip_tmp))
*                                                                      *
************************************************************************
*                                                                      *
*     Step 2. LO S1 ->S2
*
      call dcopy_(nBas**2,Work(ip_S),1,Work(ip_Save),1)
      call dcopy_(nBas**2,[Zero],0,Work(ip_T2),1)
      call dcopy_(nBas,[One],0,Work(ip_T2),nBas+1)
*
      Call Step2(iCenter,Work(ip_S),nBas,Work(ip_T2),iType,
     &           Work(ip_Save),Work(ip_tmp))
*                                                                      *
************************************************************************
*                                                                      *
*     Now set up things for the GS orthogonalization
*
*
*     Step 3. GS S2 ->S3
*
      Call Step3(iCenter,Work(ip_S), nBas,Work(ip_T3),iType)

*                                                                      *
************************************************************************
*                                                                      *
*     Now do a final Lowdin to remove virtual-virtual
*
*
*     Step 4. LO S3 ->S4
*
      Call Step4(Work(ip_S),nBas,Work(ip_T4),iType)

*                                                                      *
************************************************************************
*                                                                      *
*     Now I have LO orthog twice and GS orthog
*     The corresponding transformation matrices are T1 T2 and TGS
*     Now I generate the total transformation matrix TTOT=T1*T2*TGS
*
*        ...now T is T1*T2*T3*T4
*
      Call Ttotal(Work(ip_T1),Work(ip_T2),Work(ip_T3), Work(ip_T4),
     &            Ttot,Ttot_Inv,nBas)
*                                                                      *
************************************************************************
*                                                                      *
*      Check that all this works by applaying TTot to original overlap
*      matrix and you should get unit matrix
*
C     Call RecPrt('Old S',' ',SMatrix,nBas,nBas)
C     Call xxDGeMul(SMatrix,nBas,'N',
C    &            Ttot,nBas,'N',
C    &            Work(ip_tmp),nBas,
C    &            nBas,nBas,nBas)
C     Call xxDGeMul(Ttot,nBas,'T',
C    &            Work(ip_tmp),nBas,'N',
C    &            SMatrix,nBas,
C    &            nBas,nBas,nBas)
C     Call RecPrt('New S',' ',SMatrix,nBas,nBas)
*                                                                      *
************************************************************************
*                                                                      *
*     Dealloctate memory
*
      Call Free_Work(ip_Save)
      Call Free_Work(ip_S)
      Call Free_Work(ip_tmp)
      Call Free_Work(ip_T4)
      Call Free_Work(ip_T3)
      Call Free_Work(ip_T2)
      Call Free_Work(ip_T1)
*                                                                      *
************************************************************************
*                                                                      *
*     Write out the transformation matrix as molecular orbitals.
*
      nOcc=0
      Do iBas = 1, nBas
         If (iType(iBas).eq.Occ) nOcc = nOcc + 1
      End Do
*
      IndType(1)=0
      IndType(2)=nOcc
      IndType(3)=0
      IndType(4)=0
      IndType(5)=0
      IndType(6)=nBas-nOcc
      IndType(7)=0
*
      Call Allocate_Work(ipE,nBas)
      Call FZero(Work(ipE),nBas)
      OrbName='LPRORB'
      LuOut=20
      iUHF=0
      nSym=1
      Note='LoProp localized orbitals'
      Call WrVec_(OrbName,LuOut,'COEI',iUHF,nSym,[nBas],[nBas],
     &            TTot,[0.0d0],
     &            Work(ipE),[0.0d0],
     &            Work(ipE),[0.0d0],IndType,Note,0)
      Call Free_Work(ipE)
*                                                                      *
************************************************************************
*                                                                      *
*     Not implemented for symmetry!
*
      Filename='MD_LPR'
      Call Get_iScalar('nSym',nSym)
      If (nSym.eq.1) Call Molden_Interface(iUHF,OrbName,Filename,
     &                                     .False.)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
