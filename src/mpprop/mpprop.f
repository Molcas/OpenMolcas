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
      Subroutine MpProp(iReturn)
      Implicit Real*8 (a-h,o-z)
*                                                                      *
************************************************************************
*                                                                      *
*     - Include files
#include "MpParam.fh"
#include "WrkSpc.fh"
#include "Address.fh"
#include "MolProp.fh"
C
C Variable definition
C
C     - Molcas part
      Integer nPrim(8), nBas(8)
      Character*8 Label
      Character*80 VTitle
      Character*6 FName
      Character*8 MemLabel
      Logical LNearestAtom
      Logical LAllCenters
      Logical AveOrb, Diffuse(3)
      Logical LFirstRun
      Logical LLumOrb
      Logical Exist
      Dimension dLimmo(2),iDum(1)
*                                                                      *
************************************************************************
*                                                                      *
*     Set some zeroes
      Do i=1,mxAtomMP
       iAtomPar(i)=1
       nub(i)=0
       Do j=1,mxAtomMP
        nbi(i,j)=0
       EndDo
      EndDo
      Do i=1,mxCen
         iBondPar(i)=1
      EndDo
      Do i=1,8
         nPrim(i)=0
         nBas(i)=0
      EndDo
*                                                                      *
************************************************************************
*                                                                      *
*     Set some defaults
!      iPol=1
      LNearestAtom=.True.
      LAllCenters=.False.
      AveOrb=.False.
      LLumOrb=.False.
      Diffuse(1)=.False.
      Diffuse(2)=.False.
      Diffuse(3)=.False.
      dLimmo(1)=0.65d0
      dLimmo(2)=2.0d0
      Thrs1=1d-5
      Thrs2=1d-4
      nThrs=3
      ThrsMul=1d-2
      iPrint=1
*                                                                      *
************************************************************************
*                                                                      *
*---- Get some coords, nuc. charges and labels
*
! Memory debugging
!      Call Setmem('TRACE=ON')
*
* IO check
!      call fastio('TRACE=ON')
*
C     Call Bnnr
      Call Get_cArray('Relax Method',Method,8)
      If(Method.eq.'RHF-SCF') Then
         iPol=1
      ElseIf(Method.eq.'UHF-SCF') Then
         iPol=1
      ElseIf(Method.eq.'MBPT2') Then
         iPol=0
      EndIf
!     Runfile update
!      Call Get_nAtoms(nAtoms)
      Call Get_iScalar('Unique atoms',nAtoms)
      If (nAtoms.gt.mxAtomMP) Then
         Write (6,'(A)')'MPProp: Too many atoms'
         Call Abend()
      End If
!
      Call GetMem('Coord','Allo','Real',ip_Coor,3*nAtoms)
      Call GetMem('Atype','Allo','Real',iAtype,nAtoms)
!     Runfile update
!      Call Get_Coord(Work(ip_Coor),nAtoms)
      Call Get_dArray('Unique Coordinates',Work(ip_Coor),nAtoms*3)

!     Runfile update
!      Call Get_AtomLabel(Labe,nAtoms)
      Call Get_cArray('Unique Atom Names',Labe,LENIN*nAtoms)
!

!     Runfile update
!      Call Get_Charge(Work(iAtype),nAtoms)
      Call Get_dArray('Nuclear charge',Work(iAtype),nAtoms)
!
      Call GetMem('Qnuc','Allo','Real',iQnuc,nAtoms)
!     Runfile update
!      Call Get_Charge_Eff(Work(iQnuc),nAtoms)
      Call Get_dArray('Effective nuclear Charge',Work(iQnuc),nAtoms)
!
      Do i=1,nAtoms
         iAtomType(i) = Int(Work(iAtype+i-1))
         COR(1,i,i) = Work(ip_Coor+(i-1)*3)
         COR(2,i,i) = Work(ip_Coor+(i-1)*3+1)
         COR(3,i,i) = Work(ip_Coor+(i-1)*3+2)
      End Do
      Call Wr_Cord(nAtoms)
      Call GetMem('Atype','Free','Real',iAtype,nAtoms)
      nCenters  = nAtoms*(nAtoms+1)/2
      nSum=nAtoms*6
*                                                                      *
************************************************************************
*
*---- Get information from input

      Call Get_Mpprop_input(nAtoms,iPol,LNearestAtom,LAllCenters,AveOrb,
     &                      LLumOrb,Diffuse,dLimmo,Thrs1,Thrs2,nThrs,
     &                      ThrsMul,iPrint)
*                                                                      *
************************************************************************
*                                                                      *
*---- Read first the size of the contracted basis using the default file
*     ONEINT and COMFILE
*
!     Runfile update
!      Call Get_nSym(nSym)
      Call Get_iScalar('nSym',nSym)
!
!     Runfile update
!      Call Get_nBas(nBas)
      nIrrep = 1
      Call Get_iArray('nBas',nBas,1)
!
      nVec=0
      nOcc=0
      Do iSym = 1, nSym
         nVec = nVec + nBas(iSym)**2
         nOcc = nOcc + nBas(iSym)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Read first the size of the primitiv basis using ONEREL and COMREL
*
!     Runfile update
!      Call Get_nBas(nPrim)
      Call Get_iArray('nBas_Prim',nPrim,1)

      nSize=0
      Do iSym = 1, nSym
         nSize = nSize + nPrim(iSym)*(nPrim(iSym)+1)/2
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Do a dirty trick since the one-intergral file is not explicitly
*     opened.
      Call Put_iArray('nBas',nPrim,1)
C     Call OneBas('PRIM')
*                                                                      *
************************************************************************
*                                                                      *
*---- Read overlap, dipole moment and quadrupole moment integrals
*
*
* Do it general
*
      Do iMltpl = 0, mxMltPl
         Write (label,'(a,i2)') 'PLTPL ',iMltpl
         nComp = (iMltpl+1)*(iMltpl+2)/2
         Write(MemLabel,'(A5,i3.3)') 'MltPl',iMltpl
         Call GetMem(MemLabel,'Allo','Inte',iMltPlAd(iMltpl),nComp)
         nSum=nSum+nComp
         Do iComp = 1, nComp
            irc=-1
            iopt=1
!EB            Call RdOne (irc,iopt,label,iComp,nInt,iSmLbl)
            Call iRdOne (irc,iopt,label,iComp,iDum,iSmLbl)
            If (irc.eq.0) nInt=iDum(1)
            if (irc.ne.0) Then
               If(iComp.ne.1) Then
                  Write (6,'(2A)')'MPProp: Error reading iComp.ne.0 lab'
     &            ,'el=',label
                  Call Abend()
               Else
                  Call GetMem(MemLabel,'Free','Inte',iMltPlAd(iMltpl),
     &            nComp)
                  nMltPl=iMltPl-1
                  nSum=nSum-nComp
                  go to 100
               EndIf
            EndIf
            If (nInt.ne.0) Then
               Write(MemLabel,'(i3.3,i5.5)') iMltpl, iComp
               Call GetMem(MemLabel,'Allo','Real',
     &         iWork(iMltPlAd(iMltpl)+iComp-1),nInt+4)
               nSum=nSum+nInt+4
               irc=-1
               iopt=0
               Call RdOne (irc,iopt,label,iComp,
     &         Work(iWork(iMltPlAd(iMltpl)+iComp-1)),iSmLbl)
            Else
               Write (6,'(2A)')'MPProp: Error reading nInt=0 label='
     &         ,label
               Call Abend()
            End If
            If (irc.ne.0) Then
               Write (6,'(2A)') '2 MPProp: Error reading ',label
               Call Abend()
            End If
!???????????????????????
            If (nInt.ne.0)
     &         Call CmpInt(Work(iWork(iMltPlAd(iMltpl)+iComp-1)),nInt,
     &         nPrim,nIrrep,iSmLbl)
            Do i=1,3
               CordMltPl(i,iMltpl)=
     &         Work(iWork(iMltPlAd(iMltpl))+nInt+i-1)
            End Do
         End Do
      End Do

100   Continue

*                                                                      *
************************************************************************
*                                                                      *
*---- Allocate Memory For Multipoles on Atoms + Atoms and Bonds
*
      Do iMltpl = 0, nMltPl
         nComp = (iMltpl+1)*(iMltpl+2)/2
         Write(MemLabel,'(A5,i3.3)') 'AMtPl',iMltpl
         Call GetMem(MemLabel,'Allo','Real',iAtMltPlAd(iMltpl),
     &   nComp*nAtoms)
         Write(MemLabel,'(A5,i3.3)') 'ABMtP',iMltpl
         Call GetMem(MemLabel,'Allo','Real',iAtBoMltPlAd(iMltpl),
     &   nComp*nCenters)
         Write(MemLabel,'(A5,i3.3)') 'MtPCp',iMltpl
         Call GetMem(MemLabel,'Allo','Real',iAtBoMltPlAdCopy(iMltpl),
     &   nComp*nCenters)
         nSum=nSum+nComp*(nAtoms+nCenters)
         Do i=1,nComp*nAtoms
            Work(iAtMltPlAd(iMltpl)+i-1)=0.0D0
         EndDo
         Do i=1,nComp*nCenters
            Work(iAtBoMltPlAd(iMltpl)+i-1)=0.0D0
         EndDo
      EndDo
*                                                                      *
************************************************************************
*                                                                      *
*---- Read the P-matrix
*
      Label='P_matrix'
      irc=-1
      iopt=1
      iComp=1
!EB      Call RdOne (irc,iopt,label,iComp,nInt,iSmLbl)
      Call iRdOne (irc,iopt,label,iComp,iDum,iSmLbl)
      nInt=iDum(1)
      If (irc.ne.0) Then
         Write (6,'(2A)') 'MPProp: Error getting length of ',
     &                            label
         Write(6,*) 'Length of the vector', nInt, iSmLbl
         Write(6,*) 'irc=',irc
         Call Abend()
      End If
      If (nInt.ne.nSize) Then
         Write (6,*) 'MPProp: nInt.ne.nSize'
         Write (6,*) 'nInt=',nInt
         Write (6,*) 'nSize=',nSize
         Call Abend()
      End If
      irc=-1
      iopt=0

      Call GetMem('CenX','Allo','Real',iCenX,nInt+4)
      Call GetMem('CenY','Allo','Real',iCenY,nInt+4)
      Call GetMem('CenZ','Allo','Real',iCenZ,nInt+4)
      nSum=nSum+nInt*3
*
      iComp=1
      Call RdOne (irc,iopt,label,iComp,Work(iCenX),iSmLbl)
      iComp=2
      Call RdOne (irc,iopt,label,iComp,Work(iCenY),iSmLbl)
      iComp=3
      Call RdOne (irc,iopt,label,iComp,Work(iCenZ),iSmLbl)
*                                                                      *
************************************************************************
*                                                                      *
*     Restore
*
      Call Put_iArray('nBas',nBas,1)
C     Call OneBas('CONT')
*                                                                      *
************************************************************************
*                                                                      *
*     READ TRANSFORMATION MATRIX  nPrim(i)*nBas(i)
*
*     The TM vectors holds the contraction coefficients of each
*     basisset in the primitive base. This means that the columns in the
*     TM are composed of the coefficients for each atomic basefunction
*     and there are nBas(iSym) of them. The coefficients found in each
*     column are exactly the same as in the basisset shifted to the
*     right position in the matrix column. For a hydrogen molecule
*     calculated with 2s function on each atom the first column will
*     look like (c11,c21,...,0,0,...) where the zeros are to delete the
*     primitive gaussians on the "second" hydrogen atom.
*     The second column looks like (c12,c22,...,0,0,...) and the third
*     (0,0,...,c11,c21,...) where in the third column the coefficients
*     are shifted to zero out the "first" hydrogen atom.
*
       nTM=0
       Do iSym = 1, nSym
          nTM = nTM + nBas(iSym)*nPrim(iSym)
       End Do
       Call GetMem('TM','Allo','Real',ip_TM,nTM)
       nSum=nSum+nTM
!      Runfile update
!       Call Get_TPC(Work(ip_TM),nTM)
       Call Get_dArray('NEMO TPC',Work(ip_TM),nTM)
!
*                                                                      *
************************************************************************
*                                                                      *
*---- Read the CMO's and occupation numbers
*
*     CMO's stored as nBas(iSym)*nSO(iSym)
*
      NOCOB=0
      nOcOb_b=0
      If(LLumOrb) Then
        If(Method.eq.'UHF-SCF') Then
          Call GetMem('Vec','Allo','Real',ip_Vec,2*nVec)
          Call GetMem('Occ','Allo','Real',ip_Occ,2*nOcc)
          Call GetMem('Ene','Allo','Real',ip_Ene,2*nOcc)
          nSum=nSum+4*nOcc+2*nVec
        Else
          Call GetMem('Vec','Allo','Real',ip_Vec,nVec)
          Call GetMem('Occ','Allo','Real',ip_Occ,nOcc)
          Call GetMem('Ene','Allo','Real',ip_Ene,2*nOcc)
c                    ! The 2* is for technical reason
          nSum=nSum+2*nOcc+nVec
        EndIf
        Lu_=11
        FName='INPORB'
        If(Method.eq.'UHF-SCF') Then
          Call RdVec_(FName,Lu_,'COE',1,nSym,nBas,nBas,Work(ip_Vec),
     &         Work(ip_Vec+nVec),Work(ip_Occ),Work(ip_Occ+nOcc),
     &         Work(ip_Ene),Work(ip_Ene+nOcc),iDum,VTitle,iWarn,iErr,
     &         iWFtype)
        Else
          Call RdVec(FName,Lu_,'COE',nSym,nBas,nBas,Work(ip_Vec),
     &         Work(ip_Occ),Work(ip_Ene),iDum,VTitle,iWarn,iErr)
        EndIf
        If (Index(VTitle,'IVO').ne.0) Then
          Write (6,*) ' MpProp not implemented for IVO orbitals!'
          Call Abend
        EndIf
!       If(iPol.eq.2) Then
!          Call LauraPol()
!       End If
        If(Method.eq.'UHF-SCF') Then
          Do i=0,nOcc-1
            If (Work(ip_Occ+i).ne.0.0D0) Then
              nOcOb = nOcOb + 1
            End If
          End Do
          Do i=nOcc,2*nOcc-1
            If (Work(ip_Occ+i).ne.0.0D0) Then
              nOcOb_b = nOcOb_b + 1
            End If
          End Do
        Else
          Do i=0,nOcc-1
            If (Work(ip_Occ+i).ne.0.0D0) Then
              nOcOb = nOcOb + 1
            End If
          End Do
        EndIf
*                                                                      *
************************************************************************
*                                                                      *
*---- Project the MO's on the primitive basis
*
        nVec_p=0
        Do iSym = 1, nSym
          nVec_p = nVec_p + nPrim(iSym)*nBas(iSym)
        End Do
        Call GetMem('Vec_p','Allo','Real',ip_Vec_p,nVec_p)
        nSum=nSum+nVec_p
        iOff1=ip_Vec
        iOff2=ip_TM
        iOff3=ip_Vec_p
        Do iSym = 1, nSym
          If (nPrim(iSym).gt.0) Then
            Call DGEMM_('N','N',
     &                  nPrim(iSym),nBas(iSym),nBas(iSym),
     &                  1.0d0,Work(iOff2),nPrim(iSym),
     &                  Work(iOff1),nBas(iSym),
     &                  0.0d0,Work(iOff3),nPrim(iSym))
            iOff1 = iOff1 + nBas(iSym)**2
            iOff2 = iOff2 + nPrim(iSym)*nBas(iSym)
            iOff3 = iOff3 + nPrim(iSym)*nBas(iSym)
          End If
        End Do

        If(Method.eq.'UHF-SCF') Then
          Call GetMem('Vec_p_b','Allo','Real',ip_Vec_p_b,nVec_p)
          nSum=nSum+nVec_p
          iOff1=ip_Vec+nVec
          iOff2=ip_TM
          iOff3=ip_Vec_p_b
          Do iSym = 1, nSym
            If (nPrim(iSym).gt.0) Then
              Call DGEMM_('N','N',
     &                    nPrim(iSym),nBas(iSym),nBas(iSym),
     &                    1.0d0,Work(iOff2),nPrim(iSym),
     &                    Work(iOff1),nBas(iSym),
     &                    0.0d0,Work(iOff3),nPrim(iSym))
              iOff1 = iOff1 + nBas(iSym)**2
              iOff2 = iOff2 + nPrim(iSym)*nBas(iSym)
              iOff3 = iOff3 + nPrim(iSym)*nBas(iSym)
            End If
          End Do
        EndIf
*
        Call GetMem('Ocof','Allo','Real',iOcof,nVec_p)
        nSum=nSum+nVec_p
        Call Get_OCOF(nPrim(1),nBas(1),Work(ip_Vec_p),nVec_p,
     &  Work(iOcof))
        If(Method.eq.'UHF-SCF') Then
          Call GetMem('Ocofb','Allo','Real',iOcof_b,nVec_p)
          nSum=nSum+nVec_p
          Call Get_OCOF(nPrim(1),nBas(1),Work(ip_Vec_p_b),nVec_p,
     &         Work(iOcof_b))
        EndIf
*
        Call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),ip_D_p,nOcOb,
     &       Work(ip_Occ),Work(iOcof))
        If(Method.eq.'UHF-SCF') Then
          Call Gen_Prim_Density_Matrix(nBas(1),nPrim(1),ip_D_p_b,
     &         nOcOb_b,Work(ip_Occ+nOcc),Work(iOcof_b))
        EndIf
      Else
*                                                                      *
************************************************************************
*                                                                      *
*---- If the densities are used for expansion then
*---- Project the densities on the primitive basis
*
!       Get the orbital energy and allocate dubble memory due to UHF calculations
        If(Method.eq.'RHF-SCF') Then
          Call Get_OrbE_mpprop(ip_Ene,nOcc)
        Else
          Call GetMem('OrbE','Allo','Real',ip_Ene,2*nOcc)
          Do iEne=1,2*nOcc
            Work(ip_Ene+iEne-1)=0.0d0
          EndDo
        EndIf
        Call Get_Density_Matrix_mpprop(ip_D,nDens,nBas(1),nSym)
        Write(6,*) 'No polarizability will be calculated'
        iPol=0
        Call Get_Prim_Density_Matrix(ip_D,nBas(1),ip_D_p,nPrim(1),
     &      Work(ip_TM))
        Call Free_Work(ip_D)
      EndIf

*                                                                      *
************************************************************************
*                                                                      *
      Call Get_Prim_Atom_Tab(nAtoms,nPrim(1),Work(ip_Coor),
     &Work(iCenX),Work(iCenY),Work(iCenZ))
*                                                                      *
************************************************************************
*                                                                      *
* Get the SCF or RASSCF energy
!     Runfile update
!      Call Get_Energy(EneV)
      Call Get_dScalar('Last energy',EneV)
!
      Write(6,*)
      Write(6,'(a,f16.8)') ' Total SCF energy ', EneV
      Write(6,*)
      nOrbi = nBas(1)

* If multipoles and polarizability of the bonds should go to
* the atoms they belong to or to the nearest atom
      If(.not.LNearestAtom) Then
         Write(6,*)
         Write(6,*)' I will not move bonds and polarizabilities'
         Write(6,*)' to the nearest atom, just to the bonding pair'
         Write(6,*)
      Else
         Write(6,*)
         Write(6,*)' I WILL move bonds and polarizabilities'
         Write(6,*)' to the nearest atom'
         Write(6,*)
      EndIf

C     Get multipole properties
      LFirstRun=.True.

      Call Get_MpProp(nPrim(1),nBas(1),nAtoms,nCenters,!nOcOb,
!     &nMltPl,Work(ip_Occ),nOcc,Work(iOcof),
     &nMltPl,ip_D_p,
     &Work(iCenX),Work(iCenY),Work(iCenZ),LNearestAtom,
     &LFirstRun,LLumOrb)

      If(Method.eq.'UHF-SCF') Then
         LFirstRun=.False.
         Call Get_MpProp(nPrim(1),nBas(1),nAtoms,nCenters,!nOcOb_b,
!     &        nMltPl,Work(ip_Occ+nOcc),nOcc,Work(iOcof_b),
     &        nMltPl,ip_D_p_b,
     &        Work(iCenX),Work(iCenY),Work(iCenZ),LNearestAtom,
     &        LFirstRun,LLumOrb)
         LFirstRun=.True.
      EndIf

*
*-- If the esteemed user wants to obtain diffuse local distributions,
*   then proceed here. Much of the code is common with LoProp, hence
*   there is first a call to a routine to order some quantaties from
*   MpProp in the same say as in LoProp, then Mother Goose is called.
*
      If(Diffuse(1)) then
        Call StoreMpAsLop(nAtoms,ip_ANr,nBas(1),ip_Ttot,ip_Ttot_Inv
     &                   ,ipMP,nMltPl,ip_EC)
        Call GetMem('ToPoint','Allo','Real',iTP,nAtoms)
        Call CoreToPoint(nAtoms,ipMP,iTP)
        LuYou=IsFreeUnit(81)
        Call OpnFl('DIFFPR',LuYou,Exist)
        Call Diff_MotherGoose(Diffuse,nAtoms,nBas(1),ipMP,ip_Coor
     &                       ,nCenters,ip_EC
     &                       ,ip_ANr,ip_Ttot
     &                       ,ip_Ttot_Inv,nMltPl,iTP,dLimmo
     &                       ,Thrs1,Thrs2,nThrs,iPrint
     &                       ,ThrsMul,LuYou)
        Close(LuYou)
        Call Free_iWork(ip_ANr)
        Call GetMem('T','Free','Real',ip_Ttot,nBas(1)**2)
        Call GetMem('Tinv','Free','Real',ip_Ttot_Inv,nBas(1)**2)
        Call GetMem('ExpCent','Free','Real',ip_EC,3*nAtoms*(nAtoms+1)/2)
        nSize1=nAtoms*(nAtoms+1)/2
        nSize2=(nMltPl*(nMltPl**2+6*nMltPl+11)+6)/6
        Call GetMem('MultMom','Free','Real',ipMP,nSize1*nSize2)
        Call GetMem('ToPoint','Free','Real',iTP,nAtoms)
      Endif
*
*-- End of Diffuse.
*

      If(LLumOrb) Then
C     Get center of charge for each molecular orbital
       Call GetMem('Ocen','Allo','Real',iOcen,3*nOrbi)
       If(Method.eq.'UHF-SCF')
     & Call GetMem('Ocen_b','Allo','Real',iOcen_b,3*nOrbi)
       nSum=nSum+3*nOrbi*2
      EndIf
C     Get polarizabillities if iPol
      Call GetMem('AtPol','Allo','Real',iAtPolAd,nAtoms*6)
      Call GetMem('AtBoPol','Allo','Real',iAtBoPolAd,nCenters*6)
      nSum=nSum+6*(nCenters+nAtoms)
      Do i=0,nAtoms*6-1
         Work(iAtPolAd+i)=0.0D0
      EndDo
      Do i=0,nCenters*6-1
         Work(iAtBoPolAd+i)=0.0D0
      EndDo
      If(iPol.gt.0) Then
!EB         Call Get_OrbCen(nPrim(1),nBas(1),NORBI,Work(iWork(iMltPlAd(0)))
        Call Get_OrbCen(nPrim(1),NORBI,Work(iWork(iMltPlAd(0))),
     &  Work(iOcen),Work(iCenX),Work(iCenY),Work(iCenZ),Work(iOcof))
         If(Method.eq.'UHF-SCF')
     &     Call Get_OrbCen(nPrim(1),NORBI,Work(iWork(iMltPlAd(0))),
     &          Work(iOcen_b),Work(iCenX),Work(iCenY),Work(iCenZ),
     &          Work(iOcof_b))
         If(iPol.eq.1) Then
            If(nOcOb.lt.nOcc) Then
               Call Get_Polar(nPrim(1),nBas(1),nAtoms,nCenters,nOcOb,
!EB     &         Work(ip_Ene),Work(ip_Occ),nOcc,Work(iOcof),Work(iOcen),
     &         Work(ip_Ene),nOcc,Work(iOcof),Work(iOcen),LNearestAtom,
     &         LFirstRun)
               If(Method.eq.'UHF-SCF') Then
                  LFirstRun=.False.
                  Call Get_Polar(nPrim(1),nBas(1),nAtoms,nCenters,
     &            nOcOb_b,Work(ip_Ene+nOcc),nOcc,Work(iOcof_b),
     &            Work(iOcen_b),LNearestAtom,LFirstRun)
               EndIf
            Else
              Write(6,*)
              Write(6,*)'I will not do an analyze of the polarizability'
              Write(6,*)'no of occupied orb. is to large'
              Write(6,*)
              If(Method.eq.'UHF-SCF') Then
                 Write(6,*)' nOcOb nOcOb_b nOcc ',nOcOb,nOcOb_b,nOcc
              Else
                 Write(6,*)' nOcOb nOcc ',nOcOb,nOcc
              EndIf
              Write(6,*)
              iPol=0
            EndIf
         ElseIf(iPol.eq.2) Then
            Call LauraPol()
         EndIf
      End If

C     Write output, the properties
      Call Wr_Prop(nAtoms,nCenters,nBas(1),nMltPl,NOCOB,NOCOB_b,
     &     Work(ip_Ene),Work(ip_Ene+nOcc),iPol,LAllCenters)
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate Work
*
      Write(6,*)
      Write(6,*) 'Number of allocated real*8 words', nsum
      Write(6,'(a,f6.2,a)') ' That is ',dble(nsum)*8.0/1024.0/1024.0,
     &' MBytes'
      Write(6,*)
      Call Free_Work(ip_D_p)
!      If(iPol.gt.0) Then
         Call GetMem('AtBoPol','Free','Real',iAtBoPolAd,nCenters*6)
         Call GetMem('AtPol','Free','Real',iAtPolAd,nAtoms*6)
!      EndIf
      If(LLumorb) Then
        If(Method.eq.'UHF-SCF') Then
          Call GetMem('Ocen_b','Free','Real',iOcen_b,3*nOrbi)
          Call GetMem('Ocof','Free','Real',iOcof_b,nVec_p)
          Call GetMem('Vec_p','Free','Real',ip_Vec_p_b,nVec_p)
        EndIf
        Call GetMem('Ocen','Free','Real',iOcen,3*nOrbi)
        Call GetMem('Ocof','Free','Real',iOcof,nVec_p)
        Call GetMem('Vec_p','Free','Real',ip_Vec_p,nVec_p)
        Call GetMem('Occ','Free','Real',ip_Occ,nOcc)
        Call GetMem('Vec','Free','Real',ip_Vec,nVec)
      EndIf
      Call GetMem('Ene','Free','Real',ip_Ene,nOcc)
      Call GetMem('TM','Free','Real',ip_TM,nTM)
      Call GetMem('CenZ','Free','Real',iCenZ,nInt+4)
      Call GetMem('CenY','Free','Real',iCenY,nInt+4)
      Call GetMem('CenX','Free','Real',iCenX,nInt+4)
      Do iMltpl = 0, nMltPl
         nComp = (iMltpl+1)*(iMltpl+2)/2
         Write(MemLabel,'(A5,i3.3)') 'AMtPl',iMltpl
         Call GetMem(MemLabel,'Free','Real',iAtMltPlAd(iMltpl),
     &   nComp*nAtoms)
         Write(MemLabel,'(A5,i3.3)') 'ABMtP',iMltpl
         Call GetMem(MemLabel,'Free','Real',iAtBoMltPlAd(iMltpl),
     &   nComp*nCenters)
         Call GetMem(MemLabel,'Free','Real',iAtBoMltPlAdCopy(iMltpl),
     &   nComp*nCenters)
         Do iComp=1,nComp
            Write(MemLabel,'(i3.3,i5.5)') iMltpl, iComp
            Call GetMem(MemLabel,'Free','Real',
     &      iWork(iMltPlAd(iMltpl)+iComp-1),nInt+4)
         EndDo
         Write(MemLabel,'(A5,i3.3)') 'MltPl',iMltpl
         Call GetMem(MemLabel,'Free','Inte',iMltPlAd(iMltpl),nComp)
      EndDo
      Call GetMem('Qnuc','Free','Real',iQnuc,nAtoms)
      Call GetMem('Coord','Free','Real',ip_Coor,3*nAtoms)
      Call GetMem('Coord','Check','Real',ip_Coor,3*nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
      iReturn=0
      Return
      End
