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
* Copyright (C) Per-Olof Widmark                                       *
*               2017, Roland Lindh                                     *
************************************************************************
      Subroutine VecFind(OccSet,FermSet,CharSet,SpinSet)
************************************************************************
*                                                                      *
* This routine figure out which set of starting orbitals are used.     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
*          R. Lindh                                                    *
*          Uppsala University, Sweden, Feb 2017                        *
*          Remove Work                                                 *
*                                                                      *
************************************************************************
*                                                                      *
* Input vector mappings                                                *
*                                                                      *
* LstVec              InVec                                            *
* -1 Die               0 Core                                          *
*  0 Old SCF           1 NDDO (not used)                               *
*  1 Guessorb          2 Lumorb                                        *
*  2 Lumorb            3 Old density                                   *
*  3 Old Density       4 Restart (not used)                            *
*  4 Core              8 Old SCF                                       *
*  5 NDDO (use not)    9 Guessorb                                      *
*                                                                      *
************************************************************************
#ifdef _HDF5_
      Use mh5, Only: mh5_fetch_attr
#endif
      Implicit Real*8 (A-H,O-Z)
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Logical OccSet
      Logical FermSet
      Logical CharSet
      Logical SpinSet
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Character*80 cList
      Logical Found
      Integer mynSym
      Integer mynBas(8)
      Integer mynOrb(8)
      Character*10 infoLbl
      Real*8, Dimension(:), Allocatable:: EOrb
*----------------------------------------------------------------------*
* Setup                                                                *
*----------------------------------------------------------------------*
      nSqrSum=0
      Do iSym = 1, nSym
         nSqrSum=nSqrSum+nBas(iSym)*nBas(iSym)
      End Do
*----------------------------------------------------------------------*
* Check start orbital priority list                                    *
*----------------------------------------------------------------------*
*     Write(6,*) 'LstVec',LstVec
      Do i=1,nStOpt
         Found=.false.
         If(LstVec(i).eq.0) Then
*           Write(6,'(2i5,a)') i,LstVec(i),' Old scf orbitals'
            Call qpg_darray('SCF orbitals',Found,nData)
            If(Found .and. nData.eq.nSqrSum) Then
*              Write(6,*) 'Found orbitals'
               Call qpg_darray('OrbE',Found,nData)
               If(Found) Then
*                 Write(6,*) 'Found energies'
                  InVec=8
                  GoTo 101
               End If
            End If
         Else If(LstVec(i).eq.1) Then
*           Write(6,'(2i5,a)') i,LstVec(i),' Guessorb orbitals'
            Call qpg_darray('Guessorb',Found,nData)
            If(Found .and. nData.eq.nSqrSum) Then
               Call qpg_darray('Guessorb energies',Found,nData)
               If(Found) Then
                  InVec=9
                  GoTo 101
               End If
            End If
         Else If(LstVec(i).eq.2) Then
*           Write(6,'(2i5,a)') i,LstVec(i),' Lumorb orbitals'
            Call F_Inquire(SCF_FileOrb,Found)
            If(Found) Then
               If (isHDF5) Then
#ifdef _HDF5_
                 iVer=-1
                 Call mh5_fetch_attr(fileorb_id,'NSYM',mynSym)
                 Call mh5_fetch_attr(fileorb_id,'NBAS',mynBas)
#endif
               Else
                 Call ChkVec(SCF_FileOrb,
     $                  iVer,mynSym,mynBas,mynOrb,InfoLbl,iRc)
               End If
               If(iVer.eq.0) Found=.false.
               If(mynSym.ne.nSym) Found=.false.
               Do j=1,nSym
                 If(mynBas(j).ne.nBas(j)) Found=.false.
               End Do
            End If
            If(Found) Then
*              Write(6,*) 'Found INPORB'
               InVec=2
               GoTo 101
            End If
         Else If(LstVec(i).eq.3) Then
*           Write(6,'(2i5,a)') i,LstVec(i),' Density'
            InVec=0
         Else If(LstVec(i).eq.4) Then
*           Write(6,'(2i5,a)') i,LstVec(i),' Core orbitals'
            InVec=0
            Found=.true.
            GoTo 101
         Else
*           Write(6,'(2i5,a)') i,LstVec(i),' Die'
            InVec=0
            Found=.false.
            GoTo 101
         End If
*        Write(6,*) 'LstVec(i),Found:',LstVec(i),Found
      End Do
      Found=.false.
101   Continue
*     Write(6,*) 'VecFind: InVec=',InVec
*----------------------------------------------------------------------*
* Did we find the requested orbitals?                                  *
*----------------------------------------------------------------------*
      If(.not.Found) Then
         cList=''
         n2=0
         Do i=1,nStOpt
            If(LstVec(i).eq.0) Then
               n1=n2+1
               n2=n1+9
               Write(cList(n1:n2),'(a)') 'Old Scf, '
            Else If (LstVec(i).eq.1) Then
               n1=n2+1
               n2=n1+10
               Write(cList(n1:n2),'(a)') 'Guessorb, '
            Else If (LstVec(i).eq.2) Then
               n1=n2+1
               n2=n1+8
               If (isHDF5) Then
                  Write(cList(n1:n2),'(a)') 'HDF5, '
               Else
                  Write(cList(n1:n2),'(a)') 'Lumorb, '
               End If
            Else If (LstVec(i).eq.3) Then
               n1=n2+1
               n2=n1+13
               Write(cList(n1:n2),'(a)') 'Old density, '
            Else If (LstVec(i).eq.4) Then
               n1=n2+1
               n2=n1+6
               Write(cList(n1:n2),'(a)') 'Core, '
            Else If (LstVec(i).eq.5) Then
               n1=n2+1
               n2=n1+6
               Write(cList(n1:n2),'(a)') 'NDDO, '
            Else If (LstVec(i).eq.-1) Then
               GoTo 290
            Else
               n1=n2+1
               n2=n1+9
               Write(cList(n1:n2),'(a)') 'Unknown, '
            End If
         End Do
290      Continue
         Call SysAbendMsg('SCF:',
     &      'Can not find start orbitals according to list:',
     &      cList)
      End If
*----------------------------------------------------------------------*
* What are the defaults for the different cases?                       *
*----------------------------------------------------------------------*
      Call Get_dScalar('Total Nuclear Charge',Tot_Nuc_Charge)
      Tot_El_Charge=Tot_Charge-Tot_Nuc_Charge
      If(InVec.eq.0) Then
*
* We will use core diagonalization
*
*        Write(6,*) 'Using core diagonalization'
         OccSet=.false.
         FermSet=.true.
      Else If(Invec.eq.1) Then
*
* We will use NDDO orbitals, should not be used!
*
*        Write(6,*) 'Using NDDO orbitals'
      Else If(Invec.eq.2) Then
*
* We will use Lumorb orbitals
*
*        Write(6,*) 'Using Lumorb orbitals'
         Call ChkLumo(OccSet,FermSet,SpinSet)
      Else If(Invec.eq.3) Then
*
* We will use density as start, does it even work?
*
*        Write(6,*) 'Using density'
      Else If(Invec.eq.4) Then
*
* This is a restart case
*
*        Write(6,*) 'Using Restart'
      Else If(Invec.eq.8) Then
*
* We will use old SCF orbitals
*
*        Write(6,*) 'Using SCF orbitals'
*        Write(6,*) 'tot_charge',tot_charge
*        Write(6,*) 'tot_el_charge',tot_el_charge
*        Write(6,*) 'tot_nuc_charge',tot_nuc_charge
         Call qpg_iarray('SCF nOcc',Found,nData)
         idspin=0
         If(Found) Then
*           Write(6,*) 'vecfind: Allright, old scf orbitals it is'
            nEle=0
            If(iUHF.eq.0) Then
               Call Get_iarray('SCF nOcc',nOcc(1,1),nSym)
               Call Get_iScalar('SCF mode',iTmp)
               If(iTmp.eq.0) Then
*                 Write(6,*) 'Starting RHF with RHF orbitals'
                  Do iSym=1,nSym
                     nEle=nEle+2*nOcc(iSym,1)
                  End Do
               Else
*                 Write(6,*) 'Starting RHF with UHF orbitals'
                  Call Get_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
                  nEle1=0
                  nEle2=0
                  Do iSym=1,nSym
                     nEle1=nEle1+2*nOcc(iSym,1)
                     nEle2=nEle2+nOcc(iSym,1)+nOcc(iSym,2)
                     idspin=idspin+nOcc(iSym,1)-nOcc(iSym,2)
                  End Do
                  If(nEle1.ne.nEle2) Then
*                    Tot_Charge=Tot_Nuc_Charge-nEle1
*                    Tot_El_Charge=-nEle1
                     CharSet=.true.
                  End If
                  nEle=nEle2
*                 Write(6,*) 'After strange code'
*                 Write(6,*) 'tot_charge',tot_charge
*                 Write(6,*) 'tot_el_charge',tot_el_charge
*                 Write(6,*) 'tot_nuc_charge',tot_nuc_charge
               End If
            Else
*              Write(6,*) 'Starting UHF with RHF/UHF orbitals'
               Call Get_iarray('SCF nOcc',nOcc(1,1),nSym)
               Call qpg_iarray('SCF nOcc_ab',Found,nData)
               If(Found) Then
                  Call Get_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
               Else
                  Call Get_iarray('SCF nOcc',nOcc(1,2),nSym)
               End If
               Do iSym=1,nSym
                  nEle=nEle+nOcc(iSym,1)+nOcc(iSym,2)
                  idspin=idspin+nOcc(iSym,1)-nOcc(iSym,2)
               End Do
            End If
*           Write(6,*) 'idspin',idspin
*           Write(6,*) 'iAu_ab',iAu_ab
            idspin=idspin-iAu_ab
            If(Abs(Tot_El_Charge+nEle).gt.0.5d0.or.idspin.ne.0) Then
               If(Abs(Tot_El_Charge+nEle).gt.0.5d0) Then
*                 Write(6,*) 'System have changed charge!'
               End If
               If(idspin.ne.0) Then
*                 Write(6,*) 'System have changed spin!'
               End If
               If(CharSet.or.idspin.ne.0) Then
                  OccSet=.false.
                  FermSet=.true.
               Else
                  OccSet=.true.
                  FermSet=.false.
*                 Tot_Charge=Tot_Nuc_Charge-nEle
*                 Tot_El_Charge=-nEle
               End If
            Else
*              Write(6,*) 'System have same spin and charge'
               OccSet=.true.
               FermSet=.false.
            End If
*           Write(6,*) 'OccSet  ',OccSet
*           Write(6,*) 'FermSet ',FermSet
*           Write(6,*) 'CharSet ',CharSet
*           Write(6,*) 'SpinSet ',SpinSet
*           Write(6,*) 'nOcc',nOcc
         Else
            OccSet=.false.
            FermSet=.true.
         End If
      Else If(Invec.eq.9) Then
*
* We will use Guessorb orbitals
*
*        Write(6,*) 'Using Guessorb orbitals'
         If(OccSet) Then
*           Write(6,*) 'Occupation is set'
            Continue
         Else If(FermSet) Then
*           Write(6,*) 'Fermi is set'
            Continue
         Else
*           Write(6,*) 'Must decide if to use Fermi'
            If(nAufb(1).eq.-1) Then
               mtmp=Int(-Tot_El_Charge+0.1D0)
               If(iUHF.eq.0) Then
                  If(Mod(mtmp,2).ne.0) Then
                     Write(6,*) 'VecFind: Error in number of electrons'
                     Write(6,*) '         An even number of electrons ',
     &                          '         are required by RHF, use UHF'
                     Write(6,*)
                     Call Abend()
                  End If
                  nAufb(1)=mtmp/2
               Else
                  nAufb(2)=(mtmp-iAu_ab)/2
                  nAufb(1)=Int(-Tot_El_Charge-nAufb(2))
               End If
            End If
*           Write(6,*) 'nAufb',nAufb
*           Write(6,*) 'Now figure out homo-lumo gap'
            Call qpg_darray('Guessorb energies',Found,nData)
            Call mma_allocate(EOrb,nData,Label='EOrb')
            Call get_darray('Guessorb energies',Eorb,nData)
            If(iUHF.eq.0) Then
               Call GetGap(Eorb,nData,nAufb(1),Gap,Ealpha)
            Else
               Call GetGap(Eorb,nData,nAufb(1),tmp,Ealpha)
               Call GetGap(Eorb,nData,nAufb(2),Gap,Ebeta)
               Gap=Min(tmp,Gap)
            End If
            Call get_darray('Guessorb energies',Eorb,nData)
            If(Gap.ge.0.5d0) Then
               If(iUHF.eq.0) Then
                  iOff=0
                  Do iSym=1,nSym
                     n=0
                     Do iBas=1,nBas(iSym)
                        If(EOrb(iOff+iBas).lt.Ealpha) n=n+1
                     End Do
                     nOcc(iSym,1)=n
                     iOff=iOff+nBas(iSym)
                  End Do
               Else
                  iOff=0
                  Do iSym=1,nSym
                     n=0
                     Do iBas=1,nBas(iSym)
                        If(EOrb(iOff+iBas).lt.Ealpha) n=n+1
                     End Do
                     nOcc(iSym,1)=n
                     iOff=iOff+nBas(iSym)
                  End Do
                  iOff=0
                  Do iSym=1,nSym
                     n=0
                     Do iBas=1,nBas(iSym)
                        If(EOrb(iOff+iBas).lt.Ebeta) n=n+1
                     End Do
                     nOcc(iSym,2)=n
                     iOff=iOff+nBas(iSym)
                  End Do
               End If
               OccSet=.true.
               FermSet=.false.
*              Write(6,*) 'Decided on occupation list'
            Else
               OccSet=.false.
               FermSet=.true.
*              Write(6,*) 'Decided on Fermi aufbau'
            End If
            Call mma_deallocate(EOrb)
*           Write(6,*) 'Gap is',Gap
         End If
      Else
*
* This case should not appear
*
         Call SysAbendMsg('SCF:',
     &      'Internal error in VecFind!',
     &      'InVec have illegal value')
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
*     Write(6,'(a,i2)') 'VecFind: InVec=',InVec
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
