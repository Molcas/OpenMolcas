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
      Subroutine Init2()
      use Slapaf_Info, only: Cx, Gx, Gx0, NAC, Coor, Grd, GNrm, Lambda,
     &                       Energy, Energy0, MF, DipM
      Implicit Real*8 (a-h,o-z)
#include "sbs.fh"
#include "real.fh"
#include "nadc.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "info_slapaf.fh"
#include "print.fh"
      Logical Is_Roots_Set
      Integer, Allocatable:: Information(:)
      Real*8, Allocatable:: MMGrd(:,:), Tmp(:), DMs(:,:)
#include "SysDef.fh"
      Character*100 Get_SuperName, SuperName
      External Get_SuperName
* Beijing Test
      Logical Exist_2,lMMGrd, Found
      Real*8 Columbus_Energy(2)
************ columbus interface ****************************************
      Integer Columbus
*
*                                                                      *
************************************************************************
*                                                                      *
*     Dummy-add the TS and saddle constraints, so that the arrays have
*     a large enough size, these constraints will be actually added or
*     not later on.
*
      Call Merge_Constraints('UDC','TSC','purge.Udc',mLambda,iDummy)
      Call Merge_Constraints('purge.Udc','UDC.Saddle','purge.Udc',
     &                       mLambda,iDummy)
*                                                                      *
************************************************************************
*                                                                      *
*
*     Get the gradient information from the runfile.
*
      Call mma_allocate(Information,7,Label='Information')
      Call qpg_iArray('Slapaf Info 1',Exist,itmp)
      If (Exist) Call Get_iArray('Slapaf Info 1',Information,7)
*
      If (.Not.Exist.or.(Exist.and.Information(1).eq.-99)) Then
C        Write (6,*) 'Reinitiate Slapaf fields on runfile'
         Information(:)=0
         Information(3)=-99
         Call Put_iArray('Slapaf Info 1',Information,7)
      End If

*     iNew  =Information(1)
      iter  =Information(2)+1
      If (iter.ge.MaxItr+1) Then
         Write (6,*) 'Increase MaxItr in info_slapaf.fh'
         Call WarningMessage(2,'iter.ge.MaxItr+1')
         Call Abend()
      End If
      mTROld=Information(3)
      lOld_Implicit= Information(4).eq.1
      Call mma_deallocate(Information)
*
*---  Pick up information from previous iterations
*
      l1=3*nsAtom
      Lngth=
     &     +(MaxItr+1)            ! Ener
     &     +(MaxItr+1)            ! Ener0
     &     +(MaxItr+1)*3          ! DipM
     &     +(MaxItr+1)            ! GNrm
     &     +l1*(MaxItr+1)         ! Cx
     &     +l1*(MaxItr+1)         ! Gx
     &     +l1*(MaxItr+1)         ! Gx0
     &     +l1                    ! MF
     &     +mLambda*(MaxItr+1)    ! L
      Call GetMem('Relax','Allo','Real',ipRlx,Lngth)
      Call FZero(Work(ipRlx),Lngth)
      ipEner = ipRlx
      ipEner0= ipEner  +          MaxItr+1
      ipDipM = ipEner0 +          MaxItr+1
      ipGNrm = ipDipM  +         (MaxItr+1)*3
      ipCx   = ipGNrm  +          MaxItr+1
      ipGx   = ipCx    +      l1*(MaxItr+1)
      ipGx0  = ipGx    +      l1*(MaxItr+1)
      ipMF   = ipGx0   +      l1*(MaxItr+1)
      ipL    = ipMF    +      l1
      ipqInt = ip_Dummy
      ipdqInt= ip_Dummy
*
      Call mma_allocate(Energy,          MaxItr+1,Label='Energy')
      Energy(:) = Zero
      Call mma_allocate(Energy0,         MaxItr+1,Label='Energy0')
      Energy0(:) = Zero
      Call mma_allocate(DipM,   3,       MaxItr+1,Label='DipM')
      DipM(:,:) = Zero
      Call mma_allocate(GNrm,            MaxItr+1,Label='GNrm')
      GNrm(:) = Zero
      Call mma_allocate(Cx,     3,nsAtom,MaxItr+1,Label='Cx')
      Cx(:,:,:) = Zero
      Call mma_allocate(Gx,     3,nsAtom,MaxItr+1,Label='Gx')
      Gx(:,:,:) = Zero
      Call mma_allocate(Gx0,    3,nsAtom,MaxItr+1,Label='Gx0')
      Gx0(:,:,:) = Zero
      Call mma_allocate(MF,     3,nsAtom,         Label='MF')
      MF(:,:) = Zero
      If (mLambda>0) Then
        Call mma_allocate(Lambda,mLambda,MaxItr+1,Label='Lambda')
        Lambda(:,:)=Zero
      End If
*                                                                      *
************************************************************************
*                                                                      *
      If (iter.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
*...     Start iteration
         nqInt=0
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         SuperName=Get_Supername()
         If (SuperName.ne.'numerical_gradient') Then
            Call Get_dArray('Slapaf Info 2',Work(ipRlx),Lngth)
*
            Call DCopy_(   (MaxItr+1),Work(ipEner ),1,Energy ,1)
            Call DCopy_(   (MaxItr+1),Work(ipEner0),1,Energy0,1)
            Call DCopy_(3* (MaxItr+1),Work(ipDipM ),1,DipM   ,1)
            Call DCopy_(   (MaxItr+1),Work(ipGNrm ),1,GNrm   ,1)
            Call DCopy_(l1*(MaxItr+1),Work(ipCx   ),1,Cx     ,1)
            Call DCopy_(l1*(MaxItr+1),Work(ipGx   ),1,Gx     ,1)
            Call DCopy_(l1*(MaxItr+1),Work(ipGx0  ),1,Gx0    ,1)
            Call DCopy_(l1           ,Work(ipMF   ),1,MF     ,1)
            If (mLambda>0) Call DCopy_(mLambda*(MaxItr+1),Work(ipL),1,
     &                                                    Lambda,1)
         Else
            iter=1
         End If
         Call qpg_dArray('qInt',Found,nqInt)
         If (Found) Then
            Call GetMem(' qInt','Allo','Real',ipqInt, nqInt)
            Call GetMem('dqInt','Allo','Real',ipdqInt,nqInt)
            Call Get_dArray('qInt',Work(ipqInt),nqInt)
            Call Get_dArray('dqInt',Work(ipdqInt),nqInt)
         End If
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*-----Save coordinates and gradients from this iteration onto the list.
*
      ipOff = (iter-1)*3*nsAtom + ipCx
      call dcopy_(3*nsAtom,Coor,1,Work(ipOff),1)
      call dcopy_(3*nsAtom,Coor,1,Cx(:,:,iter),1)
      If (iter.gt.1) Then
        Temp=Zero
        Do i = 1, nsAtom
           Do j = 1, 3
              Temp=Max(Temp,Abs(Cx(j,i,iter)-Cx(j,i,iter-1)))
           End Do
           If (Temp.gt.Zero) Exit
        End Do
        If (Temp.eq.Zero) Then
          Call WarningMessage(2,'Error in Init2')
          Write (6,*)
          Write (6,*) '****************** ERROR *********************'
          Write (6,*) 'Coordinates did not change!'
          Write (6,*) 'Maybe SEWARD is not inside the loop?'
          Write (6,*) '**********************************************'
          Call Quit_OnUserError()
        End If
      End If
      call dcopy_(3*nsAtom,Grd ,1,Gx(1,1,iter),1)
* In case of a QM/MM geometry optimization, all the old MM gradients are
* replaced by the new one (both gradients are stored on the Runfile).
*
      Call Qpg_dArray('MM Grad',lMMGrd,6*nsAtom)
      lMMGrd = .False.
      If (lMMGrd) Then
         Call mma_allocate(MMGrd,3*nsAtom,2,Label='MMGrd')
         Call Get_dArray('MM Grad',MMGrd,3*nsAtom*2)
         Do iN = 1, iter-1
            Write(6,*) 'Grad at iteration :',iN
            Call RecPrt('Old:',' ',Gx(1,1,iN),3,nsAtom)
            Call DaXpY_(3*nsAtom,-One,MMGrd(:,1),1,Gx(1,1,iN),  1)
            Call DaXpY_(3*nsAtom, One,MMGrd(:,2),1,Gx(1,1,iN),  1)
            Call RecPrt('New:',' ',Gx(1,1,iN),3,nsAtom)
         End Do
         Call mma_deallocate(MMGrd)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up the reference structure from input or the run file.
*
      If (Iter.eq.1) Then
*
*        Check if reference structure is defined by Gateway/Seward
*        for the Saddle approach to find a TS.

         Call qpg_dArray('Ref_Geom',Found,nData)
         If (Found) Then
            If (.Not.Ref_Geom) Then
               Call GetMem('ipRef','Allo','Real',ipRef,3*nsAtom)
               Ref_Geom=.True.
            End If
            Call Get_dArray('Ref_Geom',Work(ipRef),3*nsAtom)
         Else
*
*           Not defined: default reference structure to the starting
*           structure.
*
            If (.Not.Ref_Geom) Then
               ipRef = ip_of_Work(Cx(1,1,1))
            End If
            Call Put_dArray('Ref_Geom',Work(ipRef),3*nsAtom)
         End If
      Else
*
*        Pick up the reference structure.
*
         If (.Not.Ref_Geom) Then
            Call GetMem('ipRef','Allo','Real',ipRef,3*nsAtom)
            Ref_Geom=.True.
         End If
         Call Get_dArray('Ref_Geom',Work(ipRef),3*nsAtom)
      End If
*
*     Align the reference structure to the current one, otherwise
*     measuring distances does not make much sense
*
*     (disabled for the moment, moving the reference affects the
*      computation of some vectors for MEP)
C     If (iter.gt.1) Call Align(Work(ipRef),Coor,nsAtom)
C     Call RecPrt('Ref_Geom',' ',Work(ipRef),3,nsAtom)
*                                                                      *
************************************************************************
*                                                                      *
*...  Get the energy of the last iteration
*
*     Check if we are running in C&M mode.
*
      Call Get_iScalar('Columbus',columbus)
*
      If (Columbus.eq.1) Then
         Call Get_dArray('MR-CISD energy',Columbus_Energy,2)
         Energy(iter)=Columbus_Energy(1)
      Else
         Is_Roots_Set = .False.
         Call Qpg_iScalar('Number of roots',Is_Roots_Set)
         nRoots = 1
         If (Is_Roots_Set) Then
            Call Get_iScalar('Number of roots',nRoots)
         End If
*        Write (6,*) 'Runfile'
*        Write (6,*) 'nRoots=',nRoots
         If (nRoots.ne.1) Then
            Call Get_iScalar('NumGradRoot',iRoot)
*           Write (6,*) 'iRoot=',iRoot
            Call mma_allocate(Tmp,nRoots,Label='Tmp')
            Call Get_dArray('Last energies',Tmp,nRoots)
*           Call RecPrt('Last Energies',' ',Tmp,1,nRoots)
            Energy(iter)=Tmp(iRoot)
            Call mma_deallocate(Tmp)
         Else
            Call Get_dScalar('Last energy',Energy(iter))
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*                                                                      *
************ columbus interface ****************************************
*                                                                      *
*     The dipole moment is required for numerical differentiation.
*     Currently it is not written to the RUNFILE if we are in C&M mode.
*
      If (columbus.eq.1) Then
      Else
         Is_Roots_Set = .False.
         Call Qpg_iScalar('Number of roots',Is_Roots_Set)
         nRoots = 1
         If (Is_Roots_Set) Then
            Call Get_iScalar('Number of roots',nRoots)
         End If
         If (nRoots.ne.1) Then
            Call Get_iScalar('NumGradRoot',iRoot)
*           Write (6,*) 'iRoot=',iRoot
            Call mma_allocate(DMs,3,nRoots,Label='DMs')
            DMs(:,:)=Zero
            Call Qpg_dArray('Last Dipole Moments',Found,nDip)
            If (Found.and.(nDip.eq.3*nRoots)) Then
               Call Get_dArray('Last Dipole Moments',DMs,3*nRoots)
            End If
            Call dCopy_(3,DMs(:,iRoot),1,DipM(:,iter),1)
            Call mma_deallocate(DMs)
         Else
            Call Qpg_dArray('Dipole moment',Found,nDip)
            If (Found.and.(nDip.eq.3)) Then
               Call Get_dArray('Dipole moment',DipM(:,iter),3)
            Else
               DipM(:,iter)=Zero
            End If
         End If
*        Call RecPrt('Dipole Moment',' ',DipM(:,iter),1,3)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Test if this is a case of an intersection calculation. Depending
*     on if we are running the calculation in M or C&M mode this is
*     done in a bit different way.
*
      If (Columbus.eq.1) Then
*
*        C&M mode
*
         Call Get_iScalar('ColGradMode',iMode)
         If (iMode.eq.2.or.iMode.eq.3) Then
*
            E0=Columbus_Energy(2)
            Energy0(iter)=E0
*
            Call qpg_dArray('Grad State2',Found,Length)
            If (.not.Found .or. Length.eq.0) Then
               Call SysAbendmsg('Get_Molecule','Did not find:',
     &                          'Grad State2')
            End If
            Call Get_dArray('Grad State2',Gx0(1,1,iter),Length)
            Gx0(:,:,iter) = -Gx0(:,:,iter)
*
         End If
         If (iMode.eq.3) Then
            Call mma_allocate(NAC,3,Length/3,Label='NAC')
            Call Get_dArray('NADC',NAC,Length)
         End If
*
      Else
*
*     M mode
*
         Call f_Inquire('RUNFILE2',Exist_2)
         If (Exist_2) Then
            Call NameRun('RUNFILE2')
*
            Is_Roots_Set = .False.
            Call Qpg_iScalar('Number of roots',Is_Roots_Set)
            nRoots = 1
            If (Is_Roots_Set) Then
               Call Get_iScalar('Number of roots',nRoots)
            End If

*           Write (6,*) 'Runfile2'
*           Write (6,*) 'nRoots=',nRoots
            If (nRoots.ne.1) Then
               Call Get_iScalar('NumGradRoot',iRoot)
C              Write (6,*) 'iRoot=',iRoot
               Call mma_allocate(Tmp,nRoots,Label='Tmp')
               Call Get_dArray('Last energies',Tmp,nRoots)
*              Call RecPrt('Last Energies',' ',Tmp,1,nRoots)
               E0=Tmp(iRoot)
               Call mma_deallocate(Tmp)
            Else
               Call Get_dScalar('Last energy',E0)
            End If
            Energy0(iter)=E0
*
            nGrad=3*nsAtom
            Call Get_Grad(Gx0(1,1,iter),nGrad)
            Gx0(:,:,iter) = -Gx0(:,:,iter)
*
            Call NameRun('RUNFILE')
            TwoRunFiles = .True.
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Value_l=20.0D0
      Call Qpg_dScalar('Value_l',Found)
      If (.Not.Found) Call Put_dScalar('Value_l',Value_l)

*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
