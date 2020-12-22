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
      use Slapaf_Info, only: Cx, Gx, Gx0, NAC, Coor, Grd,
     &                       Energy, Energy0, DipM, qInt, dqInt,
     &                       RefGeo, Get_Slapaf
      use Slapaf_Parameters, only: MaxItr, mTROld, lOld_Implicit,
     &                             TwoRunFiles
      Implicit Real*8 (a-h,o-z)
#include "sbs.fh"
#include "real.fh"
#include "nadc.fh"
#include "stdalloc.fh"
#include "info_slapaf.fh"
#include "print.fh"
      Logical Is_Roots_Set
      Real*8, Allocatable:: MMGrd(:,:), Tmp(:), DMs(:,:)
#include "SysDef.fh"
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
*     Get some basic information from the runfile.
*
      Call Get_Slapaf(iter, MaxItr, mTROld, lOld_Implicit,SIZE(Coor,2),
     &                mLambda)
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*---  Pick up information from previous iterations
*
*                                                                      *
************************************************************************
*                                                                      *
      If (iter/=1) Then
*                                                                      *
************************************************************************
*                                                                      *
         Call qpg_dArray('qInt',Found,nqInt)
         If (Found) Then
            nQQ=nqInt/MaxItr
            Call mma_allocate( qInt,nQQ,MaxItr,Label= 'qInt')
            Call mma_allocate(dqInt,nQQ,MaxItr,Label='dqInt')
            Call Get_dArray( 'qInt', qInt,nQQ*MaxItr)
            Call Get_dArray('dqInt',dqInt,nQQ*MaxItr)
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
      Cx(:,:,iter)=Coor(:,:)
      Gx(:,:,iter)=Grd(:,:)

      If (iter.gt.1) Then
        Temp=Zero
        Do i = 1, SIZE(Coor,2)
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

* In case of a QM/MM geometry optimization, all the old MM gradients are
* replaced by the new one (both gradients are stored on the Runfile).
*
      Call Qpg_dArray('MM Grad',lMMGrd,6*SIZE(Coor,2))
      lMMGrd = .False.
      If (lMMGrd) Then
         Call mma_allocate(MMGrd,3*SIZE(Coor,2),2,Label='MMGrd')
         Call Get_dArray('MM Grad',MMGrd,3*SIZE(Coor,2)*2)
         Do iN = 1, iter-1
            Write(6,*) 'Grad at iteration :',iN
            Call RecPrt('Old:',' ',Gx(:,:,iN),3,SIZE(Coor,2))
            Call DaXpY_(3*SIZE(Coor,2),-One,MMGrd(:,1),1,Gx(:,:,iN),  1)
            Call DaXpY_(3*SIZE(Coor,2), One,MMGrd(:,2),1,Gx(:,:,iN),  1)
            Call RecPrt('New:',' ',Gx(1,1,iN),3,SIZE(Coor,2))
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
            If (.Not.Allocated(RefGeo)) Then
               Call mma_allocate(RefGeo,3,SIZE(Coor,2),Label='RefGeo')
            End If
            Call Get_dArray('Ref_Geom',RefGeo,3*SIZE(Coor,2))
         Else
*
*           Not defined: default reference structure to the starting
*           structure.
*
            If (.Not.Allocated(RefGeo)) Then
               Call mma_allocate(RefGeo,3,SIZE(Coor,2),Label='RefGeo')
            End If
            RefGeo(:,:)=Cx(:,:,1)
            Call Put_dArray('Ref_Geom',RefGeo,3*SIZE(Coor,2))
         End If
      Else
*
*        Pick up the reference structure.
*
         If (.Not.Allocated(RefGeo)) Then
            Call mma_allocate(RefGeo,3,SIZE(Coor,2),Label='RefGeo')
         End If
         Call Get_dArray('Ref_Geom',RefGeo,3*SIZE(Coor,2))
      End If
*
*     Align the reference structure to the current one, otherwise
*     measuring distances does not make much sense
*
*     (disabled for the moment, moving the reference affects the
*      computation of some vectors for MEP)
C     If (iter.gt.1) Call Align(RefGeo,Coor,SIZE(Coor,2))
C     Call RecPrt('Ref_Geom',' ',RefGeo,3,SIZE(Coor,2))
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
            nGrad=3*SIZE(Coor,2)
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
