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
* Copyright (C) 2013-2015, Ignacio Fdez. Galvan                        *
************************************************************************
*  MEP_Dir
*
*> @brief
*>   Compute the new reference structure and initial coordinates for the next MEP point
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Compute the new reference structure (and direction) and initial geometry
*> for the next MEP point optimization (using the Gonz&aacute;lez--Schlegel or
*> M&uacute;ller--Brown method).
*> Calculate some properties of the path (length and curvature) between the previous
*> and current MEP points.
*> All calculations are done in weighted coordinates (mass-weighted by default).
*>
*> @param[in,out] Cx            Cartesian coordinates in all iterations
*> @param[in]     Gx            Cartesian gradient in all iterations
*> @param[in]     nAtom         Number of symmetry-unique atoms
*> @param[in]     iMEP          Number of this MEP point
*> @param[in]     iOff_iter     Iteration of the previous MEP point
*> @param[in]     iPrint        Print level
*> @param[in]     IRCRestart    Flag to mark the start of a backward IRC search
*> @param[out]    ResGrad       Residual gradient
*> @param[out]    BadConstraint Flag to signal constraint problems
************************************************************************
      Subroutine MEP_Dir(Cx,Gx,nAtom,iMEP,iOff_iter,iPrint,IRCRestart,
     &                   ResGrad,BadConstraint)
      use Symmetry_Info, only: nIrrep
      use Slapaf_Info, only: Weights, MF, RefGeo
      use Slapaf_Parameters, only: IRC, nLambda, rMEP, MEP, nMEP,
     &                             MEPNum, dMEPStep, MEP_Type,
     &                             MEP_Algo, iter
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "weighting.fh"
#include "info_slapaf.fh"
      Real*8 Cx(3*nAtom,iter+1),Gx(3*nAtom,iter+1)
      Logical IRCRestart,BadConstraint
      Parameter ( RadToDeg=180.0D0/Pi )
      Integer iDum(1)
      Real*8, Allocatable:: PrevDir(:,:), PostDir(:,:), Disp(:,:),
     &                      Grad(:,:), Dir(:,:), Cen(:,:),
     &                      Len(:), Cur(:), drdx(:,:)
*
*                                                                      *
************************************************************************
*                                                                      *
      nCoor=3*nAtom
      iPrev_iter=Max(iOff_iter,1)
      Call mma_allocate(PrevDir,3,nAtom,Label='PrevDir')
      Call mma_allocate(PostDir,3,nAtom,Label='PostDir')
      Call mma_allocate(Disp,3,nAtom,Label='Disp')
      Call mma_allocate(Grad,3,nAtom,Label='Grad')
      Call mma_allocate(Dir,3,nAtom,Label='Dir')
      Call mma_allocate(Cen,3,nAtom,Label='Cen')
*
*     Obtain some useful vectors:
*     PrevDir: difference between ref. structure and previous MEP point
*     PostDir: difference between current MEP point and ref. structure
*     Disp:    difference between current and previous MEP points
*     Grad:    gradient at current MEP point
*
      If (iter.gt.1) Then
        PrevDir(:,:) = RefGeo(:,:)
        Call DaXpY_(nCoor,-One,Cx(:,iPrev_iter),1,PrevDir(:,:),1)
        Call dCopy_(nCoor,Cx(:,iter),1,PostDir(:,:),1)
        PostDir(:,:) = PostDir(:,:) - RefGeo(:,:)
        Call dCopy_(nCoor,Cx(:,iter),1,Disp(:,:),1)
        Call DaXpY_(nCoor,-One,Cx(:,iPrev_iter),1,Disp(:,:),1)
      Else
        PrevDir(:,:)=Zero
        PostDir(:,:)=Zero
        Disp(:,:)   =Zero
      End If
      Call dCopy_(nCoor,Gx(:,iter),1,Grad(:,:),1)
*
*     Normalize the vectors in weighted coordinates
*     and compute some angles that provide information on the path
*     shape and quality.
*     Note that gradient and coordinates transform differently
*
      TWeight=Zero
      dPrevDir=Zero
      dPostDir=Zero
      dDisp=Zero
      dGrad=Zero
      dPostDirGrad=Zero
      dPrevDirGrad=Zero
      dPrevDirDisp=Zero
      dPrevDirPostDir=Zero
      iOff=0
      Do iAtom=1,nAtom
        Fact=Dble(iDeg(Cx(1+iOff,iter)))
        xWeight=Weights(iAtom)
        TWeight=TWeight+Fact*xWeight
        Do ixyz=1,3
          dPrevDir=dPrevDir+Fact*xWeight*PrevDir(ixyz,iAtom)**2
          dPostDir=dPostDir+Fact*xWeight*PostDir(ixyz,iAtom)**2
          dDisp=dDisp+Fact*xWeight*Disp(ixyz,iAtom)**2
          dGrad=dGrad+Fact*Grad(ixyz,iAtom)**2/xWeight
          dPostDirGrad=dPostDirGrad+
     &        Fact*PostDir(ixyz,iAtom)*Grad(ixyz,iAtom)
          dPrevDirGrad=dPrevDirGrad+
     &        Fact*PrevDir(ixyz,iAtom)*Grad(ixyz,iAtom)
          dPrevDirDisp=dPrevDirDisp+
     &        Fact*xWeight*PrevDir(ixyz,iAtom)*Disp(ixyz,iAtom)
          dPrevDirPostDir=dPrevDirPostDir+
     &        Fact*xWeight*PrevDir(ixyz,iAtom)*PostDir(ixyz,iAtom)
          iOff=iOff+1
        End Do
      End Do
      dPrevDir=Sqrt(dPrevDir)
      dPostDir=Sqrt(dPostDir)
      dDisp=Sqrt(dDisp)
      dGrad=Sqrt(dGrad)
      If (dPrevDir.gt.Zero)
     &  Call DScal_(nCoor,One/dPrevDir,PrevDir(:,:),1)
      If (dPostDir.gt.Zero)
     &  Call DScal_(nCoor,One/dPostDir,PostDir(:,:),1)
      If (dDisp.gt.Zero)
     &  Call DScal_(nCoor,One/dDisp,Disp(:,:),1)
      If (dGrad.gt.Zero)
     &  Call DScal_(nCoor,One/dGrad,Grad(:,:),1)
*     Any zero vector is assumed to be parallel to any other
      If (dPostDir*dGrad.gt.Zero) Then
        dPostDirGrad=dPostDirGrad/(dPostDir*dGrad)
      Else
        dPostDirGrad=One
      End If
      If (dPrevDir*dGrad.gt.Zero) Then
        dPrevDirGrad=dPrevDirGrad/(dPrevDir*dGrad)
      Else
        dPrevDirGrad=One
      End If
      If (dPrevDir*dDisp.gt.Zero) Then
        dPrevDirDisp=dPrevDirDisp/(dPrevDir*dDisp)
      Else
        dPrevDirDisp=One
      End If
      If (dPrevDir*dPostDir.gt.Zero) Then
        dPrevDirPostDir=dPrevDirPostDir/(dPrevDir*dPostDir)
      Else
        dPrevDirPostDir=One
      End If
*     A negative curvature means there is no appropriate value
      Curvature=-1.0D-12
      PathLength=dDisp/Sqrt(TWeight)
      If (MEP.and.MEP_Type.eq.'SPHERE    ') Then
*       The curvature is the inverse radius of the circle tangent to
*       both the current and previous MEP points
*       The path length is the arc length between these points
        If (One-dPrevDirPostDir.gt.Zero)
     &    Curvature=(One-dPrevDirPostDir)/Sqrt(One-dPrevDirPostDir**2)/
     &              Sqrt(dPrevDir*dPostDir/TWeight)
        If (Curvature.gt.Zero)
     &    PathLength=aCos(dPrevDirPostDir)/Curvature
      End If
*
*     Store the length and curvature values, and print results
*
      Call mma_allocate(Len,nMEP+1,Label='Len')
      Call mma_allocate(Cur,nMEP+1,Label='Cur')
      If (iMEP.ge.1) Then
        Call Get_dArray('MEP-Lengths   ',Len,nMEP+1)
        Call Get_dArray('MEP-Curvatures',Cur,nMEP+1)
        If (IRC.eq.-1) Then
          Len(1+iMEP)=-PathLength
        Else
          Len(1+iMEP)=PathLength
        End If
        Cur(1+iMEP)=Curvature
        Call Put_dArray('MEP-Lengths   ',Len,nMEP+1)
        Call Put_dArray('MEP-Curvatures',Cur,nMEP+1)
        If ((iMEP.ge.1).and.(iPrint.ge.5)) Then
          If (MEP_Type.eq.'TRANSVERSE') Then
            ConstraintAngle=aCos(dPrevDirGrad)*RadToDeg
          Else
            ConstraintAngle=aCos(dPostDirGrad)*RadToDeg
          End If
          If (MEP) Then
            PathAngle=aCos(dPrevDirPostDir)*RadToDeg
          Else
            PathAngle=aCos(-dPrevDirDisp)*RadToDeg
          End If
          ResGrad=dGrad
          Write(6,*)
          Write(6,'(a)')
     &      ' Last IRC/MEP step'//
     &      ' (in weighted coordinates / sqrt(total weight))'
          Write(6,'(a)')
     &      ' --------------------------------------------------------'
          Write(6,100) 'Residual gradient size:',ResGrad,'hartree/bohr'
          Write(6,100) 'Angle with constraint surface:',
     &                 ConstraintAngle,'degrees'
          Write(6,100) 'Path angle:',PathAngle,'degrees'
          If (Curvature.ge.Zero) Write(6,100)
     &                 'Path curvature:',Curvature,'bohr^(-1)'
          Write(6,100) 'Path length:',PathLength,'bohr'
100       Format(1X,A30,1X,F12.6,1X,A)
        End If
      Else
        Len(:)=Zero
        Cur(:)=Zero
        Call Put_dArray('MEP-Lengths   ',Len,nMEP+1)
        Call Put_dArray('MEP-Curvatures',Cur,nMEP+1)
      End If
      Call mma_deallocate(Len)
      Call mma_deallocate(Cur)
*
*     Do not mess with the geometry or reference if the next iteration
*     will be the start of a reverse IRC search.
*
      If (.Not.IRCRestart) Then
*
*       The new direction for the MEP should be the gradient in weighted
*       coordinates, but this may break additional constraints if they
*       exist, and the gradient may be close to zero.
*       Instead, we will use the direction that, on convergence, should
*       be parallel to the gradient when there are no other constraints.
*       Note that we could be following the gradient uphill
*
        If (MEP_Type.eq.'TRANSVERSE') Then
*         In the TRANSVERSE case, PrevDir is the vector parallel to the
*         gradient, but following this will never change the hyperplane
*         orientation.
*         We will try using a linear combination of PrevDir and Disp
*           dp  = Disp.PrevDir
*           Dir = dp*Disp + a*(PrevDir-dp*Disp)
*               = a*PrevDir + (1-a)*dp*Disp
*         Try different values for Fact (the "a" above).
*         I believe the true MEP direction should be more slanted from
*         the plane normal (PrevDir) than the Disp vector, therefore
*         a negative value is probably better.
          Fact=-0.5D0
          Dir(:,:) = Fact * PrevDir(:,:)
          Call DaXpY_(nCoor,(One-Fact)*dPrevDirDisp,Disp(:,:),1,
     &                                             Dir(:,:),1)
        Else
*         In the SPHERE case, PostDir is the vector to use
          Dir(:,:) = PostDir(:,:)
        End If
*
*       Special cases
*
        If (iMEP.eq.0) Then
          If (IRC.eq.0) Then
*           In the initial iteration of a MEP, use the initial direction
            Call Get_dArray('Transverse',Dir(:,:),nCoor)
          Else
*           In the initial iteration of an IRC branch, use the reaction vector
            Call dCopy_(nCoor,MF,1,Dir(:,:),1)
          End If
        End If
*
*       Project any additional constraints out of the direction vector
*       The constraint vectors are read from the dRdX file
*
        LudRdX=30
        Call DaName(LudRdX,'dRdX')
        iAd=0
        Call iDaFile(LudRdX,2,iDum,1,iAd)
        nLambda_=iDum(1)
        Call iDaFile(LudRdX,2,iDum,1,iAd)
        nCoor_=iDum(1)
        Call mma_allocate(drdx,nCoor_,nLambda_,Label='drdx')
        Call dDaFile(LudRdX,2,drdx,nLambda_*nCoor_,iAd)
        Call DaClos(LudRdX)
        Do iLambda=1,nLambda
          If (iLambda.ne.MEPnum) Then
            dd=dDot_(nCoor,drdx(:,iLambda),1,drdx(:,iLambda),1)
            drd=dDot_(nCoor,drdx(:,iLambda),1,Dir(:,:),1)
            Call DaXpY_(nCoor,-drd/dd,drdx(:,iLambda),1,Dir(:,:),1)
          End If
        End Do
        Call mma_deallocate(drdx)
*
*       Compute the length of the direction vector in weighted coordinates
*
        dDir=Zero
        iOff=0
        Do iAtom=1,nAtom
          Fact=Dble(iDeg(Cx(1+iOff,iter)))
          xWeight=Weights(iAtom)
          Do ixyz=1,3
            dDir=dDir+Fact*xWeight*Dir(ixyz,iAtom)**2
            iOff=iOff+1
          End Do
        End Do
        dDir=Sqrt(dDir)
*
*       According to the Gonzalez-Schlegel method, the reference point
*       is half-step away in the search direction.
*       First set the new reference point (except for rMEP) and then
*       compute the new starting structure at the full step distance
*       For an IRC first step, keep the initial structure as reference
*
        Fact=dMEPStep*Sqrt(TWeight)/dDir
        Call dCopy_(nCoor,Cx(:,iter),1,Cen(:,:),1)
        If (MEP_Algo.eq.'GS') Then
          If ((IRC.eq.0).or.(iMEP.ne.0))
     &      Call Find_Distance(Cx(:,iter),Cen(:,:),Dir(:,:),
     &                Half*Fact,Half*dMEPStep,nAtom,BadConstraint)
          If (.Not.rMEP) Call Put_dArray('Ref_Geom',Cen(:,:),nCoor)
          Call Find_Distance(Cen(:,:),Cx(:,iter+1),Dir(:,:),
     &              Half*Fact,Half*dMEPStep,nAtom,BadConstraint)
        Else If (MEP_Algo.eq.'MB') Then
          If (.Not.rMEP) Call Put_dArray('Ref_Geom',Cx(1,iter),nCoor)
          Call Find_Distance(Cx(:,iter),Cx(:,iter+1),Dir(:,:),
     &              Fact,dMEPStep,nAtom,BadConstraint)
        End If
*
*       Randomly displace the new geometry 1/20th of the MEP distance
*       to try to break symmetry
*
        If (nIrrep.eq.1) Then
          Call Random_Vector(nCoor,Disp(:,:),.True.)
          dDir=Zero
          iOff=0
          Do iAtom=1,nAtom
            xWeight=Weights(iAtom)
            Do ixyz=1,3
              dDir=dDir+xWeight*Disp(ixyz,iAtom)**2
              iOff=iOff+1
            End Do
          End Do
          Fact=dMEPStep*Sqrt(TWeight/dDir)
          Call DaXpY_(nCoor,0.05D0*Fact,Disp(:,:),1,Cx(:,iter+1),1)
        End If
        Call Put_dArray('Transverse',Dir(:,:),nCoor)
        If (iter.eq.1) BadConstraint=.False.
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_deallocate(PrevDir)
      Call mma_deallocate(PostDir)
      Call mma_deallocate(Disp)
      Call mma_deallocate(Grad)
      Call mma_deallocate(Dir)
      Call mma_deallocate(Cen)
      Return
      End
