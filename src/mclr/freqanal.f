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
      Subroutine Freqanal(nDeg,nrvec,H,converged,
     &                    ELEC,iel,elout,ldisp,Lu_10)
      use stdalloc, only: mma_allocate, mma_deallocate
      use input_mclr, only: nSym,nDisp,nUserPT,nSRot,UserP,ChIrr,UserT
      use temperatures, only: DefTemp
      Implicit None
      Integer nDeg(*),nrvec(*)
      Real*8 H(*)
      Logical converged(8)
      Real*8 elec(*)
      Integer iel(3)
      Real*8 elout(*)
      Integer ldisp(nsym), Lu_10

*     local variables
      logical Do_Molden
      Real*8, Allocatable:: NMod(:), EVec(:), EVec2(:,:), EVal(:),
     &                      EVal2(:), Intens(:), RedMas(:), Tmp3(:),
     &                      Temp(:)
      Integer ipNx,nModes,lModes,i1,i3,j,ii,lnm_molpac,iSym,nx,iCtl,
     &        ll,kk,i,k,iT,jpNx,ix,jx,nEig,iNeg
      Integer, external:: IsFreeUnit
      Real*8 Tmp,Fact,rNorm
*
      Call mma_allocate(NMod,nDisp**2,Label='NMod')
      Call mma_allocate(EVec,nDisp**2,Label='EVec')
      Call mma_allocate(EVec2,2,nDisp**2,Label='EVec2')
      Call mma_allocate(EVal,nDisp,Label='EVal')
      Call mma_allocate(Intens,nDisp*2,Label='Intens')
      Call mma_allocate(RedMas,nDisp,Label='RedMas')
      ipNx=1
      nModes=0
      lModes=0
*
      Write (6,*)
      Write (6,*) '     ************************************'
      Write (6,*) '     *                                  *'
      Write (6,*) '     * Harmonic frequencies in cm-1     *'
      Write (6,*) '     * Intensities in km/mole           *'
      Write (6,*) '     *                                  *'
      Write (6,*) '     * No correction due to curvilinear *'
      write (6,*) '     * representations has been done    *'
      Write (6,*) '     *                                  *'
      Write (6,*) '     ************************************'
      Write (6,*)
      i1=1
      i3=1
      j=0
      ii=1
      Write(Lu_10,'(A)') '*PERTURBATIONS'
      Write(Lu_10,*) ldisp
      Write(Lu_10,'(A)') '*BEGIN NORMAL MODES'
      Write(Lu_10,'(A)') '*NOTICE THAT THEY ARE SYMMETRY ADAPTED'
      Write(Lu_10,'(A)') '*USING ORTHOGONAL TRANSFORMATIONS '
      WRITE(Lu_10,'(A)') '*AND NOT ALASKA TYPE'
*

!*    !> open normal mode file for "normal mode molpac" ! yma
      lnm_molpac=60
      lnm_molpac=isFreeUnit(lnm_molpac)
      Call Molcas_Open(lnm_molpac,'normal_modes_molpac')


      Do_Molden=.True.
      Do iSym=1,nSym
         nX=ldisp(isym)
         If (nX.ne.0) Then
            Write(6,*)
            Write(6,*) '   Symmetry ',chirr(isym)
            Write(6,*) '  =============='
            Write(6,*)
*
            If (converged(isym))  Then
               Call mma_allocate(EVal2,2*nX,Label='EVal2')
               Call mma_allocate(Tmp3,nX**2,Label='Tmp3')
               Call FREQ(nX,H(i3),nDeg(i1),nrvec(i1),
     &                   Tmp3,EVec2,EVal2,RedMas,iNeg)
               EVec(:) = EVec2(1,:)
               EVal(i1:i1+nX-1) = EVal2(1:nX)
               Call mma_deallocate(EVal2)
               Call mma_deallocate(Tmp3)
*
               iCtl=0
               ll=0
               kk=j+1
               Do i=1,3
                  If (iel(i).eq.isym) Then
                     iCtl=1
                     ll=ll+1
                     Do k=1,nx
                        j=j+1
                        tmp=0.0d0
                        Do it=0,nx-1
                           Fact=Sqrt( DBLE(nDeg(i1+it)) )
                           tmp=tmp+EVec(1+(k-1)*nx+it)*
     &                             elec(ii+it)*Fact
                        End Do
                        elout(j)=tmp
                     End Do
                     ii=ii+nx
                  End If
               End Do
               Write(Lu_10,'(A,I1)') '*NORMAL MODES SYMMETRY: ',isym
*
*------------- Save normal modes for later generation of Molden input.
*
               call dcopy_(nX**2,EVec,1,NMod(ipNx),1)
               jpNx=ipNx

               Do iX = 1, nX
*
*                 Transform from mass-weighted cartesian to cartesian for
*                 Molden.
*
                  rNorm=0.0D0
                  Do jX = 0, nX-1
                     Fact=Sqrt(DBLE(nDeg(jX+1)))
                     NMod(ipNx+jX) = NMod(ipNx+jX)/Fact
                     rNorm=rNorm+DBLE(nDeg(jX+1))*NMod(ipNx+jX)**2
                  End Do
                  Call DScal_(nX,1.0D0/Sqrt(rNorm),NMod(ipNx),1)
*
                  ipNx=ipNx+nX
                  lModes=lModes+nX
               End Do
               nModes=nModes+nX
               call dcopy_(nX**2,NMod(jpNx),1,EVec,1)
               Call GF_Print(EVal(i1),EVec,elout(kk),
     &                       ll,nX,nX,iCtl,Intens(i1),
     &                       RedMas,Lu_10,i1-1)
            Else
               Write(6,*)
               Write (6,*)'     NOT CONVERGED'
               Write(6,*)
               Do i=1,3
                  If (iel(i).eq.isym) Then
                     j=j+1
                     ii=ii+nx
                     elout(j)=-99999999D0
                  End If
               End Do
               Do_Molden=.False.
            End If
         End If
         i3=i3+nx*(nx+1)/2
         i1=i1+nx
      End Do
      nEig = i1 - 1
*
!*     !> close the normal mode file
       close(lnm_molpac) ! This is for normal_modes_molpac -- yma
!      call NM_MOPAC_SNF(nsym,ldisp,natoms) ! f90 not support ....

      If (nsym.eq.1) Then
         Call Print_Mode_Components(NMod,EVal,
     &                              nModes,lModes,lDisp)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call mma_allocate(Temp,nEig,Label='Temp')
      call dcopy_(nEig,Eval,1,Temp,1)
*
*     For verification purpose we skip frequencies close to zero.
*
      Do i = 1, nEig
         If (Abs(Temp(i)).lt.5.0D0) Temp(i)=0.0D0
      End Do
      Call Add_Info('Harm_Freq',Temp,nEig,1)
      Call mma_deallocate(Temp)
*
      Do i = 1, nEig
         If (Abs(Intens(i)).lt.1.0D0) Intens(i)=0.0D0
      End Do
      Call Add_Info('IR_Intensities',Intens,nEig,1)
*                                                                      *
************************************************************************
*                                                                      *
      Write(Lu_10,'(A)') '*END NORMAL MODES'
*
*------------- Calculate thermodynamic properties----------
*
      If (nUserPT.eq.0 .and. nsRot.eq.0) then
        UserP=1.0d0
        nUserPT=Size(DefTemp)
        Do i=1,nUserPT
          UserT(i)=DefTemp(i)
        End Do
*       Call ThermoData(EVal,nEig)
      EndIf
      Call Thermo_Driver(UserT,UserP,nUserPT,nsRot,EVal,nEig,.False.)
*
*
*---- Write stuff on Molden input file
*
      If (Do_Molden)
     &   Call Freq_Molden(EVal,nModes,NMod,lModes,nSym,
     &                    Intens,lDisp,RedMas)
*
      Call mma_deallocate(NMod)
      Call mma_deallocate(evec)
      Call mma_deallocate(evec2)
      Call mma_deallocate(eval)
      Call mma_deallocate(intens)
      Call mma_deallocate(redmas)
*
      End Subroutine Freqanal
