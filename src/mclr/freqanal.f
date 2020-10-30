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
      Implicit Real*8(a-h,o-z)
#include "Input.fh"
#include "stdalloc.fh"
      Logical converged(8)
      Real*8 H(*),elec(*),elout(*)
      logical Do_Molden
      Integer nrvec(*),nDeg(*),iel(3),ldisp(nsym)
      Real*8, Allocatable:: NMod(:), EVec(:), EVal(:), Intens(:),
     &                      RedMas(:), Tmp2(:), Tmp3(:), Temp(:)
#include "temperatures.fh"
*
      Call mma_allocate(NMod,nDisp**2,Label='NMod')
      Call mma_allocate(EVec,2*nDisp**2,Label='EVec')
      Call mma_allocate(EVal,nDisp*2,Label='EVal')
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
               naux=Max(nx*2,nX**2)
               Call mma_allocate(Tmp2,naux,Label='Tmp2')
               Call mma_allocate(Tmp3,nX**2,Label='Tmp3')
               Call FREQ(nX,H(i3),nDeg(i1),nrvec(i1),
     &                   Tmp2,Tmp3,EVec,EVal(i1),RedMas,iNeg)
               Call mma_deallocate(Tmp3)
               Call mma_deallocate(Tmp2)
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
                           tmp=tmp+EVec(1+2*(k-1)*nx+2*it)*
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
               call dcopy_(nX**2,EVec,2,NMod(ipNx),1)
               jpNx=ipNx

* =========================================================================
*                mass-weighted normal mode print out
* =========================================================================
               Call NM_MOPAC_print(EVal(i1),EVec,
     &              elout(kk),ll,nX,nX,iCtl,Intens(i1),Lu_10,
     &              i1-1,lnm_molpac)
* =========================================================================

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
               call dcopy_(nX**2,NMod(jpNx),1,EVec,2)
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
        nUserPT=NDefTemp
        Do i=1,NDefTemp
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
      Call mma_deallocate(eval)
      Call mma_deallocate(intens)
      Call mma_deallocate(redmas)
*
      Return
      End

!      Subroutine NM_MOPAC_SNF(nsym,ldisp,nAtom)
! These used be a f90 file.. (now deleted dut to not work) -- yingjin
!
!        integer :: nsym, nAtom
!        integer :: ldisp(nsym)
!
!        write(*,*)nsym,natom
!        write(*,*)(ldisp(i),i=1,nsym)
!
!      end subroutine NM_MOPAC_SNF

*     Print normal modes in MOPAC format with symm -- yingjin
*     Only for ground state (no imaginary freqs)
      Subroutine NM_MOPAC_print(EVal,EVec,dDipM,iel,nX,nDim,ictl,IRInt,
     &                          Lu_10,iOff,lut)

      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "constants2.fh"
      Real*8 EVal(nDim), EVec(2,nX,nDim),dDipM(ndim,iel),IRInt(nDim)
      Parameter(Inc=6)
      Character*80 Format, Line*120
      Character*(LENIN6) ChDisp(3*MxAtom),Label
      character*(10) char_num
      character charx

*
*      LUt=lnm_molpac
*
      Call Get_iScalar('nChDisp',nChDisp)
      If (nChDisp.lt.nX) Then
         Write(LUt,*) 'nm_MOPAC_Print: nChDisp.lt.nX!'
         Call Abend()
      End If
      Call Get_cArray('ChDisp',ChDisp,(LENIN6)*nChDisp)
*
      !> omit the eigenvalues which less than 5 cm^-1 (abs) -- yma
      idiscard=0
      imaginary=0
      Do iHarm = 1, Inc
*         write(*,*)"Eval(",iHarm,") = ",Eval(iHarm) ! for checking
         if(abs(Eval(iHarm)).lt.5.0)then
           idiscard=idiscard+1
         else if((Eval(iHarm)).lt.-5.0)then
           imaginary=imaginary+1
         end if
      end do

      write(LUt,*)
      write(LUt,*)
     &" Frequencies and mass-weighted normal coordinates " ! yma
      write(LUt,*)
     &" ================================================ "
      write(LUt,*)
      write(LUt,*)
     &"          eigenvalues  in  cm^-1                  "
      write(LUt,*)

      Do iHarm = 1+idiscard, nDim, Inc
         Jnc=Min(Inc,nDim-iHarm+1)
         Label='ROOT NO.: '
         Write(Format,'(A,I3,A)') '(1X,A10,',Jnc,'(I10))'
         Write (LUt,Format) Label,(i-idiscard,i=iHarm,iHarm+Jnc-1)
         Label='EIGVAL. : '
         Write(Format,'(A,I3,A)') '(1X,A10,2X,',Jnc,'(f10.3))'
         Write (LUt,Format) Label,(EVal(i),i=iHarm,iHarm+Jnc-1)
      Write (LUt,*)
      end do
      Write (LUt,*)
      Write (LUt,*)
      Write (LUt,*)

      char_num="0123456789"
      iIRInt=0
      Do iHarm = 1+idiscard, nDim, Inc
           Jnc=Min(Inc,nDim-iHarm+1)
           Label='root no.'
           Write(Format,'(A,I3,A)') '(4X,A8,',Jnc,'(I5,7X))'

           Write (LUt,Format) Label,(i-idiscard,i=iHarm,iHarm+Jnc-1)
           Write (LUt,*)
*
           Write(Format,'(A,I3,A)') '(8x,',Jnc,'F12.5)'
           Line=' '
           Write (LUt,Format) (EVal(i),i=iHarm,iHarm+Jnc-1)
           Write (LUt,*)

           If (ictl.ne.0) Then
           Else
              Do i=1,Jnc
                 iIRInt=iIRInt+1
                 IRInt(iIRInt)=Zero
              enddo
           End if
*
           Do iInt = 1, nX
              inum_min=99
              inum_max=0
              do i=1,10
                charx=char_num(i:i)
                inum1=index(ChDisp(iInt+iOff)(1:LENIN6),charx)
                inum2=index(ChDisp(iInt+iOff)(1:LENIN6),charx,.true.)
                if(inum_min.gt.inum1.and.inum1.ne.0)inum_min=inum1
                if(inum_max.lt.inum2)inum_max=inum2
              end do
            Write(Format,'(A,I3,A)') '(1X,A,2x,A,A,A,',Jnc,'(F10.5,2x))'

              Write (LUt,Format) ChDisp(iInt+iOff)(LENIN6:LENIN6),
     &               ChDisp(iInt+iOff)(1:inum_min-1),
     &               ChDisp(iInt+iOff)(inum_max+1:4+inum_min),
     &               ChDisp(iInt+iOff)(inum_min:inum_max),
     &               (EVec(1,iInt,i),
     &               i=iHarm,iHarm+Jnc-1)

           End Do
           Write (LUt,*)
           Write (LUt,*)
           Write (LUt,*)
           Write (LUt,*)
           Write (LUt,*)
      End Do
           write(LUt,*)
           write(LUt,*)


      Return

c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(dDipM)
         Call Unused_integer(Lu_10)
      End If
      End Subroutine NM_MOPAC_print


