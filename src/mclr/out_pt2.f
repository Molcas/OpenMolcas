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
       SubRoutine Out_Pt2(iKapDisp,iCIDisp)
********************************************************************
*                                                                  *
********************************************************************
       Implicit Real*8 (a-h,o-z)
#include "detdim.fh"

#include "Input.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "disp_mclr.fh"
#include "cicisp_mclr.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "sa.fh"
#include "dmrginfo_mclr.fh"
#include "SysDef.fh"
       Logical CI, Is_Roots_Set
       Character*80 Note
! Added for DMRG calculation
       real*8,allocatable::tmpDe(:,:),tmpP(:),tmpDeM(:,:),tmpPM(:,:,:,:)
       Integer iKapDisp(nDisp),iCiDisp(nDisp)
       Character(Len=16) mstate
       Dimension rdum(1),idum(7,8)
*                                                                      *
************************************************************************
*                                                                      *
       itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*                                                                      *
************************************************************************
*                                                                      *
       Call QEnter('Out_PT2')
       isym=1
       CI=.true.
       Call Setup_MCLR(iSym)
       nbas_tot=0
       ntot1=0
       nDLMO=0
       nLCMO=0
       Do is=1,nsym
          nbas_tot=nbas_tot+nbas(is)
          ntot1=ntot1+nbas(is)*(nbas(is)+1)/2
          nDLMO=nDLMO+nash(is)
          nLCMO=nLCMO+nbas(is)*nbas(is)
       End Do
       nDLMO=nDLMO*(nDLMO+1)/2
       nPLMO=nDLMO*(nDLMO+1)/2
*
       Call GetMem('kappa1','Allo','Real',ipK1,  nDens2)
       Call Getmem('kappa2','ALLO','Real',ipk2,  nDens2)
       Call GetMem('ONEDEN','Allo','Real',ipDAO,nDens2)
       Call GetMem('ONEDEN','Allo','Real',ipD_CI,n1Dens)
       Call GetMem('ONEDEN','Allo','Real',ipD1,n1Dens)
       Call Getmem('TWODEN', 'ALLO','Real',ipP_CI,n2Dens)
       Call Getmem('TWODEN', 'ALLO','Real',ipP1,n2Dens)
       Call Getmem('Conn', 'ALLO','Real',ipF,nDens2)
       Call Getmem('OCCU ', 'ALLO','Real',ipO,nbas_tot)
       Call Getmem('CMO', 'ALLO','Real',ipCMON,ndens2)
*      OBS nBuf might not be def.
       Call Getmem('TMP', 'MAX','Real',ipT,nBuf)
       Call GetMem('TMPDEN','Allo','Real',ipDtmp,nDens2)

*
*
*
* 1)   CI Part
*
*      All multipliers are introduced as densities
*

       If (CI) Then
         nconf1=ncsf(State_sym)
         ilen=nconf1*nroots ! nroot = # of roots in SA
         ipcip=ipget(nconf1*nroots)
         iDisk=iCIDisp(1)
         Call dDaFile(LuTemp,2,Work(ipin(ipCIp)),iLen,iDisk)
*
*-------Calculate the densities that correct the nonvariational CI stuff
*
         Call CIDens_sa(.true.,ipCIp,ipCI,
     &                  State_sym,State_sym,
     &                  Work(ipP_CI),Work(ipD_CI)) ! \bar{d} and \bar{D}

! ===========================================================================
         if(doDMRG)then  ! yma
           call dmrg_dim_change_mclr(LRras2(1:8),ntash,0)
           call dmrg_dim_change_mclr(RGras2(1:8),ndim,0)

           ! use mma_allocate later - yma
           allocate(tmpDe(ndim,ndim))
           allocate(tmpP(ndim**2*(ndim**2+1)/2))
           allocate(tmpDeM(ntash,ntash))
           allocate(tmpPM(ntash,ntash,ntash,ntash))
           tmpDe=0.0d0
           tmpP=0.0d0
           tmpDeM=0.0d0
           tmpPM=0.0d0

           ij=0
           do i=1,ntash
             do j=1,ntash
               ij=ij+1
               if(abs(Work(ipD_CI+ij-1)).lt.1.0e-12)then
                 Work(ipD_CI+ij-1)=0.0d0
               end if
               tmpDeM(i,j)=Work(ipD_CI+ij-1)
             end do
           end do

           ij=0
           do i=1,ndim
             do j=1,ndim
               ij=ij+1
               if(i.gt.ntash.or.j.gt.ntash)then
                 tmpDe(i,j)=0.0d0
               else
                 tmpDe(i,j)=tmpDeM(i,j)
               end if
             end do
           end do

           Do i=1,ntash
             Do j=1,ntash
               Do k=1,ntash
                 Do l=1,ntash
                   ij1=ntash*(i-1)+j
                   ij2=ntash*(j-1)+i
                   kl1=ntash*(k-1)+l
                   kl2=ntash*(l-1)+k
                   if(ij1.ge.kl1)then
                   if(abs(Work(ipP_CI+itri(ij1,kl1)-1)).lt.1.0e-12)then
                     Work(ipP_CI+itri(ij1,kl1)-1)=0.0d0
                   end if
                   tmpPM(i,j,k,l)=Work(ipP_CI+itri(ij1,kl1)-1)
                   end if
                 End Do
               End Do
             End Do
           End Do

           do i=1,ndim
             do j=1,ndim
               do k=1,ndim
                 do l=1,ndim
                   ij1=ndim*(i-1)+j
                   ij2=ndim*(j-1)+i
                   kl1=ndim*(k-1)+l
                   kl2=ndim*(l-1)+k
                   if(ij1.ge.kl1)then
         if(i.gt.ntash.or.j.gt.ntash.or.k.gt.ntash.or.l.gt.ntash)then
                       tmpP(itri(ij1,kl1))=0.0d0
                     else
                       tmpP(itri(ij1,kl1))=tmpPM(i,j,k,l)
                     end if
                   end if
                 end do
               end do
             end do
           end do

           ij=0
           do i1=1,ndim
             do j1=1,ndim
               ij=ij+1
               Work(ipD_CI+ij-1)=tmpDe(i1,j1)
             end do
           end do
           do i=1,n2dens
             Work(ipP_CI+i-1)=tmpP(i)
           end do
           ! mma_deallocate later
           deallocate(tmpDe)
           deallocate(tmpDeM)
           deallocate(tmpP)
           deallocate(tmpPM)
           call dmrg_dim_change_mclr(RGras2(1:8),
     &                               ntash,0)
         end if
! ===================================================================

*
*-------Some administrative shit
*
*       Store densities in triangular form
*
         Do i=1,ntAsh
          Do j=1,i
           Work(ipD1+itri(i,j)-1)=Work(ipD_CI+(i-1)*ntash+j-1)
          End Do
         End Do

         Do i=1,ntAsh
          Do j=1,i
           ij=itri(i,j)
           ij2=i+(j-1)*ntash
           ji2=j+(i-1)*ntash
           Do k=1,i
            Do l=1,k
             kl=itri(k,l)
             kl2=k+(l-1)*ntash
             lk2=l+(k-1)*ntash
             ijkl=itri(ij2,kl2)
             jikl=itri(ji2,kl2)
             ijlk=itri(ij2,lk2)
             jilk=itri(ji2,lk2)
             Work(ipP1+itri(ij,kl)-1)=Quart*
     &             (Work(ipP_CI+ijkl-1)+Work(ipP_CI+jikl-1)+
     &              Work(ipP_CI+ijlk-1)+Work(ipP_CI+jilk-1))
            End Do
           End Do
          End Do
         End Do

        DO  K=1,NTASH
         DO L = 1, K
          KL = K*(K-1)/2 + L
          KLROW = KL*(KL-1)/2
          IF( L .EQ. K ) THEN
           IMAX = K
          ELSE
           IMAX = K-1
          END IF
          DO I = 1,IMAX
           II= I*(I+1)/2
           IIKL= KLROW + II
           Work(ipP1+IIKL-1) = Work(ipP1+IIKL-1)*Half
          End Do
         End Do
         End Do


C         Do i=1,ntAsh
C         Do j=1,i
C         ij=itri(i,j)
C         ij2=i+(j-1)*ntash
C         ji2=j+(i-1)*ntash
C         Do k=1,ntAsh
C         Do l=1,k
C          kl=itri(k,l)
C          kl2=k+(l-1)*ntash
C          ijkl=itri(ij2,kl2)
C          jikl=itri(ji2,kl2)
C          fact=Half
C          if(ij.ge.kl .and. k.eq.l) fact=Quart
C          if(ij.lt.kl .and. i.eq.j) fact=Quart
C          Work(ipP1+itri(ij,kl)-1)=
C     &        fact*(Work(ipP_CI+ijkl-1)+Work(ipP_CI+jikl-1))
C         End Do
C         End Do
C         End Do
C         End Do
C         If (debug) Call triprt('LP',' ',Work(ipp1),(ntash**2+ntash)/2)
c
c Write the 'bar' densities to disk,  not symmetry blocked.
c

!         Call Put_DLMO(Work(ipD1),ndim1) ! \bar{D} triangular  ! yma
!         Call Put_PLMO(Work(ipP1),ndim2) ! \bar{d} triangular  ! yma

         Call Put_DLMO(Work(ipD1),nDLMO) ! \bar{D} triangular ! original
         Call Put_PLMO(Work(ipP1),nPLMO) ! \bar{d} triangular ! original
*
       End If
*
*      2) Orbital response
*         ================
*
*       Read in from disk
*
       iDisk=iKapDisp(1)
       Call dDaFile(LuTemp,2,Work(ipK1),nDensC,iDisk) ! Read \bar{kappa}
       Call Uncompress(work(ipK1),Work(ipK2),1)
c
c If we want to estimate the error
c
       If (esterr) Then
*        Do iestate=1,lroots
*          Call calcerr(Work(ipK2),iestate)
*        End do
         Call calcerr(Work(ipK2),istate)
       End If
*
*----- First we fix the renormalization contribution
*
       Call Get_Fock_Occ(ipD_K,Length)
*      Calculates the effective Fock matrix
       Call Make_Conn(Work(ipF),Work(ipK2),
     &                Work(ipP_CI),work(ipD_CI))   !ipD_CI not changed
       Call DaxPy_(ndens2,One,Work(ipD_K),1,Work(ipF),1)
*      call dcopy_(ndens2,Work(ipD_K),1,Work(ipF),1)
       Call Put_Fock_Occ(Work(ipF),nTot1)
*
*      Transposed one index transformation of the density
*      (only the inactive density to store it separately)
*
       Call OITD(Work(ipK2),1,Work(ipDAO),Work(ipDtmp),.False.)
*
*      Transformation to AO basis (covariant)
*
c
c Transforms to AO differently dep on last arg.
c
       Call TCMO(Work(ipDAO),1,-2)
*
*      Fold AO density and write to disk
c Mult all terms that are not diag by 2
*
       Call FOLD2(nsym,nbas,Work(ipDAO),Work(ipK1))
*
       Call Put_DLAO(Work(ipk1),ntot1)
*
*      Now with active density too, to form the variational density
*
!      gives \tilde{D}
       Call OITD(Work(ipK2),1,Work(ipD_K),Work(ipDtmp),.True.)
*
       Do iS=1,nsym
c
c C*\tilde{\kappa} --> ipDAO
c
          If (nBas(is).ge.1)
     &       CALL DGEMM_('N','N',
     &                   NBAS(is),NBAS(is),NBAS(is),
     &                   One,Work(ipCMO+ipCM(is)-1),NBAS(is),
     &                   Work(ipK2+ipmat(is,is)-1),NBAS(is),
     &                   Zero,Work(ipDAO+ipCM(is)-1),NBAS(is))
       End Do
*
       Call Put_LCMO(Work(ipDAO),nLCMO)
*
       if(doDMRG)then  ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),ntash,0)
         call dmrg_spc_change_mclr(RGras2(1:8),nash)
       end if
*
       If (isNAC) Then
         Call Get_D1MO(ipG1q,ng1)
         iR = 0 ! set to dummy value.
       Else
         iR=iroot(istate)
         jdisk=itoc(3)
         ng1=itri(ntash,ntash)
         ng2=itri(ng1,ng1)
         Call Getmem('TMP', 'ALLO','Real',ipG1q,n1dens)
c
c Read active one el dens for state j from JOBIPH and store in ipG1q
c
         Do i=1,iR-1  ! Dummy read until state j
           Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
           Call dDaFile(LUJOB ,0,rdum,ng1,jDisk)
           Call dDaFile(LUJOB ,0,rdum,ng2,jDisk)
           Call dDaFile(LUJOB ,0,rdum,ng2,jDisk)
         End Do
         Call dDaFile(LUJOB ,2,Work(ipG1q),ng1,jDisk)
       EndIf
*
*    Construct a variationally stable density matrix. In MO
c
c D_eff = D^j + \tilde{D} +\bar{D}
c ipD_K = (ipG1q + inact) + ipD_K + ipD_CI
*
C
C       call dcopy_(ndens2, [Zero], 0, Work(ipD_K), 1)  !DEBUG
C
       If (isNAC) Then
*
** For NAC, first build DAO and then DAO_var
*
         Do is=1,nSym
c Note: no inactive part for transition densities
          Do iA=1,nash(is)
           Do jA=1,nash(is)
            i=iA+nish(is)
            j=jA+nish(is)
            iAA=iA+na(is)
            jAA=jA+na(is)
            Work(ipD_K+ipmat(is,is)-1+i-1+(j-1)*nbas(is))=
     &       Work(ipD_K+ipmat(is,is)-1+i-1+(j-1)*nbas(is))
     &      +Work(ipD_CI+iAA-1+(jAA-1)*ntash)
     &      +Work(ipG1q-1+itri(iAA,jAA))
           End Do
          End Do
         End Do
         Call Getmem('TMP', 'ALLO','Real',ipT,nBuf/2)
         Call NatOrb(Work(ipD_K),Work(ipCMO),Work(ipCMON),Work(ipO))
         Call dmat_MCLR(Work(ipCMON),Work(ipO),Work(ipT))
         Call Put_D1ao_var(Work(ipT),nTot1)
         Call Getmem('TMP', 'FREE','Real',ipT,nBuf/2)
*
** Transform the antisymmetric transition density matrix to AO
**  (there is no guarantee the symmetry will work here)
*
         iDisk=0
         LuDens=20
         Call DaName(LuDens,'MCLRDENS')
         Call dDaFile(LuDens,2,Work(ipG1q),ng1,iDisk)
         Call DaClos(LuDens)
         Call Getmem('D1ao-','ALLO','Real',ipG1m,ndens2)
         Call DCopy_(ndens2,[Zero],0,Work(ipG1m),1)
* Reconstruct the square matrix
         Do is=1,nSym
          Do iA=1,nash(is)
           i=iA+nish(is)
           iAA=iA+na(is)
           Do jA=1,iA-1
            j=jA+nish(is)
            jAA=jA+na(is)
            Work(ipG1m+ipmat(is,is)-1+i-1+(j-1)*nbas(is))=
     &           Work(ipG1q-1+itri(iAA,jAA))
            Work(ipG1m+ipmat(is,is)-1+j-1+(i-1)*nbas(is))=
     &          -Work(ipG1q-1+itri(iAA,jAA))
           End Do
           Work(ipG1m+ipmat(is,is)-1+i-1+(i-1)*nbas(is))=Zero
          End Do
         End Do
* Transform
         Call TCMO(Work(ipG1m),1,-2)
* Save the triangular form
         iOff=ipG1m
         Do is=1,nSym
          ibas=nbas(is)
          Do i=1,ibas
           Do j=1,i
            Work(iOff-1+itri(i,j))=
     &        Work(ipG1m+ipmat(is,is)-1+j-1+(i-1)*nbas(is))
           End Do
          End Do
          iOff=iOff+(ibas*ibas+ibas)/2
         End Do
         Call Put_dArray('D1ao-',Work(ipG1m),nTot1)
         Call Getmem('D1ao-','FREE','Real',ipG1m,nTot1)
*
       Else
*
** Normal SA gradient (no NAC)
*
         Do is=1,nSym
          Do i=1,nish(is)
c
c The inactive density
c
           Work(ipD_K+ipmat(is,is)-1+i-1+(i-1)*nbas(is))=
     &     Work(ipD_K+ipmat(is,is)-1+i-1+(i-1)*nbas(is))+Two
          End DO
          Do iA=1,nash(is)
           Do jA=1,nash(is)
            i=iA+nish(is)
            j=jA+nish(is)
            iAA=iA+na(is)
            jAA=jA+na(is)
c
c The active density ipG1q and \bar{D}
c
            Work(ipD_K+ipmat(is,is)-1+i-1+(j-1)*nbas(is))=
     &       Work(ipD_K+ipmat(is,is)-1+i-1+(j-1)*nbas(is))
     &      +Work(ipD_CI+iAA-1+(jAA-1)*ntash)
     &      +Work(ipG1q-1+itri(iAA,jAA))
           End Do
          End Do
         End Do
c
c Diagonalize the effective density to be able to use Prpt
c ipO eigenvalues of eff dens
c ipCMON eigenvectors (new orb coef)
c
         Call Getmem('TMP', 'ALLO','Real',ipT,nBuf/2)
         Call NatOrb(Work(ipD_K),Work(ipCMO),Work(ipCMON),Work(ipO))
         Call dmat_MCLR(Work(ipCMON),Work(ipO),Work(ipT))
         Call Put_D1ao_Var(Work(ipT),nTot1)
         Call Getmem('TMP', 'FREE','Real',ipT,nBuf/2)
         Call get_D1MO(ipT,nTot1)
         Call get_DLMO(ipTt,nTot1)
         Call DaxPy_(nTot1,1.0d0,Work(ipTt),1,Work(ipT),1)
         Call Getmem('DENS','FREE','Real',ipT,nTot1)
         Call Getmem('DENS','FREE','Real',ipTt,nTot1)
         Note='var'
         LuTmp=50
         LuTmp=IsFreeUnit(LuTmp)
         Call WrVec('TMPORB',LuTmp,'O',nSym,nBas,nBas,
     &            rDum,Work(ipO),rDum,iDum,Note)
         Call Prpt()
*                                                                      *
************************************************************************
*        There should now be dipole moments on the runfile which
*        corresponds to the gradient of the energy w.r.t. the
*        electric field. Let's update the list of values stored
*        on the runfile.
*
         Is_Roots_Set = .False.
         Call Qpg_iScalar('Number of roots',Is_Roots_Set)
         nRoots = 1
         If (Is_Roots_Set) Then
            Call Get_iScalar('Number of roots',nRoots)
         End If
*
         If (nRoots.ne.1) Then
*           Write (*,*) 'iR=',iR
            Call GetMem('DIPM', 'Allo','Real',ipDM,3)
            Call GetMem('DIPMs','Allo','Real',ipDMs,3*nROOTS)
            Call Get_dArray('Last Dipole Moments',Work(ipDMs),3*nRoots)
*           Call RecPrt('Last Dipole Moments',' ',Work(ipDMS),3,nRoots)
            Call Get_dArray('Dipole Moment',Work(ipDM),3)
*           Call RecPrt('Dipole Moment',' ',Work(ipDM),1,3)
            Call DCopy_(3,Work(ipDM),1,
     &                    Work(ipDMS+(iR-1)*3),1)
*           Call RecPrt('Last Dipole Moments',' ',Work(ipDMS),3,nRoots)
            Call Put_dArray('Last Dipole Moments',Work(ipDMs),3*nRoots)
            Call Free_Work(ipDMs)
            Call Free_Work(ipDM)
         End If
************************************************************************
*                                                                      *
       End If
       Call Getmem('TMP', 'FREE','Real',ipG1q,ng1)
C
c--------------------------  debug -----
c
       if(doDMRG)then ! yma
         call dmrg_dim_change_mclr(LRras2(1:8),ntash,0)
         call dmrg_spc_change_mclr(LRras2(1:8),nash)
       end if

c
c  Write the effective active one el density to disk in the same format as ipg1q
c
c       Call Getmem('TEMP1','ALLO','REAL',ipDeff_act,ndens2)
c       call dcopy_(nDens2,Work(ipD_K),1,Work(ipDeff_act),1)
c       Do is=1,nSym
c        Do i=1,nish(is)
c
c Subtract the inactive density
c
c         Work(ipDeff_act+ipmat(is,is)-1+i-1+(i-1)*nbas(is))=
c     &   Work(ipD_K+ipmat(is,is)-1+i-1+(i-1)*nbas(is))-Two
c        End Do
c       End Do
c
c      Call Put_DEff(Work(ipDeff_act),ndens2)
c
c       Call Getmem('TEMP1','FREE','REAL',ipDeff_act,ndens2)
c
c--------------------------------------------------
c
c Diagonalize the effective density to be able to use Prpt
c ipO eigenvalues of eff dens
c ipCMON eigenvectors (new orb coef)
c
c      Call NatOrb(Work(ipD_K),Work(ipCMO),Work(ipCMON),Work(ipO))
c      Call Getmem('TMP', 'ALLO','Real',ipT,nBuf/2)
c      Call dmat_MCLR(Work(ipCMON),Work(ipO),Work(ipT))
c      Call Put_D1ao_Var(Work(ipT),nTot1)
c      Note='var'
c      LuTmp=50
c      LuTmp=IsFreeUnit(LuTmp)
c      Call WrVec('TMPORB',LuTmp,'O',nSym,nBas,nBas,
c    &            Dum,Work(ipO),Dum,iDum,Note)
c      Call Prpt()

c
c Standard routine, ipT effective dens in AO
c
*       Call dmat_MCLR(Work(ipCMON),Work(ipO),Work(ipT))
c
*       Call Put_D1ao_Var(Work(ipT),nTot1)
c      Call Getmem('TMP', 'FREE','Real',ipT,nBuf/2)

       Call Put_iScalar('SA ready',1)
       If (isNAC) Then
         Write(mstate,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
       Else
         Write(mstate,'(I16)') irlxroot
       End If
       If (override) mstate(1:1)='+'
       Call Put_cArray('MCLR Root',mstate,16)
*
       Call GetMem('kappa1','FREE','Real',ipK1,  nDens2)
       Call Getmem('kappa2','FREE','Real',ipk2,  nDens2)
         Call GetMem('Dens','FREE','Real',ipD_K,Length)
       Call GetMem('ONEDEN','FREE','Real',ipDAO,nDens2)
       Call GetMem('ONEDEN','FREE','Real',ipD_CI,n1Dens)
       Call GetMem('ONEDEN','FREE','Real',ipD1,n1Dens)
       Call Getmem('TWODEN', 'FREE','Real',ipP_CI,n2Dens)
       Call Getmem('TWODEN', 'FREE','Real',ipP1,n2Dens)
       Call Getmem('OCCU ', 'FREE','Real',ipO,nbas_tot)
       Call Getmem('CMO', 'FREE','Real',ipCMON,ndens2)
       Call Getmem('conn', 'FREE','Real',ipF,ndens2)
       Call GetMem('TMPDEN','FREE','REAL',ipDtmp,ndens2)

       irc=ipclose(-1)
*
       Call QExit('Out_PT2')

       Return
       End

c --------------------------------------------------------------------------
c
      Subroutine OITD(rK,isym,D,Dtmp,act)
*
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "real.fh"
#include "sa.fh"
      Real*8 rK(*),D(*),Dtmp(*)
      Logical act
      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)
*
      Call QEnter('OITD')
      call dcopy_(ndens2,[Zero],0,Dtmp,1)
*
*     Note: even with NAC we set the inactive block,
*     because this is the SA density, not the transition density
      Do iS=1,nSym
        Do iB=1,nIsh(iS)
          Dtmp(1+(ipCM(iS)+(ib-1)*nOrb(iS)+ib-1)-1) = Two
        End Do
      End Do
      If (act) Then
       Do iS=1,nSym
        Do iB=1,nAsh(iS)
         Do jB=1,nAsh(iS)
          Dtmp(1+(ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nOrb(is)-1)-1)=
     &    Work(ipG1t+(itri((nA(is)+ib),(nA(is)+jb)))-1)
         End Do
        End Do
       End Do
      End If
*
      Do iS=1,nsym
         jS=ieor(iS-1,isym-1)+1
         If (nOrb(iS)*nOrb(jS).ge.1) Then
            Call DGEMM_('N','T',nOrb(iS),nOrb(jS),nOrb(iS),One,
     &                 Dtmp(1+ipCM(iS)-1),nOrb(iS),
     &                 rK(ipMat(jS,iS)),nOrb(jS),
     &                 Zero,D(ipMat(iS,jS)),nOrb(iS))
            Call DGEMM_('T','N',nOrb(iS),nOrb(jS),nOrb(jS),-One,
     &                 rK(ipMat(jS,iS)),nOrb(jS),
     &                 Dtmp(1+ipCM(jS)-1),nOrb(jS),
     &                 One,D(ipMat(iS,jS)),nOrb(iS))
         End If
      End Do
      Call QExit('OITD')
      Return
      End
      Subroutine NatOrb(Dens,CMOO,CMON,OCCN)
      Implicit Real*8(a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
#include "real.fh"
      Real*8 Dens(*),CMOO(*),CMON(*),OCCN(*)
      Call GetMem('TMP','ALLO','REAL',ips,ndens2)
      Call GetMem('TMP','ALLO','REAL',ips2,ndens2)
C
C         Diagonalize the density matrix and transform orbitals
C
      If (iAnd(kprint,8).eq.8) Then
         Write(6,*)
         Write(6,*) '           Effective natural population '
         Write(6,*) '           ============================ '
         Write(6,*)
      End If
      io=0
      Do is=1,nsym
         ij=0
         Do i=0,nbas(is)-1
            Do j=0,i
               Work(ipS2+ij)=Dens(ipMat(is,is)+i+j*nbas(is))
               ij=ij+1
            End DO
         End DO
         Call dCopy_(nBas(iS)**2,[Zero],0,Work(ipS),1)
         Call dCopy_(nBas(is),[One],0,Work(ipS),nbas(is)+1)
         CALL JACOB(Work(ipS2),Work(ipS),nbas(is),nbas(is))
         ii=0
         DO i=1,nbas(is)
            ii=ii+i
            OCCN(io+i)=Work(ips2-1+ii)
         END DO
         IST=IO+1
         IEND=IO+NBAS(is)
         If (iAnd(kprint,2).eq.2)
     &      Write (6,'(6X,A3,I2,A1,10F11.6,/,(12X,10F11.6))')
     &             'sym',iS,':',(OCCN(I),I=IST,IEND)
         If (nBas(is).ge.1)
     &      CALL DGEMM_('N','N',
     &                  NBAS(is),NBAS(is),NBAS(is),
     &                  One,CMOO(ipCM(is)),NBAS(is),
     &                  Work(ipS),NBAS(is),
     &                  Zero,CMON(ipCM(is)),NBAS(is))
         io=io+nbas(is)
      End DO
      Call GetMem('TMP','FREE','REAL',ips,ndens2)
      Call GetMem('TMP','FREE','REAL',ips2,ndens2)
      Return
      End
