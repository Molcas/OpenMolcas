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
      Subroutine test_caspt2_setup(ntimes)
      Implicit Real*8 (a-h,o-z)
      Integer nSm,nIs(8),nAs(8),nSs(8)
      Integer NumC(8)
      External frandom_molcas

* Just getting test runs of the cho_caspt2_setup routine:
      write(6,*)
      write(6,*)' Exercise runs of cho_caspt2_setup with random inputs:'
      write(6,*)' Number of runs:',ntimes
      iseed=24691357
      do itimes=1,ntimes
       write(6,*)
       write(6,*)
       write(6,*)' Exercise run nr',iTimes
       nSm=2**(nint(4.0d0*frandom_molcas(iseed)-0.5d0))
       do iSm=1,nSm
        nIs(iSm)=nint(8.0d0*frandom_molcas(iseed)-0.5d0)
        x=16.0D0/nSm
        nAs(iSm)=nint(x*frandom_molcas(iseed)-0.5d0)
        nSs(iSm)=nint(15.0d0*frandom_molcas(iseed)-0.5d0)
        NumC(iSm)=nint(15.0d0*frandom_molcas(iseed)-0.5d0)
       end do
       call cho_caspt2_setup(nSm,nIs,nAs,nSs,NumC,'Allo')
       call cho_caspt2_setup(nSm,nIs,nAs,nSs,NumC,'Free')
      end do
      write(6,*)
      write(6,*)' End of exercises for today.'
      write(6,*)
      return
      end

      Subroutine cho_caspt2_setup(nSym,nIsh,nAsh,nSsh,NumCho,mode)
* -------------------------
* This subroutine uses the input variables to compute
* ipsp,nisplit, nasplit,lsplit,ipnp,ipip each dimensioned (1:nsym),
* in common /chocaspt2/,
* and allocates fields at iwork() addresses ipsp(1:nsym),ipnp(1:nsym),
* ipunit_f(1:nsym), each with sizes lsplit(1:nsym), and
* ipip(1:nsym), each with sizes nsym*lsplit(1:nsym)
* -------------------------
      Implicit Real*8 (a-h,o-z)
      Integer nSym,nIsh(8),nAsh(8),nSsh(8)
      Integer NumCho(8)
      Integer iAorb(8),iIorb(8),iKorb(8)
      Character*(*) mode
      Character*4   modecopy

#include "WrkSpc.fh"
#include "chocaspt2.fh"
      Integer  cho_irange
      External cho_irange

C *********************************************************************
      nIc(i,j) = iWork(ipSp(j)+i-1)
******
      nAc(i,j) = iWork(ipSp(j)+nIsplit(j)+i-1)
******
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C *********************************************************************

      modecopy=mode(1:4)
      Call UpCase(modecopy)
      if (modecopy.eq.'FREE') then
         Do jSym=1,nSym
        If (NumCho(jSym).gt.0) then
            Call GetMem('Unit_F','Free','Inte',ipUnit_F(jSym),
     &                                         Lsplit(jSym))
            Call GetMem('iPorb','Free','Inte',ipiP(jSym),
     &                                        nSym*Lsplit(jSym))
            Call GetMem('nPorb','Free','Inte',ipnP(jSym),Lsplit(jSym))
            Call GetMem('Split','Free','Inte',ipSp(jSym),Lsplit(jSym))
        End If
         End Do
       Return
      end if

*      write(6,*)' Cho_CasPT2_SetUp entered.'
      do i=1,8
       ipunit_f(i)=0
       ipsp(i)=0
       lsplit(i)=0
       nisplit(i)=0
       nasplit(i)=0
       ipnp(i)=0
       ipip(i)=0
      end do

* Local arrays:
* iIorb(iSym) = Nr of inactive orbitals in earlier symmetries.
* iAorb(iSym) = Nr of   active orbitals in earlier symmetries.
* iKorb(iSym) = iIorb(iSym) + iAorb(iSym)
         iIorb(1)=0
         iAorb(1)=0
         iKorb(1)=0
         Do iSym=2,nSym
            iIorb(iSym) = iIorb(iSym-1) + nIsh(iSym-1)
            iAorb(iSym) = iAorb(iSym-1) + nAsh(iSym-1)
            iKorb(iSym) = iIorb(iSym) + iAorb(iSym)
         End Do

* nIO   = total number of inactive orbitals.
* nAO   = total number of active orbitals.
* nOkrb = total number of inactive+active orbitals.
* nOrb  = total number of orbitals
         nOkrb=0
         nIO=0
         nAO=0
         nO=0
         Do iSym=1,nSym
            nOkrb = nOkrb + nIsh(iSym) + nAsh(iSym)
            nOrb = nIsh(iSym) + nAsh(iSym) + nSsh(iSym)
            nO = nO + nOrb*(nOrb+1)/2
            nIO = nIO + nIsh(iSym)
            nAO = nAO + nAsh(iSym)
         End Do

         xO = dble(nO)

* Task: set up the necessary administration of transformed Cholesky
* vectors. Each such vector is characterized by the composite symmetry
* label JSYM. Vectors with common JSYM are stored on the same unit.
* For a given JSYM, there are NumCho(jSym) vectors. These are separated
* into subsets -- batches -- There will be nISplit(jSym) batches with
* vectors with a fixed inactive orbital index (Inactive vectors), and
* nASplit(jSym) batches with vectors with a fixed active index
* (Active vectors). Each batch has vectors with fixed index within a
* range of size iWork(ipsp(jSym)+i-1) where i is in 1..nISplit(jSym) or
*  in nISplit(jSym)+1..nISplit(jSym)+nASplit(jSym).
* The vectors themselves have a size of iWork(ipnp(jSym)+i-1)
* The vectors are generally referred to as e.g. p#k, where #k stands
* for the fixed inactive or active orbital, while p is all the orbitals
* with a symmetry such that the compound symmetry label is JSYM.
* There will be a (possibly large) number of such vectors, each indexed
* with 'J' gennerally. This is the summation index in the definition of
* the cholesky vectors.
         Do jSym=1,nSym
*      write(6,*)' NumCho(jSym)=',NumCho(jSym)
            If (NumCho(jSym).lt.1) goto 99


         Call GetMem('MaxMem','Max','Real',Kdummy,MemMx)
            xMemMx = dble(MemMx)
* xMemMx = largest allocatable field.

  5      Continue
         jfrac=1

 10         Continue
         Mem1 = 2*NumCho(jSym)/jfrac  ! hold 2 vectors in memory
            If (Mem1.eq.0 .and. jfrac.gt.1) then
            Call Cho_x_Quit('cho_caspt2_setup',
     &                          ': Sorry! Too little memory!!',' ')
            EndIf

         kfrac=0
* Go here to try using larger kfrac.
 20      Continue
         kfrac=kfrac+1
* Subdivide the set of orbitals into kfrac pieces of size .le. nKsp
         nKsp = nOkrb/kfrac
* If kfrac has grown too large, start again using kfrac=1 but larger jfrac
            If (nKsp .eq. 0) Then
                jfrac = jfrac + 1
                goto 10

            EndIf

* Loop over pieces of the orbital set:
         nPmax=0
         Do iK1=1,nOkrb,nKsp
          iK2=max(iK1-1+nKsp,nOkrb)
* For this piece, compute the number of (p,k) orbital pairs
* with p active or secondary, and with joint symmetry jSym:
          nP=0
          Do iK=iK1,iK2
             kS = cho_irange(iK,iKorb,nSym,.false.)
             jS = MulD2h(kS,jSym)
             nP = nP + nAsh(jS) + nSsh(jS)
          End Do
* nPmax will be the largest number of such pairs in any piece.
          nPmax = Max(nPmax,nP)
         End Do

* PAM:The following looks strange, since the max is over the same set
* of numbers no matter what jSym is.
            mRHS=0
            Do iS=1,nSym
               jS = MulD2h(iS,jSym)
               mRHS = Max( mRHS, Max(nAsh(jS),nSsh(jS)) )
            End Do

C --- Conversion to real*8 to avoid integer overflow on 32-bit machines

* PAM:Why would this be 'mem for right-hand side?'
            xRHS = dble(mRHS**2)           ! mem. for right hand side
            xLpk = dble(Mem1*nPmax*nKsp)   ! store Cholesky MO vectors
            xPIQK= dble((nPmax*nKsp)**2)   ! store integrals
            xmNeed= xO + xPIQK + Max(xLpk,2*xRHS) ! Fmat+integrals+rhs

         If (xmNeed .gt. xMemMx) GoTo 20

* This also looks strange -- nIAc=all the inact+act orbitals no matter what.?
                nIAc = nOkrb

            jIAc = nIAc
            LSplit(jSym)=1
            Do while (jIAc .lt. nOkrb)
               LSplit(jSym) = LSplit(jSym) + 1
               jIAc = jIAc + nIAc

            End Do

* What does this mean? Simply that lsplit(jSym) will be 'the next'
* of something?
            LSplit(jSym) = LSplit(jSym) + 1 ! separate out Ina from Act

* Need xmNeedNow units of memory
            xmNeedNow = xmNeed + dble((3+nSym)*LSplit(jSym))

            If (xmNeedNow.gt.xMemMx) Then
* Does it help to start all over again??
                xMemMx = xMemMx - dble((3+nSym)*LSplit(jSym))
                goto 5

            EndIf

*         write(6,*)' Determined splitting. lsplit(jsym)=', lsplit(jsym)
* Allocate arrays, in all (3+nSym)*lsplit(jSym) elements:
            Call GetMem('Split','Allo','Inte',ipSp(jSym),Lsplit(jSym))
            Call GetMem('nPorb','Allo','Inte',ipnP(jSym),Lsplit(jSym))
            Call GetMem('iPorb','Allo','Inte',ipiP(jSym),
     &                                        nSym*Lsplit(jSym))
            Call GetMem('Unit_F','Allo','Inte',ipUnit_F(jSym),
     &                                         Lsplit(jSym))

         nmin=Min(nIO,nIAc)
            jS=0
         if(nmin.gt.0) jS=nIO/nmin
         Do i=0,jS-1
           iWork(ipsp(jSym)+i) = nmin
            End Do
         mDiff = jS*nmin - nIO
            If (mDiff .gt. 0) Then
               iWork(ipSp(jSym)+jS) = mDiff
               jS = jS + 1
            Endif
            nIsplit(jSym) = jS

         nmin=Min(nAO,nIAc)
            jS=0
         if(nmin.gt.0) jS=nAO/nmin
         Do i=0,jS-1
           iWork(ipsp(jSym)+nisplit(jSym)+i) = nmin
            End Do
         mDiff = jS*nmin - nAO
            If (mDiff .gt. 0) Then
               iWork(ipSp(jSym)+nIsplit(jSym)+jS) = mDiff
               jS = jS + 1
            Endif
            nAsplit(jSym) = jS

* Note: The inactive orbitals will be partitioned.
* The isp=1..nisplit(jSym) is counter of the partitions.
* ioff will be the number of inactive orbitals in earlier partitions.
            ioff=0
*      write(6,*)' Partitions for inactive fixed indices:'
            Do isp=1,nIsplit(jSym)
*      write(6,*)' Partition nr isp=',isp
* ioff+1 is the first orbital of partition isp.
               lS = cho_irange(ioff+1,iIorb,nSym,.false.)
*      write(6,*)' First inactive orbital:',ioff+1
*      write(6,*)'           Its symmetry:',lS
* lS is its symmetry.
               iK = nIc(isp,jSym) + ioff
* Note: nIc(isp,jSym)=iWork(ipsp(jSym)+isp-1)
* iK is the last orbital of the partition, and kS its symmetry
               kS = cho_irange(iK,iIorb,nSym,.false.)
*      write(6,*)' Last  inactive orbital:',iK
*      write(6,*)'           Its symmetry:',kS

* Compute offsets into array at iWork(ipiP(jSym)..), formally dimensioned
*  as (1:nSym,1:nsp) where nsp is nisplit(jsym)+nasplit(jsym)
* Initialize with offsets=0, which is appropriate for non-used symmetries:
               nPorb=0
               Do jS=1,nSym
                 iWork(ipiP(jSym)+nSym*(isp-1)+jS-1)=nPorb
               End Do
* and then loop over the actual values.
* Loop over the symmetry range of inactive orbitals in the partition
               Do iS=lS,kS
*      write(6,*)'    If #k has symmetry:',iS
* jS is then the symmetry of its companion orbital in the pair.
                  jS=MulD2h(iS,jSym)
*      write(6,*)'   then p has symmetry:',jS
*      write(6,*)' Number of p orbitals :',nAsh(jS) + nSsh(jS)
*      write(6,*)' Offset to them       :',nPOrb
* ipip(jSym) is pointer to an array dimensioned iP(nSym,nisplit(jSym))
* which is used for offsets. In some other array, after the position
* iP(jS,isp) follows space for nAsh(jS)+nSsh(jS) items.
                  iWork(ipiP(jSym)+nSym*(isp-1)+jS-1) = nPorb
* Update nPorb.
                  nPorb = nPorb + nAsh(jS) + nSsh(jS)
               End Do

* ipnp(jSym) is pointer to an array dimensioned nP(nisplit(jSym))
* It gives the total size for the items mentioned above.
               iWork(ipnP(jSym)+isp-1) = nPorb
               ioff = ioff + nIc(isp,jSym)

            End Do

* Here follows similar arrays for the active orbitals.
* Addressing is the same, except we use iP(jS,nisplit(jSym)+isp)
* and nP(nisplit(jSym)+isp)
            ioff=0
            Do isp=1,nAsplit(jSym)
               lS = cho_irange(ioff+1,iAorb,nSym,.false.)

               iw = nAc(isp,jSym) + ioff
               kS = cho_irange(iw,iAorb,nSym,.false.)

* Initialize with offsets=0, which is appropriate for non-used symmetries:
               nPorb=0
               Do jS=1,nSym
                 iWork(ipiP(jSym)+nSym*(nIsplit(jSym)+isp-1)+jS-1)=nPorb
               End Do
* and then loop over the actual values.
               Do iS=lS,kS
                 jS=MulD2h(iS,jSym)
                 iWork(ipiP(jSym)+nSym*(nIsplit(jSym)+isp-1)+jS-1)=nPorb
                 nPorb = nPorb + nAsh(jS) + nSsh(jS)
               End Do

               iWork(ipnP(jSym)+nIsplit(jSym)+isp-1) = nPorb
               ioff = ioff + nAc(isp,jSym)

            End Do



99          Continue

         End Do

      write(6,*)
      write(6,*)' cho_caspt2_setup report:'
      write(6,*)
      write(6,'(1x,a,8i4)')' Inactive :',(nIsh(isym),isym=1,nSym)
      write(6,'(1x,a,8i4)')' Active   :',(nAsh(isym),isym=1,nSym)
      write(6,'(1x,a,8i4)')' Secondary:',(nSsh(isym),isym=1,nSym)
      write(6,*)
      write(6,'(1x,a,8i4)')' NumCho   :',(NumCho(isym),isym=1,nSym)
      write(6,*)
      do jsym=1,nsym
       write(6,*)' Symm:',jSym
       write(6,*)' Partition  Fixed orbitals       nPorb space'
       kend=0
       do isp=1,nisplit(jsym)
        ksta=kend+1
        kend=kend+iWork(ipsp(jSym)+isp-1)
        kstasym=cho_irange(ksta,iIorb,nSym,.false.)
        kendsym=cho_irange(kend,iIorb,nSym,.false.)
        nPorb=iWork(ipnp(jSym)+isp-1)
        write(6,'(1x,i4,5x,i4,a4,i4,2x,i1,a4,i1,5x,i4)')
     &        isp,ksta,' -- ',kend,kstasym,' -- ',kendsym,nPorb
        write(6,'(1x,a,8i8)')' iP Offsets:',(iWork(ipip(jSym)+
     &                 nSym*(isp-1)+MulD2h(jSym,iS)-1),iS=1,nSym)
       write(6,*)
       end do
       kend=0
       do isp=1,nasplit(jsym)
        ksta=kend+1
        kend=kend+iWork(ipsp(jSym)+nisplit(jSym)+isp-1)
        kstasym=cho_irange(ksta,iAorb,nSym,.false.)
        kendsym=cho_irange(kend,iAorb,nSym,.false.)
        nPorb=iWork(ipnp(jSym)+nisplit(jSym)+isp-1)
        write(6,'(1x,i4,5x,i4,a4,i4,2x,i1,a4,i1,5x,i4)')
     &        isp,ksta,' -- ',kend,kstasym,' -- ',kendsym,nPorb
        write(6,'(1x,a,8i8)')' iP Offsets:',
     &         (iWork(ipip(jSym)+nSym*(nisplit(jSym)+isp-1)+
     &                             MulD2h(jSym,iS)-1),iS=1,nSym)
       write(6,*)
       end do
      End Do


      Return
      END
