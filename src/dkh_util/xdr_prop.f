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
*Copyright (C) 2021, Rulin Feng                                        *
************************************************************************
C
C----------------------------------------------------------------------|
C
      subroutine XDR_Prop(nbas,isize,jsize,imethod,paratyp,dkhorder,
     &                    xorder,inS,inK,inV,inpVp,inX,inpXp,
     &                    inUL,inUS,clight,Label,iComp,iSizec)
C
C Driver for relativistic transformation of property integrals
C
C   also called the "picture change correction" since the physical operator
C   X is defined in four-component picture, one need to transform them
C   to two/one- component picture as well as the Hamiltonian in the two/one-
C   component relativistic calculations
C
      implicit none
#include "WrkSpc.fh"
C Input variables
      integer nbas,isize,jsize,imethod,paratyp,dkhorder,xorder,iComp
      Real*8 clight
      Real*8 inS(isize),inK(isize),inpVp(isize),inpXp(isize)
      Real*8 inV(isize),inUL(jsize),inUS(jsize)
C Input/Output variables
      Real*8 inX(isize)
C Local variables
      integer nn,i,j,k
      integer jS,jK,jV,jpVp,jX,jpXp,itmp
      character*8 Label
      integer iSizec
      integer idbg
      integer nSym ! equals one
      integer iOpt, iRC, Lu_One, lOper, nInt
      integer imagaPX, imagaPXs, imagaXP, imagaXPs
      integer imagbPX, imagbPXs, imagbXP, imagbXPs
      integer iPSO, iPSOt
      integer jComp, iPSOComp
      integer ip_Ppso
      integer, dimension(1) :: IDUM
      character*8 magLabel
      character*8 PSOLabel
C
C Convert to square matrices
C
      nn = nbas*nbas+4
      call getmem('skin ','ALLOC','REAL',jK   ,nn)
      call getmem('sSS  ','ALLOC','REAL',jS   ,nn)
      call getmem('sV   ','ALLOC','REAL',jV   ,nn)
      call getmem('spVp ','ALLOC','REAL',jpVp ,nn)
      call getmem('sX   ','ALLOC','REAL',jX   ,nn)
      call getmem('spXp ','ALLOC','REAL',jpXp ,nn)
#ifdef MOLPRO
      call square(inK  ,Work(jK),nbas,nbas)
      call square(inS  ,Work(jS),nbas,nbas)
      call square(inV  ,Work(jV),nbas,nbas)
      call square(inpVp,Work(jpVp),nbas,nbas)
      call square(inX  ,Work(jX),nbas,nbas)
      call square(inpXp,Work(jpXp),nbas,nbas)
#else
      call square(inK  ,Work(jK),nbas,1,nbas)
      call square(inS  ,Work(jS),nbas,1,nbas)
      call square(inV  ,Work(jV),nbas,1,nbas)
      call square(inpVp,Work(jpVp),nbas,1,nbas)
      call square(inX  ,Work(jX),nbas,1,nbas)
      call square(inpXp,Work(jpXp),nbas,1,nbas)
#endif
C
C Calculate the relativistic transformed property integrals
C
      if(imethod.eq.2.or.imethod.eq.3.or.
     &   (imethod.eq.1.and.xorder.ge.15) )then
C
C Handle magnetic integrals if provided by Gen1Int
C
*********************************************************************
c
c          PSO integrals with i, pre-existing PSO integrals in ONEREL
c          required. The h_UL{\dag} should have a minus sign, for
c          which a tranpose is not enough, thus add this manually.
c *         MAG:
c *          1: x_k*d/dx  4: x_k*d/dy  7: x_k*d/dz
c *          2: y_k*d/dx  5: y_k*d/dy  8: y_k*d/dz
c *          3: z_k*d/dx  6: z_k*d/dy  9: z_k*d/dz
c *         PSO1: operator (y_k*d/dz - z_k*d/dy)/r_k^3
c *          = MAG 8 - MAG 6
c *         PSO2: operator (z_k*d/dx - x_k*d/dz)/r_k^3
c *          = MAG 3 - MAG 7
c *         PSO3: operator (x_k*d/dy - y_k*d/dx)/r_k^3
c *          = MAG 4 - MAG 2
c *         PSO n = MAG a - MAG b

        If ((Label(1:3).eq.'MAG').and.
     &           ((iComp.eq.3).or.(iComp.eq.4).or.
     &            (iComp.eq.8))) then
            nSym = 1 !always C1 symmetry
            inUS = inUS * clight
            iOpt = 0
            iRC = -1
            Lu_One = 2
            call OpnOne(iRC,iOpt,'ONEREL',Lu_One)
            If (iRC.ne.0) Go To 9999
            !switch basis to primitive functions
            call OneBas('PRIM')
            ! do MAG a
            write(magLabel,'(A,A3)') 'MAGXP',Label(6:8)
            iOpt = 1
            iRC = -1
            lOper = -1
            call iRdOne(iRC,iOpt,magLabel,iComp,idum,lOper)
            nInt = IDUM(1)
            If (iRC.ne.0) Go To 9999
            call getmem('MAGaXP','ALLOC','REAL',imagaXP,nInt+4)
            iOpt = 0
            iRC = -1
            call RdOne(iRC,iOpt,magLabel,iComp,Work(imagaXP),lOper)
            If (iRC.ne.0) Go To 9999
            !call CmpInt(Work(imagaXP),nInt,nbas,nSym,lOper)
            call getmem('MAGaPX','ALLOC','REAL',imagaPX,nInt+4)
            magLabel(1:5) = 'MAGPX'
            call RdOne(iRC,iOpt,magLabel,iComp,Work(imagaPX),lOper)
            !If (iRC.ne.0) Go To 9999
            call getmem('MAGaPXs','ALLOC','REAL',imagaPXs,nn)
            call getmem('MAGaXPs','ALLOC','REAL',imagaXPs,nn)
            call square(Work(imagaPX),Work(imagaPXs),nbas,1,nbas)
            call square(Work(imagaXP),Work(imagaXPs),nbas,1,nbas)
            call getmem('MAGaXP','FREE','REAL',imagaXP,nInt+4)
            call getmem('MAGaPX','FREE','REAL',imagaPX,nInt+4)
            call merge_mag_ints(nbas, jsize, Work(imagaXPs),
     &                          Work(imagaPXs),.true.)
            ! put Work(imagaPXs) into -Work(imagaPXs)
            ! or put Work(imagaXPs) into -Work(imagaXPs)
            ! test which
            do k = 0, nn-5
                Work(imagaPXs+k) = -Work(imagaPXs+k)
            end do
            call getmem('TMP ','ALLOC','REAL',itmp,nn)
            ! rulin- disable or reenable dmxma for mag integral X2C
            ! transformations
            call dmxma(nbas,'N','N',Work(imagaPXs),inUS,Work(itmp),1.d0)
            call dmxma(nbas,'T','N',inUL,Work(itmp),Work(imagaPXs),1.d0)
            call dmxma(nbas,'N','N',Work(imagaXPs),inUL,Work(itmp),1.d0)
            call dmxma(nbas,'T','N',inUS,Work(itmp),Work(imagaXPs),1.d0)
            do k = 0, nn-5
                Work(imagaXPs+k) = Work(imagaXPs+k) + Work(imagaPXs+k)
            end do
            call getmem('MAGaPXs','FREE','REAL',imagaPXs,nn)
            ! define MAG b
            if (iComp.eq.3) then
               jComp = 7
               iPSOComp = 2
            else if (iComp.eq.4) then
               jComp = 2
               iPSOComp = 3
            else if (iComp.eq.8) then
               jComp = 6
               iPSOComp = 1
            end if
            ! do MAG b
            write(magLabel,'(A,A3)') 'MAGXP',Label(6:8)
            iOpt = 1
            iRC = -1
            lOper = -1
            call iRdOne(iRC,iOpt,magLabel,jComp,idum,lOper)
            nInt = IDUM(1)
            If (iRC.ne.0) Go To 9999
            call getmem('MAGbXP','ALLOC','REAL',imagbXP,nInt+4)
            iOpt = 0
            iRC = -1
            call RdOne(iRC,iOpt,magLabel,jComp,Work(imagbXP),lOper)
            call getmem('MAGbPX','ALLOC','REAL',imagbPX,nInt+4)
            magLabel(1:5) = 'MAGPX'
            call RdOne(iRC,iOpt,magLabel,jComp,Work(imagbPX),lOper)
            call ClsOne(iRC,iOpt)
            call getmem('MAGbPXs','ALLOC','REAL',imagbPXs,nn)
            call getmem('MAGbXPs','ALLOC','REAL',imagbXPs,nn)
            call square(Work(imagbPX),Work(imagbPXs),nbas,1,nbas)
            call square(Work(imagbXP),Work(imagbXPs),nbas,1,nbas)
            call getmem('MAGbXP','FREE','REAL',imagbXP,nInt+4)
            call getmem('MAGbPX','FREE','REAL',imagbPX,nInt+4)
            call merge_mag_ints(nbas, jsize, Work(imagbXPs),
     &                          Work(imagbPXs),.true.)
            ! put Work(imagbPXs) into -Work(imagbPXs)
            ! or put Work(imagbXP) into -Work(imagbXP)
            ! test which
            do k = 0, nn-5
                Work(imagbPXs+k) = -Work(imagbPXs+k)
            end do
            call dmxma(nbas,'N','N',Work(imagbPXs),inUS,Work(itmp),1.d0)
            call dmxma(nbas,'T','N',inUL,Work(itmp),Work(imagbPXs),1.d0)
            call dmxma(nbas,'N','N',Work(imagbXPs),inUL,Work(itmp),1.d0)
            call dmxma(nbas,'T','N',inUS,Work(itmp),Work(imagbXPs),1.d0)
            do k = 0, nn-5
                Work(imagbXPs+k) = Work(imagbXPs+k) + Work(imagbPXs+k)
            end do
            call getmem('MAGbPXs','FREE','REAL',imagbPXs,nn)
            call getmem('TMP ','FREE','REAL',itmp,nn)
            ! do PSO combinations
            call getmem('PSO ','ALLOC','REAL',iPSO,nn)
            do k = 0, nn-5
                Work(iPSO+k) = Work(imagaXPs+k) - Work(imagbXPs+k)
            end do
            call getmem('MAGaXPs','FREE','REAL',imagaXPs,nn)
            call getmem('MAGbXPs','FREE','REAL',imagbXPs,nn)
            write(PSOLabel,'(A,A3)') 'PSOI ',Label(6:8)
            !test for off-diagonal elements
            IDUM(1) = nbas
            call CmpInt(Work(iPSO),nInt,idum(1),nSym,lOper)
            ! store PSO integrals to ONEINT
            call getmem('PSOt','ALLOC','REAL',iPSOt,nInt+4)
            k = 0
            do i = 1, nbas
                do j = 1, i
                    k = k + 1
                    Work(iPSOt+k-1) = Work(iPSO+j-1+(i-1)*nbas)
                end do
            end do
            call getmem('PSO ','FREE','REAL',iPSO,nn)
            call Allocate_Work(ip_Ppso,iSizec+4)
            idbg = -1
            call repmat(idbg,Work(iPSOt),Work(ip_Ppso),.FALSE.)
            call getmem('PSOt','FREE','REAL',iPSOt,nInt+4)
            iOpt = 0
            iRC = -1
            Lu_one = 2
            Call OpnOne(iRC,iOpt,'ONEINT',Lu_one)
            If (iRC.ne.0) Go to 9999
            iRC = -1
            lOper = 255
            call WrOne(iRC,iOpt,PSOLabel,iPSOComp,Work(ip_Ppso),lOper)
            If (iRC.ne.0) Go to 9999
            iOpt = 0
            Call ClsOne(iRC,iOpt)
            call Free_Work(ip_Ppso)
            inUS = inUS / clight
        End if
*********************************************************************
        if (Label(1:3).eq.'MAG') then
          inUS = inUS * clight
          ! Put together lower and upper triangular matrices
          call merge_mag_ints(nbas, jsize, Work(jX), Work(jpXp),.true.)
          call getmem('TMP ','ALLOC','REAL',itmp ,nn)
          ! Eval U_L^{\dag} X U_S
          ! Eval U_S^{\dag} X U_L
          call dmxma(nbas,'N','N',Work(jpXp),inUS,Work(itmp),1.d0)
          call dmxma(nbas,'T','N',inUL,Work(itmp),Work(jpXp),1.d0)
          call dmxma(nbas,'N','N',Work(jX),inUL,Work(itmp),1.d0)
          call dmxma(nbas,'T','N',inUS,Work(itmp),Work(jX),1.d0)
          ! Sum
          do k=0,nbas*nbas-1
            Work(jX+k) = Work(jX+k) + Work(jpXp+k)
          end do
          ! copy jpXp back so its not a half computed matrix
          do k=0,nbas*nbas-1
            Work(jpXp+k) = Work(jX+k)
          end do
          call getmem('TMP ','FREE','REAL',itmp ,nn)
          inUS = inUS / clight
        else
C
C X2C/BSS transformation
C
C   because the transformation matrix in non-orthogonal basis picture has
C   obtained via the Hamiltonian drivers, here we just need to simply apply
C   it to property integrals ( X, pXp in four-component picture )
C
C   high order DKH can also employ this formulation, only negligible
C   contribution from higher orders is included
C
          call getmem('TMP ','ALLOC','REAL',itmp ,nn)
C         ! eval U_L^{\dag} X U_L
          call dmxma(nbas,'C','N',inUL,Work(jX),Work(itmp),1.d0)
          call dmxma(nbas,'N','N',Work(itmp),inUL,Work(jX),1.d0)
C         ! eval U_S^{\dag}pXp U_S
          call dmxma(nbas,'C','N',inUS,Work(jpXp),Work(itmp),1.d0)
          call dmxma(nbas,'N','N',Work(itmp),inUS,Work(jpXp),1.d0)
C         ! sum
          do k=0,nbas*nbas-1
            Work(jX+k) = Work(jX+k) + Work(jpXp+k)
          end do
          call getmem('TMP ','FREE','REAL',itmp ,nn)
        End if !MAG
C
      else if(imethod.eq.1)then
C
C Arbitrary order DKH transformation
C
        call dkh_prop(nbas,Work(jS),Work(jK),Work(jV),Work(jpVp),
     &              Work(jX),Work(jpXp),clight,dkhorder,xorder,paratyp)
      endif
C
C Copy transformed property integral back to inX
C
      k=0
      do i=1,nbas
        do j=1,i
          k=k+1
          inX(k)=Work(jX+j-1+(i-1)*nbas)
        enddo
      enddo
C
C Free temp memories
C
      call getmem('skin ','FREE','REAL',jK   ,nn)
      call getmem('sSS  ','FREE','REAL',jS   ,nn)
      call getmem('sV   ','FREE','REAL',jV   ,nn)
      call getmem('spVp ','FREE','REAL',jpVp ,nn)
      call getmem('sX   ','FREE','REAL',jX   ,nn)
      call getmem('spXp ','FREE','REAL',jpXp ,nn)
C
      return
 9999 Continue
      Write (6,*) ' *** Error in subroutine XDR_Prop ***'
      Write (6,*) '     Abend in subroutine OpnOne'
      Call Abend
      end
