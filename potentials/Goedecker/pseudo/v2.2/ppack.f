      subroutine ppack (rloc,gpot,hsep,r_l,pp,
     :     lpx,lpmx,nspin,nsmx,maxdim,nfit,nstring,
     :     avgl1,avgl2,avgl3,ortprj,litprj)
c     pack the fittingparameters into one array for amoeba
c     and vice versa
c     if nstring = 'init'   -> initialisation + packing
c     if nstring = 'pack'   -> packing
c     if nstring = 'unpack' -> unpacking

      implicit none

      integer lpx,lpmx,nspin,nsmx,maxdim,nfit,increm,lpxfit
      real*8 rloc,gpot(4),hsep(6,lpmx,nsmx),r_l(lpmx)
      real*8 pp(maxdim)
      real*8 h11,h12,h13,h22,h23,h33,hh11,hh12,hh13,hh22,hh23,hh33
      character*(*) nstring
      
      integer maxpar
      parameter (maxpar = 33)
      character*10 spack(maxpar)
      logical lpack(maxpar),ortprj,litprj
      integer i,j,ncount,nhsep1,nhsep2,nhsep3,ll,nn,nnmin
      logical filepr,avgl1,avgl2,avgl3

      integer llpmx,nnsmx
      parameter(llpmx=4, nnsmx=2)

      real*8 pih,xshift
      real*8 orloc,ogpot(4),orl(0:llpmx-1),ohsep(6,llpmx,nnsmx)
      save orloc,ogpot,orl,ohsep,lpack

      data spack /'rloc','gpot(1)','gpot(2)','gpot(3)','gpot(4)',
     :     'r_l(1)','hsep(1,1)','hsep(1,2)','hsep(1,3)',
     :               'hsep(1,4)','hsep(1,5)','hsep(1,6)',
     :     'r_l(2)','hsep(2,1)','hsep(2,2)','hsep(2,3)',
     :               'hsep(2,4)','hsep(2,5)','hsep(2,6)',
     :     'r_l(3)','hsep(3,1)','hsep(3,2)','hsep(3,3)',
     :               'hsep(3,4)','hsep(3,5)','hsep(3,6)',
     :     'r_l(4)','hsep(4,1)','hsep(4,2)','hsep(4,3)',
     :               'hsep(4,4)','hsep(4,5)','hsep(4,6)'/

c
c------------------------------------------------------------------
c     initialisation: set up the list of parameters to be (un)packed
c     lpack() = 1 -> (un)pack, lpack() = 0 -> keep aktual value
c     if a file "FITPAR" exist in the current directory lpack() will be
c     read from there

c      print*,'entered ppack, nstring=',nstring
c      print*,'nfit=',nfit
      if ( llpmx.lt.lpmx .or. nnsmx.lt.nsmx ) then
         write(6,*) 'array dimension problem in ppack'
         stop
      endif
      pih=2.d0*atan(1.d0)
      if (nstring.eq.'init') then
         do i=1,maxdim
            lpack(i)=.false.
            pp(i) = 0.0d0
         enddo
         INQUIRE ( FILE = 'FITPAR', EXIST = FILEPR )
         IF ( FILEPR ) THEN
            open(20,file='FITPAR')
            read(20,*)
            read(20,*) (lpack(i),i=1,5)
            write(6,*) (lpack(i),i=1,5)
            read(20,*) lpxfit
            if ( lpxfit .ne. lpx ) then
               write(6,*) 'wrong number of projectors in FITPAR'
               stop
            endif

            do i=0, lpx
               read(20,*)  (lpack(i*7+j),j=6,12)
               write(6,*)  (lpack(i*7+j),j=6,12)
            enddo
            close(20)
            write(6,*) 'fitting parameter determined by ''FITPAR'' '
         endif
c
c     if the projectors have to be orthogonalized do not count
c     h(l+1,2),h(l+1,4),h(l+1,5)  as fitting parameter 
c
         if (ortprj.or.litprj) then
            do i=0,lpx
               if (lpack( (i*7+6) +2) ) then
                  write(6,*) ' l=',i
                  write(6,*) 'warning! transformation of hsep(): ',
     :                 'no fit of h(1,2)'
                  lpack((i*7+6) +2)=.false.
               endif
               if (lpack( (i*7+6) +4) ) then
                  write(6,*) ' l=',i
                  write(6,*) 'warning! transformation of hsep(): ',
     :                 'no fit of h(1,3)'
                  lpack((i*7+6) +4)=.false.
               endif
               if (lpack( (i*7+6) +5) ) then
                  write(6,*) ' l=',i
                  write(6,*) 'warning! transformation of hsep(): ',
     :                 'no fit of h(2,3)'
                  lpack((i*7+6) +5)=.false.
               endif
            enddo
         endif
         nfit = 0
         write(6,*) 'fitting parameter:'
         do i=1,maxpar
            if (lpack(i)) then 
               write(6,*) spack(i)
               nfit = nfit +1
            endif
         enddo
c
c     double count projectors for nspin=2 and l>0
c
c         print*,'nfit before double count:',nfit
         if (nspin.eq.2) then
c     l=1   hsep(1,1)  at position 14 
            do i=14,19
               if (lpack(i)) nfit =nfit + 1
            enddo
c     l=2   hsep(1,1)  at position 21
            do i=21,26
               if (lpack(i)) nfit =nfit + 1
            enddo
c     l=3   hsep(1,1)  at position 28
            do i=28,33
               if (lpack(i)) nfit =nfit + 1
            enddo
         endif
c         print*,'nfit after double count:',nfit
         if (nfit.gt.maxdim) stop 'nfit > maxdim'
         write(6,*) 'maxdim,nfit:',maxdim,nfit
         if (maxdim.lt.nfit) stop
c
c     save initial parameter 
c
         orloc = rloc
         do i=1,4
            ogpot(i)=gpot(i)
         enddo
         do ll=0,lpx
            orl(ll)=r_l(ll+1)
            do i=1,min(2*ll+1,nspin)
               do j=1,6
                  ohsep(j,ll+1,i)= hsep(j,ll+1,i)
               enddo
c     if the projectors are othonormalized we transform ohsep to the orthognomal basis
c     and save the diagonal terms of the new ohsep() matrix
               if (ortprj) then
                  h11=ohsep(1,ll+1,i)
                  h12=ohsep(2,ll+1,i)
                  h22=ohsep(3,ll+1,i)
                  h13=ohsep(4,ll+1,i)
                  h23=ohsep(5,ll+1,i)
                  h33=ohsep(6,ll+1,i)
                  if (ll.eq.0) then

      HH11=h11 + 1.549193338482967d0*h12 + 0.975900072948533d0*h13 
     :     + 0.6d0*h22 + 0.7559289460184545d0*h23 
     :     + 0.2380952380952381d0*h33
      HH12=0.6324555320336759d0*h12 + 0.7968190728895957d0*h13 + 
     :     0.4898979485566356d0*h22 + 0.925820099772551d0*h23 
     :     + 0.3888078956798695d0*h33
      HH13=0.3563483225498992d0*h13 + 0.2760262237369417d0*h23 + 
     :     0.173880176985767d0*h33
      HH22=0.4d0*h22+1.007905261357939d0*h23 + 0.6349206349206349d0*h33
      HH23=0.02839451399999733d0*(7.937253933193773d0*h23 + 10.d0*h33)
      HH33=0.126984126984127d0*h33

                  elseif (ll.eq.1) then

      HH11=h11 + 1.690308509457033d0*h12 + 1.189176780021126d0*h13
     :     + 0.7142857142857144d0*h22 + 1.005037815259212d0*h23 
     :     + 0.3535353535353535d0*h33
      HH12=0.5345224838248488d0*h12 + 0.7521014330903549d0*h13 
     :     + 0.4517539514526256d0*h22 + 0.953462589245592d0*h23 
     :     + 0.4471907802258314d0*h33
      HH13=0.2842676218074805d0*h13 + 0.240249990052149d0*h23 
     :     + 0.1690222275826415d0*h33
      HH22=0.2857142857142857d0*h22 + 0.80403025220737d0*h23 + 
     :     0.5656565656565657d0*h33
      HH23=0.01527129183875666d0*(9.9498743710662d0*h23 + 14.d0*h33)
      HH33=0.0808080808080808d0*h33

                 elseif (ll.eq.2) then


      HH11=h11 + 1.763834207376394d0*h12 + 1.327493036606129d0*h13 + 
     :     0.7777777777777778d0*h22 + 1.170738814009927d0*h23 
     :     + 0.4405594405594406d0*h33
      HH12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + 
     :     0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + 
     :     0.4709778528806361d0*h33
      HH13=0.236524958395633d0*h13 + 0.2085954062582479d0*h23 + 
     :     0.1569926176268787d0*h33
      HH12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + 
     :     0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + 
     :     0.4709778528806361d0*h33
      HH22=0.2222222222222222d0*h22 + 0.6689936080056727d0*h23 + 
     :     0.5034965034965035d0*h33
      HH23=0.00932400932400932d0*(11.9582607431014d0*h23 + 18.d0*h33)
      HH33=0.05594405594405595d0*h33

                 elseif (ll.eq.3) then

      HH11=h11 + 1.809068067466582d0*h12 + 1.425050606388851d0*h13 + 
     :     0.818181818181818d0*h22 + 1.289006773270979d0*h23 + 
     :     0.5076923076923077d0*h33
      HH12=0.0006593070220853591d0*(646.741834119303d0*h12 + 
     :     1018.911183568028d0*h13 + 585.d0*h22 + 
     :     1382.459764333125d0*h23 + 726.d0*h33)
      HH13=0.2025478734167333d0*h13 + 0.1832114449657378d0*h23 + 
     :     0.144320484917644d0*h33
      HH22=0.1818181818181818d0*h22 + 0.5728918992315464d0*h23 + 
     :     0.4512820512820513d0*h33
      HH23=0.006184848093902844d0*(13.96424004376894d0*h23+22.d0*h33)
      HH33=0.04102564102564103d0*h33

                 endif
                 ohsep(1,ll+1,i)=HH11
                 ohsep(2,ll+1,i)=HH12
                 ohsep(3,ll+1,i)=HH22
                 ohsep(4,ll+1,i)=HH13
                 ohsep(5,ll+1,i)=HH23
                 ohsep(6,ll+1,i)=HH33
              endif
            enddo
         enddo
c         write(6,*)'saved parameter:'
c         write(6,'(6f10.5)') orloc, ogpot
c         do ll=0,lpx
c            write(6,'(f10.5)') orl(ll)
c            do i=1,min(2*ll+1,nspin)
c               write(6,'(6f10.5)') ( ohsep(j,ll+1,i),  j=1,6)
c            enddo
c         enddo
      endif
c
c------------------------------------------------------------------
c
c     pack the parameters into the array pp()
c
      if (nstring.eq.'pack'.or.nstring.eq.'init') then
         print*,'nfit=',nfit
         do i=1,nfit
            pp(i)=0.0d0
         enddo
c         write(6,*) 'packed array pp:'
c         write(6,'(20e8.2)') (pp(i),i=1,nfit)
      endif
c
c------------------------------------------------------------------
c
c     unpack array pp()
c
      if (nstring.eq.'unpack') then
c        write(6,*) 'unpacking array pp:'
c        write(6,'(20e8.2)') (pp(i),i=1,nfit)
         ncount = 1
c     rloc
         if (lpack(1)) then 
            rloc = orloc+.10d0*orloc*atan(pp(ncount))/pih 
            ncount = ncount + 1
         endif  
c     gpot(1-4)
         do i=1,4
            if (lpack(1+i)) then 
               gpot(i) = ogpot(i)+pp(ncount)
               ncount = ncount + 1
            endif  
         enddo
c     projectors for l=0: r_l(1),hsep(1,1-6)
         increm=6
         if ( lpack(increm) ) then 
            r_l(1) = orl(0)+.10d0*orl(0)*atan(pp(ncount))/pih 
            ncount = ncount + 1
         endif  
         do i=1,6
            if (lpack(increm+i)) then 
               hsep(i,1,1)= ohsep(i,1,1)+pp(ncount)
               ncount = ncount + 1
            else
               hsep(i,1,1)= ohsep(i,1,1)
            endif  
         enddo
         if (ortprj) then
c     do back transformation from orthonormal projectors to 
c     unnormalized projectors:
            h11=hsep(1,1,1)
            h22=hsep(3,1,1)
            h33=hsep(6,1,1)

            HH11=H11 + 1.5d0*H22 + 1.875d0*H33
            HH12=-1.936491673103709d0*H22 - 4.841229182759272d0*H33
            HH13=3.842606537234849d0*H33
            HH22=2.5d0*H22 + 12.5d0*H33
            HH23=-9.92156741649221d0*H33
            HH33=7.875d0*H33
            
            hsep(1,1,1)=HH11
            hsep(2,1,1)=HH12
            hsep(3,1,1)=HH22
            hsep(4,1,1)=HH13
            hsep(5,1,1)=HH23
            hsep(6,1,1)=HH33
         else if (litprj) then
            hsep(2,1,1)= -0.5d0 * sqrt(3.d0/5.d0)  *hsep(3,1,1)
            hsep(4,1,1)= +0.5d0 * sqrt(5.d0/21.d0) *hsep(6,1,1)
            hsep(5,1,1)= -0.5d0 * sqrt(100.d0/63.0d0) *hsep(6,1,1)
         endif
c     projectors for l=1: r_l(1),hsep(1,1-6,1-2)
         increm=13
         if (lpack(increm)) then 
            r_l(2) = orl(1)+.10d0*orl(1)*atan(pp(ncount))/pih 
            ncount = ncount + 1
         endif  
         do i=1,6
            if (lpack(increm+i)) then 
               hsep(i,2,1)= ohsep(i,2,1) + pp(ncount)
               if (nspin.eq.2) hsep(i,2,2)= ohsep(i,2,2)+pp(ncount+1)
               ncount = ncount + nspin
            else
               hsep(i,2,1)= ohsep(i,2,1)
               if (nspin.eq.2) hsep(i,2,2)= ohsep(i,2,2)
            endif  
         enddo
         if (ortprj) then
            do i=1,nspin
c     
c     do back transformation from orthonormal projectors to 
c     unnormalized projectors:
               h11=hsep(1,2,i)
               h22=hsep(3,2,i)
               h33=hsep(6,2,i)

               HH11=H11 + 2.5d0*H22 + 4.375d0*H33
               HH12=-2.958039891549808d0*H22 - 10.35313962042433d0*H33
               HH13=7.358031326380719d0*H33
               HH22=3.5d0*H22 + 24.5d0*H33
               HH23=-17.41228014936585d0*H33
               HH33=12.375d0*H33

               hsep(1,2,i)=HH11
               hsep(2,2,i)=HH12
               hsep(3,2,i)=HH22
               hsep(4,2,i)=HH13
               hsep(5,2,i)=HH23
               hsep(6,2,i)=HH33
                  
            enddo
         endif
         if (litprj) then
            do i=1,nspin
               hsep(2,2,i)= -0.5d0*sqrt(5.d0/7.d0)       *hsep(3,2,i)
               hsep(4,2,i)=  0.5d0*sqrt(35.d0/11.d0)/3.d0*hsep(6,2,i)
               hsep(5,2,i)= -0.5d0*14.d0/sqrt(11.d0)/3.d0*hsep(6,2,i)
            enddo
         endif
c     projectors for l=2: r_l(1),hsep(1,1-6,1-2)
         increm=20
         if (lpack(increm)) then 
            r_l(3) = orl(2)+.10d0*orl(2)*atan(pp(ncount))/pih 
            ncount = ncount + 1
         endif  
         do i=1,6
            if (lpack(increm+i)) then 
               hsep(i,3,1)= ohsep(i,3,1)+pp(ncount)
               if (nspin.eq.2) hsep(i,3,2)= ohsep(i,3,2)+pp(ncount+1)
               ncount = ncount + nspin
            else
               hsep(i,3,1)= ohsep(i,3,1)
               if (nspin.eq.2) hsep(i,3,2)= ohsep(i,3,2)
            endif  
         enddo
         if (ortprj) then
            do i=1,nspin
c     
c     do back transformation from orthonormal projectors to 
c     unnormalized projectors:
               h11=hsep(1,3,i)
               h22=hsep(3,3,i)
               h33=hsep(6,3,i)


               HH11=H11 + 3.5d0*H22 + 7.875d0*H33
               HH12=-3.968626966596886d0*H22-17.85882134968598d0*H33
               HH13=11.86446901466728d0*H33
               HH22=4.5d0*H22 + 40.5d0*H33
               HH23=-26.90608667197814d0*H33
               HH33=17.875d0*H33

               hsep(1,3,i)=HH11
               hsep(2,3,i)=HH12
               hsep(3,3,i)=HH22
               hsep(4,3,i)=HH13
               hsep(5,3,i)=HH23
               hsep(6,3,i)=HH33
                  
            enddo
         endif
         if (litprj) then
            do i=1,nspin
               hsep(2,3,i)= -0.5d0*sqrt(7.d0/9.d0)   *hsep(3,3,i)
               hsep(4,3,i)=  0.5d0*3.d0*sqrt(7.0d0/143.d0)*hsep(6,3,i)
               hsep(5,3,i)= -0.5d0*18.d0*sqrt(1/143.0d0)*hsep(6,3,i)
            enddo
         endif
c     projectors for l=3: r_l(1),hsep(1,1-6,1-2)
         increm=27
         if (lpack(increm)) then 
            r_l(4) = orl(3)+.10d0*orl(3)*atan(pp(ncount))/pih 
            ncount = ncount + 1
         endif  
         do i=1,6
            if (lpack(increm+i)) then 
               hsep(i,4,1)= ohsep(i,4,1)+pp(ncount)
               if (nspin.eq.2) hsep(i,4,2)= ohsep(i,4,2)+pp(ncount+1)
               ncount = ncount + nspin
            else
               hsep(i,4,1)= ohsep(i,4,1)
               if (nspin.eq.2) hsep(i,4,2)= ohsep(i,4,2)
            endif  
         enddo
         if (ortprj) then
            do i=1,nspin
c     
c     do back transformation from orthonormal projectors to 
c     unnormalized projectors:
               h11=hsep(1,4,i)
               h22=hsep(3,4,i)
               h33=hsep(6,4,i)

               HH11=H11 + 4.5d0*H22 + 12.375d0*H33
               HH12=-4.9749371855331d0*H22 - 27.36215452043205d0*H33
               HH13=17.36780426536412d0*H33
               HH22=5.5d0*H22 + 60.5d0*H33
               HH23=-38.40166012036458d0*H33
               HH33=24.375d0*H33

               hsep(2,4,i)=HH12
               hsep(3,4,i)=HH22
               hsep(4,4,i)=HH13
               hsep(5,4,i)=HH23
               hsep(6,4,i)=HH33
                  
            enddo
         endif
         if (litprj) then
            do i=1,nspin
               hsep(2,4,i)= -0.5d0*sqrt(9.d0/11.d0)  *hsep(3,4,i)
               hsep(4,4,i)= +0.5d0*sqrt(33.d0/65.d0) *hsep(6,4,i)
               hsep(5,4,i)= -0.5d0*22.d0/sqrt(195.d0)*hsep(6,4,i)
            enddo
         endif
c
c     if avgl is set: modify hsep() so that average potential for the highest
c     projector is zero; if the projectors have to be orthogonalized modify
c     also the corresponding offdialgonal elements
c     (only for relativistic calculations)
         if (nspin.eq.2) then
            nhsep1=0
            nhsep2=0
            nhsep3=0
            do i=1,6
               if (hsep(i,2,1).ne.0) nhsep1=i
               if (hsep(i,3,1).ne.0) nhsep2=i
               if (hsep(i,4,1).ne.0) nhsep3=i
            enddo
            if (avgl1) then
               nnmin=1
               if (nhsep1.eq.3 .and. litprj) nnmin=2
               if (nhsep1.eq.6 .and. litprj) nnmin=4
               if (.not.litprj) nnmin=nhsep1
               if (ortprj) then
c     for orthogonal projectors remove average-part form the highest 
c     orthonormal projector 
                  do i=1,nspin
                     h11=hsep(1,2,i)
                     h12=hsep(2,2,i)
                     h22=hsep(3,2,i)
                     h13=hsep(4,2,i)
                     h23=hsep(5,2,i)
                     h33=hsep(6,2,i)
      HH11=h11 + 1.690308509457033d0*h12 + 1.189176780021126d0*h13
     :     + 0.7142857142857144d0*h22 + 1.005037815259212d0*h23 
     :     + 0.3535353535353535d0*h33
      HH22=0.2857142857142857d0*h22 + 0.80403025220737d0*h23 + 
     :     0.5656565656565657d0*h33
      HH33=0.0808080808080808d0*h33
                     hsep(1,2,i)=hh11
                     hsep(3,2,i)=hh22
                     hsep(6,2,i)=hh33
                  enddo
               endif
               do nn=nnmin,nhsep1
                  xshift= hsep(nn,2,2)+2.0d0*hsep(nn,2,1)
                  hsep(nn,2,1) =  hsep(nn,2,1) - xshift*1.d0/3.d0
                  hsep(nn,2,2) =  hsep(nn,2,2) - xshift*1.d0/3.d0
               enddo
               if (ortprj) then
                  do i=1,nspin
c     do back transformation of hsep()
                     h11=hsep(1,2,i)
                     h22=hsep(3,2,i)
                     h33=hsep(6,2,i)

               HH11=H11 + 2.5d0*H22 + 4.375d0*H33
               HH12=-2.958039891549808d0*H22 - 10.35313962042433d0*H33
               HH13=7.358031326380719d0*H33
               HH22=3.5d0*H22 + 24.5d0*H33
               HH23=-17.41228014936585d0*H33
               HH33=12.375d0*H33

                     hsep(1,2,i)=HH11
                     hsep(2,2,i)=HH12
                     hsep(3,2,i)=HH22
                     hsep(4,2,i)=HH13
                     hsep(5,2,i)=HH23
                     hsep(6,2,i)=HH33

                  enddo
               endif
            endif
            if (avgl2) then
               nnmin=1
               if (nhsep1.eq.3 .and. litprj) nnmin=2
               if (nhsep1.eq.6 .and. litprj) nnmin=4
               if (.not.litprj) nnmin=nhsep2
               if (ortprj) then
c     for orthogonal projectors remove average-part form the highest 
c     orthonormal projector 
                  do i=1,nspin
                     h11=hsep(1,3,i)
                     h12=hsep(2,3,i)
                     h22=hsep(3,3,i)
                     h13=hsep(4,3,i)
                     h23=hsep(5,3,i)
                     h33=hsep(6,3,i)
      HH11=h11 + 1.763834207376394d0*h12 + 1.327493036606129d0*h13 + 
     :     0.7777777777777778d0*h22 + 1.170738814009927d0*h23 
     :     + 0.4405594405594406d0*h33
      HH12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + 
     :     0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + 
     :     0.4709778528806361d0*h33
      HH13=0.236524958395633d0*h13 + 0.2085954062582479d0*h23 + 
     :     0.1569926176268787d0*h33
      HH12=0.4714045207910317d0*h12 + 0.7095748751868991d0*h13 + 
     :     0.4157397096415491d0*h22 + 0.938679328162116d0*h23 + 
     :     0.4709778528806361d0*h33
      HH22=0.2222222222222222d0*h22 + 0.6689936080056727d0*h23 + 
     :     0.5034965034965035d0*h33
      HH23=0.00932400932400932d0*(11.9582607431014d0*h23 + 18.d0*h33)
      HH33=0.05594405594405595d0*h33


                     hsep(1,3,i)=hh11
                     hsep(3,3,i)=hh22
                     hsep(6,3,i)=hh33
                  enddo
               endif

               do nn=nnmin,nhsep2
                  xshift= 2*hsep(nn,3,2)+3.0d0*hsep(nn,3,1)
                  hsep(nn,3,1) =  hsep(nn,3,1) - xshift*1.d0/5.d0
                  hsep(nn,3,2) =  hsep(nn,3,2) - xshift*1.d0/5.d0
               enddo
               if (ortprj) then
                  do i=1,nspin
c     do back transformation of hsep()
                     h11=hsep(1,3,i)
                     h22=hsep(3,3,i)
                     h33=hsep(6,3,i)
 
               HH11=H11 + 3.5d0*H22 + 7.875d0*H33
               HH12=-3.968626966596886d0*H22-17.85882134968598d0*H33
               HH13=11.86446901466728d0*H33
               HH22=4.5d0*H22 + 40.5d0*H33
               HH23=-26.90608667197814d0*H33
               HH33=17.875d0*H33

                     hsep(1,3,i)=HH11
                     hsep(2,3,i)=HH12
                     hsep(3,3,i)=HH22
                     hsep(4,3,i)=HH13
                     hsep(5,3,i)=HH23
                     hsep(6,3,i)=HH33
                  enddo
               endif
            endif
            if (avgl3) then
               nnmin=1
               if (nhsep1.eq.3 .and. litprj) nnmin=2
               if (nhsep1.eq.6 .and. litprj) nnmin=4
               if (.not.litprj) nnmin=nhsep3
               if (ortprj) then
c     for orthogonal projectors remove average-part form the highest 
c     orthonormal projector 
                  do i=1,nspin
                     h11=hsep(1,4,i)
                     h12=hsep(2,4,i)
                     h22=hsep(3,4,i)
                     h13=hsep(4,4,i)
                     h23=hsep(5,4,i)
                     h33=hsep(6,4,i)

      HH11=h11 + 1.809068067466582d0*h12 + 1.425050606388851d0*h13 + 
     :     0.818181818181818d0*h22 + 1.289006773270979d0*h23 + 
     :     0.5076923076923077d0*h33
      HH12=0.0006593070220853591d0*(646.741834119303d0*h12 + 
     :     1018.911183568028d0*h13 + 585.d0*h22 + 
     :     1382.459764333125d0*h23 + 726.d0*h33)
      HH13=0.2025478734167333d0*h13 + 0.1832114449657378d0*h23 + 
     :     0.144320484917644d0*h33
      HH22=0.1818181818181818d0*h22 + 0.5728918992315464d0*h23 + 
     :     0.4512820512820513d0*h33
      HH23=0.006184848093902844d0*(13.96424004376894d0*h23+22.d0*h33)
      HH33=0.04102564102564103d0*h33

                     hsep(1,4,i)=hh11
                     hsep(3,4,i)=hh22
                     hsep(6,4,i)=hh33
                  enddo
               endif
               do nn=nnmin,nn
                  xshift= 5*hsep(nn,4,2)+7.0d0*hsep(nn,4,1)
                  hsep(nn,4,1) =  hsep(nn,4,1) - xshift*1.d0/7.d0
                  hsep(nn,4,2) =  hsep(nn,4,2) - xshift*1.d0/7.d0
               enddo
               if (ortprj) then
                  do i=1,nspin
c     do back transformation of hsep()
                     h11=hsep(1,4,i)
                     h22=hsep(3,4,i)
                     h33=hsep(6,4,i)

               HH11=H11 + 4.5d0*H22 + 12.375d0*H33
               HH12=-4.9749371855331d0*H22 - 27.36215452043205d0*H33
               HH13=17.36780426536412d0*H33
               HH22=5.5d0*H22 + 60.5d0*H33
               HH23=-38.40166012036458d0*H33
               HH33=24.375d0*H33

                     hsep(1,4,i)=HH11
                     hsep(2,4,i)=HH12
                     hsep(3,4,i)=HH22
                     hsep(4,4,i)=HH13
                     hsep(5,4,i)=HH23
                     hsep(6,4,i)=HH33
                  enddo
               endif
            endif

         endif
      endif  
c------------------------------------------------------------------
c      print*,'leave ppack with nfit=',nfit
      return
      end












