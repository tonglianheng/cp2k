
	subroutine crtvh(ng,lcx,lmax,xp,vh,nint,rmt,rmtg,
     1    ud,rr)
	implicit real*8 (a-h,o-z)
	dimension xp(0:ng),
     1  vh(((ng+1)*(ng+2))/2,0:lcx,((ng+1)*(ng+2))/2,0:lmax),
     1  rmt(nint,((ng+1)*(ng+2))/2,lmax+1),
     1  rmtg(nint,((ng+1)*(ng+2))/2,lmax+1),
     1  ud(nint,((ng+1)*(ng+2))/2,lcx+1),
     1  rr(nint)
c	print*,'entered crtvh'
	if (lcx.gt.3) stop 'crtvh:lcx>3'
	if (lmax.gt.3) stop 'crtvh:lmax>3'

c wavefunction products on grid
	do 8049,l=0,lmax
	do 8049,k=1,nint
	 r=rr(k)  
	ij=0
c lower triangle
	do 8049,j=0,ng
	do 8049,i=j,ng
	ij=ij+1
	rmt(k,ij,l+1)=(r**2)**l*exp(-(xp(i)+xp(j))*r**2)
	a=xp(j)+xp(i)
	tt1=2.d0*(l-a*r**2)*(1.d0/r)
	rmtg(k,ij,l+1)=tt1*rmt(k,ij,l+1)
8049	continue 


c   Hartree potential on grid 
	do 9873,k=1,nint
	r=rr(k)
c lower triangle
	ij=0
        do 2487,j=0,ng
        do 2487,i=j,ng
	ij=ij+1
        spi=1.772453850905516d0
        d=xp(i)+xp(j)
        sd=sqrt(d)
        tx=exp(-d*r**2)
        tt=spi*Derf(sd*r)
	if (lcx.ge.0) 
     1     ud(k,ij,0+1)=tt/(4.d0*sd**3*r)
	if (lcx.ge.1) 
     &     ud(k,ij,0+2)=-tx/(4.d0*d**2) + 3.d0*tt/(8.d0*sd**5*r)
	if (lcx.ge.2) 
     &     ud(k,ij,0+3)=-tx*(7.d0 + 2.d0*d*r**2)/(8.d0*d**3) +
     &     15.d0*tt/(16.d0*sd**7*r)
        if (lcx.ge.3) 
     &     ud(k,ij,0+4)=-tx*(57.d0+22.d0*d*r**2+4.d0*d**2*r**4)/
     &     (16.d0*d**4) + 105.d0*tt/(32.d0*sd**9*r)
2487    continue
9873	continue

c Coulombic integrals
c lower triangle
	ij=0
	do 100,j=0,ng
	do 100,i=j,ng
	ij=ij+1
	c=xp(i)+xp(j)
c lower triangle
	ipjp=0
	do 100,jp=0,ng
	do 100,ip=jp,ng
	ipjp=ipjp+1
	d=xp(ip)+xp(jp)
	scpd=sqrt(c+d)
        vh(ipjp,0,ij,0)=0.2215567313631895d0/(c*d*scpd)
	if (lcx.eq.0) goto 11
        vh(ipjp,1,ij,0)=
     &  0.1107783656815948d0*(2.d0*c+3.d0*d)/(c*d**2*scpd**3)
	if (lcx.eq.1) goto 11
        vh(ipjp,2,ij,0)=
     &  0.05538918284079739d0*(8.d0*c**2+20.d0*c*d+15.d0*d**2)
     &     /(c*d**3*scpd**5)
	if (lcx.eq.2) goto 11
        vh(ipjp,3,ij,0)=
     &  0.0830837742611961d0*
     &  (16.d0*c**3+56.d0*c**2*d+70.d0*c*d**2+35.d0*d**3)/
     &   (c*d**4*scpd**7)
11	continue
	if (lmax.eq.0) goto 200

        vh(ipjp,0,ij,1)=
     &  0.1107783656815948d0*(3.d0*c+2.d0*d)/(c**2*d*scpd**3)
	if (lcx.eq.0) goto 22
        vh(ipjp,1,ij,1)=
     &  0.05538918284079739d0*(6.d0*c**2+15.d0*c*d+6.d0*d**2)/
     &   (c**2*d**2*scpd**5)
	if (lcx.eq.1) goto 22
        vh(ipjp,2,ij,1)=
     &  0.02769459142039869d0*
     &   (24.d0*c**3+84.d0*c**2*d+105.d0*c*d**2+30.d0*d**3)/
     &   (c**2*d**3*scpd**7)
	if (lcx.eq.2) goto 22
        vh(ipjp,3,ij,1)=
     &  0.04154188713059803d0*
     &     (48.d0*c**4+216.d0*c**3*d+378.d0*c**2*d**2+
     &      315.d0*c*d**3+70.d0*d**4)/(c**2*d**4*scpd**9)
22	continue
	if (lmax.eq.1) goto 200

        vh(ipjp,0,ij,2)=
     &  0.05538918284079739d0*(15.d0*c**2+20.d0*c*d+8.d0*d**2)/
     &   (c**3*d*scpd**5)
	if (lcx.eq.0) goto 33
        vh(ipjp,1,ij,2)=
     &  0.02769459142039869d0*
     &   (30.d0*c**3+105.d0*c**2*d+84.d0*c*d**2+24.d0*d**3)/
     &   (c**3*d**2*scpd**7)
	if (lcx.eq.1) goto 33
        vh(ipjp,2,ij,2)=
     &  0.2077094356529901d0*
     &   (8.d0*c**4+36.d0*c**3*d+63.d0*c**2*d**2+
     &      36.d0*c*d**3+8.d0*d**4)/(c**3*d**3*scpd**9)
	if (lcx.eq.2) goto 33
        vh(ipjp,3,ij,2)=
     &  0.1038547178264951d0*
     &      (48.d0*c**5+264.d0*c**4*d+594.d0*c**3*d**2+
     &      693.d0*c**2*d**3+308.d0*c*d**4+56.d0*d**5)/
     &   (c**3*d**4*scpd**11)
33	continue
	if (lmax.eq.2) goto 200

        vh(ipjp,0,ij,3)=
     &  0.0830837742611961d0*
     &   (35.d0*c**3+70.d0*c**2*d+56.d0*c*d**2+16.d0*d**3)/
     &   (c**4*d*scpd**7)
	if (lcx.eq.0) goto 44
        vh(ipjp,1,ij,3)=
     &  0.04154188713059803d0*
     &      (70.d0*c**4+315.d0*c**3*d+378.d0*c**2*d**2+
     &      216.d0*c*d**3+48.d0*d**4)/(c**4*d**2*scpd**9)
	if (lcx.eq.1) goto 44
        vh(ipjp,2,ij,3)=
     &  0.1038547178264951d0*
     &      (56.d0*c**5+308.d0*c**4*d+693.d0*c**3*d**2+
     &      594.d0*c**2*d**3+264.d0*c*d**4+48.d0*d**5)/
     &   (c**4*d**3*scpd**11)
	if (lcx.eq.2) goto 44
        vh(ipjp,3,ij,3)=
     &  1.090474537178198d0*(16.d0*c**6+104.d0*c**5*d+286.d0*c**4*d**2+
     &    429.d0*c**3*d**3+286.d0*c**2*d**4+104.d0*c*d**5+16.d0*d**6)/
     &   (c**4*d**4*scpd**13)
44	continue
200	continue

100	continue
	return
	end

