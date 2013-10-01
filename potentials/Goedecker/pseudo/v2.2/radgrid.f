        subroutine radgrid(nrad,rr,rw,rd,a_grd,b_grd,rmax)
        implicit real*8 (a-h,o-z)
        dimension rr(nrad),rw(nrad),rd(nrad)
c generates logarithmic grid
c rr() radial grid: rr(i)=a_grd*b_grd**(i-1)-c_grd
c rw() weights for radial integration (dr/di)
c rd() di/dr

        fourpi=16.d0*atan(1.d0)
        do 100,i=1,nrad
        rr(i)=a_grd*exp(b_grd*(i-1))
        rw(i)=b_grd*rr(i)
        rd(i)=1.d0/rw(i)
        rw(i)=rw(i)*fourpi*rr(i)**2
        if (rr(i).gt.rmax) goto 200
100     continue
        write(6,*)'rmax too large, stopped in rradgrid'
        stop
 200    nrad=i-1
c modify weights at en point for improved accuracy
        rw(1)=rw(1)*17.d0/48.d0
        rw(2)=rw(2)*59.d0/48.d0
        rw(3)=rw(3)*43.d0/48.d0
        rw(4)=rw(4)*49.d0/48.d0

        return
        end
