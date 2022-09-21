program molecular_dynamics
use omp_lib
integer, parameter :: n=4000
real ::rx(n),ry(n),rz(n),vx(n),vy(n),vz(n)
real ::potentialenergy,Kineticenergy,TotalEnergy,T,xbox,ybox,zbox,rcut,Tinst,Th,Tc
real :: rxprev(n), ryprev(n), rzprev(n)
real :: fx(n),fy(n),fz(n)
real :: Wx(n),Wy(n),Wz(n)
real :: eps(n), sigma(n)
integer ::step,i,j,max_bin,k,zbin,neqstep,NVEeq, eq_Avg,neqrun_start
real :: delr,dens,delz , zpos, xpos , ypos , dt , totcfg 
real, dimension (:), allocatable ::Tzprof,Tinsprof,Tzsteady
real, dimension (:), allocatable ::hist,bi,bieq
integer :: iseed
real :: parth , partc !counter of hot and cold particle
real :: partvh , partvc !velocity counter of hot and cold particle
real :: gamah,gamac
integer :: tconfig , thermostep 
integer :: state ! 0 - means Nonequilibrium followed by Equilibrium , 1 - nonequilibrium start from last state
real :: junk
real :: Pdesire, Pdamp, Pinst,Pneq
state = 0
thermostep = 100 !calculate temperature profile after this many steps
tconfig = 500 !skip this step and save coordinate
T = 1.2
Th=1.2 
Tc = 1.0        ![enter the value of temperature]
step= 4000000	![calculate time step]
neqstep = 150000

rcut= 2.5		![enter rcut]
dt = 0.001		!time_step
eq_Avg = 100000
neq_Avg = 3950000  ! average calculate after this many steps
neqrun_start = 0
totcfg = 100.0
!Pdesire = 1.2
!Pneq = 30.0
Pdamp = 1000*dt

call initialization (n,rx,ry,rz,vx,vy,vz,xbox,ybox,zbox)
!xbox= 10.0		![enter the box length]
!ybox= 10.0		![enter the box length]
!zbox = 57.1
write(*,*) xbox,ybox,zbox
zbin = nint(zbox) 
!if (zbin > zbox) then
!	zbin = zbin-1
!endif	
delz = real(zbox/zbin)  
write(*,*) zbin , delz
dens=n/(xbox*ybox*zbox)
max_bin=100
delr=xbox/2/(max_bin)

gamah = 100.0
gamac = 100.0


allocate (hist (max_bin),bi(max_bin),bieq(max_bin))  
allocate (Tzprof (zbin),Tinsprof(zbin),Tzsteady(zbin))  

!**********************initalizing bin ***********************

	do i =1,max_bin 
		hist(i) =0.0
		bi(i)=0.0
		bieq(i)=0.0
	end do

	



iseed = 5774232 !seed for random number generator
call ranint(iseed)



!call init(T,n,rx,ry,rz,vx,vy,vz,rxprev,ryprev,rzprev)




if (state .ne. 0) then
	open(109,file='final-eq-sigma-eps.dat')
	do k =1,n
	      read(109,*)eps(k),sigma(k)
	end do 
	close(109)
else
	do k =1,zbin 
		Tzprof(k) =0.0
		Tinsprof(k)=0.0
		Tzsteady(k) = 0.0
	end do	
endif

!mood loop Equilibration
if (state .eq. 0) then
	open (unit = 135,file = 'thermo-eq.dat')
	open (unit= 33, file = 'pair_correlation-eq.dat')
	open (unit = 115,file = 'final-eq.xyz')
	open (unit = 26,file = 'final-vel-eq.dat')
	open(19,file='Equilibrium-Coordinates1.xyz')
	open(39,file='Equilibrium-velocity.xyz')
	open(109,file='stress-eq.xyz')
	open(unit = 119 , file = 'final-eq-sigma-eps.dat')
	
	do j = 1, n
		eps(j) = 1.0 !initialization of epsilon for each particle
		sigma(j) = 1.0 !initialization of diameter for each particle
	enddo
	
	
	call forcecal(n,xbox,ybox,zbox,rx,ry,rz,rcut,eps,sigma,potentialenergy,fx,fy,fz,Wx,Wy,Wz)
	do i=1,neqstep
	
		call fix_langevin(i,n,xbox,ybox,zbox,dt,rx,ry,rz,rcut,eps,sigma,vx,vy,vz,fx,fy,fz,Wx,Wy,Wz,KineticEnergy,PotentialEnergy,gamah,Th)
		do j =1, n
			Wx(j) = -vx(j)*vx(j) - Wx(j)/2.0
			Wy(j) = -vy(j)*vy(j) - Wy(j)/2.0
			Wz(j) = -vz(j)*vz(j) - Wz(j)/2.0
		end do
		call fix_thermo(i,n,xbox,ybox,zbox,rx,ry,rz,Wx,Wy,Wz,Pinst,Pxx,Pyy,Pzz)
		!if (i > 50000) then
			!call fix_brensden_iso(i,n,xbox,ybox,zbox,dt,rx,ry,rz,Pdesire,Pdamp,Pinst,Pxx,Pyy,Pzz)
		!endif
		TotalEnergy = KineticEnergy + PotentialEnergy
		if (mod(i,10) .eq. 0) then
			write (135,*) i,(TotalEnergy/n),KineticEnergy/n,PotentialEnergy/n, 2*Kineticenergy/(3.0*n-3),Pinst,xbox,ybox,zbox
		endif
		
		if (i.gt.eq_Avg) then
			if (mod(i,tconfig) .eq. 0) then
				call pair_correlation (n,rx,ry,rz,xbox,ybox,zbox,delr,max_bin,hist,dens)
		 		do k=1,max_bin
					bieq(k) = bieq(k) + hist(k)
				enddo
			
			endif
		
		end if
		
		
		if (mod(i,tconfig) .eq. 0) then !virial stress tensor store
			write(109,*) n
                        write(109,*) 'Timestep' , i , 'box length' , xbox , ybox, zbox 
			do j =1, n
				write(109,*) 'A', Wx(j) , Wy(j) , Wz(j) 
			end do
		endif	
		
		if (mod(i,tconfig) .eq. 0) then
			write(19,*) n
                        write(19,*) 'Timestep' , i , 'box length' , xbox , ybox, zbox 
			write(39,*) n
                        write(39,*) 'Timestep' , i 
			do j = 1 , n	
				write(19,*) 'A', rx(j) , ry(j) , rz(j) 
				write(39,*) 'A',  vx(j) , vy(j), vz(j)
			end do
		endif
	
		if (i .eq. neqstep) then
			write(115,*) n
			write(115,*) 'Atom'
			do j = 1, n
				xpos = rx(j) - xbox*nint(rx(j)/xbox)
				ypos = ry(j) - ybox*nint(ry(j)/ybox)
				zpos = rz(j) - zbox*nint(rz(j)/zbox)
				if (abs(zpos) < zbox/4.0 ) then
					write(115,*) 'H', xpos , ypos , zpos  
				else
					write(115,*) 'O', xpos , ypos , zpos 
				endif
				write(26,*) vx(j) , vy(j) , vz(j) , j
				write(119,*) sigma(j) , eps(j) , j
			enddo
		endif
	
	end do

	do k=1,max_bin
		write (33,*) k*delr,bieq(k)/totcfg
	end do

	 close(33)
	 close(115)
	 close(19)
	 close(39)
	 close(135)
	 close(119)
	 close(26)
	 state = state + 1
endif	 

!mood loop Temperature gradient

if (state .eq. 1) then
	open (unit = 35,file = 'thermo-neq.dat')
	open (unit = 32, file = 'pair_correlation-neq.dat')
	open (unit = 11,file = 'Temp-prof-z-100.dat')
	open (unit = 101,file = 'Temp-prof-z-steady.dat')
	open (unit = 10,file = 'particle.dat')
	open (unit = 15,file = 'final-neq.xyz')
	open (unit = 16,file = 'final-vel.dat')
	open (unit = 9,file='Steadystate-Coordinates1.xyz')
	open (unit = 59,file='Steadystate-velocity.xyz')
	open (unit = 209,file='stress-neq.xyz')
	open (unit = 219 , file = 'final-neq-sigma-eps.dat')
	i =neqrun_start+0
	
	do k =1,zbin
	      Tzprof(k) =0.0
	      Tinsprof(k)=0.0
	      Tzsteady(k) = 0.0
	end do

	call forcecal(n,xbox,ybox,zbox,rx,ry,rz,rcut,eps,sigma,potentialenergy,fx,fy,fz,Wx,Wy,Wz)
	TotalEnergy = KineticEnergy + PotentialEnergy
	write (35,*) neqrun_start+0,(TotalEnergy/n),KineticEnergy/n,PotentialEnergy/n,2*Kineticenergy/(3.0*n-3)


	do i=neqrun_start+1,step

	
	
		call fix_langevin_rgn(i,n,xbox,ybox,zbox,dt,rx,ry,rz,rcut,eps,sigma,vx,vy,vz,fx,fy,fz,Wx,Wy,Wz,&
		KineticEnergy,PotentialEnergy,gamah,gamac,Th,Tc)
		
		 do j =1, n
                        Wx(j) = -vx(j)*vx(j) - Wx(j)/2.0
                        Wy(j) = -vy(j)*vy(j) - Wy(j)/2.0
                        Wz(j) = -vz(j)*vz(j) - Wz(j)/2.0
                end do
                call fix_thermo(i,n,xbox,ybox,zbox,rx,ry,rz,Wx,Wy,Wz,Pinst,Pxx,Pyy,Pzz)
                !if (i > 1) then
                !	call fix_brensden_iso(i,n,xbox,ybox,zbox,dt,rx,ry,rz,Pneq,Pdamp,Pinst,Pxx,Pyy,Pzz)
                !endif
                TotalEnergy = KineticEnergy + PotentialEnergy
                if (mod(i,10) .eq. 0) then
                  write (35,*) i,(TotalEnergy/n),KineticEnergy/n,PotentialEnergy/n, 2*Kineticenergy/(3.0*n-3),Pinst,xbox,ybox,zbox
                endif
		
		if (mod(i,10) .eq. 0) then
			call modify_Eps2(i,n,xbox,ybox,zbox,rx,ry,rz,eps,sigma,delz,zbin,Tzprof,Th,Tc)
		endif
		if (i > 0) then
		if (mod(i,thermostep) .eq. 0) then
			write(11,*) 'timestep :' , i
			do k=1,zbin
				Tzprof(k)=Tzprof(k)/(thermostep)
				write (11,*) (k-1)*delz-zbox/2.0+delz/2.0,Tzprof(k)
			end do
			!call modify_Eps2(i,n,xbox,ybox,zbox,rx,ry,rz,eps,sigma,delz,zbin,Tzprof,Th,Tc)
		
			if (i.gt.neq_Avg) then
				do k=1,zbin
					Tzsteady(k)=Tzsteady(k) + Tzprof(k) 
				end do
			endif	
			
			
			do k =1,zbin 
				Tzprof(k) =0.0
			end do
		else if (mod(i,1) .eq. 0) then
			call temp_prof_z(i,n,xbox,ybox,zbox,rx,ry,rz,vx,vy,vz,delz,zbin,Tinsprof)
			do k=1,zbin
				Tzprof(k) = Tzprof(k) + Tinsprof(k)
			enddo
		endif
		endif
		if (i.gt.neq_Avg) then
			if (mod(i,tconfig) .eq. 0) then
				call pair_correlation (n,rx,ry,rz,xbox,ybox,zbox,delr,max_bin,hist,dens)
		 		do k=1,max_bin
					bi(k) = bi(k) + hist(k)
				enddo
			endif
		
		end if

		
	
		if (mod(i,tconfig) .eq. 0) then !virial stress tensor store
			write(209,*) n
                        write(209,*) 'Timestep' , i , 'box length' , xbox , ybox, zbox 
			do j =1, n
				write(209,*) 'A', Wx(j) , Wy(j) , Wz(j) 
			end do
		endif	

		if (mod(i,tconfig) .eq. 0) then
			write(9,*) n
                        write(9,*) 'Timestep' , i , 'box length' , xbox , ybox, zbox
			write(59,*) n
                        write(59,*) 'Timestep' , i 
			do j = 1 , n	
				write(9,*) 'A',rx(j) , ry(j) , rz(j) 
				write(59,*) 'A', vx(j) , vy(j), vz(j)
			enddo
		endif
		parth = 0.0
		partc = 0.0
		partvh = 0.0
		partvc = 0.0
		do j = 1 , n !check particular in particular region in the system
			zpos = rz(j) - zbox*anint(rz(j)/zbox)
			if (abs(zpos) < zbox/4.0) then
				parth = parth + 1.0
				partvh = partvh + sqrt(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))
			else
				partc = partc + 1.0
				partvc = partvc + sqrt(vx(j)*vx(j)+vy(j)*vy(j)+vz(j)*vz(j))
			endif
		enddo
		write(10,*) i , parth , partc , partvh , partvc
	
	
		if (i .eq. step) then
			write(15,*) n
			write(15,*) 'Atom'
			do j = 1, n
				xpos = rx(j) - xbox*nint(rx(j)/xbox)
				ypos = ry(j) - ybox*nint(ry(j)/ybox)
				zpos = rz(j) - zbox*nint(rz(j)/zbox)
				if (abs(zpos) < zbox/4.0 ) then
					write(15,*) 'H', xpos , ypos , zpos  
				else
					write(15,*) 'O', xpos , ypos , zpos 
				endif
				write(16,*) vx(j) , vy(j) , vz(j) , j
				write(219,*) sigma(j) , eps(j) , j
			enddo
		
		endif
	 
	end do
	

	do k=1,max_bin
	write (32,*) k*delr,bi(k)/(2*totcfg)
	end do

	do k=1,zbin-1
		write (101,*) (k-1)*delz-zbox/2.0+delz/2.0,Tzsteady(k)/(10*totcfg)
	end do
	 close (35)
	 close(32)
	 close(11)
	 close(101)
	 close(10)
	 close(15)
	 close(16)
	 close(9)
	 close(59)
	 close(219)
endif
end program

subroutine fix_thermo(im,natom,xboxd,yboxd,zboxd,rxd,ryd,rzd,Wxd,Wyd,Wzd,Pinst,Pxx,Pyy,Pzz)
integer :: im,natom
real :: xboxd,yboxd,zboxd, vol
real :: rxd(natom),ryd(natom),rzd(natom)
real :: Wxd(natom),Wyd(natom),Wzd(natom)
real :: Pxx, Pyy, Pzz, Pinst
integer :: i ,j 
vol = xboxd*yboxd*zboxd
Pxx = 0.0
Pyy = 0.0
Pzz = 0.0
do i=1,natom
	Pxx = Pxx + Wxd(i)
	Pyy = Pyy + Wyd(i)
	Pzz = Pzz + Wzd(i)
end do

Pxx = -Pxx/vol
Pyy = -Pyy/vol
Pzz = -Pzz/vol	

Pinst = (pxx+pyy+pzz)/3.0


return

end

subroutine fix_brensden_iso(im,natom,xboxd,yboxd,zboxd,dt,rxd,ryd,rzd,Pdesired,Pdamp,Pinst,Pxx,Pyy,Pzz)
integer :: im,natom
real :: xboxd,yboxd,dt,zboxd
real :: rxd(natom),ryd(natom),rzd(natom)
real :: Wxd(natom),Wyd(natom),Wzd(natom)
real :: Pxx, Pyy, Pzz, eta, compress , Pdesired, Pdamp, Pinst
integer :: i ,j
vol = xboxd*yboxd*zboxd
compress = 1.0


eta = (1 - dt*compress*(Pdesired-Pinst)/Pdamp)

xboxd = eta*xboxd
yboxd = eta*yboxd
zboxd = eta*zboxd

do i = 1 ,natom
        rxd(i) = eta*rxd(i)
        ryd(i) = eta*ryd(i)
        rzd(i) = eta*rzd(i)
end do

return

end


subroutine fix_brensden(im,natom,xboxd,yboxd,zboxd,dt,rxd,ryd,rzd,Pdesire,Pdamp,Pinst,Pxx,Pyy,Pzz)
integer :: im,natom
real :: xboxd,yboxd,dt,zboxd 
real :: rxd(natom),ryd(natom),rzd(natom)
real :: Wxd(natom),Wyd(natom),Wzd(natom)
real :: Pxx, Pyy, Pzz, etax, etay, etaz , compress , Pdesirev, Pdamp, Pinst
integer :: i ,j 
vol = xboxd*yboxd*zboxd
compress = 1.0


etax = (1 - dt*compress*(Pdesire-Pxx)/Pdamp)
etay = (1 - dt*compress*(Pdesire-Pyy)/Pdamp)
etaz = (1 - dt*compress*(Pdesire-Pzz)/Pdamp)  

xboxd = etax*xboxd
yboxd = etay*yboxd
zboxd = etaz*zboxd

do i = 1 ,natom
	rxd(i) = etax*rxd(i)
	ryd(i) = etay*ryd(i)
	rzd(i) = etaz*rzd(i)
end do

return

end


subroutine fix_langevin(im,natom,xboxd,yboxd,zboxd,dt,rxd,ryd,rzd,rcut,eps,sigma,vxd,vyd,vzd,fxd,fyd,fzd,Wxd,Wyd,Wzd,KE,PE,gamah,Th)
integer :: im,natom
real :: xboxd,yboxd,rcut,dt,zboxd
real :: rxd(natom),ryd(natom),rzd(natom)
real :: vxd(natom),vyd(natom),vzd(natom)
real :: fxd(natom),fyd(natom),fzd(natom)
real :: Wxd(natom),Wyd(natom),Wzd(natom)
real :: eps(natom),sigma(natom)
integer :: i ,j 
real KE,PE
real :: Th,Tc,gamah,gamac,gfrih,gfric,noiseh,noisec,random
real :: m !mass
REAL, PARAMETER :: c1 = 2.0, c2 = -2.0, c3 = 4.0/3.0, c4 = -2.0/3.0 ! Taylor series coefficients
real x,c
x = gamah * dt
IF ( x > 0.0001 ) THEN ! Use formula
	c = 1-EXP(-2.0*x)
ELSE ! Use Taylor expansion for low x
	c = x * ( c1 + x * ( c2 + x * ( c3 + x * c4 ) ) ) 
END IF
    !c = SQRT ( c )

m = 1.0


noiseh = sqrt(c*Th/m)

do j = 1 , natom !1st half kick velocity
	vxd(j) = vxd(j) + 0.5*dt*fxd(j) 
	vyd(j) = vyd(j) + 0.5*dt*fyd(j)
	vzd(j) = vzd(j) + 0.5*dt*fzd(j)
enddo

do j = 1 , natom !1st half kick position
		
		rxd(j) = rxd(j) + 0.5*dt*vxd(j)
		ryd(j) = ryd(j) + 0.5*dt*vyd(j) 
		rzd(j) = rzd(j) + 0.5*dt*vzd(j) 
enddo

do j = 1 , natom ! random velocities and friction step
	call Gaussrandom(im, random)
	vxd(j) = vxd(j)*EXP(-x) + noiseh*random
	call Gaussrandom(im, random)
	vyd(j) = vyd(j)*EXP(-x) + noiseh*random
	call Gaussrandom(im, random)
	vzd(j) = vzd(j)*EXP(-x) + noiseh*random
enddo


do j = 1 , natom !2nd half kick position
		
		rxd(j) = rxd(j) + 0.5*dt*vxd(j)
		ryd(j) = ryd(j) + 0.5*dt*vyd(j) 
		rzd(j) = rzd(j) + 0.5*dt*vzd(j) 
enddo


    	
call forcecal(natom,xboxd,yboxd,zboxd,rxd,ryd,rzd,rcut,eps,sigma,PE,fxd,fyd,fzd,Wxd,Wyd,Wzd)	


KE = 0.0
	
do j=1,natom
	vxd(j) = vxd(j) + 0.5*dt*fxd(j)/m 
	vyd(j) = vyd(j) + 0.5*dt*fyd(j)/m
	vzd(j) = vzd(j) + 0.5*dt*fzd(j)/m

	KE = KE + 0.5*(vxd(j)*vxd(j) + vyd(j)*vyd(j) + vzd(j)*vzd(j))

enddo

return
end


subroutine fix_langevin_rgn(im,natom,xboxd,yboxd,zboxd,dt,rxd,ryd,rzd,rcut,eps,sigma,vxd,vyd,vzd,fxd,fyd,fzd,Wxd,Wyd,Wzd,&
KE,PE,gamah,gamac,Th,Tc)
integer :: im,natom
real :: xboxd,yboxd,rcut,dt,zboxd
real :: rxd(natom),ryd(natom),rzd(natom)
real :: vxd(natom),vyd(natom),vzd(natom)
real :: fxd(natom),fyd(natom),fzd(natom)
real :: Wxd(natom),Wyd(natom),Wzd(natom)
real :: eps(natom),sigma(natom)
integer :: i ,j 
real KE,PE
real :: Th,Tc,gamah,gamac,gfrih,gfric,noiseh,noisec,random,zpos
real :: m !mass
real xhot,xcold,chot,ccold
REAL, PARAMETER :: c1 = 2.0, c2 = -2.0, c3 = 4.0/3.0, c4 = -2.0/3.0 ! Taylor series coefficients

m = 1.0


!noiseh = sqrt(2*gamah*Th*dt/m)
!noisec = sqrt(2*gamac*Tc*dt/m)

do j = 1 , natom !1st half kick velocity
	vxd(j) = vxd(j) + 0.5*dt*fxd(j)/m 
	vyd(j) = vyd(j) + 0.5*dt*fyd(j)/m
	vzd(j) = vzd(j) + 0.5*dt*fzd(j)/m
enddo

do j = 1 , natom !1st half kick position
		
		rxd(j) = rxd(j) + 0.5*dt*vxd(j)
		ryd(j) = ryd(j) + 0.5*dt*vyd(j) 
		rzd(j) = rzd(j) + 0.5*dt*vzd(j) 
enddo


xhot = gamah * dt
IF ( xhot > 0.0001 ) THEN ! Use formula
	chot = 1-EXP(-2.0*xhot)
ELSE ! Use Taylor expansion for low x
	chot = xhot * ( c1 + xhot * ( c2 + xhot * ( c3 + xhot * c4 ) ) ) 
END IF
!chot = SQRT ( chot )

xcold = gamac * dt
IF ( xcold > 0.0001 ) THEN ! Use formula
	ccold = 1-EXP(-2.0*xcold)
ELSE ! Use Taylor expansion for low x
	ccold = xcold * ( c1 + xcold * ( c2 + xcold * ( c3 + xcold * c4 ) ) ) 
END IF
!chot = SQRT ( ccold )

noiseh = sqrt(chot*Th/m)
noisec = sqrt(ccold*Tc/m)

do j = 1 , natom ! random velocities and friction step
	zpos = rzd(j) - zboxd*anint(rzd(j)/zboxd)
	if (abs(zpos) < zboxd/4.0) then
		call Gaussrandom(im, random)
		vxd(j) = vxd(j)*EXP(-xhot) + noiseh*random
		call Gaussrandom(im, random)
		vyd(j) = vyd(j)*EXP(-xhot) + noiseh*random
		call Gaussrandom(im, random)
		vzd(j) = vzd(j)*EXP(-xhot) + noiseh*random
	else if	(abs(zpos) > zboxd/4.0) then
		call Gaussrandom(im, random)
		vxd(j) = vxd(j)*EXP(-xcold) + noisec*random
		call Gaussrandom(im, random)
		vyd(j) = vyd(j)*EXP(-xcold) + noisec*random
		call Gaussrandom(im, random)
		vzd(j) = vzd(j)*EXP(-xcold) + noisec*random
	endif	
enddo


do j = 1 , natom !2nd half kick position
		
		rxd(j) = rxd(j) + 0.5*dt*vxd(j)
		ryd(j) = ryd(j) + 0.5*dt*vyd(j) 
		rzd(j) = rzd(j) + 0.5*dt*vzd(j) 
enddo

call forcecal(natom,xboxd,yboxd,zboxd,rxd,ryd,rzd,rcut,eps,sigma,PE,fxd,fyd,fzd,Wxd,Wyd,Wzd)
	


KE = 0.0
	
do j=1,natom
	vxd(j) = vxd(j) + 0.5*dt*fxd(j)/m 
	vyd(j) = vyd(j) + 0.5*dt*fyd(j)/m
	vzd(j) = vzd(j) + 0.5*dt*fzd(j)/m

	KE = KE + 0.5*(vxd(j)*vxd(j) + vyd(j)*vyd(j) + vzd(j)*vzd(j))

enddo

return
end

subroutine initialization (m,rxd,ryd,rzd,vxd,vyd,vzd,xboxd,yboxd,zboxd)
integer :: m,i
real :: xboxd,yboxd,zboxd
real :: rxd(m),ryd(m),rzd(m)
real :: vxd(m),vyd(m),vzd(m)
character :: a
open (unit = 2, file = 'fcc-lattice4000_disorder-0.7-2.dat' )
open (unit = 3, file = 'init_vel-4000.dat' )
read(2,*)
read(2,*) a,a ,a ,a ,xboxd,yboxd,zboxd
write(*,*) xboxd,yboxd,zboxd
do i=1,m
read (2 , *)rxd(i),ryd(i),rzd(i)
end do
do i=1,m
read (3 , *)  vxd(i),vyd(i),vzd(i)
end do
return
end


subroutine init (T,m,rxd,ryd,rzd,vxd,vyd,vzd,rxprevd,ryprevd,rzprevd)
integer :: m,i	
real :: rxd(m), ryd(m), rzd(m)
real :: vxd(m),vyd(m),vzd(m)
real :: sumvx,sumvy,sumvz,sumv2
real::rxprevd(m),ryprevd(m),rzprevd(m)
real :: lambda,dt,T
dt=.001
sumvx=0.0
sumvy=0.0
sumvz=0.0

do i=1,m
	sumvx=sumvx+vxd(i)
	sumvy=sumvy+vyd(i)
	sumvz=sumvz+vzd(i)
	sumv2=sumv2+vxd(i)**2+vyd(i)**2+vzd(i)**2
end do
	sumvx = sumvx/m
	sumvy = sumvy/m
	sumvz = sumvz/m
 !lambda= sqrt (3*m*T/sumv2)
 !lambda = 1
sumv2=0.0
do i=1,m
	vxd(i)= (vxd(i)-sumvx)
	vyd(i)= (vyd(i)-sumvy)
	vzd(i)= (vzd(i)-sumvz)
	
	sumv2=sumv2+vxd(i)**2+vyd(i)**2+vzd(i)**2
end do
 lambda= sqrt ((3*m-3)*T/sumv2)
 !lambda = 1
do i = 1,m
         vxd(i) = lambda*vxd(i)
         vyd(i) = lambda*vyd(i)
         vzd(i) = lambda*vzd(i)

	rxprevd(i)=rxd(i)-vxd(i)*dt
	ryprevd(i)=ryd(i)-vyd(i)*dt
	rzprevd(i)=rzd(i)-vzd(i)*dt
enddo

	

return
end 

!********re initialization of epsilon for each particle
subroutine modify_Eps(im,m,xboxd,yboxd,zboxd,rxd,ryd,rzd,eps,sigma,delz,zbin,Tzprof,Thot,Tcold)
integer ::i,j,m,im,zbin , bin
real :: rxd(m), ryd(m), rzd(m)
real :: eps(m),sigma(m)
real :: xboxd,yboxd,zboxd,boxhalf
real :: zpos , delz
real :: Tzprof(zbin)
real :: sigmaH, sigmaC, Thot,Tcold , mtan
boxhalf = zboxd/2.0
sigmaH = 1.0
sigmaC = 1.3
mtan = (sigmaH-sigmaC)/(Thot-Tcold)
do j = 1, m
	zpos = rzd(j) - zboxd*anint(rzd(j)/zboxd)
	bin = int((zpos+boxhalf)/delz) + 1
	
	if (Tzprof(bin) .ge. Thot) then
		sigma(j) = sigmaH
	else if (Tzprof(bin) .le. Tcold) then
		sigma(j) = sigmaC	
	else
		sigma(j) = sigmaH + mtan*(Tzprof(bin) - Thot)	
	endif
	eps(j) = 1.0
	!write(*,*) sigma(j) , j , im , Tzprof(bin) , mtan , Thot, Tcold
enddo

end

subroutine modify_Eps3(im,m,xboxd,yboxd,zboxd,rxd,ryd,rzd,eps,sigma,delz,zbin,Tzprof,Thot,Tcold)
integer ::i,j,m,im
real :: rxd(m), ryd(m), rzd(m)
integer :: zbin,bin
real :: eps(m),sigma(m),Tzprof(zbin)
real :: xboxd,yboxd,zboxd,boxhalf,delz
real :: zpos
real :: sigmaH, sigmaC, Thot,Tcold , mtan
boxhalf = zboxd/2.0
sigmaH = 1.0
sigmaC = 1.5
mtan = (sigmaH-sigmaC)/(Thot-Tcold)
do j = 1, m
	zpos = rzd(j) - zboxd*anint(rzd(j)/zboxd)
	bin = int((zpos+boxhalf)/delz) + 1
	if (Tzprof(bin) .ge. Thot) then
		sigma(j) = sigmaH
	else if (Tzprof(bin) .le. Tcold) then
		sigma(j) = sigmaC	
	else
		sigma(j) = sigmaH + mtan*(Tzprof(bin) - Thot)	
	endif
	eps(j) = 1.0
	!write(*,*) sigma(j) , j , im , Tzprof(bin) !, mtan , Thot, Tcold
	
enddo

end

!if (Tzprof(bin) .gt. Tcold .and. Tzprof(bin) .lt. Thot) then
subroutine modify_Eps2(im,m,xboxd,yboxd,zboxd,rxd,ryd,rzd,eps,sigma,delz,zbin,Tzprof,Thot,Tcold)
integer ::i,j,m,im,zbin , bin
real :: rxd(m), ryd(m), rzd(m)
real :: eps(m),sigma(m)
real :: boxd,zboxd,boxhalf
real :: zpos , delz
real :: Tzprof(zbin)
real :: sigmaH, sigmaC, Thot,Tcold , mtan
boxhalf = zboxd/2.0
mtan = (sigmaH-simgaC)/(Thot-Tcold)
sigmaH = 1.0
sigmaC = 1.7
do j = 1, m
	zpos = rzd(j) - zboxd*anint(rzd(j)/zboxd)
	if (abs(zpos) < zboxd/4.0) then
		eps(j) = 1.0
		sigma(j) = sigmaH
	else if  (abs(zpos) > zboxd/4.0) then
		eps(j) = 1.0
		sigma(j) = sigmaC
	
	endif
enddo

end




subroutine forcecal (m,xboxd,yboxd,zboxd,rxd,ryd,rzd,rcutd,eps,sigma,pot,fxd,fyd,fzd,Wxd,Wyd,Wzd)
integer ::i,j,m
real :: rxdist,rydist,rzdist
real ::rxd(m),ryd(m),rzd(m)
real :: fxd(m),fyd(m),fzd(m)
real :: Wxd(m),Wyd(m),Wzd(m)
real :: xboxd,yboxd,zboxd,pot,rijsq,rcutd,rcutdsq,ff
real :: eps(m),sigma(m)
real :: sigmaeff,epseff
do i=1,m
	fxd(i)=0.0
	fyd(i)=0.0
	fzd(i)=0.0
	Wxd(i)=0.0
	Wyd(i)=0.0
	Wzd(i)=0.0
end do
pot=0.0
!rcutdsq=rcutd**2

!$omp parallel &
!!$omp   shared ( h, n ) &
!$omp shared (rxd,ryd,rzd,sigma,eps,xboxd,yboxd,zboxd,m,rcutd) &
!$omp private (i,j,rxdist,rydist,rzdist,sigmaeff,epseff,ff,sr2,sr6,rijsq,rcutdsq)

!$omp do reduction(+:pot,fxd,fyd,fzd,Wxd,Wyd,Wzd)
do i=1,m-1
	do j=i+1,m
		rxdist=rxd(i)-rxd(j)
		rydist=ryd(i)-ryd(j)
		rzdist=rzd(i)-rzd(j)
	
		rxdist=rxdist-xboxd*nint(rxdist/xboxd)
		rydist=rydist-yboxd*nint(rydist/yboxd)
		rzdist=rzdist-zboxd*nint(rzdist/zboxd)

	 	rijsq=rxdist**2+rydist**2+rzdist**2
		sigmaeff = (sigma(i)+sigma(j))/2.0
		rcutdsq  = (rcutd*sigmaeff)**2
		if (rijsq.lt.rcutdsq) then
			!sigmaeff = (sigma(i)+sigma(j))/2.0
			sigmaeff = sigmaeff*sigmaeff
			epseff = sqrt(eps(i)*eps(j))
			sr2=sigmaeff/rijsq
			sr6=sr2**3
			ff=48*epseff*sr2*sr6*(sr6-0.5)/sigmaeff
			fxd(i)=fxd(i)+ff*rxdist
			fyd(i)=fyd(i)+ff*rydist
			fzd(i)=fzd(i)+ff*rzdist
			fxd(j)=fxd(j)-ff*rxdist
			fyd(j)=fyd(j)-ff*rydist
			fzd(j)=fzd(j)-ff*rzdist
			
			Wxd(i)=Wxd(i)+ ff*rxdist*rxdist
			Wyd(i)=Wyd(i)+ ff*rydist*rydist
			Wzd(i)=Wzd(i)+ ff*rzdist*rzdist
			Wxd(j)=Wxd(j)+ ff*rxdist*rxdist
			Wyd(j)=Wyd(j)+ ff*rydist*rydist
			Wzd(j)=Wzd(j)+ ff*rzdist*rzdist

			pot = pot+4*epseff*sr6*(sr6-1.0)
		
		end if
		
	end do
end do
!$omp end do 
!$omp end parallel 
return
end






 



subroutine temp_prof(im,m,xboxd,yboxd,zboxd,rxd,ryd,rzd,vxd,vyd,vzd,delz,zbin,Tprof)
integer :: i , j, k, m,bin,zbin,im
real :: rxd(m), ryd(m), rzd(m)
real :: vxd(m), vyd(m), vzd(m)
real :: delz,zpos,kesum
real :: Tprof(zbin),count1(zbin)
real :: xboxd,yboxd,zboxd,boxhalf,KE,track,dev
real :: vxcom , vycom, vzcom
boxhalf = zboxd/2.0
KE = 0.0
track = 0.0
dev = 0.0
do i=1,zbin
	Tprof(i)=0.0
	count1(i) = 0.0
end do

vxcom = 0.0
vycom = 0.0
vzcom = 0.0
do i=1,m
	zpos = rzd(i) - zboxd*anint(rzd(i)/zboxd)
	
	if (zpos >-1.35 .and. zpos < -0.35) then
		vxcom = vxcom + vxd(i) 
		vycom = vycom + vyd(i) 
		vzcom = vzcom + vzd(i)
		track =  track + 1.0
	endif
enddo

vxcom = vxcom/track
vycom = vycom/track
vzcom = vzcom/track

KE = 0.0
track = 0.0
dev = 0.0

do i=1,m
	zpos = rzd(i) - zboxd*anint(rzd(i)/zboxd)
	bin = int((zpos+boxhalf)/delz) + 1
	kesum = 0.5*(vxd(i)*vxd(i)+vyd(i)*vyd(i)+vzd(i)*vzd(i))
	!Tprof(bin) = Tprof(bin) + kesum
	count1(bin) = count1(bin) + 1.0
	if (zpos >-1.35 .and. zpos < -0.35) then
		!KE =  KE + 0.5*((vxd(i)-vxcom)*(vxd(i)-vxcom)+(vyd(i)-vycom)*(vyd(i)-vycom)+(vzd(i)-vzcom)*(vzd(i)-vzcom))
		KE =  KE + 0.5*(vxd(i)*vxd(i)+vyd(i)*vyd(i)+vzd(i)*vzd(i))
		track =  track + 1.0
		dev = dev + (1.77 - 0.5*(vxd(i)*vxd(i)+vyd(i)*vyd(i)+vzd(i)*vzd(i)))
	endif
enddo

!write(*,*) im , 2*KE/(3*track-3) , KE/track, dev/track

do i = 1,zbin
	if (count1(i) > 1) then
		!Tprof(i) = 2*Tprof(i)/(3*count1(i)-3)
		!write(*,*) Tprof(i), count1(i) , i*delz
	else 	
		write(*,*) count1(i) 
		!stop
	endif
enddo

return
end

subroutine temp_prof_z(im,m,xboxd,yboxd,zboxd,rxd,ryd,rzd,vxd,vyd,vzd,delz,zbin,Tprofz)
integer :: i , j, k, m,bin,zbin,im
real :: rxd(m), ryd(m), rzd(m)
real :: vxd(m), vyd(m), vzd(m)
real :: delz,zpos,kesum
real :: Tprofz(zbin),count1(zbin)
real :: xboxd,yboxd,zboxd,boxhalf,KE,track,dev
real :: vxcom(zbin) , vycom(zbin), vzcom(zbin)
boxhalf = zboxd/2.0
KE = 0.0
track = 0.0
dev = 0.0
do i=1,zbin
	Tprofz(i)=0.0
	count1(i) = 0.0
	vxcom(i) = 0.0
	vycom(i) = 0.0
	vzcom(i) = 0.0
end do


do i=1,m
	zpos = rzd(i) - zboxd*anint(rzd(i)/zboxd)
	bin = int((zpos+boxhalf)/delz) + 1
	!write(*,*)bin , zpos , im , i
	vxcom(bin) = vxcom(bin) + vxd(i) 
	vycom(bin) = vycom(bin) + vyd(i) 
	vzcom(bin) = vzcom(bin) + vzd(i)
	count1(bin) = count1(bin) + 1.0
	!write(*,*) rzd(i) , bin , i , im
enddo

vxcom(:) = vxcom(:)/count1(:)
vycom(:) = vycom(:)/count1(:)
vzcom(:) = vzcom(:)/count1(:)

KE = 0.0
track = 0.0
dev = 0.0

do i=1,m
	zpos = rzd(i) - zboxd*anint(rzd(i)/zboxd)
	bin = int((zpos+boxhalf)/delz) + 1
	!kesum = 0.5*(vxd(i)*vxd(i)+vyd(i)*vyd(i)+vzd(i)*vzd(i))
kesum = 0.5*((vxd(i)-vxcom(bin))**2+(vyd(i)-vycom(bin))**2+(vzd(i)-vzcom(bin))**2)
	Tprofz(bin) = Tprofz(bin) + kesum
	!count1(bin) = count1(bin) + 1.0
	
	!KE =  KE + 0.5*((vxd(i)-vxcom)*(vxd(i)-vxcom)+(vyd(i)-vycom)*(vyd(i)-vycom)+(vzd(i)-vzcom)*(vzd(i)-vzcom))
	!KE =  KE + 0.5*(vxd(i)*vxd(i)+vyd(i)*vyd(i)+vzd(i)*vzd(i))
	
	
enddo

!write(*,*) im , 2*KE/(3*track-3) , KE/track, dev/track

do i = 1,zbin
	if (count1(i) > 1) then
		Tprofz(i) = 2*Tprofz(i)/(3*count1(i)-3)
		!write(*,*) Tprof(i), count1(i) , i*delz
	else 	
		!write(*,*) count1(i) , (i-1)*delz-zboxd/2.0+delz/2.0, im
		!stop
	endif
enddo

return
end



subroutine pair_correlation (m,rxd,ryd,rzd,xboxd,yboxd,zboxd,delr,max_bin,hist,densd)
integer :: i,j,k,m,bin,max_bin
real :: rxd(m), ryd(m), rzd(m)
real :: Rxij,Ryij,Rzij,xboxd,yboxd,zboxd,Rijsq,rijd
real :: delr
real :: hist (max_bin)
real :: densd,Rmax,Rmin,nideal
real,parameter :: pi=3.14

  do i=1,max_bin
         hist(i)=0.0
      end do

	do i=1,m-1
   do j = i+1,m
   Rxij = rxd(i) - rxd(j)
   Ryij = ryd(i) - ryd(j)
   Rzij = rzd(i) - rzd(j)
   Rxij = Rxij - xboxd*nint(Rxij/xboxd)
   Ryij = Ryij - yboxd*nint(Ryij/yboxd)
   Rzij = Rzij - zboxd*nint(Rzij/zboxd)
   Rijsq = Rxij*Rxij+Ryij*Ryij+Rzij*Rzij
   rijd = sqrt(Rijsq)
	 bin = int (rijd/delr)+1

	if (bin.le.max_bin) then 
  		hist (bin)= hist(bin)+2
	end if

	end do

end do
	


 do i = 1,max_bin
        Rmin = (i-1)*delr
				Rmax = Rmin + delr
				nideal = densd*4*pi*(Rmax**3 - Rmin**3)/3
        !rij = delg*(i+0.5)
        !vb = (real(i+1)**3-real(i)**3)*delg**3
        !nid = (4.0/3.0)*pi*densd*vb
        !write(15,*)(3*g(i))/(Rmax**3 - Rmi**3)*0.7/real(natom)/5000.0/4/pi
        hist(i) = hist (i)/nideal/m
        
        

 end do     
 


end subroutine

!--------------------------------------------------------------------------
subroutine Gaussrandom(tm, randomd)
real :: randomd , u1,u2, ran1 , ran2,toss
integer tm
data pi/3.1415927/ 
u1 = ranf()
u2 = ranf()
ran1 = cos(2*pi*u1)*sqrt(-2.*log(u2));
ran2 = cos(2*pi*u2)*sqrt(-2.*log(u1));

toss = ranf()
if (toss > 0.5) then
   randomd = ran1
else
   randomd = ran2
endif
if (ran1-1 .eq. ran1) then 
   randomd = ran2
   !write(*,*) 'ran1 is a inf' , tm
endif
if (ran2-1 .eq. ran2) then 
   randomd = ran1
   !write(*,*) 'ran2 is a inf' , tm
endif
return
end  

!****************************************************************
      REAL function ranf()
      common /rjran/ i3,i2,i1,i0

!     berkeley random number generator
!     range changed to 0 < 1

      INTEGER  I0, I1, I2, I3, J0, J1, J2, J3, K0, K1, K2, K3
      INTEGER  M0, M1, M2, M3, MM
      REAL     T1, T2, T3, T4     

      parameter (m3=647,m2=1442,m1=3707,m0= 373)
      parameter (t4=2.0**48,t3=2.0**36,t2=2.0**24,t1=2.0**12)
      parameter (mm=4096)

!     mm = 2 ** 12

!     the random number is:
!     (i3*t3+i2*t2+i1*t1+i0)/2.0**48
!     the multiplier is:
!     (m3*t3+m2*t2+m1*t1+m0)

      ranf = float(i3)/t1 + float(i2)/t2 + float(i1)/t3 + float(i0)/t4
      if(ranf.ge.0.9999999)  ranf = 0.0

!     multiply i's and m's

      j0 = m0 * i0
      j1 = m0 * i1 + m1 * i0
      j2 = m0 * i2 + m1 * i1 + m2 * i0
      j3 = m0 * i3 + m1 * i2 + m2 * i1 + m3 * i0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      i0 = mod(k0,mm)
      i1 = mod(k1,mm)
      i2 = mod(k2,mm)
      i3 = mod(k3,mm)
      return
      end      


      subroutine ranint(istart)

!     initialize random number generator

      INTEGER  IRAN(4), ISTART

      iran(1)=12345 + istart*1000
      iran(2)=12345 + istart*1000
      iran(3)=12345 + istart*1000
      iran(4)=12345 + istart*1000
      call ranset(iran)
      call ranget(iran)

      return
      end


      SUBROUTINE RANSET(IRAN)
      COMMON /RJRAN/ II3,II2,II1,II0

      INTEGER    IRAN(4), MM, NN
      PARAMETER (MM = 4096, NN = 100000)
      INTEGER    II0, II1, II2, II3
      INTEGER    I0, I1, I2, I3, J0, J1, J2, J3

      i3 = iran(1)
      i2 = iran(2)
      i1 = iran(3)
      i0 = iran(4)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii0)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii1)
      call divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,mm,ii2)
      call divide(j3,j2,j1,j0,i3,i2,i1,i0,nn,mm,ii3)
      return
      end


      subroutine ranget(iran)
      common /rjran/ ii3,ii2,ii1,ii0

      INTEGER    mm, m10, m21, m20, m32, m31, m30

      parameter (mm = 100000)
      parameter (m10 = 4096)
      parameter (m21 =  167, m20 = 77216)
      parameter (m32 =    6, m31 = 87194, m30 = 76736)

      INTEGER    iran(4)
      INTEGER    ii3, ii2, ii1, ii0  
      INTEGER    J0, J1, J2, J3, K0, K1, K2, K3

      j0 = ii0 + m10 * ii1 + m20 * ii2 + m30 * ii3
      j1 =                   m21 * ii2 + m31 * ii3
      j2 =                               m32 * ii3
      j3 =                                       0
      k0 = j0
      k1 = j1 + k0 / mm
      k2 = j2 + k1 / mm
      k3 = j3 + k2 / mm
      iran(4) = mod(k0,mm)
      iran(3) = mod(k1,mm)
      iran(2) = mod(k2,mm)
      iran(1) = mod(k3,mm)
      return
      end


      subroutine divide(i3,i2,i1,i0,j3,j2,j1,j0,nn,id,ir)

      INTEGER     I3,I2,I1,I0,J3,J2,J1,J0,ID,IR,NN,K0,K1,K2,K3

!     given the integer i = i0 + nn * (i1 + nn * (i2 + nn * i3))
!     this routine calculates j = i / id and ir = mod(i, id)
!     j is expressed as i, ir is just an integer

      j3 = i3 / id
      k3 = mod(i3, id)
      k2 = k3 * nn + i2
      j2 = k2 / id
      k2 = mod(k2, id)
      k1 = k2 * nn + i1
      j1 = k1 / id
      k1 = mod(k1, id)
      k0 = k1 * nn + i0
      j0 = k0 / id
      ir = mod(k0, id)
      return
      end
!*********************************************************************************








