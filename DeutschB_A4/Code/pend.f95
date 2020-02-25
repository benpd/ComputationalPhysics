module setup
	use numtype
	implicit none
	
	integer, parameter :: n_eq = 3
	real(dp):: q, g, omega_d

end module setup

program pendulum
	use setup
	implicit none
	real(dp) :: t, dt, dg, gmax, gmin, tmax
	real(dp), dimension(n_eq) :: y
	integer :: i, imax
!'--Constant--'	
	q = 2._dp
	omega_d = 2._dp/3
	
	y(1) = 0._dp 	!omega
	y(2) = 0.2 		!theta
	y(3) = 0._dp 	!phi - driving phase
!'-----T-----'	
	t = 0._dp
	dt = 0.5_dp
	tmax = 4000*2*pi
!'---Step----'
	imax = 5000
	gmin = .9_dp 
	gmax = 1.6_dp 
	dg = (gmax-gmin)/ imax
!'----Loop----'
	do i=0, imax
	
		g = gmin + i * dg   
		
		t = 0._dp
		do while (t<tmax) 
 	
 	 		if( pi < y(2) ) then 
   				y(2) = y(2) - 2*pi 
   			else if (y(2) < -pi )  then 
   				y(2) = y(2) + 2*pi 
   			end if 
!'---Cutoff---'	
		if (t > 9.0_dp * tmax/10) then 
		
			if (mod(y(3), 2*pi) <= 0.01_dp)   write (50,*) g, y(1)
		
			end if
!'----RK4-----'	
			call rk4step(t, y, dt)
		end do 
	end do
	
end program pendulum