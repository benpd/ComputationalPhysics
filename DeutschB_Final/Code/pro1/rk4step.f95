subroutine rk4step(t, y, dt)

	use setup
	implicit none
	
	real(dp), intent(inout) :: t 
	real(dp), intent(in) :: dt 
	real(dp), dimension(n_eq), intent(inout) :: y 
	real(dp), dimension(n_eq) :: k1, k2, k3, k4, dy
	
	call deriv(t, y, dt, k1) 
	call deriv(t + dt/2, y + k1/2, dt, k2) 
	call deriv(t + dt/2, y + k2/2, dt, k3)
	call deriv(t + dt, y + k3, dt, k4)
	
	dy = (k1 + 2*k2 + 2*k3 + k4)/6
	
	t = t + dt 
	y = y + dy 
	
	contains
	
		subroutine deriv(t, y, dt, k)
		
			implicit none
			
			real(dp), intent(in) :: t, dt 
			real(dp), dimension(n_eq), intent(in) :: y
			real(dp), dimension(n_eq), intent(out) :: k
			real(dp), dimension(n_eq) :: f
			
			f(1) = -(y(2)+y(3))
			f(2) = (y(1)+a*(y(2)))
			f(3) = (b+y(3)*(y(1)-c))
			
			k = f*dt
			
		end subroutine deriv
	
end subroutine rk4step