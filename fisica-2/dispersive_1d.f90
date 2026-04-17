program dispersive_1d
  implicit none

  ! Parámetros
  integer, parameter :: N = 200
  integer, parameter :: steps = 500
  real(8), parameter :: dx = 1.0d0
  real(8), parameter :: dt = 0.01d0

  ! Parámetros físicos
  real(8), parameter :: eps_inf = 1.0d0
  real(8), parameter :: mu0 = 1.0d0
  real(8), parameter :: wp = 1.0d0
  real(8), parameter :: w0 = 1.0d0
  real(8), parameter :: gamma = 0.0d0   ! sin pérdidas

  ! Campos
  real(8) :: E(N), H(N)
  real(8) :: P(N), V(N)

  integer :: i, t

  ! Inicialización
  E = 0.0d0
  H = 0.0d0
  P = 0.0d0
  V = 0.0d0

  ! Pulso inicial
  do i = 80, 120
     E(i) = exp(-((i-100.0d0)/10.0d0)**2)
  end do

  ! Loop temporal
  do t = 1, steps

     ! --- Actualizar H ---
     do i = 1, N-1
        H(i) = H(i) - (dt/mu0)*(E(i+1) - E(i))/dx
     end do

     ! --- Actualizar V (oscilador) ---
     do i = 1, N
        V(i) = V(i) + dt*(wp**2 * eps_inf * E(i) - w0**2 * P(i) - gamma * V(i))
     end do

     ! --- Actualizar P ---
     do i = 1, N
        P(i) = P(i) + dt * V(i)
     end do

     ! --- Actualizar E ---
     do i = 2, N
        E(i) = E(i) + (dt/eps_inf)*((H(i) - H(i-1))/dx - V(i))
     end do

     ! Condiciones de borde simples
     E(1) = 0.0d0
     E(N) = 0.0d0

     ! Imprimir cada cierto tiempo
     if (mod(t,50) == 0) then
        print *, "Paso:", t, "E centro:", E(N/2)
     end if

  end do

end program dispersive_1d
