!===============================================================================
! cfd.f90 -- 2D incompressible Navier-Stokes boilerplate (single file)
!
!   Nondimensional equations, rho = 1:
!     du/dt + (u.grad)u = -grad p + (1/Re) lap u,     div u = 0
!
!   Default case: lid-driven square cavity (top wall moves at u_lid, others
!   no-slip). Change the parameters and apply_velocity_bc for other flows.
!
!   Discretisation
!     - staggered MAC grid: p at cell centres, u on x-faces, v on y-faces
!     - explicit Euler in time, 2nd-order central differences in space
!     - Chorin projection: predictor -> pressure Poisson (SOR) -> corrector
!
!   Build / run:
!     gfortran -O2 -o cfd cfd.f90 && ./cfd
!   Output:
!     centerline_u.dat  u(y) along x = lx/2   (compare with Ghia et al. 1982)
!     centerline_v.dat  v(x) along y = ly/2
!     field.vtk         cell-centred u, v, p  (open in ParaView / VisIt)
!===============================================================================
program cfd
  implicit none
  integer, parameter :: dp = kind(1.0d0)

  ! ---------------------------------------------------------------- parameters
  integer,  parameter :: nx = 64, ny = 64          ! interior cells
  real(dp), parameter :: lx = 1.0_dp, ly = 1.0_dp  ! domain size
  real(dp), parameter :: re = 100.0_dp             ! Reynolds number, u_lid*lx/nu
  real(dp), parameter :: u_lid = 1.0_dp            ! lid speed
  real(dp), parameter :: t_end = 30.0_dp           ! end time
  real(dp), parameter :: safety = 0.5_dp           ! fraction of stable dt
  real(dp), parameter :: omega = 1.7_dp            ! SOR relaxation, 1 < omega < 2
  real(dp), parameter :: tol_p = 1.0e-6_dp         ! Poisson residual (RMS) tolerance
  integer,  parameter :: itmax_p = 10000           ! max SOR sweeps per step
  integer,  parameter :: n_print = 500             ! diagnostics every n steps

  ! -------------------------------------------------------------------- fields
  ! u(i,j): x-face between cells i and i+1, i = 0..nx;   j = 0..ny+1 (ghosts)
  ! v(i,j): y-face between cells j and j+1, j = 0..ny;   i = 0..nx+1 (ghosts)
  ! p(i,j): cell centre,  i = 0..nx+1, j = 0..ny+1 (ghost ring)
  real(dp) :: u(0:nx, 0:ny+1), v(0:nx+1, 0:ny), p(0:nx+1, 0:ny+1)
  real(dp) :: f(0:nx, 0:ny+1), g(0:nx+1, 0:ny), rhs(nx, ny)

  real(dp) :: dx, dy, nu, dt, t, res, du2dx, duvdy, duvdx, dv2dy, lap
  real(dp) :: umax, vmax, divmax
  integer  :: i, j, it, step

  dx = lx / nx
  dy = ly / ny
  nu = 1.0_dp / re
  u = 0.0_dp;  v = 0.0_dp;  p = 0.0_dp
  f = 0.0_dp;  g = 0.0_dp
  t = 0.0_dp;  step = 0
  call apply_velocity_bc(u, v)

  do while (t < t_end)
    ! ── stable time step, explicit Euler + central differences ──
    ! diffusive limit dt <= 1/(2 nu (1/dx^2 + 1/dy^2)), convective (CFL) limit
    ! dt <= min(dx/|u|, dy/|v|); cell Reynolds number |u|dx/nu should stay < 2.
    umax = max(maxval(abs(u)), 1.0e-12_dp)
    vmax = max(maxval(abs(v)), 1.0e-12_dp)
    dt = safety * min(0.5_dp / (nu * (1.0_dp/dx**2 + 1.0_dp/dy**2)), dx/umax, dy/vmax)
    if (t + dt > t_end) dt = t_end - t

    ! ── predictor: F = u + dt (nu lap u - d(uu)/dx - d(uv)/dy), G likewise ──
    ! Harlow & Welch 1965, Phys. Fluids 8, 2182 -- MAC staggered grid
    ! convective terms in conservative form, fluxes interpolated to cell corners
    do j = 1, ny
      do i = 1, nx-1
        du2dx = ((0.5_dp*(u(i,j)+u(i+1,j)))**2 - (0.5_dp*(u(i-1,j)+u(i,j)))**2) / dx
        duvdy = (0.25_dp*(v(i,j)+v(i+1,j))*(u(i,j)+u(i,j+1)) &
               - 0.25_dp*(v(i,j-1)+v(i+1,j-1))*(u(i,j-1)+u(i,j))) / dy
        lap   = (u(i+1,j) - 2.0_dp*u(i,j) + u(i-1,j)) / dx**2 &
              + (u(i,j+1) - 2.0_dp*u(i,j) + u(i,j-1)) / dy**2
        f(i,j) = u(i,j) + dt * (nu*lap - du2dx - duvdy)
      end do
    end do
    do j = 1, ny-1
      do i = 1, nx
        duvdx = (0.25_dp*(u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j)) &
               - 0.25_dp*(u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j))) / dx
        dv2dy = ((0.5_dp*(v(i,j)+v(i,j+1)))**2 - (0.5_dp*(v(i,j-1)+v(i,j)))**2) / dy
        lap   = (v(i+1,j) - 2.0_dp*v(i,j) + v(i-1,j)) / dx**2 &
              + (v(i,j+1) - 2.0_dp*v(i,j) + v(i,j-1)) / dy**2
        g(i,j) = v(i,j) + dt * (nu*lap - duvdx - dv2dy)
      end do
    end do
    ! wall-normal velocities are prescribed, so F, G there equal the BC values;
    ! this makes the Poisson problem's Neumann condition dp/dn = 0 consistent.
    f(0,:) = u(0,:);  f(nx,:) = u(nx,:)
    g(:,0) = v(:,0);  g(:,ny) = v(:,ny)

    ! ── pressure Poisson: lap p = div(F,G)/dt, SOR with dp/dn = 0 on walls ──
    ! Chorin 1968, Math. Comp. 22, 745 -- projection method
    do j = 1, ny
      do i = 1, nx
        rhs(i,j) = ((f(i,j) - f(i-1,j))/dx + (g(i,j) - g(i,j-1))/dy) / dt
      end do
    end do
    do it = 1, itmax_p
      call apply_pressure_bc(p)
      do j = 1, ny
        do i = 1, nx
          p(i,j) = (1.0_dp - omega)*p(i,j) + omega / (2.0_dp/dx**2 + 2.0_dp/dy**2) &
                 * ((p(i+1,j) + p(i-1,j))/dx**2 + (p(i,j+1) + p(i,j-1))/dy**2 - rhs(i,j))
        end do
      end do
      call apply_pressure_bc(p)
      res = 0.0_dp
      do j = 1, ny
        do i = 1, nx
          res = res + ((p(i+1,j) - 2.0_dp*p(i,j) + p(i-1,j))/dx**2 &
                     + (p(i,j+1) - 2.0_dp*p(i,j) + p(i,j-1))/dy**2 - rhs(i,j))**2
        end do
      end do
      res = sqrt(res / (nx*ny))
      if (res < tol_p) exit
    end do
    ! pure-Neumann problem: p is defined up to a constant, pin its mean to zero
    p(1:nx,1:ny) = p(1:nx,1:ny) - sum(p(1:nx,1:ny)) / (nx*ny)
    call apply_pressure_bc(p)

    ! ── corrector: u = F - dt dp/dx, v = G - dt dp/dy ──
    ! Chorin 1968, Math. Comp. 22, 745
    do j = 1, ny
      do i = 1, nx-1
        u(i,j) = f(i,j) - dt * (p(i+1,j) - p(i,j)) / dx
      end do
    end do
    do j = 1, ny-1
      do i = 1, nx
        v(i,j) = g(i,j) - dt * (p(i,j+1) - p(i,j)) / dy
      end do
    end do
    call apply_velocity_bc(u, v)

    t = t + dt
    step = step + 1
    if (mod(step, n_print) == 0 .or. t >= t_end) then
      divmax = 0.0_dp
      do j = 1, ny
        do i = 1, nx
          divmax = max(divmax, abs((u(i,j)-u(i-1,j))/dx + (v(i,j)-v(i,j-1))/dy))
        end do
      end do
      write(*,'(a,i8,a,f10.4,a,es10.3,a,i6,a,es10.3,a,es10.3)') 'step', step, '  t', t, &
           '  dt', dt, '  sor', it, '  res', res, '  max|div u|', divmax
    end if
  end do

  call write_output()

contains

  ! Ghost-cell velocity BCs. Tangential velocity is set by reflection so the
  ! wall value (average of interior and ghost) equals the wall speed.
  ! Edit here for inflow/outflow, periodic, slip walls, etc.
  subroutine apply_velocity_bc(u, v)
    real(dp), intent(inout) :: u(0:nx, 0:ny+1), v(0:nx+1, 0:ny)
    u(0, :)  = 0.0_dp                          ! left wall,   normal
    u(nx, :) = 0.0_dp                          ! right wall,  normal
    u(:, 0)    = -u(:, 1)                      ! bottom wall, tangential = 0
    u(:, ny+1) = 2.0_dp*u_lid - u(:, ny)       ! lid,         tangential = u_lid
    v(:, 0)  = 0.0_dp                          ! bottom wall, normal
    v(:, ny) = 0.0_dp                          ! lid,         normal
    v(0, :)    = -v(1, :)                      ! left wall,   tangential = 0
    v(nx+1, :) = -v(nx, :)                     ! right wall,  tangential = 0
  end subroutine apply_velocity_bc

  ! Homogeneous Neumann dp/dn = 0 on all walls via ghost copy.
  subroutine apply_pressure_bc(p)
    real(dp), intent(inout) :: p(0:nx+1, 0:ny+1)
    p(0, 1:ny)    = p(1, 1:ny)
    p(nx+1, 1:ny) = p(nx, 1:ny)
    p(1:nx, 0)    = p(1:nx, 1)
    p(1:nx, ny+1) = p(1:nx, ny)
  end subroutine apply_pressure_bc

  subroutine write_output()
    integer :: ic, jc, un
    ! u(y) on the vertical centreline (x-face index nx/2 sits at x = lx/2 for even nx)
    open(newunit=un, file='centerline_u.dat', status='replace')
    write(un,'(a)') '# y  u(x=lx/2, y)'
    write(un,'(2es16.8)') 0.0_dp, 0.0_dp
    do jc = 1, ny
      write(un,'(2es16.8)') (jc-0.5_dp)*dy, u(nx/2, jc)
    end do
    write(un,'(2es16.8)') ly, u_lid
    close(un)
    ! v(x) on the horizontal centreline
    open(newunit=un, file='centerline_v.dat', status='replace')
    write(un,'(a)') '# x  v(x, y=ly/2)'
    write(un,'(2es16.8)') 0.0_dp, 0.0_dp
    do ic = 1, nx
      write(un,'(2es16.8)') (ic-0.5_dp)*dx, v(ic, ny/2)
    end do
    write(un,'(2es16.8)') lx, 0.0_dp
    close(un)
    ! legacy ASCII VTK, cell-centred values
    open(newunit=un, file='field.vtk', status='replace')
    write(un,'(a)') '# vtk DataFile Version 3.0'
    write(un,'(a)') 'cfd.f90 output'
    write(un,'(a)') 'ASCII'
    write(un,'(a)') 'DATASET STRUCTURED_POINTS'
    write(un,'(a,3i6)')     'DIMENSIONS ', nx, ny, 1
    write(un,'(a,3es16.8)') 'ORIGIN ', 0.5_dp*dx, 0.5_dp*dy, 0.0_dp
    write(un,'(a,3es16.8)') 'SPACING ', dx, dy, 1.0_dp
    write(un,'(a,i10)')     'POINT_DATA ', nx*ny
    write(un,'(a)') 'SCALARS pressure double 1'
    write(un,'(a)') 'LOOKUP_TABLE default'
    do jc = 1, ny
      do ic = 1, nx
        write(un,'(es16.8)') p(ic, jc)
      end do
    end do
    write(un,'(a)') 'VECTORS velocity double'
    do jc = 1, ny
      do ic = 1, nx
        write(un,'(3es16.8)') 0.5_dp*(u(ic-1,jc) + u(ic,jc)), 0.5_dp*(v(ic,jc-1) + v(ic,jc)), 0.0_dp
      end do
    end do
    close(un)
  end subroutine write_output

end program cfd
