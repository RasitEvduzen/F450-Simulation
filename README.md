# F450 Quadcopter Flight Simulator

A modular MATLAB simulator for an F450 quadcopter. It combines a nonlinear
flight dynamics model, a multi-rate cascaded flight controller that mirrors a
real flight stack, and trajectory generators for waypoint missions and
figure-8 paths. The model is built on the equations of motion from Stephan,
*Multicopter Flight Control* (Springer, 2025).

The simulator is intended as a validated testbed for developing and testing
flight control algorithms and trajectory generation.

<!-- Add a hero image / GIF of the animation here -->

## Repository layout

| File | Purpose |
|------|---------|
| `quadParams.m` | All configuration: physical model, controller gains, simulation rates, visual settings |
| `quadDynamics.m` | Nonlinear flight dynamics model + RK4 integrator |
| `quadControl.m` | Multi-rate cascaded controller and control allocation (mixer) |
| `quadTrajectory.m` | Trajectory generators (waypoint mission, lemniscate) |
| `quadVisualize.m` | STL animation and analysis figures |
| `quatUtils.m` | Quaternion helper functions |
| `runMission.m` | Top-level script: run a trajectory (mission or lemniscate) |

## Running

Put the STL files (`Base.STL`, `MotorPropCW.STL`, `MotorPropCCW.STL`) on the
MATLAB path, then run `runMission.m`. Select the trajectory at the top of the
script:

```matlab
mode = 'mission';      % 'mission' (square) or 'lemniscate'
```

Everything tunable lives in `quadParams.m`.

## Dynamics

The plant is a nonlinear rigid-body model integrated at 1000 Hz with a
fixed-step RK4 scheme. The state vector has 17 elements:

```
x = [Omega(4); q(4); omega(3); v(3); r(3)]
```

where `Omega` are the four rotor speeds, `q` is the body-to-geodetic attitude
quaternion (scalar first), `omega` is the body angular rate, and `v`, `r` are
the velocity and position in a North-East-Down (NED) geodetic frame. The input
is the four per-rotor PWM commands in the range [0, 1].

### Equations of motion

The model follows Eq. 3.101 of the reference text:

```
Motor:        JR*Omdot_i = kT*(Ehat*delta_i - Omega_i) - Q_i - MF
Kinematics:   qdot       = 1/2 * Theta(q) * omega
Rotational:   Jb*omdot   = tau_Sigma - omega x ( Jb*omega + h_R )
Translation:  m*vdot     = R(q)*F_Sigma + m*g
Position:     rdot       = v
```

The per-rotor thrust and drag torque are quadratic in rotor speed, and the
total body force and torque are summed over the four rotors:

```
T_i = cT*Omega_i^2,   Q_i = cQ*Omega_i^2
tau_Sigma = sum_i  p_i x (n_i*T_i) - th_i*n_i*Q_i
F_Sigma   = sum_i  n_i*T_i
h_R       = sum_i  th_i*JR*n_i*Omega_i      (rotor angular momentum)
```

Here `n_i` is the rotor thrust axis, `p_i` its position, and `th_i = +/-1` its
spin direction. Gravity is `g = [0; 0; +g]` in NED (down is positive z).

### Modelling assumptions

1. Rigid body with fixed mass and inertia.
2. Near-hover / low-speed flight: rotor inflow is taken as zero, so the
   in-plane force vanishes and thrust and drag are quadratic in rotor speed.
3. Out of ground effect.
4. No airframe or rotor aerodynamic drag, no wind, no external disturbances.
5. Rotors horizontal, thrust axis vertical in the body frame.
6. Gyroscopic momentum coupling is kept; the rotor spin-up reaction torque is
   omitted (consistent with the control-design model, since JR is small).
7. Battery voltage held constant (no discharge dynamics).

Assumptions 2 and 4 are the first to break at high speed, where rotor inflow
becomes nonzero and aerodynamic drag becomes significant.

### Parameters (F450)

| Symbol | Value | Description |
|--------|-------|-------------|
| m | 1.5 kg | mass |
| Jb | diag(0.02, 0.02, 0.04) kg m^2 | body inertia |
| arm | 0.225 m | rotor radius from center of gravity (450 mm span) |
| JR | 6e-5 kg m^2 | rotor + motor inertia |
| kT | 0.003 | motor torque constant (tau_mot = JR/kT = 20 ms) |
| cT | 1.5e-5 | thrust coefficient |
| cQ | 2.5e-7 | drag coefficient |

The X-configuration uses PX4 motor numbering: M1 front-right (CCW),
M2 rear-left (CCW), M3 front-left (CW), M4 rear-right (CW).

<!-- Add a diagram of the airframe / motor layout here -->

## Control architecture

The controller is a multi-rate cascaded P / PID controller whose loop rates
mirror a real PX4 flight stack. Each loop runs at its own rate and holds its
output zero-order between updates, exactly as a real flight stack schedules its
tasks.

| Loop | Rate | Input | Output |
|------|------|-------|--------|
| Rate (inner) | 1000 Hz | body rate error | body torque |
| Attitude | 250 Hz | attitude error | body rate setpoint |
| Position / velocity (outer) | 50 Hz | position error | acceleration, attitude setpoint |

The cascade flows from outside in:

```
position (P) -> velocity (PID) -> acceleration
             -> attitude bridge (acceleration to tilt + collective)
             -> attitude (quaternion P) -> rate setpoint
             -> rate (PID) -> body torque
             -> control allocation (mixer) -> per-rotor PWM
```

The position loop produces a velocity setpoint with velocity feedforward. The
velocity loop produces an acceleration setpoint. The acceleration-to-attitude
bridge converts the desired acceleration into a tilt direction (desired
attitude quaternion) and a collective thrust, since a quadcopter can only
produce thrust along its body z axis. The attitude loop uses a yaw-weighted
reduced-attitude formulation to prioritise tilt over heading. The rate loop is
a PID with anti-windup on the integrator. The mixer maps the desired collective
thrust and body torques to four rotor commands using the airframe geometry.

The control laws match the core PX4 v1.15 implementation
(`rate_control.cpp`, `AttitudeControl.cpp`, `PositionControl.cpp`).

The controller is pluggable: any controller that implements the common
interface (`mixer`, `initState`, `step`) can be registered in `quadControl.m`,
while the allocation and rotor-command mapping stay shared.

<!-- Add tracking / states / motor analysis figures here -->

## Trajectories

### Waypoint mission (square)

Fly straight, then turn in place. The vehicle translates to each waypoint with
its heading held constant (an LSPB position profile, so it flies straight with
no mid-flight yaw), then on arrival holds position and rotates to the next
heading with a cosine-eased yaw. The default mission takes off, flies a square
loop stopping and turning at each corner, closes the loop, and lands.

<!-- Add square mission animation / trajectory plot here -->

### Lemniscate (figure-8)

A lemniscate of Gerono in the x-y plane with a slow altitude ramp:

```
x = A sin(w t)
y = (B/2) sin(2 w t)
z = z0 + (z1 - z0) t / Ttot
```

The nose is kept tangent to the path (yaw = atan2(vy, vx)), like a coordinated
turn. The trajectory yaw rate is fed forward to the attitude loop so the
controller does not lag behind the continuously turning heading.

<!-- Add figure-8 animation / trajectory plot here -->

## To do

- [ ] System identification: closed-loop chirp excitation with batch-LS ARX
      estimation (separate `runSysId.m` entry point)
- [ ] INDI controller (incremental nonlinear dynamic inversion)
- [ ] LQR controller (state-space optimal control)
- [ ] MPC controller (model predictive control with constraints)
- [ ] Add linear aerodynamic drag to relax the near-hover assumption for
      high-speed flight

## Reference

Stephan, J. *Multicopter Flight Control*. Springer, 2025.
