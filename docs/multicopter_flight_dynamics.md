# Multicopter Flight Dynamics

A complete nonlinear flight-mechanical model of a multicopter, covering the
propulsion system (rotor aerodynamics, motor dynamics, battery discharge), the
rigid-body equations of motion (rotational and translational), and the
aerodynamic effects acting on the airframe. The formulation follows Stephan,
*Multicopter Flight Control* (Springer, 2025), Chapter 3. This is the full
verification model: it is presented as written in the reference, without
introducing the simplifications used later for control design.

Notation: vectors are written in bold or with explicit components. The
body-fixed frame is denoted $b$ and the geodetic (North-East-Down) frame is
denoted $g$. $R_{gb}$ is the rotation matrix from body to geodetic frame. A
rotor index $i = 1, \dots, N$ runs over the $N$ propulsion units.

---

## 1. Overview and State

The complete model couples three sub-systems:

1. The **propulsion system**: for each rotor, the aerodynamic loads (thrust,
   in-plane force, drag torque) from blade-element momentum theory, the
   brushless motor dynamics, and the shared battery discharge.
2. The **rigid-body equations of motion**: rotational kinematics, rotational
   dynamics (with gyroscopic coupling), and translational motion.
3. The **airframe aerodynamics**: drag forces and torques not produced by the
   rotors.

The system state comprises the rotor speeds, the attitude quaternion, the body
angular rate, the translational velocity, and the position:

$$
\mathbf{x} = \big[\, \Omega_1, \dots, \Omega_N,\;\; \mathbf{q},\;\; \boldsymbol{\omega},\;\; \mathbf{v}_g,\;\; \mathbf{r}_g \,\big]
$$

where

- $\Omega_i$ are the rotor speeds [rad/s],
- $\mathbf{q} = (q_0, \mathbf{q}_v)^\top \in \mathbb{S}^3$ is the attitude quaternion (body to geodetic, scalar first),
- $\boldsymbol{\omega} = (p, q, r)^\top$ is the body angular rate [rad/s],
- $\mathbf{v}_g = (v_N, v_E, v_D)^\top$ is the velocity in the geodetic frame [m/s],
- $\mathbf{r}_g = (r_N, r_E, r_D)^\top$ is the position in the geodetic frame [m].

The control input is the PWM command vector $\boldsymbol{\delta} = (\delta_1, \dots, \delta_N)^\top$, with $\delta_i \in [0, 1]$, supplied by the flight control system. External disturbances enter as a torque $\boldsymbol{\tau}_d$ and a force $\mathbf{F}_d$.

---

## 2. Propulsion System

### 2.1 Rotor Aerodynamics (Blade-Element Momentum Theory)

For each rotor, the thrust $T$, the in-plane force $H \ge 0$, and the drag
torque $Q$ are computed from blade-element momentum theory (BEMT) as functions
of the rotor speed $\Omega$, the unperturbed inflow $\mathbf{V}^\infty = (V_h^\infty, V_z^\infty)^\top$, and the air density $\rho$:

$$
\begin{aligned}
    T_i = c_{GE}\, c_T(\rho, \lambda_i, \mu_i)\, \Omega_i^2, \\
    H_i = c_H(\rho, \lambda_i, \mu_i)\, \Omega_i^2, \\
    Q_i = c_Q(\rho, \lambda_i, \mu_i)\, \Omega_i^2
\end{aligned}
$$

The thrust, in-plane force, and drag coefficients depend on the inflow ratio
$\lambda$ and the advance ratio $\mu$:

$$
\begin{aligned}
    \lambda_i = \frac{V_{z,i}^\infty + \nu_{\text{ind},i}}{\Omega_i R}, \\
    \mu_i = \frac{V_{h,i}^\infty}{\Omega_i R}
\end{aligned}
$$

where $R$ is the rotor radius and $\nu_{\text{ind}}$ is the induced velocity,
itself the solution of

$$
\begin{aligned}
    \nu_{\text{ind},i} = \frac{T_i}{2 \rho A \, \lVert \mathbf{V}_{R,i} \rVert}, \\
    \mathbf{V}_{R,i} = \mathbf{V}_i^\infty + (0,\; \nu_{\text{ind},i})^\top
\end{aligned}
$$

This is a coupled, implicit system, typically solved with Newton's method for
each rotor. The aerodynamic coefficients themselves follow from the airfoil and
blade geometry:

$$
c_T = \frac{\sigma \pi}{4}\, \rho\, C_{l\alpha}\, R^4 (\theta_0 - \alpha_0)\!\left(1 + \tfrac{3}{2}\mu^2 - \tfrac{2}{3}\lambda\right)
$$

$$
c_H = \frac{\sigma \pi}{4}\, \rho\, C_{l\alpha}\, R^4 \mu \!\left(\tfrac{C_d}{C_{l\alpha}} + (\theta_0 - \alpha_0)\lambda\right)
$$

$$
c_Q = \frac{\sigma \pi}{8}\, \rho\, C_d\, R^5 \!\left(1 + 3\mu^2 + \tfrac{R}{C_d}(\lambda c_T - \mu c_H)\right)
$$

with solidity ratio $\sigma = N_b c / (\pi R)$, lift slope $C_{l\alpha}$, zero-lift angle $\alpha_0$, blade twist $\theta_0$, and airfoil drag coefficient $C_d$. The ground-effect coefficient $c_{GE} \ge 1$ accounts for the increase in thrust near the ground and depends on the rotor height above ground.

The unperturbed inflow itself depends on the vehicle motion. With $p_i$ the rotor position relative to the centre of gravity:

$$
\mathbf{v}_{A,i} = R_{bg}\big(\mathbf{v} - \mathbf{v}_w\big) + (\boldsymbol{\omega} - \boldsymbol{\omega}_w) \times \mathbf{p}_i
$$

$$
\begin{aligned}
    V_{h,i}^\infty = \big\lVert (I_3 - \mathbf{n}_i \mathbf{n}_i^\top)\, \mathbf{v}_{A,i} \big\rVert, \\
    V_{z,i}^\infty = \mathbf{n}_i^\top \mathbf{v}_{A,i}
\end{aligned}
$$

so the rotor loads are coupled to the body rate, the translational velocity, and the wind. This is the coupling that disappears in the near-hover simplification.

### 2.2 Motor Dynamics

Each rotor is driven by a brushless motor through an electronic speed controller
(ESC). With the electrical dynamics assumed fast, the motor current is

$$
I_{\text{mot}} = \frac{E\, \delta - k_e\, \Omega}{R_{\text{mot}}}
$$

where $E$ is the ESC terminal voltage, $\delta$ the PWM signal, $k_e$ the
electrical motor constant, and $R_{\text{mot}}$ the motor resistance. The motor
torque is $M_R = k_T I_{\text{mot}}$, and with $k_T = k_e$ (a common BLDC
approximation) the motor dynamics become

$$
J_R\, \dot{\Omega}_i = k_T\big(E\,\delta_i - k_T\,\Omega_i\big)\big/R_{\text{mot}} - M_F - Q_i
$$

where $J_R$ is the combined rotor and motor inertia, $M_F$ is motor friction,
and $Q_i$ is the aerodynamic drag torque loading the rotor. Equivalently, when
the speed command $\Omega_c$ is mapped to PWM, the closed-loop motor behaves as
a first-order lag:

$$
\begin{aligned}
    T_{\text{ESC}}\, \dot{\Omega}_i = \Omega_{c,i} - \Omega_i, \\
    T_{\text{ESC}} = \frac{J_R R_{\text{mot}}}{k_T^2}
\end{aligned}
$$

with $T_{\text{ESC}}$ the speed-tracking time constant.

### 2.3 Battery Dynamics

All propulsion units share one battery. Its open-circuit voltage depends on the
state of charge $\xi \in (0, \xi_0]$:

$$
\hat{E}(\xi) = E_0 - E_{\text{pol}}\frac{\xi_0 - \xi}{\xi} + E_{\text{exp}} \exp\!\left(-\frac{\xi_0 - \xi}{\xi_{\text{exp}}}\right)
$$

where $E_0$ is the nominal voltage, $E_{\text{pol}}$ the polarization voltage,
$E_{\text{exp}}$ the exponential-zone overshoot, and $\xi_{\text{exp}}$ the
related decay constant. The state of charge is depleted by the total current
drawn across all rotors:

$$
\dot{\xi} = -\sum_{i=1}^{N} \delta_i\, \frac{\hat{E}(\xi)\,\delta_i - k_T\,\Omega_i}{R_{\text{mot}}}
$$

The full discharge curve has a steep initial drop (set by $E_{\text{exp}}$,
$\xi_{\text{exp}}$), a slow central plateau where $\hat{E} \approx E_0$, and a
rapid collapse near deep discharge as the polarization resistance rises.

---

## 3. Rigid-Body Equations of Motion

### 3.1 Rotational Kinematics

The attitude quaternion evolves with the body angular rate through the $4 \times 3$ matrix $\Theta(\mathbf{q})$:

$$
\dot{\mathbf{q}} = \tfrac{1}{2}\, \Theta(\mathbf{q})\, \boldsymbol{\omega}
$$

Quaternions are used to avoid gimbal lock. Numerical integration can let the
norm $\mathbf{q}^\top \mathbf{q}$ drift, so a norm-restoring variant is often
used in practice:

$$
\dot{\mathbf{q}} = \tfrac{1}{2}\, \Theta(\mathbf{q})\, \boldsymbol{\omega} + K\big(1 - \mathbf{q}^\top \mathbf{q}\big)\, \mathbf{q}
$$

with a small gain $K > 0$. At the nominal condition $\mathbf{q}^\top \mathbf{q} = 1$ this reduces to the original kinematics. The quaternion maps to the rotation matrix by

$$
R_{gb} = I_3 + 2 q_0 [\mathbf{q}_v \times] + 2 [\mathbf{q}_v \times]^2
$$

### 3.2 Angular Momentum and Gyroscopic Coupling

The total angular momentum of the vehicle is the rigid-body momentum plus the
momentum stored in the spinning rotors:

$$
\begin{aligned}
    \mathbf{h} = J_b\, \boldsymbol{\omega} + \mathbf{h}_R, \\
    \mathbf{h}_R = \sum_{i=1}^{N} \vartheta_i\, J_R\, \mathbf{n}_i\, \Omega_i
\end{aligned}
$$

where $J_b = J_b^\top \succeq 0$ is the inertia tensor of the complete vehicle,
$\mathbf{n}_i$ is the rotor axis, and $\vartheta_i \in \{-1, +1\}$ is the spin
direction ($+1$ counterclockwise, $-1$ clockwise, viewed from above). The rotor
momentum has a time derivative driven by the rotor angular acceleration:

$$
\dot{\mathbf{h}}_R = \sum_{i=1}^{N} \vartheta_i\, J_R\, \mathbf{n}_i\, \dot{\Omega}_i
$$

### 3.3 Rotational Dynamics

Applying the principle of angular momentum gives the rotational equation of
motion:

$$
J_b\, \dot{\boldsymbol{\omega}} + \dot{\mathbf{h}}_R + \boldsymbol{\omega} \times \big(J_b\, \boldsymbol{\omega} + \mathbf{h}_R\big) = \boldsymbol{\tau}_\Sigma + \boldsymbol{\tau}_D
$$

The terms are:

- $\boldsymbol{\tau}_\Sigma = (L_\Sigma, M_\Sigma, N_\Sigma)^\top$ is the
  aerodynamic propulsion torque from the rotors (roll, pitch, yaw).
- $\boldsymbol{\tau}_D$ is the airframe drag torque, capturing aerodynamic
  torque not produced by the rotors.
- $\dot{\mathbf{h}}_R$ is the rotor spin-up reaction torque: changing rotor
  speed changes rotor momentum, and the airframe feels the reaction.
- $\boldsymbol{\omega} \times (J_b \boldsymbol{\omega} + \mathbf{h}_R)$ is the
  gyroscopic coupling. The part $\boldsymbol{\omega} \times J_b \boldsymbol{\omega}$
  is the rigid-body cross-coupling; the part $\boldsymbol{\omega} \times \mathbf{h}_R$
  is the gyroscopic effect of the spinning rotors, which resists tilting and
  couples the axes.

The spin-up term $\dot{\mathbf{h}}_R$ vanishes only if the rotors are not
accelerating ($\dot{\Omega}_i = 0$); the gyroscopic term $\boldsymbol{\omega} \times \mathbf{h}_R$ vanishes only at rest ($\boldsymbol{\omega} = 0$) or with zero rotor momentum.

### 3.4 Propulsion Torque and Force

The propulsion torque sums, over all rotors, a force-induced component (thrust
and in-plane force acting through the lever arm $\mathbf{p}_i$) and a drag
component (rotor drag torque):

$$
\boldsymbol{\tau}_\Sigma = \sum_{i=1}^{N} \mathbf{p}_i \times \big(\mathbf{n}_i T_i + \mathbf{m}_i H_i\big) - \vartheta_i\, \mathbf{n}_i\, Q_i
$$

The propulsion force sums the thrust and in-plane force of each rotor:

$$
\mathbf{F}_\Sigma = \sum_{i=1}^{N} \mathbf{n}_i T_i + \mathbf{m}_i H_i
$$

where $\mathbf{m}_i$ is the direction of the in-plane force, parallel to the
horizontal inflow. The roll and pitch torque arise from thrust through the arm;
the yaw torque arises from rotor drag (and, through $\dot{\mathbf{h}}_R$, from
spin-up).

### 3.5 Translational Motion

Newton's second law in the geodetic frame, with the propulsion force rotated
from body to geodetic frame:

$$
m\, \dot{\mathbf{v}} = R_{gb}\, \mathbf{F}_\Sigma + \mathbf{F}_D + m\, \mathbf{g}, \qquad \mathbf{g} = (0, 0, g)^\top
$$

$$
\dot{\mathbf{r}} = \mathbf{v}
$$

where $m$ is the vehicle mass, $\mathbf{F}_D$ the airframe drag force, and
gravity points in the positive (downward) $z$ direction of the NED frame.

---

## 4. Airframe Aerodynamics

Separate from the rotor loads, the airframe itself produces aerodynamic drag.
The air-relative velocity and angular rate are

$$
\begin{aligned}
    \mathbf{v}_A = \mathbf{v} - \mathbf{v}_w, \\
    \boldsymbol{\omega}_A = \boldsymbol{\omega} - \boldsymbol{\omega}_w
\end{aligned}
$$

with $\mathbf{v}_w$ and $\boldsymbol{\omega}_w$ the wind. The scalar airframe
drag combines a viscous (linear) term, dominant at low airspeed, and an inertial
(quadratic) term, dominant at higher airspeed:

$$
D = \eta\, c_D\, \beta\, V_A + \tfrac{1}{2}\, \rho\, c_D'\, \beta^2\, V_A^2
$$

with airspeed $V_A = \lVert \mathbf{v}_A \rVert$. As a force and torque on the
airframe, this is written

$$
\mathbf{F}_D = -\big(c_v + c_v'\, \lVert \mathbf{v}_A \rVert\big)\, \mathbf{v}_A + \mathbf{F}_d
$$

$$
\boldsymbol{\tau}_D = -\big(C_\omega + C_\omega'\, \lVert \boldsymbol{\omega}_A \rVert\big)\, \boldsymbol{\omega}_A + \boldsymbol{\tau}_d
$$

where $c_v, c_v' \ge 0$ are the translational drag coefficients, $C_\omega, C_\omega' \succeq 0$ the rotational drag coefficients, and $\mathbf{F}_d$, $\boldsymbol{\tau}_d$ external disturbances or modelling error. If the drag application point is offset from the centre of gravity by $\Delta\mathbf{p}$, the drag torque gains an airspeed-induced term $\Delta\mathbf{p} \times R_{bg}\, \mathbf{F}_D$.

---

## 5. Complete Model (Eq. 3.101)

Collecting all of the above gives the complete flight-mechanical model:

$$
J_R\, \dot{\Omega}_i = k_T\big(\hat{E}\,\delta_i - k_T\,\Omega_i\big)\big/R_{\text{mot}} - Q_i - M_F, \qquad i = 1, \dots, N
$$

$$
\dot{\mathbf{q}} = \tfrac{1}{2}\, \Theta(\mathbf{q})\, \boldsymbol{\omega}
$$

$$
\begin{aligned}
J_b\, \dot{\boldsymbol{\omega}} = {}& \boldsymbol{\tau}_\Sigma
- \big(C_\omega + C_\omega'\, \lVert \boldsymbol{\omega}_A \rVert\big)\, \boldsymbol{\omega}_A + \boldsymbol{\tau}_d \\
& - \sum_{i=1}^{N} \vartheta_i\, J_R\, \mathbf{n}_i\, \dot{\Omega}_i
- \boldsymbol{\omega} \times \Big(J_b\, \boldsymbol{\omega} + \sum_{i=1}^{N} \vartheta_i\, J_R\, \mathbf{n}_i\, \Omega_i\Big)
\end{aligned}
$$

$$
\begin{aligned}
m\, \dot{\mathbf{v}} = {}& R_{gb}\, \mathbf{F}_\Sigma
- \big(c_v + c_v'\, \lVert \mathbf{v}_A \rVert\big)\, \mathbf{v}_A + \mathbf{F}_d
+ m\, \mathbf{g}
\end{aligned}
$$

$$
\dot{\mathbf{r}} = \mathbf{v}
$$

with the propulsion torque and force

$$
\begin{aligned}
    \boldsymbol{\tau}_\Sigma = \sum_{i=1}^{N} \mathbf{p}_i \times \big(\mathbf{n}_i T_i + \mathbf{m}_i H_i\big) - \vartheta_i\, \mathbf{n}_i\, Q_i, \\
    \mathbf{F}_\Sigma = \sum_{i=1}^{N} \mathbf{n}_i T_i + \mathbf{m}_i H_i
\end{aligned}
$$

$$
\begin{aligned}
    \boldsymbol{\omega}_A = \boldsymbol{\omega} - \boldsymbol{\omega}_w, \\
    \mathbf{v}_A = \mathbf{v} - \mathbf{v}_w, \\
    R_{gb} = I_3 + 2 q_0 [\mathbf{q}_v \times] + 2 [\mathbf{q}_v \times]^2
\end{aligned}
$$

The thrust $T_i$, in-plane force $H_i$, and rotor drag torque $Q_i$ follow from
the BEMT model of Section 2.1 as functions of the rotor speed $\Omega_i$ and the
inflow. The battery voltage $\hat{E}$ follows from the state of charge $\xi$ and
its discharge dynamics.

The inputs to the system are the PWM command vector $\boldsymbol{\delta}$, the
wind $\mathbf{v}_w$ and $\boldsymbol{\omega}_w$, and the disturbances
$\boldsymbol{\tau}_d$ and $\mathbf{F}_d$.

---

## 6. Model Parameters

| Symbol | Meaning |
|--------|---------|
| $m$ | Vehicle mass |
| $J_b$ | Vehicle inertia tensor (including rotating parts) |
| $J_R$ | Combined rotor and motor inertia of one unit |
| $k_T$ | Motor torque constant |
| $R_{\text{mot}}$ | Motor ohmic resistance |
| $M_F$ | Motor friction |
| $\vartheta_i$ | Rotor spin direction $(\pm 1)$ |
| $\mathbf{p}_i$ | Rotor position relative to centre of gravity |
| $\mathbf{n}_i$ | Rotor thrust axis |
| $c_T, c_H, c_Q$ | Thrust, in-plane force, drag coefficients (BEMT) |
| $c_{GE}$ | Ground-effect coefficient |
| $c_v, c_v'$ | Translational airframe drag coefficients |
| $C_\omega, C_\omega'$ | Rotational airframe drag coefficients |
| $E_0, E_{\text{pol}}, E_{\text{exp}}, \xi_{\text{exp}}$ | Battery discharge parameters |

---

## Reference

Stephan, J. *Multicopter Flight Control*. Springer, 2025. Chapter 3 (Flight
Mechanical Model); complete model summarised in Eq. 3.101.
