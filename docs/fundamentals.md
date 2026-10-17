# Theory and foundations

## Model foundations

This chapter follows a soil column from rainfall input to pressure changes and slope stability. Each spatial slope unit uses the same column model with its own geometry and material properties.

The flow calculation is based on the one-dimensional Richards equation, with exponential Gardner laws for conductivity and water content. These laws allow an analytical solution: formulas describe the response to changes in infiltration. The model then estimates water-table rise from retained water and evaluates an infinite-slope factor of safety.

The analytical framework follows Srivastava and Yeh (1991) and the TRIGRS publications by Baum and colleagues (2002, 2008, 2010, with the 2013 correction). The pressure, retention, and storage choices specific to unsatstab are explained beside the relevant equations. Full source details are in [equation provenance](#equation-provenance) and [references](#references).

The default Gardner/Srivastava–Yeh formulation preserves the intrinsic Gardner relations. The alternative Baum et al. formulation uses the TRIGRS conversions for comparison with published results. The tables below distinguish mathematical consistency with the stated governing laws from reproduction of a published formulation; they are not an official erratum to TRIGRS. Practical configuration is described in [inputs and outputs](inputs-outputs.md#single-column-parameter-conventions).

## Notation and equation sources

Depth is positive downward from the ground surface; time is measured from the start of the simulated rainfall history. The equations below use SI seconds when conductivity and flux are in m/s. Public model time inputs and outputs use hours, with the required factor of 3600 applied internally; see [time convention](inputs-outputs.md#public-time-convention). The main quantities are:

| Symbol | Meaning | Unit |
| --- | --- | --- |
| $z$ | Vertical depth below the ground surface | m |
| $z_{\mathrm{w},0}$, $z_{\mathrm{w}}(t)$ | Initial and current water-table depths | m |
| $z_{\mathrm{max}}$ | Maximum vertical depth of the reported soil column; finite-base depth when that boundary is selected | m |
| $z_{\mathrm{f}}$ | Depth selected for pressure and factor-of-safety time series | m |
| $\delta$ | Ground-surface slope angle | degrees or radians, as required by the formula |
| $P(t)$, $I(t)$ | Rainfall intensity and effective infiltration | m/s |
| $I_{\mathrm{ZLT}}$ | Background infiltration maintaining the initial steady state | m/s |
| $q(z,t)$, $q_{\mathrm{b}}(t)$ | Downward flux and its value at the initial water-table depth | m/s |
| $q_{\mathrm{c}}$ | Effective basal drainage capacity | m/s |
| $K(z,t)$, $K_{\mathrm{s}}$ | Hydraulic conductivity and its saturated value | m/s |
| $\psi(z,t)$ | Pressure head: pressure divided by water unit weight | m |
| $\theta(z,t)$, $\theta_{\mathrm{s}}$, $\theta_{\mathrm{r}}$ | Volumetric water content, saturated value, and residual value | m³/m³ |
| $\alpha$, $a$ | Intrinsic Gardner parameter and slope-adjusted coefficient | m⁻¹ |
| $\Delta h(t)$ | Estimated water-table rise | m |
| $S(t)$, $V_{\mathrm{f}}(z,t)$ | Retained water and pore volume available to fill, each per unit area | m of water |
| $\mathrm{FS}$ | Factor of safety: available resistance relative to driving stress | dimensionless |

A field such as $K(z,t)$ varies with depth and time. A parameter such as $K_{\mathrm{s}}$ is constant within a homogeneous column. In the governing-equation derivations, $K(z,t)=K(\psi(z,t))$ and $\theta(z,t)=\theta(\psi(z,t))$ denote the Gardner constitutive relations evaluated at the local pressure head; their arguments are sometimes omitted to keep the algebra readable. The implemented profile reconstruction is described separately below.

Pressure head is negative where water is under suction, zero at the water table, and normally positive below it. It is expressed in metres, not pascals; pore pressure is $p_{\mathrm{w}}=\gamma_{\mathrm{w}}\psi$, where $\gamma_{\mathrm{w}}$ is the water unit weight.

*Source* notes identify published relationships. *Adaptation* notes identify implementation choices or algebraic rewrites. Descriptive subscripts such as $\mathrm{s}$ (saturated), $\mathrm{r}$ (residual), $\mathrm{b}$ (basal), and $\mathrm{a}$ (antecedent) are labels. The index $i$ denotes an integration step; $r_m$ denotes the $m$th positive characteristic root, a dimensionless variable introduced in the analytical response. Original equation numbers are collected in [equation provenance](#equation-provenance).

## Model overview

unsatstab uses analytical formulas to describe how infiltration changes conductivity in the unsaturated portion of the soil between the ground surface and the initial water table. It then estimates water-table rise from accumulated water, updates pressure and water content, and calculates the infinite-slope factor of safety.

The flow formulas keep the original domain even when the estimated water-table rises. The updates to pressure and water content use the additional assumptions described in this chapter. Together, these steps provide a simplified model of wetting and slope stability; the program does not solve a new Richards flow problem with a moving boundary at each time step.

$$
P(t) \longrightarrow I(t) \longrightarrow q_{\mathrm{b}}(t) \longrightarrow S(t) \longrightarrow \Delta h(t) \longrightarrow \psi(z,t) \longrightarrow \mathrm{FS}(z,t)
$$

The sequence separates water entering the surface, water reaching the base, water retained there, and the resulting pressure and stability response. Water entering at the surface does not immediately arrive at the water table.

Here, $P$ is gross rainfall intensity and $I=\min(P,K_{\mathrm{s}})$ is the effective infiltration applied to each column; rainfall above $K_{\mathrm{s}}$ is discarded locally. The quantity $q_{\mathrm{b}}$ is flux at the initial water-table depth, $\Delta h$ is estimated water-table rise, $\psi$ is pressure head, and $\mathrm{FS}$ is factor of safety. The analytical background is described in the [model references](#references), including the TRIGRS publications and the Gardner and Srivastava–Yeh formulations.[^note-1]

## Geometry, column base, and initial water table {#geometry-and-initial-water-table}

The ground surface is $z=0$. The complete reported column extends to $z_{\mathrm{max}}$. Pressure, water content, conductivity, and factor of safety are evaluated throughout this depth range. An observation plane at $z_{\mathrm{f}}$ selects a pressure/FS time series; it does not truncate the column or alter its hydraulic boundary. If no plane is selected, the series refer to $z_{\mathrm{f}}=z_{\mathrm{max}}$.

Both the initial water-table depth and maximum depth are vertical distances measured downward. A smaller depth means a higher position. Their treatment depends on the basal assumption:

| Initial configuration | Finite restrictive base | Infinite basal domain |
| --- | --- | --- |
| Requested water table above the base, $z_{\mathrm{w},0}<z_{\mathrm{max}}$ | Retain that water table. The column includes an initial saturated layer down to the base. | Retain that water table and report the profile down to $z_{\mathrm{max}}$. |
| water table at the base, $z_{\mathrm{w},0}=z_{\mathrm{max}}$ | Zero initial water-column height above the base; pressure head is zero there. | Pressure head is zero at the reporting limit. |
| Requested water table below the base, $z_{\mathrm{w},0}>z_{\mathrm{max}}$ | Set the effective initial water table at the base, with zero initial water-column height. Retained water can then build a perched saturated region upward. | Retain the deeper water table. The analytical domain continues below the reported column, assuming the same homogeneous material. |

For a finite base, the effective initial depth is

$$
z_{\mathrm{w},0}=\min\left(z_{\mathrm{w},0}^{\mathrm{requested}},z_{\mathrm{max}}\right).
$$

Here $z_{\mathrm{w},0}^{\mathrm{requested}}$ denotes the input depth, while $z_{\mathrm{w},0}$ denotes the effective reference used in subsequent equations. Requested and effective depths are recorded separately. Moving a deeper requested water table to the restrictive base is a modelling assumption, not a prediction of its measured regional position.[^finite-base]

The current flow solution ends at this effective initial water table, where conductivity is fixed at its saturated value. The basal capacity represents removal of water reaching that level through the lower domain. The saturated region between the water table and a finite base is reconstructed using the linear pressure assumption below; transient saturated diffusion is not solved. Thus the finite base restricts drainage through a capacity, rather than imposing a new exact zero-flux Richards boundary.

As retained water changes, the estimated depth is

$$
z_{\mathrm{w}}(t)=z_{\mathrm{w},0}-\Delta h(t),
$$

where $\Delta h$ is rise above the effective initial water table. The observed plane remains fixed. Above the current water table it uses the unsaturated reconstruction; at the water table pressure head is zero, and below it the saturated reconstruction applies. A plane can become saturated during rainfall, but saturation alone does not imply failure: strength and geometry also enter FS.

## Unsaturated flow response

<span id="governing-equation-and-constitutive-relations"></span>

### Richards governing equation {#richards-equation}

With vertical depth $z$ positive downward, the one-dimensional Richards equation for a slope of angle $\delta$ is[^note-2]

$$
\frac{\partial\theta}{\partial t}
=\frac{\partial}{\partial z}\left[
K(\psi)\left(\frac{1}{\cos^2\delta}\frac{\partial\psi}{\partial z}-1\right)
\right].
$$

The publications use $Z$ for vertical depth and $z$ for the slope-normal coordinate; here the vertical coordinate is renamed $z$. The $1/\cos^2\delta$ factor accounts for that coordinate transformation. For horizontal ground, it reduces to the downward-positive form of Eq. (4a) in the 2010 article. The derivative of water content represents local storage, and the right-hand side represents the divergence of the pressure- and gravity-driven flow.

### Gardner laws and slope coefficient {#gardner-constitutive-laws-and-parameters}

The Gardner constitutive laws used to derive the analytical solution are[^note-3]

$$
K(\psi)=K_{\mathrm{s}}\exp(\alpha\psi^*),
\qquad
\theta(\psi)=\theta_{\mathrm{r}}+
(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})\exp(\alpha\psi^*),
\qquad \psi^*=\psi-\psi_0.
$$

These algebraic laws apply in the unsaturated domain; conductivity and water content are capped at $K_{\mathrm{s}}$ and $\theta_{\mathrm{s}}$ in the saturated domain. The reference allows a shifted pressure origin $\psi_0$ at the top of the capillary fringe. The current implementation uses the default value $\psi_0=0$, with the fixed analytical domain ending at the initial water table $z_{\mathrm{w},0}$. Support for a nonzero capillary-fringe pressure offset and the associated capillary-fringe treatment is planned for a future extension.

Here $\alpha$ is the intrinsic Gardner parameter of the soil, with units of m⁻¹. It controls how rapidly conductivity and water content vary with pressure head in the exponential constitutive laws. It is treated as constant within each homogeneous column.

For flow expressed in vertical coordinates on a slope, we define the derived coefficient $a$:

$$
a = \alpha \cos^2\delta.
$$

The angle $\delta$ is the local ground-surface slope. The coefficient $a$ also has units of m⁻¹ and combines the soil parameter with the slope geometry. It is not an additional fitted soil parameter or an independent input. For horizontal ground, $\delta=0$ and $a=\alpha$; for an inclined surface, $a<\alpha$. Both are constants for a given column, whereas $\psi(z,t)$, $K(z,t)$, and $\theta(z,t)$ are evolving fields. Introducing $a$ makes the $1/(\alpha\cos^2\delta)$ coefficient in the transformed flow equation read simply $1/a$.[^note-4]

### Flux and water conservation

The flux $q(z,t)$ represents volumetric flow per unit area, with units of m/s, and positive values denote downward flow. It is not total discharge in m³/s or the velocity of water within the pores. Water conservation requires

$$
\frac{\partial\theta(z,t)}{\partial t}
=-\frac{\partial q(z,t)}{\partial z}.
$$

This balance relates local changes in water content to the difference between incoming and outgoing flux. Recall the vertical-coordinate Richards equation introduced above:

$$
\frac{\partial\theta(z,t)}{\partial t}
=\frac{\partial}{\partial z}\left[
K\bigl(\psi(z,t)\bigr)
\left(\frac{1}{\cos^2\delta}\frac{\partial\psi(z,t)}{\partial z}-1\right)
\right].
$$

Comparing its right-hand side with $-\partial q/\partial z$, the expression inside the brackets is the negative of the downward-positive Darcy flux:

$$
-q(z,t)=K\bigl(\psi(z,t)\bigr)
\left[\frac{1}{\cos^2\delta}\frac{\partial\psi(z,t)}{\partial z}-1\right].
$$

The flux is identified from Darcy's law as used in Richards' equation; conservation alone specifies its spatial derivative. Multiplying the bracketed expression by $-1$ therefore gives

$$
q(z,t)=K\bigl(\psi(z,t)\bigr)
\left[1-\frac{1}{\cos^2\delta}
\frac{\partial\psi(z,t)}{\partial z}\right].
$$

### Gravity in Darcy's law {#gravitational-contribution}

The origin of the term $1$ is easiest to see for horizontal ground ($\delta=0$). With depth positive downward, elevation head $e(z)$ decreases by one metre for every metre of depth:

$$
e(z)=e_0-z,
\qquad \frac{\partial e}{\partial z}=-1.
$$

Here $e_0$ is the elevation head at the ground surface, measured relative to an arbitrary datum. Total hydraulic head is pressure head plus elevation head, $H(z,t)=\psi(z,t)+e(z)$. Darcy's law therefore gives

$$
q(z,t)=-K\bigl(\psi(z,t)\bigr)\frac{\partial H(z,t)}{\partial z}
=-K\bigl(\psi(z,t)\bigr)
\left[\frac{\partial\psi(z,t)}{\partial z}-1\right]
=K\bigl(\psi(z,t)\bigr)
\left[1-\frac{\partial\psi(z,t)}{\partial z}\right].
$$

Thus the positive $1$ comes from the negative elevation gradient combined with the negative sign in Darcy's law. It represents the gravitational contribution, while the pressure-head gradient can reinforce or oppose it. Gravity does not require a nonzero downward flux: for horizontal hydrostatic equilibrium, $\partial\psi/\partial z=1$ and $q=0$.

For inclined ground, the transformed vertical-coordinate flux above retains the term $1$, and the pressure-gradient term includes $1/\cos^2\delta$. This factor accounts for the slope-coordinate transformation in Eqs. (4b)–(4c) of [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321).

### Flux in conductivity form

Differentiating the Gardner conductivity law gives

$$
\frac{\partial K(z,t)}{\partial z}
=\alpha K(z,t)\frac{\partial\psi(z,t)}{\partial z},
\qquad
K(z,t)\frac{\partial\psi(z,t)}{\partial z}
=\frac{1}{\alpha}\frac{\partial K(z,t)}{\partial z}.
$$

Substituting this identity into the slope-adjusted vertical-coordinate Darcy flux derived under [flux and water conservation](#flux-and-water-conservation), $q=K[1-(1/\cos^2\delta)\,\partial\psi/\partial z]$, and then using $a=\alpha\cos^2\delta$ yields

$$
q(z,t)=K(z,t)-\frac{1}{\alpha\cos^2\delta}\frac{\partial K(z,t)}{\partial z}
=K(z,t)-\frac{1}{a}\frac{\partial K(z,t)}{\partial z}.
$$

Thus flux is not generally equal to conductivity: it also depends on the conductivity gradient. As a consistency check, the [steady initial conductivity profile $K_0(z)$](#initial-condition-and-reconstructed-profiles), defined by the background infiltration before the simulated rainfall history, satisfies $\partial K_0/\partial z=a(K_0-I_{\mathrm{ZLT}})$. Substitution into the flux expression gives $q(z,0)=I_{\mathrm{ZLT}}$ at every depth.

### Linearization in conductivity

To derive the conductivity equation, assume that $K_{\mathrm{s}}$, $\alpha$, $\theta_{\mathrm{s}}$, $\theta_{\mathrm{r}}$, $\delta$, and $\psi_0$ are constant. Pressure head is a field $\psi=\psi(z,t)$; conductivity and water content depend on it through the Gardner constitutive laws:

$$
K(z,t)\equiv K\bigl(\psi(z,t)\bigr),
\qquad
\theta(z,t)\equiv\theta\bigl(\psi(z,t)\bigr).
$$

In the derivation below, $K$, $\theta$, and $\psi$ abbreviate these fields at $(z,t)$; their arguments are omitted only for readability. Thus $\partial K/\partial z$ and $\partial K/\partial t$ are derivatives of the composite field $K(\psi(z,t))$, not of a constant. Both Gardner laws contain the same exponential, so eliminating it gives

$$
\exp[\alpha(\psi-\psi_0)]=\frac{K}{K_{\mathrm{s}}},
\qquad
\theta=\theta_{\mathrm{r}}
+\frac{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}{K_{\mathrm{s}}}K.
$$

Differentiating with respect to time yields

$$
\frac{\partial\theta}{\partial t}
=\frac{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}{K_{\mathrm{s}}}
\frac{\partial K}{\partial t}.
$$

The same spatial derivative used in the flux derivation is

$$
\frac{\partial K}{\partial z}
=\alpha K\frac{\partial\psi}{\partial z},
\qquad
K\frac{\partial\psi}{\partial z}
=\frac{1}{\alpha}\frac{\partial K}{\partial z}.
$$

Substitution into Richards' equation therefore gives

$$
\frac{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}{K_{\mathrm{s}}}
\frac{\partial K}{\partial t}
=\frac{\partial}{\partial z}\left[
\frac{1}{\alpha\cos^2\delta}\frac{\partial K}{\partial z}-K
\right]
=\frac{1}{a}\frac{\partial^2K}{\partial z^2}
-\frac{\partial K}{\partial z}.
$$

Multiplying by $a=\alpha\cos^2\delta$ gives the linear conductivity equation[^note-5]

$$
\frac{a(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})}{K_{\mathrm{s}}}
\frac{\partial K}{\partial t}
=\frac{\partial^2K}{\partial z^2}-a\frac{\partial K}{\partial z},
\qquad 0<z<z_{\mathrm{w},0}.
$$

The coefficient of $\partial^2K/\partial z^2$ when isolated on the right is $K_{\mathrm{s}}/[a(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})]$, with units of diffusivity (m²/s).

This equation is linear in $K(z,t)$: the field and its derivatives occur to the first power, with constant coefficients, and are not multiplied by one another. The algebraic transformation is exact under the adopted Gardner laws; it is not a small-pressure-change approximation. The constitutive laws themselves approximate the soil's hydraulic behavior. Linearity of this flow equation does not imply that the subsequent storage, moving-water-table, and stability calculations are linear.

### Boundary conditions

The flow solution needs a condition at each end of the soil column. At the ground surface, the incoming flux is the effective infiltration $I(t)$. At the initial water-table depth, conductivity is fixed at its saturated value:

$$
q(0,t)=I(t),
\qquad K(z_{\mathrm{w},0},t)=K_{\mathrm{s}}.
$$

The surface condition sets how much water enters the column. Here $I(t)=\min[P(t),K_{\mathrm{s}}]$, so any rainfall above the local conductivity limit is discarded. Rainfall intensity can change from one interval to the next. The column starts from the [steady initial conductivity profile](#initial-condition-and-reconstructed-profiles).

The lower condition is imposed at the effective initial water table, which need not coincide with the maximum reported depth. Water can still flow through that boundary. The flow solution calculates how much water reaches it; the selected basal drainage model determines how much can drain away. Water retained after the storage calculation can raise the water table, as explained under [water-table rise](#water-table-rise).[^note-6]

### Steady initial conductivity profile {#initial-condition-and-reconstructed-profiles}

Before the simulated rainfall history begins, the column is in a steady state maintained by $I_{\mathrm{ZLT}}$. The initial conductivity profile follows by imposing steady state and the boundary conditions on the flow equation derived above.

In steady state, water content does not change with time. By conservation,

$$
\frac{\partial\theta}{\partial t}
=-\frac{\partial q}{\partial z}=0.
$$

Flux is therefore constant with depth. Since the background infiltration $I_{\mathrm{ZLT}}$ enters at the surface,

$$
q(z,0)=I_{\mathrm{ZLT}}.
$$

Each part of the column receives and releases the same amount of water, without accumulation.

Substituting this flux into Darcy's expression in terms of $K$ gives a first-order differential equation for the initial conductivity profile $K_0(z)$. Subtracting the constant background infiltration and solving leaves an exponential dependence on depth, with one integration constant.

The lower boundary condition sets conductivity equal to $K_{\mathrm{s}}$ at the initial water table. This determines the integration constant and gives[^note-7]

$$
K_0(z) = I_{\mathrm{ZLT}} - (I_{\mathrm{ZLT}}-K_{\mathrm{s}})\exp[-a(z_{\mathrm{w},0}-z)], \qquad 0 \le z \le z_{\mathrm{w},0}.
$$

Steady state means that the profile remains unchanged with time, not that conductivity is uniform with depth. Its depth variation allows it to transmit the same background flux throughout the column.

### Pressure and water content: equation comparison {#pressure-and-water-content-equation-sets}

After calculating $K(z,t)$, the default set obtains pressure head by directly inverting the intrinsic Gardner law:[^note-8]

$$
\psi(z,t)=\frac{1}{\alpha}\ln\left(\frac{K(z,t)}{K_{\mathrm{s}}}\right)
=\frac{\cos^2\delta}{a}\ln\left(\frac{K(z,t)}{K_{\mathrm{s}}}\right).
$$

Srivastava and Yeh use an upward-positive vertical coordinate. Reversing its direction changes spatial derivative signs, not the scalar pressure-to-conductivity relation. Slope enters the flow equation through $a=\alpha\cos^2\delta$; intrinsic $\alpha$ remains in the material law.

Water content then follows the same constitutive relation:[^note-9]

$$
\theta(z,t)=\theta_{\mathrm{r}}+(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})\exp[\alpha\psi(z,t)]
=\theta_{\mathrm{r}}+(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})\frac{K(z,t)}{K_{\mathrm{s}}}.
$$

The identical exponential in both laws gives this linear relation between water content and conductivity. In the table below, effective saturation means $(\theta-\theta_{\mathrm{r}})/(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})$: zero at residual water content and one at saturation. Subscripts $\mathrm{B}$ and $\mathrm{SY}$ identify the Baum et al. and Srivastava–Yeh sets.

The published TRIGRS pressure conversion and its Fortran retention reconstruction differ from these identities when there is a slope. For compact notation, define the relative hydraulic conductivity $\rho(z,t)=K(z,t)/K_{\mathrm{s}}$. It is dimensionless: $\rho=1$ means saturated conductivity, while $\rho=0.1$ means conductivity is one tenth of the saturated value. This ratio is calculated from $K$; it is not an additional soil parameter. The table omits $(z,t)$ from $\rho$ for readability:

| Quantity | TRIGRS / Baum et al. set | Gardner / Srivastava–Yeh set |
| --- | --- | --- |
| Pressure head | $\psi_{\mathrm{B}}=(\cos\delta/a)\ln\rho$ | $\psi_{\mathrm{SY}}=(1/\alpha)\ln\rho=(\cos^2\delta/a)\ln\rho$ |
| Reconstructed water content | $\theta_{\mathrm{B}}=\theta_{\mathrm{r}}+(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})\exp(a\psi_{\mathrm{B}})$ | $\theta_{\mathrm{SY}}=\theta_{\mathrm{r}}+(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})\exp(\alpha\psi_{\mathrm{SY}})$ |
| Effective saturation expressed through $\rho$ | $\rho^{\cos\delta}$ | $\rho$ |

The first TRIGRS entry is Eq. (9) in Baum et al. (2008, 2010); the second reproduces the slope-adjusted retention exponent $a\psi_{\mathrm{B}}$ used by the released TRIGRS calculation. The latter is different from the intrinsic $\alpha$ written in the article's Eq. (4e). The 2013 correction does not amend Eq. (9).

For the same conductivity ratio and soil parameter, substituting the TRIGRS conversion back into the Gardner conductivity law $K(\psi)=K_{\mathrm{s}}\exp(\alpha\psi)$, with the intrinsic soil parameter $\alpha$ and $\psi_0=0$, gives $K_{\mathrm{s}}\rho^{1/\cos\delta}$ instead of the prescribed $K_{\mathrm{s}}\rho$. This is a constitutive inconsistency under the definitions used in the governing-equation derivation. A change of spatial coordinates alone does not rescale scalar pressure head. The default set preserves the constitutive identities; the Baum et al. set uses the conversions from the TRIGRS comparisons.

At $\delta=0$, both pressure and retention expressions coincide. For $\delta=40^\circ$, $\alpha=3$ m⁻¹, and $\rho=0.1$, pressure is about −1.002 m in the TRIGRS set and −0.768 m in the default set. With $\theta_{\mathrm{s}}=0.41$ and $\theta_{\mathrm{r}}=0.06$, water content is about 0.120 and 0.095, respectively. These differences propagate into fillable pore volume, water-table rise, and FS; their combined effect on stability is not necessarily monotonic.

A [controlled comparison notebook](verification/01_baum_etal_vs_srivastava_yeh.ipynb) isolates these conversions and shows their combined effects on reconstructed water-table rise and stability. It also examines the separate short-time response approximation.

### Analytical response to changing rainfall {#response-to-changing-rainfall}

Rainfall is supplied as successive intervals of constant intensity. Each column first converts rainfall to effective infiltration, $I(t)=\min[P(t),K_{\mathrm{s}}]$. The surface input can change abruptly at an interval boundary, but the soil's conductivity, moisture, and pressure respond progressively. Water entering at the surface does not immediately reach the water table.

To see how the history is represented, consider a column with zero background infiltration and a conductivity limit of at least 10 mm/h:

| Time | Effective infiltration |
| --- | ---: |
| From 0 to 2 h | 6 mm/h |
| From 2 to 5 h | 10 mm/h |
| After 5 h | 0 mm/h |

This is equivalent to three changes: $\Delta I_0=+6$ mm/h at 0 h, $\Delta I_1=+4$ mm/h at 2 h, and $\Delta I_2=-10$ mm/h at 5 h. The negative change stops the surface input; it does not instantly remove water from the column.

The first change is measured relative to the background rate, $\Delta I_0=I(0)-I_{\mathrm{ZLT}}$; subsequent changes are differences between successive imposed rates. Let $R_K(z,\tau)$ be the conductivity response per unit change in infiltration, after an elapsed time $\tau=t-t_j$. It is zero before the change occurs. Because the conductivity equation is linear, the response to the whole history is

$$
K(z,t)=K_0(z)+\sum_j\Delta I_j R_K(z,t-t_j),
$$

where $t_j$ is the time of each change. For the example above, this becomes

$$
K(z,t)=K_0(z)
+\Delta I_0 R_K(z,t)
+\Delta I_1 R_K(z,t-2\,\mathrm{h})
+\Delta I_2 R_K(z,t-5\,\mathrm{h}).
$$

Use consistent units for the rates, conductivity, and elapsed times. At 6 h, the three changes have been acting for 6, 4, and 1 h, respectively. Their effects therefore do not immediately cancel even though infiltration is now zero. The analytical column retains the effects of previous rainfall and can continue transferring water to its base.

This addition of individual responses is called *superposition*. It applies to the analytical conductivity and flux calculations. Pressure and factor of safety are calculated afterward using their own relationships; their separate responses cannot simply be added in the same way.

The initial profile $K_0$ is included once, and subsequent changes are added to it. This preserves the fixed lower condition $K(z_{\mathrm{w},0},t)=K_{\mathrm{s}}$ when rainfall ends. When using the literal finite-interval sums in the published Eq. (8), the forcing history must be continued with a dry interval before calculating later dry-period results. Extrapolating a truncated pulse sum can cancel the lower conductivity. This is a forcing-history convention, not evidence of an error in the paper's sustained-rate examples.

### Characteristic roots and response time

The response to one infiltration change is built from a series of spatial patterns. Adding these patterns describes how conductivity changes over the whole depth. Each pattern has a different shape and decays at a different rate. These are mathematical components of the solution, not separate portions of water or physical waves.

The two boundary conditions determine which patterns are allowed: imposed flux at the surface and fixed conductivity at the initial water table. Each allowed pattern is associated with a positive characteristic root $r_m$ satisfying

$$
\sin(az_{\mathrm{w},0}r_m)+2r_m\cos(az_{\mathrm{w},0}r_m)=0.
$$

The roots $r_m$ are dimensionless; the product $az_{\mathrm{w},0}$ is also dimensionless. The index $m$ numbers the roots. Their values are calculated from $a$ and the domain depth; they are not additional soil properties supplied by the user.[^note-10]

For each infiltration change, elapsed time is expressed in dimensionless form as

$$
t^*=\frac{aK_{\mathrm{s}}\tau}{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}},
\qquad \tau=t-t_j.
$$

This measures elapsed time relative to the response scale of the soil and slope. The transient components contain decay factors of the form

$$
\exp\left[-\left(r_m^2+\frac14\right)t^*\right].
$$

Components with larger roots decay quickly; those with smaller roots persist longer. If the imposed infiltration remains constant, the transient components gradually disappear and the profile approaches the new steady state. Each later rainfall change starts its own response with its own elapsed time.

### Conductivity response to an infiltration change

The roots enter directly into the series that calculates $R_K(z,\tau)$. They do not specify depths or times at which the model is evaluated; each root identifies a spatial pattern and its decay rate. Define an auxiliary series function that groups the transient terms:

$$
\mathcal{S}_K(z,\tau)=\sum_{m=1}^{\infty}
\frac{\sin[ar_m(z_{\mathrm{w},0}-z)]\sin(ar_mz_{\mathrm{w},0})}
{1+az_{\mathrm{w},0}/2+2az_{\mathrm{w},0}r_m^2}
\exp\left[-\left(r_m^2+\frac14\right)t^*\right],
\qquad \tau>0,
$$

where $t^*=aK_{\mathrm{s}}\tau/(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})$. The response per unit change in infiltration is then

$$
R_K(z,\tau)=1-\exp[-a(z_{\mathrm{w},0}-z)]
-4\exp\left(\frac{az}{2}\right)\mathcal{S}_K(z,\tau),
\qquad \tau>0.
$$

$\mathcal{S}_K$ is an auxiliary dimensionless function, not an additional soil parameter or the retained-water storage $S(t)$. The response $R_K$ is defined as zero for $\tau\le0$.

The structure of $R_K$ is a steady response plus a transient correction. The term $1-\exp[-a(z_{\mathrm{w},0}-z)]$ gives the final conductivity change per unit infiltration change. The term $-4\exp(az/2)\mathcal{S}_K(z,\tau)$ adjusts that final response so that the column begins at its original profile and evolves gradually. At the instant of the change, the two parts cancel; at long times, the transient correction disappears and only the steady response remains.[^note-11][^note-12]

The factor $\sin[ar_m(z_{\mathrm{w},0}-z)]$ describes variation with depth. The exponential in $t^*$ describes decay with elapsed time, while the remaining factors determine each component's amplitude. Contributions can be positive or negative: their sum represents the conductivity response, not separate amounts of water. The term $1-\exp[-a(z_{\mathrm{w},0}-z)]$ is the response that remains after the transient components have decayed.

The calculation first obtains the roots from $a$ and $z_{\mathrm{w},0}$. It then evaluates $\mathcal{S}_K$ at the desired depth and elapsed time to calculate $R_K$. For a single infiltration change $\Delta I$, the resulting conductivity is

$$
K(z,t)=K_0(z)+\Delta I\,R_K(z,\tau).
$$

For several changes, the program adds these responses using the respective start times, as described above. $R_K$ is dimensionless because it relates a conductivity change to a flux change; multiplying it by $\Delta I$ gives a conductivity increment in the same units as $K$.

At $z=z_{\mathrm{w},0}$, both the steady response term and every sine term vanish, giving $R_K(z_{\mathrm{w},0},\tau)=0$. The infiltration changes therefore preserve $K(z_{\mathrm{w},0},t)=K_{\mathrm{s}}$. At long times the transient terms vanish, giving the new steady profile for the imposed flux. At very short times, the infinite series must cancel the steady response term to approach zero; this is why numerical evaluation needs particular care near a change in rainfall.

### Evaluating the response at short times

The mathematical series has infinitely many terms; the program evaluates a finite number. Near a sudden infiltration change, more terms can be needed to represent the response accurately. Later, many components have already decayed and contribute little. The number of retained series terms controls this approximation; additional terms are evaluated when needed for convergence in the default formulation. Depth sampling and storage integration also affect accuracy; see [numerical controls](inputs-outputs.md#numerical-controls).

At extremely short elapsed times, the change is concentrated near the surface and the distant lower boundary has negligible influence on that change. The default implementation then evaluates a *half-space* response: a short-time solution in which the lower boundary is treated as sufficiently far away. The physical column and its initial water table remain finite.

The surface still imposes infiltration through

$$
K(0,t)-\frac{1}{a}\left.\frac{\partial K}{\partial z}\right|_{z=0}=I(t).
$$

This is called a *Robin boundary condition* because it combines the value of conductivity and its gradient. It specifies surface flux rather than prescribing conductivity alone.

The short-time comparison concerns $\Delta K(0,t)$, the surface-conductivity change caused by an infiltration change $\Delta I$. For a fixed column, the leading behavior as elapsed time tends to zero is

| Early surface response | TRIGRS / Baum et al. approximation | Consistent linear flow problem |
| --- | --- | --- |
| $\Delta K(0,t)/\Delta I$ | $\tfrac12\exp(az_{\mathrm{w},0}/2)\sqrt{t^*/\pi}$ | $2\sqrt{t^*/\pi}$ |

In both expressions the conductivity change tends to zero at the instant the surface flux changes. The consistent expression describes a gradual initial change proportional to the square root of elapsed time. The older approximation includes an additional length-dependent factor, so it can overestimate or underestimate that initial response.

The Baum et al. formulation retains the empirical factor-of-four reduction used in the released TRIGRS early-time approximation. The default set instead uses the resolved series or the short-time half-space response. Tests check its imposed surface flux and agreement with an independent numerical solution. This numerical issue is separate from the pressure/retention conversion differences discussed above.[^note-13]

## Water-table rise

### Basal transfer and drainage capacity

The analytical response gives the flux arriving at the initial water-table depth, $q_{\mathrm{b}}(t)=q(z_{\mathrm{w},0},t)$. The selected basal boundary provides an effective capacity $q_{\mathrm{c}}$ for draining this water. The current model defines

$$
q_{\mathrm{c}}=\eta K_{\mathrm{s}}\sin^2\delta.
$$

Here $\eta\in[0,1]$ is a parameter defining the drainage efficiency at the basal boundary. Efficiency multiplies a rate capacity, not a percentage of rainfall. This external capacity assumption is the same in both equation sets, so changing the constitutive conversion does not silently change the drainage setting. For example, on a 40° slope with $K_{\mathrm{s}}=36$ mm/h, efficiency 1 gives a capacity of about 14.87 mm/h; efficiency 0.015 gives about 0.2231 mm/h.

The efficient-drainage end member uses $\eta=1$, while the restricted-drainage end member uses $\eta=0.1$. Either can accumulate water when basal transfer is sufficient to give a positive storage rate, as defined below; efficiency 1 does not impose a fixed water table. A completely blocked drainage capacity corresponds to $\eta=0$.[^note-14]

### Retained water

The part of basal transfer above the drainage capacity is $\max(q_{\mathrm{b}}-q_{\mathrm{c}},0)$. The implemented storage rate also subtracts the background contribution:

$$
g(t)=\max(q_{\mathrm{b}}(t)-q_{\mathrm{c}},0)-I_{\mathrm{ZLT}}.
$$

Here $g(t)$ is the storage rate before applying the zero floor, in m/s. With zero background infiltration, storage grows whenever basal transfer exceeds capacity. With nonzero background, positive growth requires $q_{\mathrm{b}}>q_{\mathrm{c}}+I_{\mathrm{ZLT}}$. A negative $g$ reduces existing storage until it reaches zero; it does not create negative storage. Background flow alone, $q_{\mathrm{b}}=I_{\mathrm{ZLT}}$, therefore cannot raise the water table in this storage convention, even when it exceeds $q_{\mathrm{c}}$. The background subtraction treats storage as a change from the reference state, rather than a complete inlet-minus-outlet groundwater balance.

The calculation integrates this rate over successive time intervals using Simpson quadrature. Retained water $S_i$ is updated as

$$
S_i=\max\left[S_{i-1}+\int_{t_{i-1}}^{t_i}g(\tau)\,d\tau,0\right],
\qquad S_0=0.
$$

In this integral, $\tau$ is a time integration variable, not the elapsed time since a particular rainfall change used in $R_K$. The zero floor prevents a dry period from leaving a negative storage debt that would absorb the next storm. Storage starts at zero at the beginning of the complete history; an antecedent simulation can carry positive storage into the triggering event.[^note-15]

### From retained water to water-table depth

At a depth $z$, the pore volume still available to fill between that depth and the initial water table is calculated from the water-content profile before applying the new saturated region. Thus $\theta$ in the following integral is the analytical water-content profile, not the final profile after the water-table reconstruction:

$$
V_{\mathrm{f}}(z,t)=\int_z^{z_{\mathrm{w},0}}
[\theta_{\mathrm{s}}-\theta(\xi,t)]\,d\xi.
$$

Here $\xi$ is the integration variable for vertical depth, in metres. It runs from the candidate depth $z$ to the initial water-table depth $z_{\mathrm{w},0}$, while time $t$ is held fixed. It denotes the same depth coordinate as $z$, not a separate physical coordinate or soil parameter; $d\xi$ is an infinitesimal soil thickness.

The difference $\theta_{\mathrm{s}}-\theta$ is the local water deficit. Integrating it over depth gives a water volume per unit area, expressed in metres of water. This is why $S$ and $V_{\mathrm{f}}$ can be compared even though neither is a geometric water-table depth.

The calculation finds the depth where the available pore volume matches retained water, $V_{\mathrm{f}}(z_{\mathrm{w}}(t),t)=S(t)$, by interpolation between sampled depths. It then reports $\Delta h(t)=z_{\mathrm{w},0}-z_{\mathrm{w}}(t)$. If no water is retained, the water table remains at its initial depth. If retained water fills the whole sampled column, the estimated water table is capped at the ground surface.[^note-16]

This is a storage-based approximation, and the retention conversion described above affects it.

If retained water exceeds the total fillable pore volume, the surface limit constrains the estimated water-table position but does not automatically remove that excess storage. This differs from rainfall rejected at the surface. The model has no separate saturated-flow or aquifer-recession calculation to dispose of all retained water; strong, prolonged input with restricted drainage can exceed the approximation's physical range. These limitations are examined in the [antecedent water-balance checks](verification.md#antecedent-water-balance-checks).

## Basal and pressure assumptions

After estimating $z_{\mathrm{w}}(t)$, the package reconstructs the column profile. Above that depth it keeps the analytical unsaturated pressure and water-content profiles. At and below it, water content is set to $\theta_{\mathrm{s}}$, conductivity to $K_{\mathrm{s}}$, and pressure head to

$$
\psi(z,t)=\beta[z-z_{\mathrm{w}}(t)],
\qquad z\ge z_{\mathrm{w}}(t).
$$

Here $\beta$ is the slope of the reconstructed pressure-head profile below the water table: it gives the increase in pressure head per metre of additional depth. It is dimensionless (metres of pressure head per metre of soil depth) and is calculated from the ground-surface slope and background infiltration, not supplied as an independent soil property. It sets the pressure gradient, not the amount of water-table rise. The selected background gradient differs between sets. Let $r_{\mathrm{b}}=I_{\mathrm{ZLT}}/K_{\mathrm{s}}$ be the dimensionless background infiltration relative to saturated conductivity:

| Background quantity | TRIGRS / Baum et al. set | Gradient consistent with the stated vertical Darcy law |
| --- | --- | --- |
| $\beta$ | $\cos^2\delta-r_{\mathrm{b}}$ | $\cos^2\delta(1-r_{\mathrm{b}})$ |
| Darcy flux for saturated $K$ | $I_{\mathrm{ZLT}}/\cos^2\delta$ | $I_{\mathrm{ZLT}}$ |

The right-hand gradient is derived by setting $q=I_{\mathrm{ZLT}}$ in $q=K_{\mathrm{s}}[1-(\partial\psi/\partial z)/\cos^2\delta]$. It also matches the derivative of the intrinsic-Gardner steady pressure immediately above the initial water table. The left-hand gradient is the Iverson/TRIGRS saturated-pressure assumption from Eq. (1), not a new parameter from the 1991 paper. The mismatch appears when combining that assumption with the unsaturated transformed Darcy equation; the comparison does not establish that every original saturated TRIGRS formulation is erroneous.

At the water table the reconstruction gives zero pressure head. For a positive gradient, head increases with depth below it; a rising water table increases pressure at a fixed failure plane. In the Baum et al. set, some background inputs can make the gradient negative; the default set remains nonnegative for the accepted range $0\le I_{\mathrm{ZLT}}\le K_{\mathrm{s}}$.

Recall that $q_{\mathrm{c}}=\eta K_{\mathrm{s}}\sin^2\delta$ is prescribed separately. The selected gradient determines the reconstructed pressure, while drainage efficiency controls the amount of basal water retained. The complete output profile includes the saturated region from the current water table down to $z_{\mathrm{max}}$, whether or not an observation plane is selected.[^note-17]

The pressure reconstruction does not include transient saturated-pressure diffusion. It uses the linear approximation described above.

The flow calculation continues on the initial analytical domain. Updating the water table and replacing the saturated part of the profile does not recompute a fully coupled moving-boundary flow solution. These assumptions can create abrupt joins in profile plots; they should not be interpreted as a resolved wetting front.

## Infinite-slope stability

The infinite-slope model compares shear resistance with the downslope driving stress on a plane parallel to the ground surface. Its factor of safety is

$$
\mathrm{FS}(z,t)=\frac{\tan\phi'}{\tan\delta}+
\frac{c'-\chi\gamma_{\mathrm{w}}\psi(z,t)\tan\phi'}
{\gamma z\sin\delta\cos\delta}.
$$

Here $c'$ is effective cohesion (Pa), $\phi'$ is the effective friction angle, $\gamma$ is soil unit weight (N/m³), and $\gamma_{\mathrm{w}}$ is water unit weight (N/m³). Multiplying pressure head by $\gamma_{\mathrm{w}}$ gives pore pressure in Pa. The first term is the friction contribution; the second includes cohesion and the effect of pore pressure.

Positive pressure reduces the effective stress available for frictional resistance. Negative pressure represents suction and contributes to resistance through the dimensionless factor $\chi$. The implemented suction factor is

$$
\chi=\begin{cases}
\min\left[1,\max\left(0,\dfrac{\theta-\theta_{\mathrm{r}}}{\theta_{\mathrm{s}}-\theta_{\mathrm{r}}}\right)\right], & \psi<0,\\
1, & \psi\ge0.
\end{cases}
$$

The water-content ratio measures how wet the soil is between its residual and saturated limits. It bounds the contribution of suction; it is not a probability of failure.[^note-18]

The complete stability profile is evaluated over $0\le z\le z_{\mathrm{max}}$. The selected-plane result is $\mathrm{FS}(z_{\mathrm{f}},t)$, with $0<z_{\mathrm{f}}\le z_{\mathrm{max}}$; by default the plane is the column base. A value above 1 indicates resistance greater than driving stress in this model; a value at or below 1 reaches its deterministic failure threshold. This is distinct from the minimum FS over all sampled depths. It does not establish that a landslide will occur or quantify its probability. Large values are capped for output, with a special near-surface convention; see [numerical controls](inputs-outputs.md#numerical-controls).

## Finite antecedent infiltration history

### Background state and previous rainfall

Antecedent rain is rainfall before the triggering event. The package simulates that previous history, then continues into the event without resetting moisture, conductivity, or retained water. The same analytical response is used for both stages.

The background infiltration $I_{\mathrm{ZLT}}$ defines the steady state **before** the antecedent window. The antecedent rate is a separate forcing. If only amount and duration are known, a uniform-rainfall assumption gives constant intensity: 500 mm over 90 days gives $500/(90\times24)=0.23148$ mm/h. Different time distributions of that same amount can produce different responses. Every interval is capped at the local $K_{\mathrm{s}}$; excess rainfall is discarded. Evapotranspiration and interception are outside this model.

### Drainage during the two stages

The antecedent drainage efficiency $\eta_{\mathrm{a}}$ determines the capacity during the previous rainfall history. If drainage conditions are assumed unchanged, the same efficiency applies during the antecedent period and the triggering event. Different efficiencies may be assigned to the two stages as an explicit boundary assumption. The antecedent capacity is

$$
q_{\mathrm{c,a}}=\eta_{\mathrm{a}}K_{\mathrm{s}}\sin^2\delta.
$$

This is the same [capacity formula](#basal-transfer-and-drainage-capacity) with the antecedent setting. Assigning different efficiencies is a boundary assumption, not a property inferred from rainfall amount. Neither high efficiency nor efficiency 1 guarantees a fixed water table: accumulation depends on basal transfer, capacity, and the background correction.

### Event-relative results

The start of the triggering event may be used as the time origin for interpreting results. The reference water-table depth $z_{\mathrm{w},0}$ still describes the beginning of the complete history, so retained storage and water-table rise can already be positive when the event begins. Initial FS can change because previous rain changes suction, water-table depth, or both. The response time factor remains $aK_{\mathrm{s}}/(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})$; there is no fitted multiplier that accelerates time according to wetness.

A measured rainfall history should retain its changes in intensity and its dry intervals. This distinguishes a known hyetograph from a uniform approximation based only on amount and duration. The [antecedent rainfall tutorial](tutorials/02_antecedent_wetness.ipynb) compares drainage assumptions, while [inputs and outputs](inputs-outputs.md#antecedent-rainfall) describes how to supply the histories.

## Equation provenance

The references below identify verified equation numbers. Letter case follows the respective publication: for example, the retention equation is **(4C)** in the 2008 report and **(4e)** in the 2010 article. A source number indicates the relevant original relationship, not automatic equivalence of every package assumption. Formula notes state the adaptations beside the equations. References to the 2010 formulas use [Baum and Godt (2013)](https://doi.org/10.1002/jgrf.20100): Eqs. (5c), (5d), (8), (12a), and (12b) replace the extraneous $K_{\mathrm{s}}\exp(a\psi_0)$ term with $K_{\mathrm{s}}$; Eq. (10a) requires $q_{Z\mathrm{max}}=-K_{\mathrm{s}}(\beta-1)$. The correction states that the original program was coded correctly and that the reported results were unaffected.

| Relationship used in unsatstab | Baum et al. (2010) | Baum et al. (2008) | Correspondence and adaptation |
| --- | --- | --- | --- |
| Local infiltration cap $I=\min(P,K_{\mathrm{s}})$ | Eqs. (16a)–(16b) | Eqs. (16A)–(16B) | Same local cap with upslope runoff set to zero; package discards excess rather than routing it. |
| Richards equation in vertical coordinates | Eq. (4c) | Eq. (4A) | Same governing form with the source's vertical coordinate $Z$ renamed $z$; the analytical conductivity solution is evaluated on the initial domain. |
| Gardner conductivity $K(\psi)$ | Eq. (4d) | Eq. (4B) | Intrinsic exponential law; source uses shifted head $\psi^*=\psi-\psi_0$. |
| Intrinsic Gardner retention $\theta(\psi)$ | Eq. (4e) | Eq. (4C) | Standalone retention function uses intrinsic $\alpha$ and $\psi_0=0$. The default set preserves intrinsic $\alpha$; the TRIGRS code reconstruction uses $a$ instead, as shown in the equation-set comparison. |
| Slope transform $a=\alpha\cos^2\delta$ | Definition after Eq. (4f) | Definition after Eq. (4C), before Eq. (5A) | Same parameter transform, with $a$ renamed from $\alpha_1$; the definition itself has no separate equation number. |
| Linearized conductivity equation | Eq. (4f) | Linearization described after Eq. (4C) | Same analytical flow equation with $a=\alpha_1$; reconstructed retention and moving water table are additional approximations. |
| Decay/time scale $aK_{\mathrm{s}}/(\theta_{\mathrm{s}}-\theta_{\mathrm{r}})$ | Eq. (5b) | Eq. (5B) | Same decay parameter with $a=\alpha_1$. |
| Initial conductivity $K_0(z)$ | Eq. (5c), corrected in 2013 | Eq. (5C) | Renamed depth and slope parameter, with $\psi_0=0$; corrected source removes the spurious exponential multiplier. |
| Constant-flux conductivity response | Eq. (5d), corrected in 2013 | Eq. (5D) | Corrected $K_{\mathrm{s}}$ term is used in the analytical response; flux increments are superposed for variable forcing. |
| Characteristic roots | Eq. (6) | Eq. (6) | Pole-free sine/cosine rewrite of the tangent equation; root variables and indices remain italic. |
| Time-dependent basal flux/conductivity | Eqs. (7)–(8); Eq. (8) corrected in 2013 | Eqs. (8A)–(8B) | Same analytical response family; package handles flux increments and reports total basal flux, with its specified post-event forcing. |
| Unsaturated pressure from conductivity | Eq. (9), TRIGRS conversion | Eq. (9), TRIGRS conversion | The default set instead inverts Srivastava–Yeh Eq. (2a), giving $1/\alpha=\cos^2\delta/a$; see the parallel comparison. |
| Background pressure gradient | Definition with Eq. (1) | Definition with Eq. (1) | The Baum et al. gradient is $\cos^2\delta-I_{\mathrm{ZLT}}/K_{\mathrm{s}}$; the default gradient is derived from the stated Darcy law and background flux. |
| Drainage capacity and basal excess | Eqs. (10a)–(10c); Eq. (10a) corrected in 2013 | Eqs. (10A)–(10C) | Package $\eta$ corresponds to source $c_{\mathrm{d}}$; total flux/background correction and stage-specific overrides are implementation adaptations. The external capacity remains the same in both sets; the reconstructed pressure gradient is selected separately. |
| Cumulative storage with a zero floor | Section 4.2.4, paragraph 29, following Eqs. (10b)–(10c) | Water Table Rise section following Eqs. (10B)–(10C) | Same numerical volume-integration concept; the interval-wise floor is a package addition, not a numbered source formula. |
| Fillable pore volume $V_{\mathrm{f}}$ | Eq. (11a); conductivity primitives (12a)–(12b) corrected in 2013 | Eq. (11A) | Same pore-volume concept, integrated numerically from the package retention profile rather than the published conductivity primitives. |
| Linear pressure below the moving water table | Eq. (3b); steady term of Eq. (1) | Eq. (3B); steady term of Eq. (1) | Linear reconstruction retained with a selectable background gradient; it remains an extension beyond the source's initially upward-flow case. Does not implement saturated diffusion Eqs. (13)–(14). |
| Infinite-slope $\mathrm{FS}$ and unsaturated $\chi$ | Eq. (18); section 6, paragraph 42 | Eq. (18); Slope Stability Computations discussion | Effective-stress adaptation with $\chi$ for suction. Package uses constant soil unit weight and numerical clipping. |
| Antecedent rate, matching efficiency, onset comparison, budget residual | Derived from the preceding relationships | No directly corresponding numbered equation | Package definitions and scenario choices; do not attribute a new equation number to Baum. |

## References

- Alvioli, M., & Baum, R. L. (2016). Parallelization of the TRIGRS model for rainfall-induced landslides using the message passing interface. *Environmental Modelling & Software*, 81, 122–135. [https://doi.org/10.1016/j.envsoft.2016.04.002](https://doi.org/10.1016/j.envsoft.2016.04.002)
- Baum, R. L., Savage, W. Z., & Godt, J. W. (2002). *TRIGRS—A Fortran Program for Transient Rainfall Infiltration and Grid-Based Regional Slope-Stability Analysis*. U.S. Geological Survey Open-File Report 02-424.
- Baum, R. L., Savage, W. Z., & Godt, J. W. (2008). *TRIGRS—A Fortran Program for Transient Rainfall Infiltration and Grid-Based Regional Slope-Stability Analysis, Version 2.0*. U.S. Geological Survey Open-File Report 2008-1159. [https://doi.org/10.3133/ofr20081159](https://doi.org/10.3133/ofr20081159)
- Baum, R. L., Godt, J. W., & Savage, W. Z. (2010). Estimating the timing and location of shallow rainfall-induced landslides using a model for transient, unsaturated infiltration. *Journal of Geophysical Research: Earth Surface*, 115, F03013. [https://doi.org/10.1029/2009JF001321](https://doi.org/10.1029/2009JF001321)
- Baum, R. L., & Godt, J. W. (2013). Correction to “Estimating the timing and location of shallow rainfall-induced landslides using a model for transient, unsaturated infiltration.” *Journal of Geophysical Research: Earth Surface*, 118, 1999. [https://doi.org/10.1002/jgrf.20100](https://doi.org/10.1002/jgrf.20100)
- Srivastava, R., & Yeh, T.-C. J. (1991). Analytical solutions for one-dimensional, transient infiltration toward the water table in homogeneous and layered soils. *Water Resources Research*, 27(5), 753–762. [https://doi.org/10.1029/90WR02772](https://doi.org/10.1029/90WR02772)

[^note-1]: **Source/adaptation:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eqs. (16a)–(16b), with upslope runoff set to zero; unlike TRIGRS, this package discards locally rejected rain without surface routing.

[^note-2]: **Source:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (4c), and [Baum et al. (2008)](https://doi.org/10.3133/ofr20081159), Eq. (4A).

[^note-3]: **Source:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eqs. (4d)–(4e), and [Baum et al. (2008)](https://doi.org/10.3133/ofr20081159), Eqs. (4B)–(4C).

[^note-4]: **Source:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), the definition of $\alpha_1$ immediately after Eq. (4f); here $a=\alpha_1$. This is a parameter definition, not Eq. (4f) itself.

[^note-5]: **Source:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (4f), with $a=\alpha_1$ and renamed vertical depth. The 2008 report describes this linearization after Eq. (4C), without assigning a separate equation number.

[^note-6]: **Source:** the surface-flux and lower-head conditions are described with Eqs. (5c)–(5d) in [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), paragraph 24, using the corrected formulas from [Baum and Godt (2013)](https://doi.org/10.1002/jgrf.20100).

[^note-7]: **Source:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (5c), using the correction in [Baum and Godt (2013)](https://doi.org/10.1002/jgrf.20100). With $\psi_0=0$, $z_{\mathrm{w},0}=d_{\mathrm{u}}$ and $a=\alpha_1$, this is the same steady-conductivity expression with renamed symbols.

[^note-8]: **Derivation:** inversion of [Srivastava and Yeh (1991)](https://doi.org/10.1029/90WR02772), Eq. (2a), with $\psi_0=0$.

[^note-9]: **Source:** [Srivastava and Yeh (1991)](https://doi.org/10.1029/90WR02772), Eqs. (2a)–(2b).

[^note-10]: **Source:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (6), corresponding to [Srivastava and Yeh (1991)](https://doi.org/10.1029/90WR02772), Eq. (10), after transforming the coordinates and parameters. The sine/cosine form avoids the poles of the equivalent tangent equation $\tan(az_{\mathrm{w},0}r_m)+2r_m=0$.

[^note-11]: **Derived unit response:** subtract the initial steady profile from the constant-flux solution in [Srivastava and Yeh (1991)](https://doi.org/10.1029/90WR02772), Eq. (12a), or the corrected Eq. (5d) of [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), and divide by the change in surface flux. Here the coordinate increases downward, $a$ is the slope-adjusted coefficient, and $\psi_0=0$. The corrected conductivity expression is given by [Baum and Godt (2013)](https://doi.org/10.1002/jgrf.20100).

[^note-12]: **Solution method:** Srivastava and Yeh derive the homogeneous-soil solution using a Laplace transform in time and invert it by calculating residues at the transform's poles, leading to Eq. (12a). The positive characteristic roots identify the transient poles and the corresponding decaying components. Describing the response as a sum of spatial patterns gives a physical reading of the series; its terms are mathematical components, not separate water volumes.

[^note-13]: **Numerical implementation:** [Alvioli and Baum (2016)](https://doi.org/10.1016/j.envsoft.2016.04.002), sections 3.1–3.3, describe corrections to early-time basal-flux calculations in TRIGRS 2.1. The package also compares basal transfer with the released Fortran routines, with explicit handling of their different early-time cutoff; see [verification](verification.md).

[^note-14]: **Adaptation:** drainage follows Eqs. (10a)–(10c) of [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), with the sign correction in [Baum and Godt (2013)](https://doi.org/10.1002/jgrf.20100). The formula above is the source-based capacity obtained with $\beta_{\mathrm{TR}}=\cos^2\delta-I_{\mathrm{ZLT}}/K_{\mathrm{s}}$. It is retained as a separate boundary assumption even when the default pressure reconstruction uses a different, Darcy-consistent gradient. The multiplier $\eta$ plays the role of the source's $c_{\mathrm{d}}$.

[^note-15]: **Implementation adaptation:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eqs. (10b)–(10c) and the numerical volume integration in section 4.2.4, paragraph 29. Applying the zero floor at every integration interval is a package choice.

[^note-16]: **Adaptation:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (11a), defines the same fillable-pore-volume concept. Here it is integrated numerically from the implemented water-content profile rather than evaluated through the source's conductivity primitive.

[^note-17]: **Pressure adaptation:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (3b), and [Baum et al. (2008)](https://doi.org/10.3133/ofr20081159), Eq. (3B), supply the same linear form. The package extends it beyond the source's initially upward-flow case using the selected background gradient. It does not implement the transient saturated diffusion of Eqs. (13)–(14), or the saturated surface/later-time solutions of [Alvioli and Baum (2016)](https://doi.org/10.1016/j.envsoft.2016.04.002), Eqs. (1)–(4).

[^note-18]: **Adaptation:** [Baum et al. (2010)](https://doi.org/10.1029/2009JF001321), Eq. (18), with the unsaturated effective-stress modification in section 6, paragraph 42; also [Baum et al. (2008)](https://doi.org/10.3133/ofr20081159), Eq. (18). The package uses constant soil unit weight and numerical clipping.

[^finite-base]: **Geometry and adaptation:** the [TRIGRS 2.0 manual](https://doi.org/10.3133/ofr20081159), Table 3, distinguishes maximum evaluation depth from initial water-table depth and requires the latter not to exceed the former. This package instead accepts a deeper requested water table and moves it to the finite restrictive base, with zero initial water-column height. The capacity-based reconstruction remains an approximation, without the manual's saturated finite-depth diffusion solution.
