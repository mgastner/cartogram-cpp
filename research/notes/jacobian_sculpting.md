# Sculpting the Jacobian: smoothing ln g instead of the density

Status: idea for after the ICC 2027 paper (recorded 9 Nov 2026). Nothing here
is implemented. Measured statements cite the research plan
(`~/.claude/plans/analyze-the-research-plan-imperative-donut.md`, sections in
brackets). Everything else is argument or hypothesis and is marked as such.

## 1. The question

The conformality postprocessor composes a cartogram T0 with an
area-preserving map A. It keeps the Jacobian g0 = det DT0 pointwise and
minimizes distortion over the maps with that Jacobian. The cartogram
constraint itself is weaker: it fixes only the region totals
∫_{U_i} g dx = T_i (plus g > 0 and the rectangle). Within each region, g is
free. cartogram-cpp chooses it implicitly, through rasterization, blur,
the number of passes and the stopping rule.

So the Jacobian is a design variable that we have not yet used. The question
is which g, with the same region totals, admits the least distorted map.

## 2. Evidence that g matters more than the optimization

- **The Jacobian outweighs the start** (China 16², α = ½; [F22]).
  - Same Jacobian, two starts 4.9 lattice units apart: the optimized maps end
    0.28 (flow Jacobian) and 0.34 (diffusion Jacobian) apart, and E½ differs
    by 0.26% and 0.08%.
  - Different Jacobians (flow vs diffusion): the optimized maps end 3.26
    apart, every province by 2.0–6.2, and E½ differs by 10%
    (0.616 vs 0.553).
  - The same direction elsewhere [F21]: E½ after A is 5% lower with the
    diffusion Jacobian for the USA (0.2349 vs 0.2233) and 3.6% lower for
    Switzerland.
- **The two Jacobians differ only inside regions** [F19]. Region areas
  agree, but per 2 × 2 field cell the ratio g0(diffusion)/g0(flow) has
  p5–p95 0.92–1.09 (USA), 0.83–1.38 (China) and 0.92–1.11 (Switzerland).
- **The distortion left after A follows Δ ln g0** [F13]. The Spearman
  correlation between |Δ ln g0| and the excess MIPS energy m after A is 0.66
  (USA 64²), 0.57 (China 32²) and 0.45 (World), all at α = 0. The 10% of
  cells with the largest |Δ ln g0| hold 33–36% of the remaining m on USA and
  China.

## 3. Why Δ ln g, and not the density's gradient

**Local obstruction (standard).** If T is conformal at x with det DT = g,
then DT = √g R with R a rotation. T is then an isometry from the metric
g |dx|² to the flat plane. A metric e^φ |dx|² has Gaussian curvature
K = −½ e^{−φ} Δφ. Hence, with φ = ln g:

- A conformal map with Jacobian g exists locally if and only if ln g is
  harmonic. ("If": on a simply connected patch, take the harmonic conjugate
  h of ½ ln g and integrate f′ = exp(½ ln g + i h).)
- The curvature of g |dx|² is what any map with Jacobian g must absorb as
  distortion.

**Global obstruction.** A conformal map of the rectangle onto itself that
maps sides to sides is a similarity, so g ≡ const. For any non-uniform g,
some distortion is unavoidable. The design question is where it sits and how
large it is.

**Minimum principle (already in the paper plan, F2).** A harmonic ln g has no
interior extrema. A dense city surrounded by sparse land is an interior
maximum of ln g, so the map must be non-conformal around it. This matches the
persistent convergence points (western China, the US Mountain West).

**Tobler's smoothness targets the wrong quantity.** Pycnophylactic
interpolation (Tobler 1979) minimizes ∫ |∇ρ|² with the region totals fixed,
so ρ is harmonic away from the multipliers. For conformality, ln ρ should be
harmonic instead. A one-dimensional example shows the difference:

| Density ramp | Tobler's choice? | Δ ln ρ | Compatible with a conformal map? |
|---|---|---|---|
| ρ = 1 + a x (linear) | yes, ρ is harmonic | −a²/(1 + a x)² ≠ 0 | no |
| ρ = e^{c x} (exponential) | no | 0 | yes, locally: f(z) = e^{cz/2} has \|f′\|² ∝ e^{cx} |

So "smooth pycnophylactic" is the right family, but the smoothing target
should be ln g, and the obstruction is its Laplacian, not the gradient of ρ.

## 4. A candidate problem (hypothesis)

Choose φ = ln g on the lattice to minimize

  F[φ] = ∫ w(x) (Δφ)² dx   subject to   ∫_{U_i} e^φ dx = T_i for every
  region i, including the ocean,

with wall conditions consistent with the cosine basis (∂_n φ = 0).

- **Euler–Lagrange equation:** Δ(w Δφ) = λ_i e^φ in U_i, with one multiplier
  per region. For w = 1 this is Δ²φ = λ_i e^φ: a fourth-order, thin-plate-type
  problem in ln g. It can be solved by Newton iteration with DCT-based
  bi-Laplacian solves, the same machinery as the flow.
- **Consequence:** F < ∞ needs φ ∈ H², so ln g is C¹ across region borders.
  The density can no longer jump at a border. The curvature then spreads
  over a band around each border instead of sitting on it.
- **Weight w (open).**
  - w = 1 measures curvature per unit original area.
  - The intrinsic choice is the Calabi energy ∫ K² dA_g = ¼ ∫ e^{−φ}(Δφ)² dx,
    per unit cartogram area.
  - Which matches E_α (α = 0, ½, 1) is open. An inverse-map symmetry argument
    like the one for α = ½ might decide it.
- **Caveat.** F is a surrogate. I know no bound that turns a small F into a
  small achievable E_α. The evidence is the correlation in §2, and the next
  section shows that a crude version of this measure fails.

## 5. Counter-evidence: the crude measure fails

The mean |Δ ln g0| over field cells is about the same for the flow and
diffusion Jacobians (within ±10%; China slightly *higher* for diffusion)
[F19]. Yet E½ after A is 10% lower with the diffusion Jacobian. So the mean
roughness does not rank the Jacobians. Possible reasons, untested:

- what matters is where the curvature sits (relative to the weights and to
  region borders), not its mean magnitude;
- the relevant norm is weaker than a pointwise mean, e.g. an H⁻¹-type norm
  of K, because A can spread distortion;
- finite differences of a piecewise-affine T0 are noisy at the cell scale.

The first experiment below tests whether F does better.

## 6. Pipeline sketch

1. Rasterize the targets, which gives a piecewise-constant ρ.
2. Solve §4 for φ, and set ρ_design = e^φ.
3. Equalize ρ_design. Its Jacobian is then ρ_design/ρ̄, up to
   discretization.
   - Problem: cartogram-cpp re-rasterizes piecewise-constant densities from
     the polygons in every pass, which would erase the design after pass 1.
   - Options: one accurate equalization of ρ_design (diffusion, or the flow
     at a fine grid), or transport the designed field with the map instead of
     re-rasterizing.
4. Run the conformality postprocessor A as now.

## 7. Experiments to run first

1. Evaluate F (w = 1 and the Calabi weight) on the flow and diffusion
   Jacobians of USA, China 16²/32² and Switzerland. Does F rank them as E½
   does, where mean |Δ ln g0| failed?
2. China 16², α = ½: design g by §4, equalize, optimize A, and compare E½
   with the flow (0.616) and diffusion (0.553) Jacobians. Use the same
   registration and distance measures as in F22.
3. Vary the weight w and the treatment of the ocean.
4. Check that polygon measures (Miaji et al. 2026) improve together with E.

## 8. Related work to check before citing

Verified (in `paper/references.bib`):
- Gastner, Miaji & Singhania (2022), *Smooth pycnophylactic interpolation
  produced by density-equalising map projections*, Kartografija i
  geoinformacije 21(37) 60–69. It obtains smooth pycnophylactic
  interpolations from density-equalizing maps. Check how its notion of
  smoothness relates to Δ ln g.

Not yet verified (from memory; fetch and check before citing):
- Tobler (1979), smooth pycnophylactic interpolation, JASA.
- Calabi (1982), extremal Kähler metrics: the Calabi energy ∫ K² dA.
- Springborn, Schröder & Pinkall (2008), conformal equivalence of triangle
  meshes: log scale factors as the variable.
- Sawhney & Crane (2017), boundary first flattening: conformal flattening
  with prescribed curvature, through the Yamabe equation for the log scale
  factor.

## 9. Relation to other ideas in the plan

- **Interleaving corrections with density equalization (F13)** does not
  change the feasible set, because g is unchanged. Sculpting changes g, so
  it is the lever that interleaving lacks.
- **Partial equalization (F3(b)):** equalizing ρ^{1/2} gives region areas
  ∝ √(A_i P_i). That is a different trade-off: it gives up area accuracy,
  whereas sculpting keeps the region totals exact.
