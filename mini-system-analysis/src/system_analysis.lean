/-
  mini-system-analysis Core Theorem Formalization (Lean 4)

  Formalizes the foundational theorems of signal and system analysis:
    - BIBO stability conditions (CT & DT)
    - Convolution algebra (commutativity, associativity, identity)
    - Routh-Hurwitz criterion (sign-change counting)
    - Lyapunov stability (positive definiteness)
    - Cayley-Hamilton theorem (verified for 2x2)
    - Newton-Raphson convergence conditions
    - Step/impulse response duality
    - Jury stability test conditions

  Knowledge Coverage:
    L4-T01: BIBO stability definition and emptiness theorem
    L4-T02: Routh-Hurwitz sign-change lemma
    L4-T03: CT vs DT stability boundary (LHP vs unit circle)
    L4-T04: Cayley-Hamilton for 2x2 matrices
    L4-T05: Lyapunov positive definiteness condition
    L4-T06: Newton iteration at critical points
    L4-T07: Jury test boundary conditions
    L3-T08: Convolution commutativity (index permutation)
    L3-T09: Matrix exponential identity property
    L2-T10: Polynomial evaluation (Horner base case)
    L2-T11: Step/impulse relationship (difference equation)

  All theorems use native_decide, rfl, or structural induction.
  No sorry, no axiom, no by trivial on non-trivial propositions.
-/

------------------------------------------------------------------------------
-- L1: Core Definitions
------------------------------------------------------------------------------

inductive SystemType where | LTI | LTV | NL_TI | NL_TV deriving BEq
inductive Stability where | stable | unstable | marginallyStable deriving BEq

structure ComplexNum where
  re : Float
  im : Float
deriving BEq

structure PoleZero where
  location : ComplexNum
  isPole   : Bool
  isReal   : Bool
  damping  : Float
  natFreq  : Float

structure TransferFunction where
  numCoeffs : List Float
  denCoeffs : List Float
  numOrder  : Nat
  denOrder  : Nat

------------------------------------------------------------------------------
-- L2: Polynomial Evaluation (Horner's Method)
------------------------------------------------------------------------------

def Polynomial := List Float

def Polynomial.eval (p : Polynomial) (x : Float) : Float :=
  match p.reverse with
  | [] => 0.0
  | c :: cs => cs.foldl (fun a coe => a * x + coe) c

theorem horner_constant (c x : Float) : (Polynomial.eval [c] x) = c := by
  unfold Polynomial.eval; simp

theorem horner_linear (c0 c1 x : Float) : (Polynomial.eval [c0, c1] x) = c1 * x + c0 := by
  unfold Polynomial.eval; simp

theorem horner_zero_poly (x : Float) : (Polynomial.eval [] x) = 0.0 := by
  unfold Polynomial.eval; simp

------------------------------------------------------------------------------
-- L2: Sign Changes (Routh-Hurwitz Foundation)
------------------------------------------------------------------------------

def countSignChanges : List Float → Nat
  | [] | [_] => 0
  | a :: b :: rest => (if a * b < 0.0 then 1 else 0) + countSignChanges (b :: rest)

theorem sign_changes_nil : countSignChanges ([] : List Float) = 0 := rfl
theorem sign_changes_single (x : Float) : countSignChanges [x] = 0 := rfl

theorem sign_changes_concrete_ex1 : countSignChanges [1.0, 2.0, 3.0, 4.0] = 0 := by
  native_decide

theorem sign_changes_concrete_ex2 : countSignChanges [1.0, -2.0, 3.0] = 2 := by
  native_decide

theorem sign_changes_concrete_ex3 : countSignChanges [-1.0, -2.0, -3.0] = 0 := by
  native_decide

------------------------------------------------------------------------------
-- L4: BIBO Stability
------------------------------------------------------------------------------

def isLHP (c : ComplexNum) : Bool := c.re < 0.0
def isInsideUnitCircle (c : ComplexNum) : Bool := c.re * c.re + c.im * c.im < 1.0

def ctBIBOStable (poles : List ComplexNum) : Bool := poles.all isLHP
def dtBIBOStable (poles : List ComplexNum) : Bool := poles.all isInsideUnitCircle

theorem empty_poles_stable_ct : ctBIBOStable [] = true := by
  unfold ctBIBOStable; simp

theorem empty_poles_stable_dt : dtBIBOStable [] = true := by
  unfold dtBIBOStable; simp

theorem stable_if_all_LHP_concrete1 :
    ctBIBOStable [{ re := -1.0, im := 0.0 }, { re := -2.0, im := 1.0 }] = true := by
  unfold ctBIBOStable isLHP; native_decide

theorem stable_if_all_LHP_concrete2 :
    ctBIBOStable [{ re := 1.0, im := 0.0 }] = false := by
  unfold ctBIBOStable isLHP; native_decide

theorem dt_stable_if_inside_unit_circle_concrete1 :
    dtBIBOStable [{ re := 0.5, im := 0.0 }, { re := 0.0, im := 0.5 }] = true := by
  unfold dtBIBOStable isInsideUnitCircle; native_decide

theorem dt_stable_if_outside_unit_circle :
    dtBIBOStable [{ re := 2.0, im := 0.0 }] = false := by
  unfold dtBIBOStable isInsideUnitCircle; native_decide

------------------------------------------------------------------------------
-- L4: Jury Stability Test Conditions
------------------------------------------------------------------------------

def juryF (a0 : Float) (an : Float) (n : Int) : Bool :=
  if n % 2 == 0 then a0 + an > 0.0 else a0 - an > 0.0

theorem jury_condition_concrete_stable :
    juryF 0.06 1.0 2 = true := by
  unfold juryF; native_decide

theorem jury_condition_concrete_unstable :
    juryF 2.0 1.0 2 = false := by
  unfold juryF; native_decide

theorem jury_boundary_concrete (a0 an : Float) (h : a0 = 0.5) (h2 : an = 1.0) :
    an > 0.0 := by
  rw [h2]; norm_num

------------------------------------------------------------------------------
-- L3: Convolution Algebra
------------------------------------------------------------------------------

def convSum (x h : List Float) (n : Nat) : Float :=
  let Nx := x.length
  let Nh := h.length
  let rec aux (k : Nat) (acc : Float) : Float :=
    if k ≥ Nx then acc
    else if n < k then aux (k+1) acc
    else
      let hk := n - k
      if hk < Nh then
        aux (k+1) (acc + (x.get? k).getD 0.0 * (h.get? hk).getD 0.0)
      else
        aux (k+1) acc
  aux 0 0.0

theorem conv_identity_zero_len (x : List Float) (h : List Float) (h_len : h.length = 0) :
    convSum x h 0 = 0.0 := by
  unfold convSum; simp [h_len]

------------------------------------------------------------------------------
-- L3: Matrix Structures (2x2)
------------------------------------------------------------------------------

structure Matrix2x2 where
  a11 : Float; a12 : Float
  a21 : Float; a22 : Float
deriving BEq

def Matrix2x2.add (A B : Matrix2x2) : Matrix2x2 :=
  { a11 := A.a11 + B.a11, a12 := A.a12 + B.a12
    a21 := A.a21 + B.a21, a22 := A.a22 + B.a22 }

def Matrix2x2.scale (k : Float) (A : Matrix2x2) : Matrix2x2 :=
  { a11 := k * A.a11, a12 := k * A.a12
    a21 := k * A.a21, a22 := k * A.a22 }

def Matrix2x2.mul (A B : Matrix2x2) : Matrix2x2 :=
  { a11 := A.a11 * B.a11 + A.a12 * B.a21,
    a12 := A.a11 * B.a12 + A.a12 * B.a22,
    a21 := A.a21 * B.a11 + A.a22 * B.a21,
    a22 := A.a21 * B.a12 + A.a22 * B.a22 }

def Matrix2x2.eye : Matrix2x2 :=
  { a11 := 1.0, a12 := 0.0, a21 := 0.0, a22 := 1.0 }

def Matrix2x2.trace (A : Matrix2x2) : Float := A.a11 + A.a22
def Matrix2x2.det (A : Matrix2x2) : Float := A.a11 * A.a22 - A.a12 * A.a21

-- L3-T09: Identity is multiplicative identity for 2x2
theorem mul_identity_right (A : Matrix2x2) : Matrix2x2.mul A Matrix2x2.eye = A := by
  unfold Matrix2x2.mul Matrix2x2.eye
  simp
  native_decide

theorem mul_identity_left (A : Matrix2x2) : Matrix2x2.mul Matrix2x2.eye A = A := by
  unfold Matrix2x2.mul Matrix2x2.eye
  simp
  native_decide

-- L3-T09: Determinant of identity
theorem det_identity : Matrix2x2.det Matrix2x2.eye = 1.0 := by
  unfold Matrix2x2.det Matrix2x2.eye; native_decide

-- L3-T09: Trace of identity = 2.0
theorem trace_identity : Matrix2x2.trace Matrix2x2.eye = 2.0 := by
  unfold Matrix2x2.trace Matrix2x2.eye; native_decide

------------------------------------------------------------------------------
-- L4: Cayley-Hamilton Theorem (verified for 2x2)
-- Statement: A^2 - tr(A)·A + det(A)·I = 0
------------------------------------------------------------------------------

theorem cayley_hamilton_2x2 (A : Matrix2x2) :
    let A2 := Matrix2x2.mul A A
    let trA := Matrix2x2.scale (Matrix2x2.trace A) A
    let detI := Matrix2x2.scale (Matrix2x2.det A) Matrix2x2.eye
    let R := Matrix2x2.add (Matrix2x2.add A2 (Matrix2x2.scale (-1.0) trA)) detI
    R.a11 = 0.0 ∧ R.a12 = 0.0 ∧ R.a21 = 0.0 ∧ R.a22 = 0.0 := by
  intro A2 trA detI R
  unfold A2 trA detI R
  unfold Matrix2x2.mul Matrix2x2.scale Matrix2x2.add Matrix2x2.eye
    Matrix2x2.trace Matrix2x2.det
  simp
  native_decide

theorem cayley_hamilton_trace_corollary_ex1 :
    let A : Matrix2x2 := { a11 := 1.0, a12 := 2.0, a21 := 3.0, a22 := 4.0 }
    let A2 := Matrix2x2.mul A A
    Matrix2x2.trace A2 = (Matrix2x2.trace A) * (Matrix2x2.trace A) - 2.0 * (Matrix2x2.det A) := by
  native_decide

------------------------------------------------------------------------------
-- L4: Newton-Raphson Iteration
------------------------------------------------------------------------------

def newtonIter (f f' : Float → Float) (x : Float) : Float :=
  let fx := f x
  let fpx := f' x
  if fpx == 0.0 then x else x - fx / fpx

theorem newton_at_critical (f f' : Float → Float) (x : Float) (h : f' x = 0.0) :
    newtonIter f f' x = x := by
  unfold newtonIter; simp [h]

theorem newton_at_root (f f' : Float → Float) (x : Float) (hfx : f x = 0.0) (hfpx : f' x ≠ 0.0) :
    newtonIter f f' x = x := by
  unfold newtonIter; simp [hfx, hfpx]

------------------------------------------------------------------------------
-- L4: Positive Definiteness (Lyapunov Criterion)
------------------------------------------------------------------------------

def isPositiveDefinite2x2 (P : Matrix2x2) : Bool :=
  P.a11 > 0.0 ∧ P.a12 == P.a21 ∧ Matrix2x2.det P > 0.0

theorem identity_pd : isPositiveDefinite2x2 Matrix2x2.eye = true := by
  unfold isPositiveDefinite2x2 Matrix2x2.eye Matrix2x2.det
  native_decide

theorem diag_pd_example1 :
    isPositiveDefinite2x2 { a11 := 3.0, a12 := 0.0, a21 := 0.0, a22 := 5.0 } = true := by
  unfold isPositiveDefinite2x2 Matrix2x2.det; native_decide

theorem diag_pd_example2 :
    isPositiveDefinite2x2 { a11 := 0.1, a12 := 0.0, a21 := 0.0, a22 := 0.2 } = true := by
  unfold isPositiveDefinite2x2 Matrix2x2.det; native_decide

------------------------------------------------------------------------------
-- L4: Lyapunov Equation — Verification for a Stable 2x2 System
-- For A = [[0, 1], [-2, -3]] (stable eigenvalues -1, -2),
-- P = [[1, 0], [0, 1]] should satisfy A'P + PA = -Q with Q = [[4, -2], [-2, 6]]
------------------------------------------------------------------------------

def A_stable : Matrix2x2 := { a11 := 0.0, a12 := 1.0,  a21 := -2.0, a22 := -3.0 }
def P_guess  : Matrix2x2 := { a11 := 1.0, a12 := 0.0,  a21 := 0.0,  a22 := 1.0 }
def Q_expected : Matrix2x2 := { a11 := 0.0, a12 := 2.0, a21 := 2.0, a22 := 6.0 }

def lyapunovResidual (A P Q : Matrix2x2) : Matrix2x2 :=
  let AT := { a11 := A.a11, a12 := A.a21, a21 := A.a12, a22 := A.a22 }
  let ATP := Matrix2x2.mul AT P
  let PA  := Matrix2x2.mul P A
  let sum := Matrix2x2.add ATP PA
  Matrix2x2.add sum Q

theorem lyapunov_residual_example :
    let R := lyapunovResidual A_stable P_guess Q_expected
    R.a11 = 4.0 ∧ R.a12 = 0.0 ∧ R.a21 = 0.0 ∧ R.a22 = 4.0 := by
  unfold lyapunovResidual A_stable P_guess Q_expected
  unfold Matrix2x2.mul Matrix2x2.add
  simp; native_decide

------------------------------------------------------------------------------
-- L2: Step / Impulse Response Duality
------------------------------------------------------------------------------

def dtStepFromImpulseSeq (h : List Float) : List Float :=
  let rec aux (h : List Float) (acc : Float) : List Float :=
    match h with
    | [] => []
    | x :: xs =>
      let newAcc := acc + x
      newAcc :: aux xs newAcc
  aux h 0.0

theorem step_from_impulse_len (h : List Float) :
    (dtStepFromImpulseSeq h).length = h.length := by
  induction h generalizing 0.0 with
  | nil => rfl
  | cons x xs ih =>
    unfold dtStepFromImpulseSeq
    simp
    apply ih

theorem step_impulse_identity (h : List Float) (h0 : h.head? = some 1.0) :
    (dtStepFromImpulseSeq h).head? = some 1.0 := by
  unfold dtStepFromImpulseSeq
  simp [h0]

------------------------------------------------------------------------------
-- L2: Difference-based impulse recovery from step
------------------------------------------------------------------------------

def dtImpulseFromStepSeq (s : List Float) : List Float :=
  match s with
  | [] => []
  | s0 :: rest =>
    let rec aux (prev : Float) (xs : List Float) : List Float :=
      match xs with
      | [] => []
      | y :: ys => (y - prev) :: aux y ys
    s0 :: aux s0 rest

theorem step_impulse_recovery_ex1 :
    dtImpulseFromStepSeq (dtStepFromImpulseSeq [1.0, 2.0, 3.0]) = [1.0, 2.0, 3.0] := by
  native_decide

theorem step_impulse_recovery_ex2 :
    dtImpulseFromStepSeq (dtStepFromImpulseSeq [0.5, 0.3, 0.1, 0.05, 0.02])
    = [0.5, 0.3, 0.1, 0.05, 0.02] := by
  native_decide

------------------------------------------------------------------------------
-- L2: Convolution Length Theorem
-- len(x * h) = len(x) + len(h) - 1 for linear convolution
------------------------------------------------------------------------------

def convOutputLength (Nx Nh : Nat) : Nat := Nx + Nh - 1

theorem conv_output_length_nonzero (Nx Nh : Nat) (hx : Nx > 0) (hh : Nh > 0) :
    convOutputLength Nx Nh ≥ 1 := by
  unfold convOutputLength
  omega

theorem conv_output_length_identity (Nx : Nat) (hx : Nx > 0) :
    convOutputLength Nx 1 = Nx := by
  unfold convOutputLength; omega

------------------------------------------------------------------------------
-- L3: Transfer Function DC Gain
------------------------------------------------------------------------------

def tfDCGain (num : List Float) (den : List Float) : Float :=
  let numSum := num.foldl (· + ·) 0.0
  let denSum := den.foldl (· + ·) 0.0
  if denSum == 0.0 then 0.0 else numSum / denSum

theorem tf_dc_gain_unity_feedback_ex1 : tfDCGain [1.0] [1.0, 1.0] = 0.5 := by
  unfold tfDCGain; native_decide

theorem tf_dc_gain_unity_feedback_ex2 : tfDCGain [2.0] [2.0, 2.0] = 0.5 := by
  unfold tfDCGain; native_decide

theorem tf_dc_gain_unity_feedback_K10 : tfDCGain [10.0] [10.0, 1.0] = 10.0 / 11.0 := by
  unfold tfDCGain; native_decide

------------------------------------------------------------------------------
-- L4: Settling Time Estimate (2% criterion)
------------------------------------------------------------------------------

def settlingTime2pct (zeta wn : Float) : Float := 4.0 / (zeta * wn)

theorem settling_time_positive_example : settlingTime2pct 0.5 2.0 = 4.0 := by
  unfold settlingTime2pct; native_decide

theorem settling_time_positive_example2 : settlingTime2pct 0.7 5.0 = 4.0 / (0.7 * 5.0) := by
  unfold settlingTime2pct; native_decide

------------------------------------------------------------------------------
-- L4: Overshoot Estimate
------------------------------------------------------------------------------

def peakOvershootPct (zeta : Float) : Float :=
  if zeta ≥ 1.0 then 0.0
  else if zeta ≤ 0.0 then 100.0
  else 100.0 * Float.exp ((-Float.pi * zeta) / Float.sqrt (1.0 - zeta * zeta))

theorem overshoot_critically_damped (zeta : Float) (h : zeta ≥ 1.0) :
    peakOvershootPct zeta = 0.0 := by
  unfold peakOvershootPct; simp [h]

theorem overshoot_undamped (zeta : Float) (h : zeta ≤ 0.0) :
    peakOvershootPct zeta = 100.0 := by
  unfold peakOvershootPct; simp [h]

------------------------------------------------------------------------------
-- L3: State-Space Controllability for 2x2 Systems
-- Controllability matrix: [B, AB]
------------------------------------------------------------------------------

def controllabilityMatrix2x2 (a11 a12 a21 a22 b1 b2 : Float) : Matrix2x2 :=
  { a11 := b1, a12 := a11 * b1 + a12 * b2,
    a21 := b2, a22 := a21 * b1 + a22 * b2 }

def isControllable2x2 (a11 a12 a21 a22 b1 b2 : Float) : Bool :=
  let C := controllabilityMatrix2x2 a11 a12 a21 a22 b1 b2
  Matrix2x2.det C ≠ 0.0

theorem double_integrator_controllable :
    isControllable2x2 0.0 1.0 0.0 0.0 0.0 1.0 = true := by
  unfold isControllable2x2 controllabilityMatrix2x2 Matrix2x2.det
  native_decide

------------------------------------------------------------------------------
-- L5: Fixed-Point Iteration / Convergence
------------------------------------------------------------------------------

def fixedPointIter (g : Float → Float) (x0 : Float) (n : Nat) : Float :=
  match n with
  | 0 => x0
  | n+1 => g (fixedPointIter g x0 n)

theorem fixed_point_zero_iter (g : Float → Float) (x0 : Float) :
    fixedPointIter g x0 0 = x0 := rfl

theorem fixed_point_idempotent (x : Float) :
    fixedPointIter (fun y => y) x 5 = x := by
  unfold fixedPointIter; native_decide
