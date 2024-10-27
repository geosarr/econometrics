// pub(crate) fn cdf_beta(x: f64, a: f64, b: f64) -> f64 {
//     if x.is_nan() | a.is_nan() | b.is_nan() {
//         return f64::NAN;
//     }
//     if (a < 0.) | (b < 0.) {
//         return f64::NAN;
//     }
//     if x <= 0. {
//         return 0.;
//     }
//     if x >= 1. {
//         return 1.;
//     }
//     if (a == 0.) | (b == 0.) || a.is_infinite() | b.is_infinite() {
//         // NB:  0 < x < 1 :
//         if (a == 0.) && (b == 0.) {
//             // point mass 1/2 at each of {0,1} :
//             return 0.5;
//         }
//         if (a == 0.) || (a / b == 0.) {
//             // point mass 1 at 0 ==> P(X <= x) = 1, all x > 0
//             return 1.;
//         }
//         if (b == 0.) || (b / a == 0.) {
//             // point mass 1 at 1 ==> P(X <= x) = 0, all x < 1
//             return 0.;
//         }
//         // else, remaining case:  a = b = Inf : point mass 1 at 1/2
//         return if x < 0.5 { 0. } else { 1. };
//     }
//     let mut x1 = 0.5 - x + 0.5;
//     let mut w = 0.;
//     let mut wc = 0.;
//     let mut ierr: isize = 0;
//     //====
//     bratio(a, b, x, x1, &mut w, &mut wc, &mut ierr); /* -> ./toms708.c */
//     //====

//     // ierr in {10,14} <==> bgrat() error code ierr-10 in 1:4; for 1 and 4,
// warned     // *there* if(ierr && ierr != 11 && ierr != 14)
//     // MATHLIB_WARNING4("pbeta_raw(%g, a=%g, b=%g, ..) -> bratio() gave error
// code     // %d", x, a,b, ierr);

//     return w;
// }

// const LNB: f64 = 0.69314718055995;
// fn fpser(a: f64, b: f64, x: f64, eps: f64) -> f64 {
//     /* -----------------------------------------------------------------------
// *

//     * EVALUATION OF I (A,B)
//     * X

//     * FOR B < MIN(EPS, EPS*A) AND X <= 0.5

//     * -----------------------------------------------------------------------
//       */
//     let (mut ans, mut c, mut s, mut t, mut an, mut tol) = (0., 0., 0., 0.,
// 0., 0.);

//     /* SET  ans := x^a : */
//     if a > eps * 0.001 {
//         t = a * x.ln();
//         if t < f64::MIN_EXP as f64 {
//             /* exp(t) would underflow */
//             return 0.;
//         }
//         ans = t.exp();
//     } else {
//         ans = 1.;
//     }

//     /* NOTE THAT 1/B(A,B) = B */
//     ans *= b / a;

//     tol = eps / a;
//     an = a + 1.;
//     t = x;
//     s = t / an;
//     while f64::abs(c) > tol {
//         an += 1.;
//         t = x * t;
//         c = t / an;
//         s += c;
//     }

//     ans *= a * s + 1.;
//     return ans;
// } /* fpser */
// fn swap(a: f64, b: f64, x: f64, y: f64, do_swap: bool) -> (f64, f64, f64,
// f64) {     let (a0, b0, x0, y0) = if do_swap { (b, a, y, x) } else { (a, b,
// x, y) };     (a0, b0, x0, y0)
// }
// fn bratio(
//     a: f64,
//     b: f64,
//     x: f64,
//     y: f64,
//     w: &mut f64,
//     w1: &mut f64,
//     ierr: &mut isize,
//     // log_p: bool,
// ) {
//     // let do_swap: bool = false;
//     let (n, ierr1) = (0, 0);
//     let (z, a0, b0, x0, y0, lambda) = (0., 0., 0., 0., 0., 0.);

//     /*  eps is a machine dependent constant: the smallest
//     * floating point number for which   1. + eps > 1.
//     * NOTE: for almost all purposes it is replaced by 1e-15 (~= 4.5 times
//       larger) below */
//     let mut eps = f64::EPSILON; /* == DBL_EPSILON (in R, Rmath) */
//     /* -----------------------------------------------------------------------
//      */
//     *w = 0.;
//     *w1 = 0.;

//     // safeguard, preventing infinite loops further down
//     if x.is_nan() | y.is_nan() | a.is_nan() | b.is_nan() {
//         *ierr = 9;
//         return;
//     }
//     if (a < 0.) | (b < 0.) {
//         *ierr = 1;
//         return;
//     }
//     if (a == 0.) && (b == 0.) {
//         *ierr = 2;
//         return;
//     }
//     if (x < 0.) | (x > 1.) {
//         *ierr = 3;
//         return;
//     }
//     if (y < 0.) | (y > 1.) {
//         *ierr = 4;
//         return;
//     }

//     /* check that  'y == 1 - x' : */
//     z = x + y - 0.5 - 0.5;

//     if (z.abs() > eps * 3.) {
//         *ierr = 5;
//         return;
//     }

//     // R_ifDEBUG_printf("bratio(a=%g, b=%g, x=%9g, y=%9g, .., log_p=%d): ",
//     //         a,b,x,y, log_p);
//     *ierr = 0;

//     let a_lt_b = (a < b); //line was moved

//     // if (x == 0.) {
//     //     // goto L200
//     // };
//     // if (y == 0.) {
//     //     // goto L210
//     // };

//     // if (a == 0.) {
//     //     // goto L211
//     // };
//     // if (b == 0.){
//     //     // goto L201
//     // };
//     eps = max(eps, 1e-15);

//     if (/* max(a,b) */a.max(b) < eps * 0.001) {
//         /* procedure for a and b <
//         0.001 * eps */
//         // L230:  -- result *independent* of x (!)
//         // *w  = a/(a+b)  and  w1 = b/(a+b) :
//         *w = b / (a + b);
//         *w1 = a / (a + b);
//         // R_ifDEBUG_printf("a & b very small -> simple ratios
// (%g,%g)\n",*w,*w1);         return;
//     }

//     // #define SET_0_noswap \
//     // a0 = a;  x0 = x; \
//     // b0 = b;  y0 = y;

//     // #define SET_0_swap   \
//     // a0 = b;  x0 = y; \
//     // b0 = a;  y0 = x;

//     if (min(a, b) <= 1.) {
//         /*------------------------ a <= 1  or  b <= 1 ----
//          */
//         let mut do_swap = (x > 0.5);
//         let (a0, b0, x0, y0) = if do_swap { (b, a, y, x) } else { (a, b, x,
// y) };         /* now have  x0 <= 1/2 <= y0  (still  x0+y0 == 1) */
//         // R_ifDEBUG_printf(" min(a,b) <= 1, do_swap=%d;", do_swap);

//         if b0 < min(eps, eps * a0) {
//             /* L80: */
//             *w = fpser(a0, b0, x0, eps);
//             *w1 = 0.5 - *w + 0.5;
//             // R_ifDEBUG_printf("  b0 small -> w := fpser(*) = %.15g\n", *w);
//             // goto L_end;
//         }

//         if (a0 < min(eps, eps * b0)) && (b0 * x0 <= 1.) {
//             /* L90: */
//             *w1 = apser(a0, b0, x0, eps);
//             // R_ifDEBUG_printf("  a0 small -> w1 := apser(*) = %.15g\n",
// *w1);             // goto L_end_from_w1;
//         }

//         let mut did_bup = false;
//         if (max(a0, b0) > 1.) {
//             /* L20:  min(a,b) <= 1 < max(a,b) */
//             // R_ifDEBUG_printf("\n L20:  min(a,b) <= 1 < max(a,b); ");
//             if (b0 <= 1.) {
//                 // goto L_w_bpser
//             };

//             if (x0 >= 0.29) { /* was 0.3, PR#13786 */
//                 // goto L_w1_bpser
//             };

//             if (x0 < 0.1 && pow(x0 * b0, a0) <= 0.7) {
//                 // goto L_w_bpser
//             };

//             if (b0 > 15.) {
//                 *w1 = 0.;
//                 // goto L131;
//             }
//         } else {
//             /* a, b <= 1 */
//             // R_ifDEBUG_printf("\n      both a,b <= 1; ");
//             if (a0 >= min(0.2, b0)) {
//                 // goto L_w_bpser
//             };

//             if (pow(x0, a0) <= 0.9) {
//                 // goto L_w_bpser
//             };

//             if (x0 >= 0.3) {
//                 // goto L_w1_bpser
//             };
//         }
//         n = 20; /* goto L130; */
//         // *w1 = bup(b0, a0, y0, x0, n, eps, false); TODO
//         did_bup = true;
//         // R_ifDEBUG_printf("  ... n=20 and *w1 := bup(*) = %.15g; ", *w1);
//         b0 += n;
//         // L131:
//         // R_ifDEBUG_printf(" L131: bgrat(*, w1=%.15g) ", *w1);
//         // bgrat(b0, a0, y0, x0, w1, 15*eps, &ierr1, false);
//         // #ifdef DEBUG_bratio
//         // REprintf(" ==> new w1=%.15g", *w1);
//         // if(ierr1) REprintf(" ERROR(code=%d)\n", ierr1) ; else
// REprintf("\n");         // #endif
//         if (*w1 == 0) || ((0 < *w1) && (*w1 < f64::MIN)) {
//             // w1=0 or very close:
//             // "almost surely" from underflow, try more: [2013-03-04]
//             // FIXME: it is even better to do this in bgrat *directly* at
// least for the case             //  !did_bup, i.e., where *w1 = (0 or -Inf) on
// entry             // R_ifDEBUG_printf(" denormalized or underflow (?) ->
// retrying: ");             if did_bup {
//                 // re-do that part on log scale:
//                 *w1 = bup(b0 - n, a0, y0, x0, n, eps, true);
//             } else {
//                 *w1 = f64::NEG_INFINITY
//             }; // = 0 on log-scale
//             bgrat(b0, a0, y0, x0, w1, 15 * eps, &ierr1, true);
//             if ierr1 {
//                 *ierr = 10 + ierr1
//             };
//             // #ifdef DEBUG_bratio
//             //     REprintf(" ==> new log(w1)=%.15g", *w1);
//             //     if(ierr1) REprintf(" Error(code=%d)\n", ierr1) ; else
//             // REprintf("\n"); #endif
//             // goto L_end_from_w1_log;
//         } else {
//             if ierr1 {
//                 *ierr = 10 + ierr1
//             };
//             // RustTODO if (*w1 < 0)
//             // MATHLIB_WARNING4("bratio(a=%g, b=%g, x=%g): bgrat() -> w1 =
// %g",             //         a,b,x, *w1);
//             // goto L_end_from_w1;
//         }
//     } else {
//         /* L30: -------------------- both  a, b > 1  {a0 > 1  &  b0 > 1}
//         ---*/
//         /* lambda := a y - b x  =  (a + b)y  =  a - (a+b)x    {using x + y ==
// 1},
//          * ------ using the numerically best version : */
//         lambda = if (a + b).is_finite() {
//             if a > b {
//                 (a + b) * y - b
//             } else {
//                 a - (a + b) * x
//             }
//         } else {
//             a * y - b * x
//         };
//         do_swap = (lambda < 0.);
//         if do_swap {
//             lambda = -lambda;
//             SET_0_swap;
//         } else {
//             SET_0_noswap;
//         }

//         //     R_ifDEBUG_printf("  L30:  both  a, b > 1; |lambda| = %#g,
// do_swap =         // %d\n",             lambda, do_swap);

//         if (b0 < 40.) {
//             // R_ifDEBUG_printf("  b0 < 40;");
//             // if (b0 * x0 <= 0.7
//             // || (log_p && lambda > 650.)) // << added 2010-03; svn r51327
//             // goto L_w_bpser;
//             // else
//             // goto L140;
//         } else if (a0 > b0) { /* ----  a0 > b0 >= 40  ---- */
//             // R_ifDEBUG_printf("  a0 > b0 >= 40;");
//             // if (b0 <= 100. || lambda > b0 * 0.03)
//             // goto L_bfrac;
//         } else if (a0 <= 100.) {
//             // R_ifDEBUG_printf("  a0 <= 100; a0 <= b0 >= 40;");
//             // goto L_bfrac;
//         } else if (lambda > a0 * 0.03) {
//             // R_ifDEBUG_printf("  b0 >= a0 > 100; lambda > a0 * 0.03 ");
//             // goto L_bfrac;
//         }

//         //     /* else if none of the above    L180: */
//         //     *w = basym(a0, b0, lambda, eps * 100., log_p);
//         //     *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
//         //     R_ifDEBUG_printf("  b0 >= a0 > 100; lambda <= a0 * 0.03: *w:=
//         // basym(*) =%.15g\n",             *w);
//         //     goto L_end;

//         //     } /* else: a, b > 1 */
//         //     /*            EVALUATION OF THE APPROPRIATE ALGORITHM */
//         //     L_w_bpser: // was L100
//         //     *w = bpser(a0, b0, x0, eps, log_p);
//         //     *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
//         //     R_ifDEBUG_printf(" L_w_bpser: *w := bpser(*) = %.15g\n", *w);
//         //     goto L_end;

//         //     L_w1_bpser:  // was L110
//         //     *w1 = bpser(b0, a0, y0, eps, log_p);
//         //     *w  = log_p ? R_Log1_Exp(*w1) : 0.5 - *w1 + 0.5;
//         //     R_ifDEBUG_printf(" L_w1_bpser: *w1 := bpser(*) = %.15g\n",
// *w1);         //     goto L_end;

//         //     L_bfrac:
//         //     *w = bfrac(a0, b0, x0, y0, lambda, eps * 15., log_p);
//         //     *w1 = log_p ? R_Log1_Exp(*w) : 0.5 - *w + 0.5;
//         //     R_ifDEBUG_printf(" L_bfrac: *w := bfrac(*) = %g\n", *w);
//         //     goto L_end;

//         //     L140:
//         //     /* b0 := fractional_part( b0 )  in (0, 1]  */
//         //     n = (int) b0;
//         //     b0 -= n;
//         //     if (b0 == 0.) {
//         //     --n; b0 = 1.;
//         //     }

//         //     *w = bup(b0, a0, y0, x0, n, eps, false);

//         //     if(*w < DBL_MIN && log_p) { /* do not believe it; try bpser()
// :         // */     R_ifDEBUG_printf(" L140: bup(b0=%g,..)=%.15g < DBL_MIN
//         // - not used; ", b0, *w);     /*revert: */ b0 += n; /* which is only
//         //   valid if b0 <= 1 || b0*x0 <= 0.7 */ goto L_w_bpser; }
//         //   R_ifDEBUG_printf(" L140: *w := bup(b0=%g,..) = %.15g; ", b0,
// *w);         //   if (x0 <= 0.7) { /* log_p :  TODO:  w = bup(.) + bpser(.)
// -- not         //   so easy to use
//         // log-scale */     *w += bpser(a0, b0, x0, eps, /* log_p = */
//         // false);     R_ifDEBUG_printf(" x0 <= 0.7: *w := *w + bpser(*)
//         // = %.15g\n", *w);     goto L_end_from_w;
//         //     }
//         //     /* L150: */
//         //     if (a0 <= 15.) {
//         //     n = 20;
//         //     *w += bup(a0, b0, x0, y0, n, eps, false);
//         //     R_ifDEBUG_printf("\n a0 <= 15: *w := *w + bup(*) = %.15g;",
// *w);         //     a0 += n;
//         //     }
//         //     R_ifDEBUG_printf(" bgrat(*, w=%.15g) ", *w);
//         //     bgrat(a0, b0, x0, y0, w, 15*eps, &ierr1, false);
//         //     if(ierr1) *ierr = 10 + ierr1;
//         //     #ifdef DEBUG_bratio
//         //     REprintf("==> new w=%.15g", *w);
//         //     if(ierr1) REprintf(" Error(code=%d)\n", ierr1) ; else
//         // REprintf("\n");     #endif
//         //     goto L_end_from_w;

//         //     /* TERMINATION OF THE PROCEDURE */
//         //     L200:
//         //     if (a == 0.) { *ierr = 6;    return; }
//         //     // else:
//         //     L201: *w  = R_D__0; *w1 = R_D__1; return;

//         //     L210:
//         //     if (b == 0.) { *ierr = 7;    return; }
//         //     // else:
//         //     L211: *w  = R_D__1; *w1 = R_D__0; return;

//         //     L_end_from_w:
//         //     if(log_p) {
//         //     *w1 = log1p(-*w);
//         //     *w  = log(*w);
//         //     } else {
//         //     *w1 = 0.5 - *w + 0.5;
//         //     }
//         //     goto L_end;

//         //     L_end_from_w1:
//         //     if(log_p) {
//         //     *w  = log1p(-*w1);
//         //     *w1 = log(*w1);
//         //     } else {
//         //     *w = 0.5 - *w1 + 0.5;
//         //     }
//         //     goto L_end;

//         //     L_end_from_w1_log:
//         //     // *w1 = log(w1) already; w = 1 - w1  ==> log(w) = log(1 - w1)
// =         // log(1 - exp(*w1))     if(log_p) {
//         //     *w = R_Log1_Exp(*w1);
//         //     } else {
//         //     *w  = /* 1 - exp(*w1) */ -expm1(*w1);
//         //     *w1 = exp(*w1);
//         //     }
//         //     goto L_end;

//         //     L_end:
//         //     if (do_swap) { /* swap */
//         //     double t = *w; *w = *w1; *w1 = t;
//         //     }
//         //     return;
//     }
// } /* bratio */
