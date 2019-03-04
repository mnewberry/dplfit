(* Generate synthetic branching datasets on-the-fly *)
let usage = Printf.sprintf 
  "Usage: %s [options] [[-t] (symmetric|asymmetric|random|randombeta|multifurcating)]\n%!" Sys.argv.(0)

let argset ref v =
  if !ref = None then ref := Some v
  else raise (Arg.Bad "Only one argument is supported")

let beta1 = ref None and beta2 = ref None and alpha = ref None 
  and dist = ref None and xmin = ref None and xinit = ref None
  and err = ref None

let argspec = [
  ("-t", Arg.String (argset dist), "Type of branching process");
  ("-x", Arg.Float (argset xmin), "x_min value") ;
  ("-i", Arg.Float (argset xinit), "largest x value") ;
  ("-a", Arg.Float (argset alpha), "alpha value");
  ("-b", Arg.Float (argset beta1), "beta value");
  ("-B", Arg.Float (argset beta2), "secondary beta value (asymmetric trees)") ;
  ("-e", Arg.Float (argset err), "percent error (randombeta trees)")
]

let argfail _ = Arg.usage argspec usage ; exit 1
let () = Arg.parse argspec (argset dist) usage

let req s z = match !z with Some x -> x
  | _ -> Printf.fprintf stderr "Error: %s\n%!" s ; argfail () 
let def d z = match !z with Some x -> x | _ -> d

(* Initialize random number generator with a known value for deterministic
   pseudo-randomness *)
let seed = Nativeint.of_int 31337
let rng = Gsl.Rng.make (Gsl.Rng.default ())
let () = Gsl.Rng.set rng seed
let pr = Printf.printf

(* prophylactic fold inclusion *)
let rec fold kons knil = function
    [] -> knil
  | h :: t -> fold kons (kons h knil) t

let rec symmetric_dichotomous beta xmin x =
  if x < xmin then () else (
  pr "%f\n%!" x ;
  symmetric_dichotomous beta xmin (x *. beta) ;
  symmetric_dichotomous beta xmin (x *. beta))

let rec asymmetric_dichotomous beta1 beta2 xmin x = 
  if x < xmin then () else (
  pr "%f\n%!" x ;
  asymmetric_dichotomous beta1 beta2 xmin (x *. beta1) ;
  asymmetric_dichotomous beta1 beta2 xmin (x *. beta2))

let rec random_asymmetry_dichotomous alpha xmin x =
  if x < xmin then () else (
  pr "%f\n%!" x ;
  let beta1 = Gsl.Rng.uniform rng in
  let beta2 = (1. -. (beta1) ** alpha) ** (1./.alpha) in
  random_asymmetry_dichotomous alpha xmin (x *. beta1) ;
  random_asymmetry_dichotomous alpha xmin (x *. beta2))

let rec random_dbeta_dichotomous alpha err xmin x =
  if x < xmin then () else (
  pr "%f\n%!" x ;
  let beta1 = Gsl.Randist.lognormal ~zeta:(log (2.**(-.1./.alpha)))
    ~sigma:(log (1. +. err /. 100.)) rng in
  (* Seems roughly equivalent: *)
  (* let beta2 = (1. -. (beta1) ** alpha) ** (1./.alpha) in *)
  let beta2 = Gsl.Randist.lognormal ~zeta:(log (2.**(-.1./.alpha)))
    ~sigma:(log (1. +. err /. 100.)) rng in
  random_dbeta_dichotomous alpha err xmin (x *. beta1) ;
  random_dbeta_dichotomous alpha err xmin (x *. beta2))

let rec rec_n n f i = if n = 0 then i else rec_n (n - 1) f (f i)

let rec random_multifurcating alpha xmin x =
  if x < xmin then () else (
  pr "%f\n%!" x ;
  let nn = 1 + Gsl.Randist.geometric rng ~p:0.5 in
  let beta = float_of_int nn ** (-. 1. /. alpha) in
  rec_n nn (fun () -> 
    random_multifurcating alpha xmin (x *. beta)) ())

let () = match req "Type must be one of symmetric or asymmetric" dist with
    "symmetric" -> symmetric_dichotomous
      (req "beta value required" beta1)
      (req "xmin value required" xmin)
      (req "initial x value required" xinit)
  | "asymmetric" -> asymmetric_dichotomous
      (req "beta value required" beta1)
      (req "secondary beta value required" beta2)
      (req "xmin value required" xmin)
      (req "initial x value required" xinit)
  | "random" -> random_asymmetry_dichotomous
      (req "alpha value required" alpha)
      (req "xmin value required" xmin)
      (req "initial x value required" xinit)
  | "randombeta" -> random_dbeta_dichotomous
      (req "alpha value required" alpha)
      (req "err value required" err)
      (req "xmin value required" xmin)
      (req "initial x value required" xinit)
  | "multifurcating" -> random_multifurcating
      (req "alpha value required" alpha)
      (req "xmin value required" xmin)
      (req "initial x value required" xinit)
  | _ -> Printf.fprintf stderr "Unknown or no argument.\n%s\n%!" usage
