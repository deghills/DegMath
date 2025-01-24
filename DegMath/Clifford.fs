namespace DegMath

open System

(*
    Main resources I learned from:
    
    Youtube:
    https://youtu.be/htYh-Tq7ZBI?feature=shared
    https://youtu.be/60z_hpEAtD8?feature=shared
    https://youtu.be/0i3ocLhbxJ4?feature=shared
    https://youtu.be/2AKt6adG_OI?feature=shared
    https://youtu.be/Idlv83CxP-8?feature=shared
    https://youtube.com/playlist?list=PLsSPBzvBkYjxrsTOr0KLDilkZaw7UE2Vc&feature=shared
    https://www.youtube.com/playlist?list=PLpzmRsG7u_gqaTo_vEseQ7U8KFvtiJY4K

    Documentation:
    https://bivector.net/PGADYN.html
    https://geometry.mrao.cam.ac.uk/2016/10/geometric-algebra-2016/ (lecture 7)
*)

(*
    Note:

    Still haven't really found an elegant basis-independent implementation of exp and ln for bivectors. It's hard to
    implement partial functions, especially as the bitvector-blade representation is essentially encoding a
    dynamic type system, so you have to have run-time type checking/errors or just have the API allow undefined behaviour.
    Even my implementation of the inverse makes me uncomfy as lots of multivector are non-invertible. Idk I kept meaning
    to get around to it but have found it pretty easy to implement these functions on a case-by-case basis, so I will
    probably keep it that way.
*)

module Clifford =

    type Blade = byte*float32
    type Multivector = Blade list

    let getBitFromIndex b i = b >>> (i-1) &&& 1uy
    let getIndeces b = 
        let rec loop n acc =
            if n = 1 
                then
                    if (getBitFromIndex b 1 = 1uy) 
                        then 1 :: acc 
                        else acc
                else 
                    if (getBitFromIndex b n = 1uy) 
                        then loop (n - 1) (n :: acc)
                        else loop (n - 1) (acc)
        loop 8 []

    let zero = 0uy, 0f

    type CliffordAlgebra =
        {
            Signature :int * int * int
            Size :int

            GradeOfMultivector  :Multivector -> int Set
            GradeOfBlade        :Blade -> int
            MagSqr              :Multivector -> float32
            Mag                 :Multivector -> float32
            Scale               :float32 -> Multivector -> Multivector
            ScaleInv            :float32 -> Multivector -> Multivector
            Normalize           :Multivector -> Multivector*float32

            Dual            :Multivector -> Multivector
            DualInv         :Multivector -> Multivector
            Negate          :Multivector -> Multivector
            Reverse         :Multivector -> Multivector
            VersorInverse   :Multivector -> Multivector
            GradeProject    :int -> Multivector -> Multivector

            Add :Multivector -> Multivector -> Multivector
            Sub :Multivector -> Multivector -> Multivector
            Mul :Multivector -> Multivector -> Multivector
            Dot :Multivector -> Multivector -> Multivector
            Wdg :Multivector -> Multivector -> Multivector
            Reg :Multivector -> Multivector -> Multivector

            Project     :Multivector -> Multivector -> Multivector 
            Sandwich    :Multivector -> Multivector -> Multivector
        }

    let Cl (p, q, n) =
        let size = p+q+n

        if (p < 0 || q < 0 || n < 0) then failwith "SIGNATURE CANNOT HAVE NEGATIVE VALUES" 
        elif (size > 8) then failwith "the maximum size for a clifford algebra signature is 8" else

        let isUnipotent i = i <= p

        let isAntiUnipotent i = p < i && i <= p + q

        let isNilpotent i = i > p+q

        let getSquareOfEi = function
            | i when (i |> isUnipotent)     -> 1f
            | i when (i |> isAntiUnipotent) -> -1f
            | i when (i |> isNilpotent)     -> 0f
            | _ -> failwith "input is outside of this clifford algebra"

        let signFromSquares (a :byte) (b :byte) = 
            a &&& b
            |> getIndeces
            |> List.map getSquareOfEi
            |> function
                | [] -> 1f
                | x -> List.reduce (*) x

        // as bitvectors are already sorted, you only need a single mergesort iteration to count inversions
        let signFromSwaps (a :byte) (b :byte) =
            let rec checkInversion (lhs :int list) (rhs :int list) (inversioncount: int) =
                match lhs, rhs with
                | [], _ | _ , []    
                    -> inversioncount

                | xh :: xt, yh :: yt when (xh > yh) 
                    -> checkInversion (xh :: xt) yt (inversioncount + xt.Length + 1)

                | _ :: xt, yh :: yt 
                    -> checkInversion xt (yh :: yt) inversioncount

            match (checkInversion (getIndeces a) (getIndeces b) 0) % 2 with
            | 0 -> 1f
            | 1 -> -1f
            | _ -> failwith "unexpected result from modulus operation"

        let sign a b = (signFromSquares a b) * (signFromSwaps a b)

        let bldGrade = getIndeces >> List.length


        //XOR with 1111... flips every basis vector, getting the orthogonal complement
        let bldDual (b, mag) =
            let rec buildRepunit acc n =
                if n = 1
                    then (acc ||| 1uy)
                    else 
                        let ndecr = n-1
                        buildRepunit (acc ||| (1uy <<< ndecr)) ndecr

            let b' = (buildRepunit 0uy size) ^^^ b
            let sign = signFromSwaps b b'
            b', sign * mag

        let bldDualInv (b, mag) =
            let rec buildRepunit acc n =
                if n = 1
                    then (acc ||| 1uy)
                    else 
                        let ndecr = n-1
                        buildRepunit (acc ||| (1uy <<< ndecr)) ndecr
            let b' = (buildRepunit 0uy size) ^^^ b
            let sign = signFromSwaps b' b
            b', sign * mag
        
        let bldReverse (b, mag) = 
            match (bldGrade b) with
            | 2 | 3 | 6 | 7 -> 
                b, -mag

            | _ -> 
                b, mag

        let bldProduct (bld1, mag1) (bld2, mag2) =
            let bld3 = bld1 ^^^ bld2
            let mag3 = mag1 * mag2 * (sign bld1 bld2)
            bld3, mag3
        
        let bldOuter (bld1, mag1) (bld2, mag2) =
            let areOrthogonal = (bld1 &&& bld2) = 0uy
            if areOrthogonal
                then bldProduct (bld1, mag1)(bld2, mag2)
                else zero
                    
        let bldInner (bld1, mag1) (bld2, mag2) =
            let isSubsetOf set potentialSubset = (set ||| potentialSubset) = set
            match bld1, bld2 with
            | x, y when (x |> isSubsetOf y) -> bldProduct (x, mag1) (y, mag2) //left contraction
            | x, y when (y |> isSubsetOf x) -> bldProduct (x, mag1) (y, mag2) //right contraction
            | _ -> zero

        let bldRegress a b = 
            bldDualInv (bldOuter (bldDual a) (bldDual b))

        let simplify m = 
            let rec loop (cache :Map<byte,float32>) (m :Multivector) =
                match m with
                | [] 
                    -> cache |> List.ofSeq |> List.map (fun kvp -> kvp.Key, kvp.Value)

                | (_, 0f) :: tail
                    -> loop cache tail

                | (bl, mag) :: tail when (cache |> Map.containsKey bl)
                    -> loop (cache |> Map.add bl (cache.[bl] + mag)) tail

                | (bl, mag) :: tail
                    -> loop (cache |> Map.add bl mag) tail
            loop Map.empty m

        let dual = List.map bldDual

        let dualInv = List.map bldDualInv

        let rev = List.map bldReverse

        let gradeProject (grade: int) = 
            List.choose (fun (b :Blade) -> if (b |> fst |> bldGrade = grade) then Some b else None)

        let getRealPart = 
            List.choose 
                (function | 0uy, (x :float32) -> Some x | _ -> None) 
                >> function | x when x.Length = 0 -> 0f | x -> List.reduce (+) x

        //returns the grade of the multivector tupled with a flag if the vector is of pure grade
        let grade m = 
            let rec addToSet (mv :Multivector) set =
                match mv with
                | []            -> set, (set.Count = 1) :int Set * bool
                | (h, _) :: t   -> addToSet t (Set.add (bldGrade h) set)
            addToSet m Set.empty
       
        //binary operators have arguments swapped so they can be used infix e.g. mul b a = a |> mul b = ab
        let add a b = simplify (a @ b)

        let neg = List.map (fun ((bl, mag) :Blade) -> bl, -mag)

        let sub b a = add a (neg b)

        //geometric product
        let mul b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldProduct blade1 blade2)) |> simplify

        //inner product
        let dot b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldInner blade1 blade2)) |> simplify

        //exterior product
        let wdg b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldOuter blade1 blade2)) |> simplify

        //regressive product
        let reg b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldRegress blade1 blade2)) |> simplify

        let magSqr m = mul m (rev m) |> getRealPart

        let mag = magSqr >> MathF.Sqrt

        let scale (s: float32) (m: Multivector) = mul [0uy, s] m

        let scaleInv s m = scale (1f/s) m

        //returns mhat and as well as |m|
        //zero divisors can't be normalized and return the input
        let normalize m = 
            match mag m with
            |0f -> m, 0f
            |s  -> scaleInv s m, s

        //technically the pseudo-inverse; only works for versors
        let inv m = m |> rev |> scaleInv (magSqr m)

        let project b a = a |> dot b |> mul (inv b)

        let sandwich a b = a |> mul b |> mul (inv a) // RVR^-1

        {
            Signature = p,q,n
            Size = size

            GradeOfMultivector = grade >> fst
            GradeOfBlade = fst >> bldGrade
            MagSqr = magSqr
            Mag = mag
            Scale = scale
            ScaleInv = scaleInv
            Normalize = normalize

            Negate  = neg
            Reverse = rev
            Dual    = dual
            DualInv = dualInv

            VersorInverse = inv

            GradeProject = gradeProject

            Add = add
            Sub = sub

            Mul = mul
            Dot = dot
            Wdg = wdg
            Reg = reg

            Project = project
            Sandwich = sandwich
        }