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

    let getBit (b:byte) index = b >>> (index) &&& 1uy
    let getIndeces (b:byte) = seq {
        for i in 0..7 do
        if (getBit b i = 1uy) then yield i
    }

    let zero = 0uy, 0f

    type Cl(p, q, n) =
        let size = p+q+n

        let getPotency b =
            let rec loop acc i =
                match getBit b i with
                |1uy when i < p ->
                    loop acc (i+1)
                |1uy when i < q ->
                    loop (acc * -1f) (i+1)
                |1uy when i < n ->
                    loop (acc * 0f) (i+1)
                |0uy when i < size ->
                    loop acc (i+1)
                |_ -> acc
            loop 1f 0

        let signFromSquares a b =
            a &&& b |> getPotency

        // as bitvectors are already sorted, you only need a single mergesort iteration to count inversions
        let signFromSwaps (a :byte) (b :byte) =
            let rec checkInversion lhs rhs inversioncount =
                match lhs, rhs with
                | [], _ | _ , []    
                    -> inversioncount

                | x :: xs, y :: ys when (x > y) 
                    -> checkInversion (x :: xs) ys (inversioncount + xs.Length + 1)

                | _ :: xs, y :: ys 
                    -> checkInversion xs (y :: ys) inversioncount

            match (checkInversion (a |> getIndeces |> List.ofSeq) (b |> getIndeces |> List.ofSeq) 0) % 2 with
            | 0 -> 1f
            | 1 -> -1f
            | _ -> failwith "unexpected result from modulus operation"

        let sign a b = (signFromSquares a b) * (signFromSwaps a b)

        let bldGrade = getIndeces >> Seq.length
                    
        //XOR with 1111... flips every basis vector, getting the orthogonal complement
        let bldDual (bld, mag) =
            let rec buildRepunit acc n =
                if n = 1
                    then (acc ||| 1uy)
                    else 
                        let ndecr = n-1
                        buildRepunit (acc ||| (1uy <<< ndecr)) ndecr

            let bld' = (buildRepunit 0uy size) ^^^ bld
            let sign = signFromSwaps bld bld'
            bld', sign * mag

        let bldDualInv (bld, mag) =
            let rec buildRepunit acc n =
                if n = 1
                    then (acc ||| 1uy)
                    else 
                        let ndecr = n-1
                        buildRepunit (acc ||| (1uy <<< ndecr)) ndecr
            let bld' = (buildRepunit 0uy size) ^^^ bld
            let sign = signFromSwaps bld' bld
            bld', sign * mag
        
        let bldReverse (bld, mag) = 
            match (bldGrade bld) with
            | 2 | 3 | 6 | 7 -> 
                bld, -mag

            | _ -> 
                bld, mag

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
                    -> cache |> Seq.map (fun kvp -> kvp.Key, kvp.Value) |> List.ofSeq

                | (_, 0f) :: tail
                    -> loop cache tail

                | (bl, mag) :: tail when (cache |> Map.containsKey bl)
                    -> loop (cache |> Map.add bl (cache.[bl] + mag)) tail

                | (bl, mag) :: tail
                    -> loop (cache |> Map.add bl mag) tail
            loop Map.empty m

        member _.Dual = List.map bldDual

        member _.DualInv = List.map bldDualInv

        member _.Reverse = List.map bldReverse

        member _.GradeProject (grade: int) = 
            List.choose (fun (b :Blade) -> if (b |> fst |> bldGrade = grade) then Some b else None)

        member _.getRealPart = 
            List.choose 
                (function | 0uy, (x :float32) -> Some x | _ -> None) 
                >> function | x when x.Length = 0 -> 0f | x -> List.reduce (+) x

        //returns the grade of the multivector tupled with a flag if the vector is of pure grade
        member _.Grade m = 
            let rec addToSet (mv :Multivector) set =
                match mv with
                | []            -> set, (set.Count = 1) :int Set * bool
                | (h, _) :: t   -> addToSet t (Set.add (bldGrade h) set)
            addToSet m Set.empty

        //binary operators have arguments swapped so they can be used infix e.g. mul b a = a |> mul b = ab

        member _.Add a b = simplify (a @ b)

        member _.Neg = List.map (fun ((bl, mag) :Blade) -> bl, -mag)

        member this.Sub b a = this.Add a (this.Neg b)

        //geometric product
        member _.Mul b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldProduct blade1 blade2)) |> simplify

        //inner product
        member _.Dot b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldInner blade1 blade2)) |> simplify

        //exterior product
        member _.Wdg b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldOuter blade1 blade2)) |> simplify

        //regressive product
        member _.Reg b a = 
            a
            |> List.collect (fun blade1 -> 
                b |> List.map(fun blade2 -> bldRegress blade1 blade2)) |> simplify

        member this.MagSqr m = this.Mul m (this.Reverse m) |> this.getRealPart

        member this.Mag = this.MagSqr >> MathF.Sqrt

        member this.Scale (s: float32) (m: Multivector) = this.Mul [0uy, s] m

        member this.ScaleInv s m = this.Scale (1f/s) m      

        //returns mhat and as well as |m|
        //zero divisors can't be normalized and return the input
        member this.Normalize m = 
            match this.Mag m with
            |0f -> m, 0f
            |s  -> this.ScaleInv s m, s

        member this.VersorInv m = m |> this.Reverse |> this.ScaleInv (this.MagSqr m)

        member this.Project b a = a |> this.Dot b |> this.Mul (this.VersorInv b)

        member this.Sandwich a b = a |> this.Mul b |> this.Mul (this.VersorInv a) // RVR^-1