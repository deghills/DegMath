namespace DegMath

open System

[<AutoOpen>]
module private BitvectorOperations =
    let getBit (b:byte) index = 
        (b >>> index) &&& 1uy

    let getIndeces (b:byte) = [
        for i in 0..7 do
        if (getBit b i = 1uy) then yield i
    ]

    let byteToBitString :byte -> string =
        fun b -> 
            getIndeces b
            |> List.map (fun i -> i+1)
            |> List.map (fun i -> i.ToString())
            |> function | [] -> ["real"] | x -> "e" :: x
            |> List.reduce (+)

module Clifford =

    type Blade = byte*float32
    type Multivector = Map<byte, float32>
    
    [<RequireQualifiedAccess>]
    module Multivector =
        let getBlade : byte -> Multivector -> float32 =
            fun b m ->
                match Map.tryFind b m with
                |Some x -> x
                |None -> 0f

        let print : Multivector -> unit =
            fun m ->
                for KeyValue(bld, mag) in m do
                    printfn $"{(byteToBitString bld, mag)}"

        let scale : float32 -> Multivector -> Multivector = 
            fun s -> Map.map (fun _ mag -> s * mag)

        let scaleInv : float32 -> Multivector -> Multivector =
            fun s -> Map.map (fun _ mag -> (1f/s) * mag)

        let zero = Multivector[]

    let multivector : Blade seq -> Multivector =
        Seq.fold (
            fun (m:Multivector) (bld, mag) -> 
                m.Add (bld, (mag + Multivector.getBlade bld m)
            )
        ) Multivector.zero
        >> Map.filter (fun _ mag -> mag <> 0f)

    type Cl(p, q, n) =
        let size = 
            match p+q+n with
            |x when x > 8 -> failwith "the maximum size for a clifford algebra signature is 8"
            |x -> x

        let bldZero = 0uy, 0f

        let getPotency b =
            let rec loop acc i =
                match getBit b i with
                |1uy when i < p ->
                    loop acc (i+1)
                |1uy when i < p+q ->
                    loop (acc * -1f) (i+1)
                |1uy when i < p+q+n ->
                    loop (acc * 0f) (i+1)
                |0uy when i < size ->
                    loop acc (i+1)
                |_ -> acc
            loop 1f 0

        let signFromSquares : byte -> byte -> float32 =
            fun a b ->
                getPotency(a &&& b)

        // as bitvectors are already sorted, you only need a single mergesort iteration to count inversions
        let signFromSwaps : byte -> byte -> float32 =
            fun a b ->
                let rec checkInversion lhs rhs inversioncount =
                    match lhs, rhs with
                    | [], _ | _ , []    
                        -> inversioncount

                    | x :: xs, y :: ys when (x > y) 
                        -> checkInversion (x :: xs) ys (inversioncount + xs.Length + 1)

                    | _ :: xs, y :: ys 
                        -> checkInversion xs (y :: ys) inversioncount

                match (checkInversion (getIndeces a) (getIndeces b) 0) % 2 with
                | 0 -> 1f
                | 1 -> -1f
                | _ -> failwith "unexpected result from modulus operation"

        let sign a b = (signFromSquares a b) * (signFromSwaps a b)

        let bldGrade = getIndeces >> List.length

        //XOR with 1111... flips every basis vector, getting the orthogonal complement
        let bldDual (bld, mag) =
            let rec buildRepunit acc n =
                if n = 1 then 
                    (acc ||| 1uy)
                else 
                    let ndecr = n-1
                    buildRepunit (acc ||| (1uy <<< ndecr)) ndecr

            let bld' = (buildRepunit 0uy size) ^^^ bld
            let sign = signFromSwaps bld bld'
            bld', sign * mag

        let bldDualInv (bld, mag) =
            let rec buildRepunit acc n =
                if n = 1 then 
                    (acc ||| 1uy)
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
                else bldZero

        let bldInner (bld1, mag1) (bld2, mag2) =
            let isSubsetOf set potentialSubset = (set ||| potentialSubset) = set
            match bld1, bld2 with
            | x, y when (x |> isSubsetOf y) -> bldProduct (x, mag1) (y, mag2) //left contraction
            | x, y when (y |> isSubsetOf x) -> bldProduct (x, mag1) (y, mag2) //right contraction
            | _ -> bldZero

        let bldRegress a b = 
            bldDualInv (bldOuter (bldDual a) (bldDual b))

        let mergeLinear f a b =
            let allKeys = Set.union (Set <| Map.keys a) (Set <| Map.keys b)
            Map [|
                for key in allKeys ->
                    let find = Map.tryFind key
                    match find a, find b with
                    |Some x, Some y -> key, f x y
                    |Some x, None -> key, f x 0f
                    |None, Some y -> key, f 0f y
                    |_ -> failwith "unexpected key"
            |] |> Map.filter (fun _ mag -> mag <> 0f)

        let mergeBilinear : (Blade -> Blade -> Blade) -> Multivector -> Multivector -> Multivector =
            fun f b a ->
                multivector [|
                    for KeyValue bladea in a do
                    for KeyValue bladeb in b do
                        yield f bladea bladeb
                |]
        
        ///all basis vectors in the clifford algebra, sorted by their potency
        member _.Vectors = Map [
            1,  [| for i in 0..(p-1) -> 1uy <<< i |]
            -1, [| for i in p..(p+q-1) -> 1uy <<< i |]
            0,  [| for i in p+q..(p+q+n-1) -> 1uy <<< i |]
        ]

        ///all basis bivectors in the clifford algebra
        member this.Bivectors = 
            let vs = [|
                for KeyValue (_, v) in this.Vectors do
                    yield! v
            |]
            let f arr = [|
                let length = Array.length arr
                for i in 0..(length-1) do
                    for j in (i+1)..(length-1) ->
                        arr[i] ^^^ arr[j]
            |]
            f vs 

        ///orthogonal complement
        member _.Dual = Map.toSeq >> Seq.map bldDual >> Map.ofSeq

        ///Dual >> DualInv = id
        member _.DualInv = Map.toSeq >> Seq.map bldDualInv >> Map.ofSeq

        ///negates all blades of grades 2, 3, 6, and 7
        member _.Reverse = Map.map (
            fun bld mag ->
                match bldGrade bld with
                |2 |3 |6| 7 -> -mag
                |_ -> mag
        )

        ///removes every blade not of the specified grade
        member _.GradeProject (grade: int) =
            Map.filter (fun bld _ -> 
                bldGrade bld = grade)                        

        ///returns the grade of the multivector
        member _.Grade m = 
            Set [
                for KeyValue(bld, _) in m -> bldGrade bld 
            ]

        ///negates every component of the multivector
        member _.Neg : Multivector -> Multivector =
            Map.map (fun _ mag -> -mag)

        //binary operators have arguments swapped so they can be used infix e.g. mul b a = a |> mul b = ab
        ///Sum
        member _.Add : Multivector -> Multivector -> Multivector =
            fun b a ->
                mergeLinear (+) a b
        
        ///Difference
        member _.Sub : Multivector -> Multivector -> Multivector =
            fun b a ->             
                mergeLinear (-) a b

        ///Geometric product
        member _.Mul : Multivector -> Multivector -> Multivector = 
            mergeBilinear bldProduct

        ///Inner product
        member _.Dot : Multivector -> Multivector -> Multivector =
            mergeBilinear bldInner

        ///Outer product
        member _.Wdg : Multivector -> Multivector -> Multivector =
            mergeBilinear bldOuter

        ///Regressive product
        member _.Reg : Multivector -> Multivector -> Multivector =
            mergeBilinear bldRegress

        ///Squared magnitude
        member this.MagSqr : Multivector -> float32 =
            fun m -> 
                this.Mul m (this.Reverse m) |> Multivector.getBlade 0uy

        ///Magnitude
        member this.Mag : Multivector -> float32 = 
            this.MagSqr >> MathF.Sqrt

        ///Returns mhat and as well as |m|;
        ///Zero divisors can't be normalized and return the input
        member this.Normalize : Multivector -> (Multivector*float32) =
            fun m -> 
                match this.Mag m with
                |0f -> m, 0f
                |s  -> Multivector.scaleInv s m, s

        ///Calculates the inverse of versors;
        ///Behaviour is undefined for arbitrary multivectors
        member this.VersorInv : Multivector -> Multivector =
            fun m -> 
                m 
                |> this.Reverse 
                |> Multivector.scaleInv (this.MagSqr m)

        ///Project a onto b (exclude mutually orthogonal parts)
        member this.Project : Multivector -> Multivector -> Multivector =
            fun b a -> 
                a 
                |> this.Dot b 
                |> this.Mul (this.VersorInv b)

        ///Sandwich product;
        ///First argument should be a versor;
        ///RVR^-1
        member this.Sandwich : Multivector -> Multivector -> Multivector =
            fun a b -> 
                a 
                |> this.Mul b 
                |> this.Mul (this.VersorInv a)