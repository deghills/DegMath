open System
open DegMath
open Clifford

let testMultivectorPrint(m:Multivector) =
    Multivector.print m

let testMultivectorOfBlade =
    Multivector.ofBlade >> Multivector.print

let testMultivectorOfBlades : Blade seq -> unit =
    Multivector.ofBlades >> Multivector.print

let testGetBlade : byte -> Multivector -> unit =
    fun b m -> Multivector.getBlade b m |> Console.WriteLine

let testVectors(p, q, n) =
    let cl = Cl(p, q, n)
    let vecs =
        cl.Vectors
        |> Seq.map (
            function KeyValue(potency, vecs) ->
                        vecs 
                        |> Seq.map (fun v -> potency, v)) 
                        |> Seq.collect id
    for (potency, v) in vecs do printfn $"{(potency, v)}"

let testBivectors(p, q, n) =
    for bv in Cl(p, q, n).Bivectors do printfn "%A" bv

let testDual : (int*int*int) -> Multivector -> unit =
    fun (p, q, n) ->
        Cl(p, q, n).Dual >> Multivector.print

let testDualInv : (int*int*int) -> Multivector -> unit =
    fun (p, q, n) ->
        Cl(p, q, n).DualInv >> Multivector.print

let testReverse : Multivector -> unit =
    Cl(3, 3, 2).Reverse >> Multivector.print

let testGradeProject : int -> Multivector -> unit =
    fun i m -> Cl(3, 3, 2).GradeProject i m |> Multivector.print

let testGrade : Multivector -> unit =
    Cl(3, 3, 2).Grade
    >> fun grades -> 
        for g in grades do Console.WriteLine g

let testNeg : Multivector -> unit =
    Cl(3, 3, 2).Neg >> Multivector.print

let testAdd : Multivector -> Multivector -> unit =
    fun a b -> Cl(3, 3, 2).Add a b |> Multivector.print

let testSub : Multivector -> Multivector -> unit =
    fun a b -> Cl(3, 3, 2).Sub a b |> Multivector.print


let input1 = Multivector.ofBlades [
    0b10uy, 10f
    0b01uy, System.MathF.PI
    0uy, 1f
    0b11uy, -99f
    0b11001011uy, 69f
]

let input2 = Multivector.ofBlades [
    0b1101uy, 111f
    0b11001011uy, 1.5f
]


[<EntryPoint>]
let main _ =
    do
        input1 |> testSub input2


    0